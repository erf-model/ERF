#include <ERFPC.H>
#include <ERFPCParticleToMesh.H>
#include <ERF_Constants.H>
#include <ERF_TerrainConversion.H>
#include <AMReX_ParticleLocator.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

namespace {
    struct ERFPCLevelGeom {
        GpuArray<Real,AMREX_SPACEDIM> plo;
        GpuArray<Real,AMREX_SPACEDIM> dxi;
        int k_max;
        int ref_ratio;  // cumulative vertical refinement ratio
    };
}

void ERFPC::massDensity ( MultiFab&        a_mf,
                          const MultiFab&  a_z_phys_nd,
                          const int&       a_lev,
                          const int&       a_comp ) const
{
    BL_PROFILE("ERFPC::massDensity()");
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const ERFPC::ParticleTileType::ConstParticleTileDataType& ptd, int i) {
            return ptd.m_rdata[ERFParticlesRealIdx::mass][i];
        });
}

/*! \brief Diagnostic: count particles per level and, for level 1, in the halo.
 *
 * Prints particle counts for all levels. For AMR with >= 3 levels, also
 * counts level-1 particles in halo cells (not covered by level 2).
 *
 * \param[in] a_finest_level Current finest AMR level (>= 0)
 */
void ERFPC::CountParticlesPerLevelAndHalo (int a_finest_level)
{
    BL_PROFILE("ERFPC::CountParticlesPerLevelAndHalo()");

    // Per-level total counts
    amrex::Vector<Long> n_per_lev(a_finest_level + 1, 0);
    for (int lev = 0; lev <= a_finest_level; lev++) {
        n_per_lev[lev] = NumberOfParticlesAtLevel(lev, true, true);
    }
    for (int lev = 0; lev <= a_finest_level; lev++) {
        ParallelDescriptor::ReduceLongSum(n_per_lev[lev]);
    }

    amrex::Print() << "[" << m_name << "] Particle counts per level (after Redistribute): ";
    for (int lev = 0; lev <= a_finest_level; lev++) {
        amrex::Print() << "L" << lev << "=" << n_per_lev[lev];
        if (lev < a_finest_level) { amrex::Print() << " "; }
    }
    amrex::Print() << "\n";

    if (a_finest_level < 2) { return; }

    // Level-1 halo cells not covered by level 2 (1=halo, 0=covered)
    const int lev1 = 1;
    iMultiFab halo_mask = amrex::makeFineMask(ParticleBoxArray(lev1), ParticleDistributionMap(lev1),
                                               IntVect(0), ParticleBoxArray(lev1 + 1), m_gdb->refRatio(lev1),
                                               m_gdb->Geom(lev1).periodicity(), 1, 0);

    const amrex::Geometry& geom1 = m_gdb->Geom(lev1);
    const auto plo = geom1.ProbLoArray();
    const auto dxi = geom1.InvCellSizeArray();
    const amrex::Box& domain1 = geom1.Domain();

    Long n_halo = 0;
    for (ParIterType pti(*this, lev1); pti.isValid(); ++pti) {
        const auto& aos = ParticlesAt(lev1, pti).GetArrayOfStructs();
        const int np = aos.numParticles();
        if (np == 0) { continue; }

        const amrex::Box& box = pti.validbox();
        amrex::Array4<int const> const& mask_arr = halo_mask.const_array(pti);
        auto* p_pbox = aos().data();

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(np, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
            {
                auto const& p = p_pbox[i];
                if (p.id() <= 0) { return {0L}; }
                int ix = static_cast<int>(amrex::Math::floor((p.pos(0) - plo[0]) * dxi[0]))
                         + domain1.smallEnd(0);
                int iy = static_cast<int>(amrex::Math::floor((p.pos(1) - plo[1]) * dxi[1]))
                         + domain1.smallEnd(1);
                int iz = static_cast<int>(amrex::Math::floor((p.pos(AMREX_SPACEDIM-1) - plo[AMREX_SPACEDIM-1])
                                                              * dxi[AMREX_SPACEDIM-1]))
                         + domain1.smallEnd(2);
                iz = amrex::max(domain1.smallEnd(2), amrex::min(iz, domain1.bigEnd(2)));
                amrex::IntVect iv(AMREX_D_DECL(ix, iy, iz));
                if (box.contains(iv) && mask_arr(iv) == 1) {
                    return {1L};
                }
                return {0L};
            });
        n_halo += amrex::get<0>(reduce_data.value(reduce_op));
    }
    ParallelDescriptor::ReduceLongSum(n_halo);

    amrex::Print() << "[" << m_name << "] Level-1 halo particle count (cells not covered by L2): " << n_halo
                   << " (L1 total=" << n_per_lev[1] << ")\n";
}

void ERFPC::ExtractAndRouteOORParticles (int a_lev)
{
    BL_PROFILE("ERFPC::ExtractAndRouteOORParticles()");
    AMREX_ALWAYS_ASSERT(a_lev > 0);

    // GPU spatial index for the fine level (used to detect OOR particles).
    ParticleLocator<DenseBins<Box>> locator_fine;
    locator_fine.build(ParticleBoxArray(a_lev), m_gdb->Geom(a_lev));
    auto ag_fine = locator_fine.getGridAssignor();

    // Find a local L0 grid index to receive routed particles.  Particles that
    // leave the fine level always fall to L0 (which covers the full domain).
    int dest_grid_L0 = -1;
    {
        const int my_proc = ParallelDescriptor::MyProc();
        const auto& dm0 = ParticleDistributionMap(0);
        for (int i = 0; i < static_cast<int>(dm0.size()); i++) {
            if (dm0[i] == my_proc) { dest_grid_L0 = i; break; }
        }
    }
    DefaultAssignor cell_assignor;

    if (dest_grid_L0 >= 0) {
      for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
        auto& src_tile = ParticlesAt(a_lev, pti);
        auto& aos = src_tile.GetArrayOfStructs();
        auto* p_pbox = aos().data();
        const int np = aos.numParticles();
        if (np == 0) { continue; }

        Gpu::DeviceVector<int> oor_mask(np);
        auto* mask_ptr = oor_mask.data();
        ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
        {
            auto& p = p_pbox[i];
            if (p.id() <= 0) { mask_ptr[i] = 0; return; }
            const int grd = ag_fine(p, 0, cell_assignor).first;
            mask_ptr[i] = (grd < 0) ? 1 : 0;
        });

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using TT = typename decltype(reduce_data)::Type;
        reduce_op.eval(np, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> TT { return {mask_ptr[i]}; });
        const int n_oor = amrex::get<0>(reduce_data.value(reduce_op));
        if (n_oor == 0) { continue; }

        ParticleTileType tmp_tile;
        tmp_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
        tmp_tile.resize(n_oor);
        [[maybe_unused]] int nc = amrex::filterParticles(tmp_tile, src_tile, mask_ptr,
                                                         int(0), int(0), np);
        AMREX_ASSERT(nc == n_oor);

        auto& dst_tile = DefineAndReturnParticleTile(0, dest_grid_L0, 0);
        const auto dst_start = static_cast<int>(dst_tile.numParticles());
        dst_tile.resize(dst_start + n_oor);
        amrex::copyParticles(dst_tile, tmp_tile, int(0), dst_start, n_oor);

        ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
            if (mask_ptr[i]) { p_pbox[i].id() = -1; }
        });
        Gpu::synchronize();
      }
    }

    Redistribute(0, 0);
    Redistribute(a_lev, a_lev);
}

void ERFPC::ConvertZetaToZ (const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd)
{
    BL_PROFILE("ERFPC::ConvertZetaToZ()");
    const int finest = finestLevel();
    for (int lev = 0; lev <= finest; lev++) {
        if (!a_z_phys_nd[lev]) { continue; }
        const Geometry& geom = m_gdb->Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int grid = pti.index();
            auto& aos = ParticlesAt(lev, pti).GetArrayOfStructs();
            const int np = aos.numParticles();
            if (np == 0) { continue; }
            auto* p_pbox = aos().data();
            auto zheight = (*a_z_phys_nd[lev])[grid].array();
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                p.pos(AMREX_SPACEDIM-1) = static_cast<ParticleReal>(
                    ERF::ParticlePos::z_from_zeta(p.pos(0), p.pos(1),
                                                  p.pos(AMREX_SPACEDIM-1),
                                                  plo, dxi, zheight));
            });
        }
    }
}

void ERFPC::ConvertZToZeta (const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd)
{
    BL_PROFILE("ERFPC::ConvertZToZeta()");
    const int finest = finestLevel();
    for (int lev = 0; lev <= finest; lev++) {
        if (!a_z_phys_nd[lev]) { continue; }
        const Geometry& geom = m_gdb->Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        const Box& dom = geom.Domain();
        const int k_max = dom.bigEnd(AMREX_SPACEDIM-1) - dom.smallEnd(AMREX_SPACEDIM-1);
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int grid = pti.index();
            auto& aos = ParticlesAt(lev, pti).GetArrayOfStructs();
            const int np = aos.numParticles();
            if (np == 0) { continue; }
            auto* p_pbox = aos().data();
            auto zheight = (*a_z_phys_nd[lev])[grid].array();
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                p.pos(AMREX_SPACEDIM-1) = static_cast<ParticleReal>(
                    ERF::ParticlePos::zeta_from_z(p.pos(0), p.pos(1),
                                                  p.pos(AMREX_SPACEDIM-1),
                                                  plo, dxi, zheight, k_max));
            });
        }
    }
}

#endif

