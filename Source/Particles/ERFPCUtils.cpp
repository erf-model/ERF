#include <ERFPC.H>
#include <ERFPCParticleToMesh.H>
#include <ERF_Constants.H>
#include <AMReX_ParticleLocator.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

// Static Assignor state — initialized to a negative sentinel so the Assignor
// falls back to using idata(k) directly until RefreshAssignorState() is called.
AMREX_GPU_MANAGED amrex::Real ERFParticlesAssignor::s_dxi_z_finest = amrex::Real(-1);

void ERFPC::RefreshAssignorState ()
{
    const Geometry& finest_geom = m_gdb->Geom(finestLevel());
    ERFParticlesAssignor::s_dxi_z_finest =
        finest_geom.InvCellSizeArray()[AMREX_SPACEDIM-1];
}

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
            return ptd.m_rdata[ERFParticlesRealIdxSoA::mass][i];
        });
}

/*! \brief Fix k-indices for all particles after AMR regrid.
 *
 *  Two-slot convention: idata(k) is the cell index at the particle's current
 *  AMR level (read by AMReX cic_/mac_interpolate_mapped_z), and idata(k_finest)
 *  is the cell index at the finest AMR level (used by ERFParticlesAssignor to
 *  scale to any queried level via dxi).  This decouples the Redistribute
 *  placement requirement from the kernel stencil requirement, so a single
 *  multi-level Redistribute places particles correctly without the
 *  reset/recompute cycles the single-slot convention required.
 *
 *  Steps:
 *   1. Refresh static Assignor state (s_dxi_z_finest).
 *   2. Set idata(k_finest) for every particle from its z-position using the
 *      finest-level geometry.
 *   3. Full multi-level Redistribute — Assignor scales idata(k_finest) to
 *      each queried level's cell, so particles land on the correct level.
 *   4. Per-level: walk terrain heights via update_location_idata to set
 *      idata(k) for the level the particle actually lives on, and refresh
 *      idata(k_finest) = level_k * ref_to_finest_at_lev.
 */
void ERFPC::FixKIndexAMR (const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd)
{
    BL_PROFILE("ERFPC::FixKIndexAMR()");

    // Step 1: refresh Assignor static state for the current finest level.
    RefreshAssignorState();

    const int finest = finestLevel();

    // Cumulative refinement ratio from level 0 to the finest level.
    int finest_ref = 1;
    for (int lev = 0; lev < finest; lev++) {
        finest_ref *= m_gdb->refRatio(lev)[AMREX_SPACEDIM-1];
    }

    // Z-levels for non-uniform vertical grids (level-0 cell interfaces)
    const Real* zlevels = m_zlevels_d.empty() ? nullptr : m_zlevels_d.data();
    const int nz_lev0 = m_zlevels_d.empty() ? 0
                       : static_cast<int>(m_zlevels_d.size()) - 1;

    // Step 2: set idata(k_finest) from each particle's z-position using
    // finest-level geometry.  This is regardless of which level the particle
    // currently lives on; Step 3's Redistribute will fix that.
    {
        const Geometry& geom_fine = m_gdb->Geom(finest);
        const auto plo_fine = geom_fine.ProbLoArray();
        const auto dxi_fine = geom_fine.InvCellSizeArray();
        const int k_max_fine = geom_fine.Domain().bigEnd(AMREX_SPACEDIM-1);

        const auto& particles = GetParticles();
        for (int lev = 0; lev <= finest; lev++) {
            if (lev >= static_cast<int>(particles.size())) { continue; }
            if (particles[lev].empty()) { continue; }

            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                auto& ptile = ParticlesAt(lev, pti);
                auto& aos = ptile.GetArrayOfStructs();
                auto* p_pbox = aos().data();
                const int np = aos.numParticles();

                ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                {
                    auto& p = p_pbox[i];
                    if (p.id() <= 0) { return; }
                    p.idata(ERFParticlesIntIdxAoS::k_finest) =
                        compute_k_from_z(Real(p.pos(AMREX_SPACEDIM-1)),
                                         plo_fine[AMREX_SPACEDIM-1],
                                         dxi_fine[AMREX_SPACEDIM-1],
                                         k_max_fine,
                                         zlevels, nz_lev0, finest_ref);
                });
            }
        }
        Gpu::synchronize();
    }

    // Step 3: full multi-level Redistribute.  Assignor uses idata(k_finest)
    // and scales by dxi to compute the correct cell at each queried level.
    Redistribute();

    // Step 4: per-level walk through terrain heights to set idata(k) for the
    // level the particle ended up on.  update_location_idata also refreshes
    // idata(k_finest) = level_k * ref_to_finest_at_lev.
    for (int lev = 0; lev <= finest; lev++) {
        const auto& particles = GetParticles();
        if (lev >= static_cast<int>(particles.size())) { continue; }
        if (particles[lev].empty()) { continue; }

        int ref_to_finest = 1;
        for (int l = lev; l < finest; l++) {
            ref_to_finest *= m_gdb->refRatio(l)[AMREX_SPACEDIM-1];
        }

        const Geometry& geom_lev = m_gdb->Geom(lev);
        const auto plo = geom_lev.ProbLoArray();
        const auto dxi = geom_lev.InvCellSizeArray();
        const int k_max_lev = geom_lev.Domain().bigEnd(AMREX_SPACEDIM-1);

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            int grid = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos = ptile.GetArrayOfStructs();
            auto* p_pbox = aos().data();
            const int np = aos.numParticles();

            auto zheight = (*a_z_phys_nd[lev])[grid].array();

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                auto& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                // Seed idata(k) with a coordinate-based level-k guess; the
                // height walk inside update_location_idata refines from there.
                int k_guess = static_cast<int>(amrex::Math::floor(
                    (Real(p.pos(AMREX_SPACEDIM-1)) - plo[AMREX_SPACEDIM-1])
                    * dxi[AMREX_SPACEDIM-1]));
                p.idata(ERFParticlesIntIdxAoS::k) =
                    amrex::max(0, amrex::min(k_guess, k_max_lev));
                update_location_idata(p, plo, dxi, zheight, ref_to_finest);
            });
        }
        Gpu::synchronize();
    }
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
                int iz = p.idata(ERFParticlesIntIdxAoS::k);
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

/*! \brief Extract OOR particles from level a_lev, route to finest covering level,
 *         then redistribute. Avoids multi-level redistribute and excessive MPI traffic.
 *
 * \param[in] a_lev       Fine level with potentially OOR particles
 * \param[in] a_z_phys_nd Terrain height data (all levels)
 *
 * idata(k_finest) is set ONCE per OOR particle from its z-position; the
 * per-level AssignGrid query then uses ERFParticlesAssignor which scales
 * idata(k_finest) to the queried level via dxi.  After routing, idata(k) is
 * recomputed from the target level's geometry.
 */
void ERFPC::ExtractAndRouteOORParticles ( int                                        a_lev,
                                          const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    BL_PROFILE("ERFPC::ExtractAndRouteOORParticles()");
    amrex::ignore_unused(a_z_phys_nd);

    AMREX_ALWAYS_ASSERT(a_lev > 0);

    // The Assignor relies on s_dxi_z_finest; ensure it matches the current
    // finest-level geometry before any AssignGrid query below.
    RefreshAssignorState();

    ERFParticlesAssignor cell_assignor;
    const int n_levels = a_lev + 1;

    // Z-levels data for k computation (level-0 cell interfaces)
    const Real* zlevels = m_zlevels_d.empty() ? nullptr
                        : m_zlevels_d.data();
    const int nz_lev0 = m_zlevels_d.empty() ? 0
                      : static_cast<int>(m_zlevels_d.size()) - 1;

    // Build per-level ParticleLocators (GPU spatial indices)
    Vector<ParticleLocator<DenseBins<Box>>> locators(n_levels);
    for (int lev = 0; lev < n_levels; lev++) {
        locators[lev].build(ParticleBoxArray(lev), m_gdb->Geom(lev));
    }

    // Collect the per-level AssignGrid objects into a device array
    using AGType = AssignGrid<DenseBinIteratorFactory<Box>>;
    Gpu::PinnedVector<AGType> h_assignors(n_levels);
    for (int lev = 0; lev < n_levels; lev++) {
        h_assignors[lev] = locators[lev].getGridAssignor();
    }
    Gpu::DeviceVector<AGType> d_assignors(n_levels);
    Gpu::copyAsync( Gpu::hostToDevice,
                    h_assignors.begin(),
                    h_assignors.end(),
                    d_assignors.begin());

    // Per-level geometry data for compute_k_from_z on the GPU
    Gpu::PinnedVector<ERFPCLevelGeom> h_lg(n_levels);
    for (int lev = 0; lev < n_levels; lev++) {
        const auto& geom = m_gdb->Geom(lev);
        h_lg[lev].plo = geom.ProbLoArray();
        h_lg[lev].dxi = geom.InvCellSizeArray();
        h_lg[lev].k_max = geom.Domain().bigEnd(AMREX_SPACEDIM-1);
        int ref = 1;
        for (int l = 0; l < lev; l++) {
            ref *= m_gdb->refRatio(l)[AMREX_SPACEDIM-1];
        }
        h_lg[lev].ref_ratio = ref;
    }
    Gpu::DeviceVector<ERFPCLevelGeom> d_lg(n_levels);
    Gpu::copyAsync(Gpu::hostToDevice, h_lg.begin(), h_lg.end(), d_lg.begin());
    Gpu::synchronize();

    auto* ag_ptr = d_assignors.data();

    // Cumulative refinement ratio from level 0 to the finest level (= a_lev here).
    int finest_ref = 1;
    for (int l = 0; l < a_lev; l++) {
        finest_ref *= m_gdb->refRatio(l)[AMREX_SPACEDIM-1];
    }
    const Geometry& geom_fine = m_gdb->Geom(a_lev);
    const Real plo_fine_z = geom_fine.ProbLoArray()[AMREX_SPACEDIM-1];
    const Real dxi_fine_z = geom_fine.InvCellSizeArray()[AMREX_SPACEDIM-1];
    const int  k_max_fine = geom_fine.Domain().bigEnd(AMREX_SPACEDIM-1);

    // Find one local grid per level for receiving extracted particles
    Gpu::PinnedVector<int> dest_grids(n_levels, 0);
    {
        const int my_proc = ParallelDescriptor::MyProc();
        for (int lev = 0; lev < n_levels; lev++) {
            const auto& dm = ParticleDistributionMap(lev);
            for (int i = 0; i < static_cast<int>(dm.size()); i++) {
                if (dm[i] == my_proc) { dest_grids[lev] = i; break; }
            }
        }
    }

    int finest = a_lev;
    Vector<int> levels_modified(n_levels, 0);
    Vector<Long> n_routed_to(n_levels, 0);

    {
        const int lev = a_lev;
        auto src_lev_grid = locators[lev].getGridAssignor();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            auto& src_tile = ParticlesAt(lev, pti);
            auto& aos = src_tile.GetArrayOfStructs();
            auto* p_pbox = aos().data();
            const int np = aos.numParticles();
            if (np == 0) { continue; }

            // Create OOR mask on GPU
            Gpu::DeviceVector<int> mask(np);
            auto* mask_ptr = mask.data();

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                auto& p = p_pbox[i];
                if (p.id() <= 0) { mask_ptr[i] = 0; return; }
                int grd = src_lev_grid(p, 0, cell_assignor).first;
                mask_ptr[i] = (grd < 0) ? 1 : 0;
            });

            // Count OOR particles
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(np, reduce_data,
                [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                { return {mask_ptr[i]}; });
            int n_oor = amrex::get<0>(reduce_data.value(reduce_op));
            if (n_oor == 0) { continue; }

            // Extract OOR particles to a temporary tile
            ParticleTileType tmp_tile;
            tmp_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
            tmp_tile.resize(n_oor);

            [[maybe_unused]] int n_copied = amrex::filterParticles(tmp_tile, src_tile, mask_ptr,
                                                                   int(0), int(0), np);
            AMREX_ASSERT(n_copied == n_oor);

            // Find the correct target level for each particle.  idata(k_finest)
            // is set ONCE from the position; the per-level Assignor scales it.
            Gpu::DeviceVector<int> target_lev_vec(n_oor);
            auto* tgt_ptr = target_lev_vec.data();
            auto* tmp_pbox = tmp_tile.GetArrayOfStructs()().data();

            ParallelFor(n_oor, [=] AMREX_GPU_DEVICE (int i)
            {
                auto& p = tmp_pbox[i];
                p.idata(ERFParticlesIntIdxAoS::k_finest) =
                    compute_k_from_z( Real(p.pos(AMREX_SPACEDIM-1)),
                                      plo_fine_z, dxi_fine_z,
                                      k_max_fine, zlevels, nz_lev0, finest_ref );

                int found_lev = -1;
                for (int tl = finest; tl >= 0; tl--) {
                    int grd = ag_ptr[tl](p, 0, cell_assignor).first;
                    if (grd >= 0) { found_lev = tl; break; }
                }
                tgt_ptr[i] = found_lev;
            });
            Gpu::synchronize();

            // Route particles to their target levels.
            for (int tl = 0; tl < n_levels; tl++) {
                if (tl == lev) { continue; }

                Gpu::DeviceVector<int> tl_mask(n_oor);
                auto* tl_mask_ptr = tl_mask.data();
                int target = tl;
                ParallelFor(n_oor, [=] AMREX_GPU_DEVICE (int i)
                {
                    tl_mask_ptr[i] = (tgt_ptr[i] == target) ? 1 : 0;
                });

                ReduceOps<ReduceOpSum> tl_reduce_op;
                ReduceData<int> tl_reduce_data(tl_reduce_op);
                using TLReduceTuple = typename decltype(tl_reduce_data)::Type;
                tl_reduce_op.eval(n_oor, tl_reduce_data,
                    [=] AMREX_GPU_DEVICE (int i) -> TLReduceTuple
                    { return {tl_mask_ptr[i]}; });
                int n_to_lev = amrex::get<0>(tl_reduce_data.value(tl_reduce_op));
                if (n_to_lev == 0) { continue; }

                ParticleTileType tl_tile;
                tl_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                tl_tile.resize(n_to_lev);
                [[maybe_unused]] int nc = amrex::filterParticles(tl_tile, tmp_tile, tl_mask_ptr,
                                                                 int(0), int(0), n_oor);
                AMREX_ASSERT(nc == n_to_lev);

                // Recompute idata(k) for the actual target level.
                auto* tl_pbox = tl_tile.GetArrayOfStructs()().data();
                const auto tl_plo = h_lg[tl].plo;
                const auto tl_dxi = h_lg[tl].dxi;
                int tl_kmax = h_lg[tl].k_max;
                int tl_ref  = h_lg[tl].ref_ratio;
                ParallelFor(n_to_lev, [=] AMREX_GPU_DEVICE (int i)
                {
                    auto& p = tl_pbox[i];
                    p.idata(ERFParticlesIntIdxAoS::k) =
                        compute_k_from_z( Real(p.pos(AMREX_SPACEDIM-1)),
                                          tl_plo[AMREX_SPACEDIM-1],
                                          tl_dxi[AMREX_SPACEDIM-1],
                                          tl_kmax, zlevels, nz_lev0, tl_ref);
                });
                Gpu::synchronize();

                auto& dst_tile = DefineAndReturnParticleTile(tl, dest_grids[tl], 0);
                auto dst_start = static_cast<int>(dst_tile.numParticles());
                dst_tile.resize(dst_start + n_to_lev);
                amrex::copyParticles(dst_tile, tl_tile, int(0), dst_start, n_to_lev);
                levels_modified[tl] = 1;
                n_routed_to[tl] += n_to_lev;
            }

            // Mark OOR particles as dead in the original tile
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                if (mask_ptr[i]) { p_pbox[i].id() = -1; }
            });
            Gpu::synchronize();
        }
    }

    // Ensure all ranks agree on which levels were modified
    ParallelDescriptor::ReduceIntMax(levels_modified.data(), n_levels);

    // Per-level redistribute: only for levels that need it.
    // Level a_lev always needs it (particles were advected, some marked dead).
    // Lower levels need it only if particles were routed there.
    for (int lev = 0; lev <= a_lev; lev++) {
        if (lev == a_lev || levels_modified[lev]) {
            Redistribute(lev, lev);
        }
    }
}


#endif

