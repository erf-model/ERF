#include <ERFPC.H>
#include <ERFPCParticleToMesh.H>
#include <ERF_Constants.H>
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
            return ptd.m_rdata[ERFParticlesRealIdxSoA::mass][i];
        });
}

/*! \brief Redistribute particles after AMR regrid.  With pos stored in
 *         computational coordinates the cell index is exact in one pass, so
 *         no terrain-correction dance is needed. */
void ERFPC::FixKIndexAMR (const Vector<std::unique_ptr<MultiFab>>& /*a_z_phys_nd*/)
{
    BL_PROFILE("ERFPC::FixKIndexAMR()");
    Redistribute();
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



#endif

