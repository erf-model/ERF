/**
 * @file ERF_DustAtmCoupling.cpp
 * @brief Implementation of dust-to-atmosphere coupling functions.
 */

#ifdef ERF_USE_DUST

#include <ERF_DustAtmCoupling.H>
#include <AMReX_MFIter.H>

using namespace amrex;

void apply_dust_tendency_to_cc_source(
    MultiFab&       cc_source,
    const MultiFab& Q_dust_atm,
    const MultiFab& z_phys_cc,
    const Geometry& geom_atm,
    int             dust_scalar_comp,
    Real            feedback,
    bool            dust_debug)
{
    if (feedback <= 0.0) return;

    const Box& domain = geom_atm.Domain();
    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);
    const auto& dx = geom_atm.CellSize();
    Real dz_avg = dx[2];

    for (MFIter mfi(cc_source, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto src_arr  = cc_source.array(mfi);
        auto q_arr    = Q_dust_atm.const_array(mfi);
        auto z_arr    = z_phys_cc.const_array(mfi);
        const int comp = dust_scalar_comp;
        const Real fb  = feedback;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Injection at k=0 only. Higher levels receive zero tendency here.
            if (k != klo) return;

            // Cell depth: use z_phys_cc if terrain is active, else dz_avg.
            Real dz = (k < khi) ? (z_arr(i,j,k+1) - z_arr(i,j,k)) : dz_avg;
            if (dz <= 1.0e-10) dz = dz_avg;

            // d(RhoDust)/dt = F_dust * feedback / dz
            src_arr(i, j, k, comp) += fb * q_arr(i, j, 0) / dz;
        });
    }

    if (dust_debug) {
        Real F_max   = Q_dust_atm.max(0);
        Real tend_max= cc_source.max(dust_scalar_comp);
        Real tend_sum= cc_source.sum(dust_scalar_comp);
        amrex::Print() << "[DUST COUPLING] F_dust_max=" << F_max
                       << " kg/m^2/s  RhoDust_tend max=" << tend_max
                       << " sum=" << tend_sum << " [kg/m^3/s]\n";
    }
}

#endif // ERF_USE_DUST
