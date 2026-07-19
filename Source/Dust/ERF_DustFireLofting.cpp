#include <ERF_DustFireLofting.H>

#ifdef ERF_USE_DUST

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

using namespace amrex;

void apply_fire_lofting_to_emission_flux(
    MultiFab&       dust_emission_flux,
    const MultiFab& fire_heat_scratch,
    int             n_bins,
    Real            k_loft,
    Real            Q_threshold,
    Real            Q_ref,
    bool            dust_debug,
    int             step)
{
    if (Q_ref <= 0.0_rt) return;

    for (MFIter mfi(dust_emission_flux, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto heat = fire_heat_scratch.const_array(mfi);
        for (int b = 0; b < n_bins; ++b) {
            auto flux = dust_emission_flux.array(mfi, b);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real Q      = heat(i, j, k);
                Real excess = amrex::max(Q - Q_threshold, 0.0_rt);
                Real factor = 1.0_rt + amrex::min(k_loft * excess / Q_ref, k_loft);
                flux(i, j, k) *= factor;
            });
        }
    }

    if (dust_debug) {
        Real flux_max = dust_emission_flux.max(0);
        Real heat_max = fire_heat_scratch.max(0);
        Print() << "[DUST DEBUG] Phase 3 fire-lofting: step=" << step
                << " fire_heat_max=" << heat_max << " W/m2"
                << " emission_flux_max_after_loft=" << flux_max << " kg/m2/s\n";
    }
}

#endif // ERF_USE_DUST
