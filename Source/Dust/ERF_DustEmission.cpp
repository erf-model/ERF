/**
 * @file ERF_DustEmission.cpp
 * @brief Implementation of dust emission flux computation kernel.
 */

#include <ERF_DustEmission.H>

#ifdef ERF_USE_DUST

#include <AMReX_MFIter.H>

using namespace amrex;

void compute_dust_emission_flux(MultiFab& dust_emission_flux,
                                 const MultiFab& dust_ustar_t,
                                 const MultiFab& dust_ustar_in,
                                 const MultiFab& dust_silt_fraction,
                                 int n_size_bins, Real rho_air)
{
    // n_size_bins must match nComp of dust_emission_flux.
    AMREX_ALWAYS_ASSERT(dust_emission_flux.nComp() == n_size_bins);

    for (MFIter mfi(dust_emission_flux, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto flux_arr  = dust_emission_flux.array(mfi);
        auto ut_arr    = dust_ustar_t.const_array(mfi);
        auto us_arr    = dust_ustar_in.const_array(mfi);
        auto silt_arr  = dust_silt_fraction.const_array(mfi);
        const Real rho = rho_air;
        const int nbins = n_size_bins;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real ustar   = us_arr(i, j, k);
            Real ustar_t = ut_arr(i, j, k);
            Real silt    = silt_arr(i, j, k);

            // Horizontal saltation flux: Owen (1964) via Marticorena & Bergametti (1995).
            Real Qs = compute_saltation_flux(ustar, ustar_t, rho);

            // Sandblasting efficiency from silt fraction.
            // Marticorena & Bergametti (1995), Eq. 9.
            Real alpha = compute_sandblasting_efficiency(silt);

            // Vertical flux per bin. In Phase 6 all bins use the same alpha.
            // Phase 11 will apply per-bin size-dependent alpha corrections.
            for (int b = 0; b < nbins; ++b) {
                flux_arr(i, j, k, b) = compute_vertical_emission_flux(alpha, silt, Qs);
            }
        });
    }
}

#endif // ERF_USE_DUST
