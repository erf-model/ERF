/**
 * @file ERF_DustEmission.H
 * @brief Saltation and vertical dust emission flux for the ERF-Dust layer.
 *
 * Implements the Marticorena & Bergametti (1995) saltation bombardment model
 * for horizontal and vertical dust emission flux per particle size bin.
 *
 * Horizontal saltation flux (Owen 1964 form):
 *   Q_s = C_s * rho_a / g * u*^3 * (1 - u*_t/u*) * (1 + u*_t/u*)^2
 *   for u* > u*_t; Q_s = 0 otherwise.
 *   Units: [kg m^-1 s^-1]
 *
 * Vertical emission flux per size bin (Marticorena & Bergametti 1995):
 *   F_i = alpha_i * silt_fraction * Q_s
 *   Units: [kg m^-2 s^-1]
 *
 * where alpha_i is the sandblasting efficiency for bin i [m^-1],
 * derived from the silt fraction via the empirical relation:
 *   log10(alpha_i) = 0.134 * clay_fraction - 6.0
 * For mine tailings, clay_fraction is approximated as 0.2 * silt_fraction.
 *
 * References:
 *   Marticorena, B., & Bergametti, G. (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Owen, P.R. (1964). Saltation of uniform grains in air. J. Fluid Mech.,
 *     20, 225-242. https://doi.org/10.1017/S0022112064001173
 *   Shao, Y., & Lu, H. (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 *   Bagnold, R.A. (1941). The Physics of Blown Sand and Desert Dunes.
 *     Methuen, London.
 *   U.S. EPA AP-42 Chapter 13.2.2 (Unpaved Roads).
 *     https://www.epa.gov/air-emissions-factors-and-quantification/
 *     ap-42-compilation-air-emission-factors
 */

#ifndef ERF_DUST_EMISSION_H
#define ERF_DUST_EMISSION_H

#ifdef ERF_USE_DUST

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelFor.H>
#include <AMReX_MFIter.H>
#include <cmath>

namespace DustEmissionConst {
    // Saltation efficiency coefficient C_s [dimensionless].
    // Reference: Marticorena & Bergametti (1995), Eq. 5.
    constexpr amrex::Real C_S             = 2.61;

    // Standard air density [kg/m^3].
    constexpr amrex::Real RHO_AIR_STD     = 1.225;

    // Gravitational acceleration [m/s^2].
    constexpr amrex::Real G_ACCEL         = 9.81;

    // Sandblasting efficiency intercept (log10 scale).
    // log10(alpha) = ALPHA_SLOPE * clay_frac + ALPHA_INTERCEPT
    // Reference: Marticorena & Bergametti (1995), Eq. 9.
    constexpr amrex::Real ALPHA_SLOPE     = 0.134;
    constexpr amrex::Real ALPHA_INTERCEPT = -6.0;

    // Proxy for clay fraction from silt fraction at mine tailings.
    // clay_frac = CLAY_FROM_SILT * silt_frac
    constexpr amrex::Real CLAY_FROM_SILT  = 0.2;

    // Minimum u* / u*_t ratio before emission is considered active.
    // Prevents floating-point issues near the threshold.
    constexpr amrex::Real USTAR_RATIO_MIN = 1.001;

    // Maximum emission flux cap [kg/m^2/s].
    // Prevents physically unreasonable values from driving instabilities.
    constexpr amrex::Real FLUX_MAX        = 1.0e-3;
}

/**
 * @brief Compute horizontal saltation flux Q_s [kg/m/s].
 *
 * Owen (1964) form as used by Marticorena & Bergametti (1995):
 *   Q_s = C_s * rho_a / g * u*^3 * (1 - u*_t/u*) * (1 + u*_t/u*)^2
 * Returns 0 when u* <= u*_t.
 *
 * @param ustar    Friction velocity u* [m/s].
 * @param ustar_t  Threshold friction velocity u*_t [m/s].
 * @param rho_a    Air density [kg/m^3].
 * @return Q_s [kg/m/s], non-negative.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_saltation_flux(amrex::Real ustar,
                                    amrex::Real ustar_t,
                                    amrex::Real rho_a) noexcept
{
    if (ustar <= ustar_t * DustEmissionConst::USTAR_RATIO_MIN) return 0.0;
    amrex::Real rho = (rho_a > 0.0) ? rho_a : DustEmissionConst::RHO_AIR_STD;
    amrex::Real ratio = ustar_t / ustar;
    amrex::Real Qs = DustEmissionConst::C_S * rho / DustEmissionConst::G_ACCEL
                   * ustar * ustar * ustar
                   * (1.0 - ratio) * (1.0 + ratio) * (1.0 + ratio);
    return amrex::max(0.0, Qs);
}

/**
 * @brief Compute sandblasting efficiency alpha [m^-1] from silt fraction.
 *
 * log10(alpha) = ALPHA_SLOPE * clay_fraction + ALPHA_INTERCEPT
 * clay_fraction = CLAY_FROM_SILT * silt_fraction
 *
 * Reference: Marticorena & Bergametti (1995), Eq. 9.
 *
 * @param silt_frac  Surface silt mass fraction [0,1].
 * @return alpha [m^-1].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_sandblasting_efficiency(amrex::Real silt_frac) noexcept
{
    amrex::Real sf = amrex::min(amrex::max(silt_frac, 0.0), 1.0);
    amrex::Real clay = DustEmissionConst::CLAY_FROM_SILT * sf;
    amrex::Real log_alpha = DustEmissionConst::ALPHA_SLOPE * clay
                           + DustEmissionConst::ALPHA_INTERCEPT;
    return std::pow(10.0, log_alpha);
}

/**
 * @brief Compute vertical emission flux for one size bin [kg/m^2/s].
 *
 *   F_i = alpha_i * silt_fraction * Q_s
 *
 * @param alpha   Sandblasting efficiency [m^-1].
 * @param silt    Silt mass fraction [0,1].
 * @param Qs      Horizontal saltation flux [kg/m/s].
 * @return F_i [kg/m^2/s], clamped to [0, FLUX_MAX].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_vertical_emission_flux(amrex::Real alpha,
                                            amrex::Real silt,
                                            amrex::Real Qs) noexcept
{
    amrex::Real F = alpha * silt * Qs;
    return amrex::min(amrex::max(F, 0.0), DustEmissionConst::FLUX_MAX);
}

/**
 * @brief Fill dust_emission_flux MultiFab on the 2D dust grid.
 *
 * For each dust grid cell, computes saltation flux and vertical emission
 * per size bin using:
 *   - ustar from dust_ustar_in (placeholder u* until Phase 9)
 *   - ustar_t from dust_ustar_t
 *   - silt_fraction from dust_silt_fraction
 *
 * GPU ParallelFor over tiles. Pattern follows ERF_FireHeatFlux.H.
 *
 * @param[out] dust_emission_flux  Emission flux MultiFab (ncomp=n_size_bins).
 * @param[in]  dust_ustar_t        Threshold friction velocity [m/s].
 * @param[in]  dust_ustar_in       Friction velocity u* [m/s].
 *                                 Placeholder uniform value until Phase 9.
 * @param[in]  dust_silt_fraction  Silt fraction [0,1].
 * @param[in]  n_size_bins         Number of size bins.
 * @param[in]  rho_air             Air density [kg/m^3].
 */
void compute_dust_emission_flux(amrex::MultiFab&       dust_emission_flux,
                                 const amrex::MultiFab& dust_ustar_t,
                                 const amrex::MultiFab& dust_ustar_in,
                                 const amrex::MultiFab& dust_silt_fraction,
                                 int                    n_size_bins,
                                 amrex::Real            rho_air);

#endif // ERF_USE_DUST

#endif // ERF_DUST_EMISSION_H
