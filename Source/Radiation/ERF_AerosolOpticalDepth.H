#ifndef ERF_AEROSOL_OPTICAL_DEPTH_H_
#define ERF_AEROSOL_OPTICAL_DEPTH_H_

#include <AMReX_GpuControl.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

/**
 * @file ERF_AerosolOpticalDepth.H
 * @brief Prescribed bulk aerosol/turbidity optical depth diagnosis.
 *
 * Implements profile-based aerosol optical depth diagnosis for use in TwoStream radiation.
 * Supports Constant, Exponential, and Table (future) profile types.
 *
 * Aerosol optical depth is added on top of existing tau contributions (tau_base + tau_cloud + tau_dynamic).
 *
 * Constant profile:
 *   tau_aerosol(k) = aerosol_tau_per_layer  [uniform at all levels]
 *
 * Exponential profile:
 *   tau_aerosol(k) = aerosol_tau_surface * exp(-z(k) / aerosol_scale_height_m)
 *   [decay with height]
 *
 * Table profile (future):
 *   tau_aerosol(k) = lookup_table[k]  [per-level prescribed values]
 *
 * References:
 * -----------
 * - See Radiation Development notes (Source/Radiation/RAD_DEVELOPMENT.md) for implementation details.
 */

/**
 * @brief Diagnose aerosol optical depth for Constant profile.
 *
 * Returns uniform aerosol tau_per_layer (independent of height/level).
 * Clamped to [0, 100] for physical reasonableness.
 *
 * @param[in] tau_aerosol_const Constant aerosol optical depth per layer [dimensionless]
 * @return Clamped aerosol optical depth [0, 100]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_aerosol_constant(amrex::Real tau_aerosol_const)
{
    // Guard against invalid input
    if (!std::isfinite(tau_aerosol_const)) {
        return 0.0;
    }
    
    // Clamp to physically reasonable range [0, 100]
    amrex::Real tau_clipped = tau_aerosol_const;
    if (tau_clipped < 0.0) tau_clipped = 0.0;
    if (tau_clipped > 100.0) tau_clipped = 100.0;
    
    return tau_clipped;
}

/**
 * @brief Diagnose aerosol optical depth for Exponential profile.
 *
 * Returns height-dependent aerosol tau via exponential decay:
 *   tau_aerosol(k) = tau_surface * exp(-z(k) / scale_height)
 *
 * Clamped to [0, 100] for physical reasonableness.
 * Handles invalid parameters safely with fallback to 0.
 *
 * @param[in] z_level Height at level k [m]
 * @param[in] tau_surface Total-column aerosol optical depth at surface [dimensionless]
 * @param[in] scale_height_m Decay scale height [m]
 * @return Clamped aerosol optical depth [0, 100]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_aerosol_exponential(amrex::Real z_level,
                                            amrex::Real dz,                // NEW: layer thickness
                                            amrex::Real tau_surface,
                                            amrex::Real scale_height_m)
{
    if (!std::isfinite(z_level) || !std::isfinite(dz) ||
        !std::isfinite(tau_surface) || !std::isfinite(scale_height_m)) {
        return 0.0;
    }
    if (scale_height_m <= 0.0 || dz <= 0.0) {
        return 0.0;
    }
    if (tau_surface < 0.0) {
        return 0.0;
    }
    if (z_level < 0.0) {
        z_level = 0.0;
    }

    amrex::Real arg = -z_level / scale_height_m;
    if (arg < -100.0) {
        return 0.0;
    }

    // Per-layer contribution = (tau_surface / scale_height) * dz * exp(-z/H)
    // so the vertical integral over all layers approaches tau_surface as domain height >> H
    amrex::Real tau_aerosol = (tau_surface / scale_height_m) * dz * std::exp(arg);

    if (tau_aerosol < 0.0)   tau_aerosol = 0.0;
    if (tau_aerosol > 100.0) tau_aerosol = 100.0;

    return tau_aerosol;
}


/**
 * @brief Diagnose aerosol optical depth for Table profile (future extension).
 *
 * Currently returns 0.0 as a placeholder. Future phases can implement
 * per-level table lookup via fixed-size array or MultiFab.
 *
 * @param[in] k Level index
 * @return Zero (not yet implemented)
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_aerosol_table(int k)
{
    // Future extension: implement per-level table lookup
    // For now, return 0 (no contribution)
    (void)k;  // Suppress unused parameter warning
    return 0.0;
}

#endif // ERF_AEROSOL_OPTICAL_DEPTH_H_
