/**
 * @file ERF_PBLScaleAwareBlending.H
 * @brief GPU-inline blending functions for PBL diffusivity at grey-zone resolution.
 *
 * Provides scale-aware reduction of PBL vertical diffusivity (K_h) when the
 * horizontal grid spacing falls between mesoscale and LES resolution.
 *
 * Designed as a standalone header callable from any PBL scheme (MRF, YSU,
 * MYNN) without code duplication. All functions are AMREX_GPU_HOST_DEVICE
 * AMREX_FORCE_INLINE free functions with no state.
 *
 * Usage in an MFIter loop after K_h is computed:
 *   K_h_eff = pbl_kh_blend_and_cap(K_h, dx, L_blend, C_s, c_max, SmnSmn, use_smag);
 *
 * Gating: pbl_kh_blend_and_cap returns K_h unchanged when L_blend <= 0.
 * Setting L_blend = 0 recovers original PBL behaviour exactly.
 *
 * References:
 *   Boutle, I.A. et al. (2014). Mon. Wea. Rev., 142, 1655.
 *     https://doi.org/10.1175/MWR-D-13-00229.1
 *   Smagorinsky, J. (1963). Mon. Wea. Rev., 91, 99.
 *     https://doi.org/10.1175/1520-0493(1963)091<0099:GCEWTP>2.3.CO;2
 *   Honnert, R. et al. (2011). J. Atmos. Sci., 68, 2742.
 *     https://doi.org/10.1175/JAS-D-11-025.1
 *   Beare, R.J. (2014). Bound.-Layer Meteor., 153, 345.
 *     https://doi.org/10.1007/s10546-013-9881-3
 *   Hong, S.-Y., & Pan, H.-L. (1996). Mon. Wea. Rev., 124, 2322.
 *     https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
 */

#ifndef ERF_PBL_SCALE_AWARE_BLENDING_H
#define ERF_PBL_SCALE_AWARE_BLENDING_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

/**
 * @brief Boutle et al. (2014) blending factor.
 *
 * f = (dx/L_blend)^2 / (1 + (dx/L_blend)^2)
 *
 * Returns 1 when L_blend <= 0 (gate: no blending).
 * Asymptotes to 1 as dx >> L_blend (mesoscale grid, full PBL K_h).
 * Asymptotes to 0 as dx << L_blend (LES grid, no PBL K_h contribution).
 *
 * This is the weight on the parameterized (subgrid) part of the turbulence, so it
 * must grow with dx: the coarser the grid, the less of the boundary layer overturning
 * is resolved and the more the PBL scheme has to supply. It was previously written as
 * 1/(1+(dx/L_blend)^2), which has both limits backwards -- it switched the PBL scheme
 * off at mesoscale resolution and left it at full strength at LES resolution, on top of
 * whatever the resolved motions were already doing (issue #3619).
 *
 * @param dx       Horizontal grid spacing [m].
 * @param L_blend  Blending length scale [m]. Typical: 750 m. 0 = disabled.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pbl_blend_factor(amrex::Real dx,
                              amrex::Real L_blend) noexcept
{
    if (L_blend <= 0.0) return 1.0;
    const amrex::Real r = dx / L_blend;
    return (r * r) / (1.0 + r * r);
}

/**
 * @brief Smagorinsky-Lilly K_h ceiling.
 *
 * K_h_max = (C_s * dx)^2 * sqrt(2 * SmnSmn)
 *
 * Used when strain rate SmnSmn is available (e.g. from SmnSmn_lev in ERF).
 *
 * @param dx      Horizontal grid spacing [m].
 * @param SmnSmn  S_ij*S_ij strain rate squared [s^-2]. Must be >= 0.
 * @param C_s     Smagorinsky coefficient. Typical: 0.17.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pbl_kh_ceiling_smag(amrex::Real dx,
                                  amrex::Real SmnSmn,
                                  amrex::Real C_s) noexcept
{
    const amrex::Real Lmix = C_s * dx;
    return Lmix * Lmix * std::sqrt(2.0 * amrex::max(SmnSmn, amrex::Real(0.0)));
}

/**
 * @brief Power-law K_h ceiling (fallback when SmnSmn not available).
 *
 * K_h_max = c_max * dx^(4/3)
 *
 * Derived from dimensional analysis of the inertial subrange.
 * Reference: Honnert et al. (2011). https://doi.org/10.1175/JAS-D-11-025.1
 *
 * @param dx     Horizontal grid spacing [m].
 * @param c_max  Coefficient [m^(2/3) s^-1]. Typical: 0.1.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pbl_kh_ceiling_powerlaw(amrex::Real dx,
                                     amrex::Real c_max) noexcept
{
    return c_max * std::pow(dx, amrex::Real(4.0) / amrex::Real(3.0));
}

/**
 * @brief Apply scale-aware blending and ceiling to a single K_h value.
 *
 * K_h_eff = min(K_h * f_blend(dx, L_blend), K_h_max)
 *
 * Gate: returns K_h unchanged when L_blend <= 0.
 *
 * Ceiling selection:
 *   use_smag = true  and SmnSmn >= 0: Smagorinsky ceiling.
 *   otherwise: power-law ceiling.
 *
 * Callable from any PBL scheme (MRF, YSU, MYNN) without modification.
 *
 * @param K_h      PBL vertical diffusivity [m^2/s] before blending.
 * @param dx       Horizontal grid spacing [m].
 * @param L_blend  Boutle blending length [m]. 0 = no blending (returns K_h).
 * @param C_s      Smagorinsky coefficient (used when use_smag=true).
 * @param c_max    Power-law coefficient [m^(2/3)/s] (used when use_smag=false).
 * @param SmnSmn   Strain rate squared [s^-2]. Use -1 when not available.
 * @param use_smag Use Smagorinsky ceiling when SmnSmn >= 0.
 * @return         Blended and capped K_h [m^2/s].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pbl_kh_blend_and_cap(amrex::Real K_h,
                                   amrex::Real dx,
                                   amrex::Real L_blend,
                                   amrex::Real C_s,
                                   amrex::Real c_max,
                                   amrex::Real SmnSmn,
                                   bool        use_smag) noexcept
{
    if (L_blend <= 0.0) return K_h;

    const amrex::Real f       = pbl_blend_factor(dx, L_blend);
    const amrex::Real K_blend = K_h * f;

    amrex::Real K_max;
    if (use_smag && SmnSmn >= 0.0) {
        K_max = pbl_kh_ceiling_smag(dx, SmnSmn, C_s);
    } else {
        K_max = pbl_kh_ceiling_powerlaw(dx, c_max);
    }
    return amrex::min(K_blend, K_max);
}

#endif // ERF_PBL_SCALE_AWARE_BLENDING_H
