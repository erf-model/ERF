#ifndef ERF_FIRE_WIND_EXTRACT_H
#define ERF_FIRE_WIND_EXTRACT_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <cmath>
#include <ERF_FireGrid.H>

/**
 * @file ERF_FireWindExtract.H
 * @brief Wind extraction and adjustment for fire simulation
 *
 * Implements:
 * 1. MOST wind reconstruction at reference height (Dyer 1974)
 * 2. Wind Adjustment Factor (WAF) for forest canopy effects
 * 3. FARSITE terrain wind corrections (ridge, shelter, valley, deflection)
 *
 * References:
 * - Rothermel (1972): USDA Forest Service Research Paper INT-115
 * - Dyer (1974): Chapter 2 in Modelling of Atmospheric Flow Fields
 * - Andrews (1986): BEHAVE: Fire Behavior Prediction and Fuel Modeling System
 * - FARSITE documentation (Finney 1998)
 */


/**
 * @brief Extract wind from direct vertical interpolation on fire grid
 *
 * Maps atmospheric wind to fire grid at reference height by interpolating
 * vertically through the atmospheric column at the fire cell location.
 * Uses cell-center wind values (averaged from face-staggered data) and
 * linearly interpolates to the target height above terrain.
 *
 * When fuel_model_mf is non-null and per-fuel-specific heights are available,
 * each fire cell looks up its fuel model ID and uses the corresponding per-fuel
 * wind height (following WRF-SFIRE fcwh convention) instead of the global z_ref.
 *
 * @param[out] fire_wind_ref   u,v at z_ref on fire grid [m/s]
 * @param[out] fire_wind_extract_z_mf  Wind extraction height z_target [m]
 * @param[in]  xvel_mf         x-velocity (face-staggered at i+1/2, j, k)
 * @param[in]  yvel_mf         y-velocity (face-staggered at i, j+1/2, k)
 * @param[in]  z_phys_cc_mf    Cell-center terrain height z_phys_cc[0] on atm grid
 * @param[in]  fg              Fire grid
 * @param[in]  z_ref           Global fallback wind reference height [m]
 * @param[in]  nz              Number of vertical levels
 * @param[in]  fuel_model_mf   Per-cell fuel model ID; nullptr = uniform fuel
 * @param[in]  d_fcwh          Device array of per-fuel wind heights [m];
 *                             size must be nfuelcats+1 with index 0 unused.
 *                             nullptr = use global z_ref
 * @param[in]  nfuelcats       Number of fuel categories (13 for Anderson FBFM13);
 *                             0 if d_fcwh is nullptr
 *
 * Algorithm per fire cell (i_f, j_f):
 * 1. Map to atmospheric column: i_a = i_f / C, j_a = j_f / C
 * 2. Get surface height: z_surf = z_phys_cc(i_a, j_a, 0)
 * 3. Lookup wind height:
 *    - If fuel_model_mf != nullptr and d_fcwh != nullptr:
 *      fuel_code = clamp(fuel_model_mf(i_f,j_f,0), 1, nfuelcats)
 *      z_ref_cell = d_fcwh[fuel_code]
 *    - Otherwise: z_ref_cell = z_ref (global)
 * 4. Compute target height: z_target = z_surf + z_ref_cell
 * 5. Walk k from 0 to nz-2: find k_lo such that
 *    z_phys_cc(i_a, j_a, k_lo) <= z_target < z_phys_cc(i_a, j_a, k_lo+1)
 *    Clamp to k_lo = nz-2 if z_target exceeds column top
 * 6. Compute interpolation weight:
 *    alpha = (z_target - z_phys_cc(i_a, j_a, k_lo)) /
 *            (z_phys_cc(i_a, j_a, k_lo+1) - z_phys_cc(i_a, j_a, k_lo))
 *    Clamp alpha to [0, 1]
 * 7. Average u/v from bracketing face values to get cell centers:
 *    u_cc(i,j,k) = 0.5 * (xvel(i,j,k) + xvel(i+1,j,k))
 *    v_cc(i,j,k) = 0.5 * (yvel(i,j,k) + yvel(i,j+1,k))
 * 8. Interpolate: fw(i_f,j_f,0,0) = u_cc(k_lo) + alpha*(u_cc(k_lo+1) - u_cc(k_lo))
 *                fw(i_f,j_f,0,1) = v_cc(k_lo) + alpha*(v_cc(k_lo+1) - v_cc(k_lo))
 * 9. Store z_target in fire_wind_extract_z_mf
 */
void fill_fire_wind_from_interpolation(
    amrex::MultiFab&       fire_wind_ref,
    amrex::MultiFab&       fire_wind_extract_z_mf,
    const amrex::MultiFab& xvel_mf,
    const amrex::MultiFab& yvel_mf,
    const amrex::MultiFab& z_phys_cc_mf,
    const FireGrid&        fg,
    amrex::Real            z_ref,
    int                    nz,
    const amrex::MultiFab* fuel_model_mf = nullptr,
    const amrex::Real*     d_fcwh         = nullptr,
    int                    nfuelcats      = 0);

/**
 * @brief Compute unsheltered Wind Adjustment Factor (Andrews/Albini & Baughman 1979)
 *
 * @param[in] delta_ft Fuel bed depth [ft]
 * @return WAF (≈ 0.3–0.4 for typical forest)
 *
 * Formula: WAF = 1.83 / ln((20 + 0.36*h) / (0.13*h))
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_unsheltered(amrex::Real delta_ft) noexcept
{
    delta_ft = amrex::max(delta_ft, 0.1);
    amrex::Real numerator = 20.0 + 0.36 * delta_ft;
    amrex::Real denominator = 0.13 * delta_ft;
    return 1.83 / std::log(numerator / denominator);
}

/**
 * @brief Compute sheltered Wind Adjustment Factor (FARSITE formula)
 *
 * @param[in] f_c      Canopy closure [fraction, 0–1]
 * @param[in] h_c_ft   Canopy height [ft]
 * @param[in] delta_ft Fuel bed depth [ft]
 * @return WAF
 *
 * Formula: WAF = 0.555 / (sqrt(f_c*h_c) * ln((20+0.36*h_c)/(0.13*h_c)))
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_sheltered(amrex::Real f_c, amrex::Real h_c_ft,
                                   amrex::Real delta_ft) noexcept
{
    f_c = amrex::max(f_c, 0.01);
    h_c_ft = amrex::max(h_c_ft, 0.1);
    delta_ft = amrex::max(delta_ft, 0.1);

    amrex::Real ln_term = std::log((20.0 + 0.36 * h_c_ft) / (0.13 * h_c_ft));
    amrex::Real denom = std::sqrt(f_c * h_c_ft) * ln_term;
    denom = amrex::max(denom, 1.0e-6);

    return 0.555 / denom;
}

/**
 * @brief Compute BehavePlus Wind Adjustment Factor (linear formula)
 *
 * @param[in] delta_ft Fuel bed depth [ft]
 * @return WAF
 *
 * Formula: WAF = 0.36 + 0.004*h where h is in feet
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_behaviorplus(amrex::Real delta_ft) noexcept
{
    return 0.36 + 0.004 * delta_ft;
}

/**
 * @brief Apply canopy Beer-Lambert extinction to WAF
 *
 * @param[in] waf_open   WAF in open (f_c = 0)
 * @param[in] f_c        Canopy closure [0–1]
 * @param[in] alpha_c    Canopy extinction coefficient
 * @return WAF_sheltered
 *
 * Formula: WAF_sheltered = waf_open * exp(-alpha_c * f_c)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_canopy_exp(amrex::Real waf_open, amrex::Real f_c,
                                    amrex::Real alpha_c) noexcept
{
    return waf_open * std::exp(-alpha_c * f_c);
}

/**
 * @brief Compute Maximum Effective Wind Speed cap (Rothermel MEWS)
 *
 * Derives the wind speed cap from Rothermel reaction intensity using
 * φ_w_max = 0.9 * I_R and the inverse of the wind factor formula.
 *
 * @param[in] I_R              Reaction intensity [BTU/ft²/min]
 * @param[in] C, B             Wind factor coefficients
 * @param[in] beta_ratio_E     (β/β_opt)^E factor
 * @return U_max [ft/min]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_max_ews_ftmin(amrex::Real I_R, amrex::Real C,
                                   amrex::Real B, amrex::Real beta_ratio_E) noexcept
{
    const amrex::Real phi_w_max = 0.9 * I_R;

    // φ_w = C * U^B * (β/β_opt)^E  =>  U = (φ_w / (C * (β/β_opt)^E))^(1/B)
    if (B > 0.0 && beta_ratio_E > 0.0 && C > 0.0) {
        amrex::Real denom = C * beta_ratio_E;
        amrex::Real base = phi_w_max / denom;
        if (base <= 0.0) return 0.0;
        return std::pow(base, 1.0 / B);
    }
    return 0.0;
}

/**
 * @brief Compute terrain curvature on fire grid
 *
 * Curvature κ = d²z/dx² + d²z/dy² (positive for convex peaks).
 * Used in FARSITE ridge/valley detection.
 *
 * @param[out] curvature       1-component MultiFab on fire grid
 * @param[in]  fire_slopes     2-component (dz/dx, dz/dy) on fire grid
 * @param[in]  geom_fire       Fire grid geometry
 */
void compute_terrain_curvature(
    amrex::MultiFab& curvature,
    const amrex::MultiFab& fire_slopes,
    const amrex::Geometry& geom_fire);

/**
 * @brief Apply FARSITE terrain wind corrections
 *
 * Modifies wind speed and direction based on terrain:
 * - Ridge speed-up: windward ridge faces get k_ridge factor
 * - Sheltering: lee-side slopes experience wind reduction
 * - Valley channeling: wind is deflected into valleys
 * - Directional deflection: wind can turn up to ±45° on slopes
 *
 * @param[in,out] fire_wind_eff  u, v on fire grid (modified in-place)
 * @param[in]     fire_slopes    dz/dx, dz/dy
 * @param[in]     fire_curvature Terrain curvature
 * @param[in]     k_ridge        Ridge speed-up factor (≈ 1.5)
 * @param[in]     k_shelter      Lee-side shelter factor (≈ 0.6)
 * @param[in]     k_valley       Valley channeling factor (≈ 0.8)
 * @param[in]     k_deflect      Wind deflection factor (≈ 0.3)
 *
 * Algorithm (per cell):
 * 1. Compute slope magnitude |∇z| and wind magnitude |U|
 * 2. Compute wind-upslope alignment: cos_upslope = (∇z · U) / (|∇z| |U|)
 * 3. Classify terrain:
 *    - Ridge (curv > 0.01, |∇z| > 0.05, cos_upslope > 0): factor = 1 + (k_ridge-1)*min(|∇z|,1)
 *    - Shelter (cos_upslope < -0.5, curv > 0.01): factor = k_shelter
 *    - Valley (curv < -0.01): factor = k_valley + (1-k_valley)*(1-min(|∇z|,1))
 *    - Flat/other: factor = 1
 * 4. Clamp: factor = clamp(factor, 0.1, 3.0)
 * 5. Apply speed scaling: U *= factor
 * 6. For valleys: apply deflection using 2D rotation by angle k_deflect*asin(sin_cross/|U||∇z|*factor)
 */
void apply_farsite_terrain_wind(
    amrex::MultiFab& fire_wind_eff,
    const amrex::MultiFab& fire_slopes,
    const amrex::MultiFab& fire_curvature,
    amrex::Real k_ridge, amrex::Real k_shelter,
    amrex::Real k_valley, amrex::Real k_deflect);

#endif
