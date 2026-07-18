#ifndef ERF_ROTHERMEL_H
#define ERF_ROTHERMEL_H

#include <AMReX_MultiFab.H>
#include <cmath>

/**
 * @file ERF_Rothermel.H
 * @brief Rothermel fire spread model implementation
 *
 * Implements the Rothermel (1972) fire spread model for rate-of-spread (ROS)
 * computation as a function of fuel type, moisture, wind, and slope.
 *
 * Key references:
 * - Rothermel, R.C. (1972): "A Mathematical Model for Predicting Fire Spread
 *   in Wildland Fuels", USDA Forest Service Research Paper INT-115
 * - Anderson, H.E. (1982): "Aids to determining fuel models for estimating fire
 *   behavior", General Technical Report INT-122, USDA Forest Service
 * - Albini, F.A. (1976): "Estimating wildfire behavior and effects",
 *   General Technical Report INT-30, USDA Forest Service
 *
 * Ported from: wildfire_levelset/src/rothermel_model.H and
 *              wildfire_levelset/src/compute_rothermel_R.H
 */

/**
 * @struct FuelModelParams
 * @brief Anderson FBFM13 fuel model parameters
 *
 * Stores parameters for one of the 13 standard fuel models used by
 * the fire behavior community (Anderson 1982).
 */
struct FuelModelParams
{
    // Fuel load [lb/ft² or kg/m², converted as needed]
    amrex::Real w_d1;       ///< 1-hour dead fuel load
    amrex::Real w_d10;      ///< 10-hour dead fuel load
    amrex::Real w_d100;     ///< 100-hour dead fuel load
    amrex::Real w_lh;       ///< Live herbaceous fuel load
    amrex::Real w_lw;       ///< Live woody fuel load

    // Surface area to volume ratio [1/ft or 1/m]
    amrex::Real sigma_d1;   ///< 1-hour SAV
    amrex::Real sigma_lh;   ///< Herbaceous SAV
    amrex::Real sigma_lw;   ///< Woody SAV

    // Fuel bed depth [ft or m]
    amrex::Real delta;

    // Fuel moisture of extinction [fraction]
    amrex::Real Mx;

    // Heat content [BTU/lb or J/kg]
    amrex::Real heat_content;

    // Fuel particle density [lb/ft³ or kg/m³]
    amrex::Real rho_p;
};

/**
 * @brief Get Anderson FBFM13 fuel model parameters
 *
 * @param[in] model_id Anderson fuel model ID (1–13)
 * @return FuelModelParams for the specified model
 *
 * All 13 Anderson fuel models are hard-coded with standard values.
 */
FuelModelParams get_anderson_fuel_params(int model_id);

/**
 * @struct RothermelComputed
 * @brief Precomputed Rothermel coefficients (fuel + moisture dependent)
 *
 * Stores intermediate results computed once per timestep to avoid
 * redundant calculations in the per-cell GPU kernel.
 */
struct RothermelComputed
{
    amrex::Real R0;              ///< No-wind, no-slope ROS [m/s]
    amrex::Real C, B;            ///< Wind factor φ_w = C * U^B
    amrex::Real beta_ratio_E;    ///< (β/β_opt)^E packing ratio factor
    amrex::Real beta;            ///< Packing ratio
    amrex::Real phi_s_const;     ///< Base slope factor (before terrain effects)
    amrex::Real U_max_ftmin;     ///< MEWS wind speed cap [ft/min]
    amrex::Real wind_conv;       ///< m/s to ft/min conversion factor
    amrex::Real ros_conv;        ///< ft/min to m/s conversion factor
    amrex::Real I_R;             ///< Reaction intensity [BTU/ft²/min]
};

/**
 * @brief Compute Rothermel fuel/moisture coefficients
 *
 * Precomputes all Rothermel intermediate quantities that depend only on
 * fuel type and moisture content (not on wind/slope).
 *
 * @param[in] fp FuelModelParams for selected fuel model
 * @return RothermelComputed with precomputed coefficients
 *
 * Implements the full Rothermel model:
 * - Reaction intensity I_R from fuel load, SAV, moisture, heat content
 * - Propagating flux ratio Ξ and damping coefficients
 * - No-wind ROS R0 from I_R and fuel characteristics
 * - Wind and slope factor coefficients C, B, E
 * - MEWS cap U_max_ftmin
 */
RothermelComputed compute_rothermel_params(const FuelModelParams& fp,
                                           amrex::Real moisture_1hr,
                                           amrex::Real moisture_10hr,
                                           amrex::Real moisture_100hr);

/**
 * @brief Compute rate-of-spread for a single fire cell (GPU kernel)
 *
 * Computes ROS = R0 * (1 + φ_w + φ_s) where:
 * - φ_w is wind factor (uses effective midflame wind after WAF and terrain)
 * - φ_s is slope factor (vectorial combination with wind direction)
 * - Includes MEWS cap on effective wind speed
 *
 * @param[in] ux_eff, uy_eff   Effective midflame wind [m/s] after WAF and terrain
 * @param[in] sx, sy           Terrain slopes [dz/dx, dz/dy]
 * @param[in] rc               RothermelComputed coefficients
 * @return ROS [m/s]
 *
 * Wind factor:
 *   φ_w = C * U_eff^B * (β/β_opt)^E   where U_eff is effective wind speed [ft/min]
 *
 * Slope factor (vectorial combination):
 *   - Compute slope angle θ = atan(|∇z|)
 *   - Project slope elevation gradient in wind direction
 *   - φ_s = tan(θ) * cos(angle_between_wind_and_slope_gradient)
 *   - If slope opposes wind (e.g., upslope): φ_s can be negative (downwind spread)
 *   - φ_s is computed locally; actual terrain forcing comes from terrain wind corrections
 *
 * MEWS cap:
 *   - If U_eff exceeds U_max_ftmin, use U_max_ftmin for wind factor only
 *   - Preserves steady-state behavior: E ≤ 0.9 * I_R always
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real rothermel_ros_cell(
    amrex::Real ux_eff, amrex::Real uy_eff,
    amrex::Real sx, amrex::Real sy,
    const RothermelComputed& rc) noexcept;

/**
 * @brief Compute rate-of-spread field on fire grid
 *
 * Applies rothermel_ros_cell() to all cells on the fire grid.
 *
 * @param[out] fire_ros        Output ROS [m/s]
 * @param[in]  fire_wind       Effective midflame wind after WAF + terrain [m/s]
 * @param[in]  fire_slopes     Terrain slopes [dz/dx, dz/dy]
 * @param[in]  rc              RothermelComputed coefficients
 */
void compute_ros_field(
    amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_wind,
    const amrex::MultiFab& fire_slopes,
    const RothermelComputed& rc);

#endif
