#ifndef ERF_BALBI_MODEL_H
#define ERF_BALBI_MODEL_H

#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <cmath>
#include <vector>
#include <ERF_Rothermel.H>
#include <ERF_FireParams.H>

/**
 * @file ERF_BalbiModel.H
 * @brief Balbi (2009) physical fire spread model implementation
 *
 * Implements the Balbi (2009) fire spread model for rate-of-spread (ROS)
 * computation as a function of fuel type, moisture, wind, and slope.
 *
 * The model derives the rate of spread from the radiant energy flux received
 * by unburned fuel from a tilted flame front:
 *
 *   R       = A × (1 + sin α − cos α)                 [m/s]
 *   A       = χ σ_m δ_m / (2 τ₀ B*)                   [m/s]
 *   χ       = r₀₀ σ_m / (1 + r₀₀ σ_m)                 [−]
 *   B*      = (C_pf (T_i − T_a) + M_f Λ) / h          [−]
 *   tan α   = U / v_b + |∇z|
 *   v_b     = sqrt(g δ_m (T_f − T_a) / T_a)           [m/s]
 *
 * where χ is the fraction of combustion heat radiated forward, B* is the
 * dimensionless energy needed to ignite the fuel, v_b is the vertical flame
 * velocity set by thermal forcing, and α is the flame tilt angle produced by
 * wind and terrain slope.
 *
 * Inputs:
 *   U:   wind speed [m/s]
 *   δ_m: fuel bed depth [m]
 *   σ_m: surface-area-to-volume ratio of the 1-hr dead fuel [1/m]
 *   h:   heat of combustion [J/kg]
 *   M_f: fuel moisture content [fraction]
 *   All US/Rothermel unit inputs (sigma_d1 [1/ft], delta [ft],
 *   heat_content [BTU/lb]) are converted to SI internally.
 *
 * Limitations:
 *   Not appropriate for fuels with very low fuel bed depth (δ_m < 0.01 m)
 *   — returns zero ROS for such cases.
 *
 * This is a parity port of wildfire_levelset/src/balbi_model.H and
 * wildfire_levelset/src/compute_balbi_R.H
 * (https://github.com/hgopalan/wildfire_levelset).
 *
 * References:
 *   Balbi, J.-H., Rossi, J.-L., Marcelli, T. & Santoni, P.-A. (2009).
 *     A physical model for wildland fires. Combustion and Flame, 156(12),
 *     2217–2230.
 *   Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling
 *     System. USDA Forest Service GTR INT-194.
 */

/**
 * @struct BalbiComputed
 * @brief Precomputed Balbi coefficients (fuel + moisture dependent)
 */
struct BalbiComputed
{
    amrex::Real A_coeff;  ///< ROS amplitude coefficient A [m/s]
    amrex::Real v_b;      ///< Buoyancy velocity [m/s]
};

// Unit conversions, US/Rothermel → SI
inline constexpr amrex::Real BALBI_FT_TO_M         = 0.3048;   ///< [ft] → [m]
inline constexpr amrex::Real BALBI_FT_INV_TO_M_INV = 3.28084;  ///< [1/ft] → [1/m]
inline constexpr amrex::Real BALBI_BTU_LB_TO_J_KG  = 2326.0;   ///< [BTU/lb] → [J/kg]

/// Gravitational acceleration [m/s²]
inline constexpr amrex::Real BALBI_GRAVITY = 9.81;

/// Minimum safe buoyancy velocity [m/s]. Prevents division by zero in the
/// tilt-angle kernel while staying orders of magnitude below physically
/// meaningful v_b values (≥ 1 m/s).
inline constexpr amrex::Real BALBI_V_B_MIN = 1.0e-3;

/// Minimum dimensionless ignition energy B* before a fuel is treated as
/// non-burnable (A_coeff = 0). Physical values are 0.01–0.10.
inline constexpr amrex::Real BALBI_B_STAR_MIN = 1.0e-6;

/// Minimum fuel bed depth [m] below which the model is not applicable.
inline constexpr amrex::Real BALBI_MIN_DEPTH_M = 0.01;

/// Extreme-wind guard [m/s]: cells above this speed (or holding a non-finite
/// wind) fall back to the domain-mean wind before the ROS is evaluated.
inline constexpr amrex::Real BALBI_MAX_WIND_SPEED_MS = 25.0;
inline constexpr amrex::Real BALBI_MAX_WIND_SPEED_SQ =
    BALBI_MAX_WIND_SPEED_MS * BALBI_MAX_WIND_SPEED_MS;

/// Balbi table size: Anderson FBFM13 codes 1–13 plus index 0 (non-burnable).
inline constexpr int BALBI_TABLE_SIZE = 14;

/**
 * @brief Compute Balbi fuel/thermal coefficients
 *
 * Precomputes A_coeff and v_b, which depend on fuel type, fuel moisture and
 * ambient conditions but not on wind or slope.
 *
 * @param[in] fp   FuelModelParams for the selected fuel model (US units)
 * @param[in] bp   BalbiParams with thermal properties and constants
 * @param[in] M_f  Fuel moisture content [fraction]
 * @return BalbiComputed with precomputed coefficients
 */
inline BalbiComputed compute_balbi_params(
    const FuelModelParams& fp,
    const FireParams::BalbiParams& bp,
    amrex::Real M_f) noexcept
{
    BalbiComputed bc;

    // Convert fuel parameters from US/Rothermel units to SI
    const amrex::Real sigma_m = fp.sigma_d1     * BALBI_FT_INV_TO_M_INV;  // [1/m]
    const amrex::Real delta_m = fp.delta        * BALBI_FT_TO_M;          // [m]
    const amrex::Real h_si    = fp.heat_content * BALBI_BTU_LB_TO_J_KG;   // [J/kg]

    // Guard against degenerate fuel beds. The depth threshold is in metres,
    // matching the documented limit of the model.
    if (delta_m < BALBI_MIN_DEPTH_M || sigma_m <= 0.0 || h_si <= 0.0) {
        bc.A_coeff = 0.0;
        bc.v_b     = 1.0;  // avoid division by zero in the ROS kernel
        return bc;
    }

    // Radiation fraction χ = χ₀ / (1 + χ₀), with χ₀ = r₀₀ σ_m dimensionless
    const amrex::Real chi_0 = bp.r_00 * sigma_m;
    const amrex::Real chi   = chi_0 / (1.0 + chi_0);

    // Dimensionless ignition energy B* = (C_pf (T_i − T_a) + M_f Λ) / h.
    // Λ (bp.delta_H) is the latent heat of vaporisation of the fuel water;
    // h (fp.heat_content) is the heat of combustion. These are distinct.
    const amrex::Real E_ig   = bp.C_pf * (bp.T_i - bp.T_a) + M_f * bp.delta_H;
    const amrex::Real B_star = E_ig / h_si;

    // Buoyancy velocity v_b = sqrt(g δ_m (T_f − T_a) / T_a)
    const amrex::Real dT_f = bp.T_f - bp.T_a;
    const amrex::Real v_b  = (dT_f > 0.0)
                           ? std::sqrt(BALBI_GRAVITY * delta_m * dT_f / bp.T_a)
                           : BALBI_V_B_MIN;

    // ROS amplitude coefficient. σ_m δ_m is dimensionless, so dividing by
    // τ₀ [s/m] leaves [m/s].
    const amrex::Real A_coeff = (B_star > BALBI_B_STAR_MIN)
                              ? chi * sigma_m * delta_m / (2.0 * bp.tau_0 * B_star)
                              : 0.0;

    bc.A_coeff = amrex::max(A_coeff, 0.0);
    bc.v_b     = amrex::max(v_b, BALBI_V_B_MIN);
    return bc;
}

/**
 * @brief Build a per-fuel-code lookup table of Balbi coefficients
 *
 * Mirrors build_fcwh_table(): the table is indexed by integer fuel code, with
 * index 0 reserved for the non-burnable code (zero spread) and 1–13 holding
 * the Anderson FBFM13 models.
 *
 * @param[in] bp   BalbiParams with thermal properties and constants
 * @param[in] M_f  Fuel moisture content [fraction]
 * @return host-side vector of size BALBI_TABLE_SIZE
 */
inline std::vector<BalbiComputed> build_fuel_balbi_table(
    const FireParams::BalbiParams& bp,
    amrex::Real M_f)
{
    std::vector<BalbiComputed> table(BALBI_TABLE_SIZE);

    // Code 0 is non-burnable: no spread.
    table[0].A_coeff = 0.0;
    table[0].v_b     = 1.0;

    for (int code = 1; code < BALBI_TABLE_SIZE; ++code) {
        table[code] = compute_balbi_params(get_anderson_fuel_params(code), bp, M_f);
    }

    return table;
}

/**
 * @brief Compute the flame tilt angle for the Balbi model
 *
 * Given wind speed U and terrain slope θ, compute the resultant angle α
 * via tan α = U/v_b + tan θ. Provided as a diagnostic helper; the ROS kernel
 * below obtains sin α and cos α from tan α directly.
 *
 * @param[in] U        Wind speed [m/s]
 * @param[in] theta    Terrain slope [radians]
 * @param[in] v_b      Buoyancy velocity [m/s]
 * @return angle α [radians]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_balbi_angle(
    amrex::Real U,
    amrex::Real theta,
    amrex::Real v_b) noexcept
{
    // Local constexpr copies so the constants are usable from device code
    constexpr amrex::Real v_b_min = BALBI_V_B_MIN;
    if (v_b < v_b_min) v_b = v_b_min;  // Avoid division by zero

    amrex::Real tan_alpha = U / v_b + std::tan(theta);
    return std::atan(tan_alpha);
}

// NOTE: macarthur_ros() is defined in ERF_MacArthurModel.H — do not duplicate here.

/**
 * @brief Fill fire_ros [m/s] using the Balbi (2009) model
 *
 * Mirrors the GPU kernel structure of compute_ros_field() in ERF_Rothermel.H.
 * When fuel_model_mf and d_balbi_table are both supplied, per-cell Balbi
 * coefficients are looked up from the device table by fuel code; codes outside
 * the table fall back to bc_default.
 *
 * @param[out] fire_ros          Rate of spread [m/s]
 * @param[in]  fire_wind_eff     Effective wind after WAF + terrain [m/s]
 * @param[in]  fire_slopes       Terrain slopes [dz/dx, dz/dy]
 * @param[in]  bc_default        Default Balbi coefficients (uniform fuel)
 * @param[in]  fuel_model_mf     Per-cell fuel model ID (nullptr for uniform)
 * @param[in]  d_balbi_table     Device array of Balbi coeffs indexed by fuel code
 * @param[in]  balbi_table_size  Number of entries in d_balbi_table
 */
inline void fill_balbi_ros(
    amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_wind_eff,
    const amrex::MultiFab& fire_slopes,
    const BalbiComputed& bc_default,
    const amrex::MultiFab* fuel_model_mf = nullptr,
    const BalbiComputed* d_balbi_table = nullptr,
    int balbi_table_size = 0)
{
    const bool use_fuel_table = (fuel_model_mf != nullptr &&
                                 d_balbi_table != nullptr &&
                                 balbi_table_size > 0);

    const amrex::Real A_default  = bc_default.A_coeff;
    const amrex::Real vb_default = bc_default.v_b;

    // Domain-mean wind, used as the fallback for cells whose wind is
    // unphysically large or non-finite.
    const amrex::Long ncells = fire_wind_eff.boxArray().numPts();
    const amrex::Real ux_mean = (ncells > 0) ? fire_wind_eff.sum(0) / amrex::Real(ncells)
                                             : amrex::Real(0.0);
    const amrex::Real uy_mean = (ncells > 0) ? fire_wind_eff.sum(1) / amrex::Real(ncells)
                                             : amrex::Real(0.0);

    for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ros_arr    = fire_ros.array(mfi);
        auto wind_arr   = fire_wind_eff.const_array(mfi);
        auto slopes_arr = fire_slopes.const_array(mfi);

        amrex::Array4<const amrex::Real> fuel_model;
        if (use_fuel_table) {
            fuel_model = fuel_model_mf->const_array(mfi);
        }

        const BalbiComputed* tbl = d_balbi_table;
        const int tbl_sz         = balbi_table_size;

        // Local constexpr copies so the constants are usable from device code
        constexpr amrex::Real v_b_min      = BALBI_V_B_MIN;
        constexpr amrex::Real max_wind_sq  = BALBI_MAX_WIND_SPEED_SQ;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // Wind components. The comparison is written so that a non-finite
            // wind also takes the domain-mean branch.
            amrex::Real u = wind_arr(i, j, k, 0);
            amrex::Real v = wind_arr(i, j, k, 1);
            if (!(u*u + v*v <= max_wind_sq)) {
                u = ux_mean;
                v = uy_mean;
            }
            const amrex::Real wind_speed = std::sqrt(u*u + v*v);

            // Per-cell or global Balbi coefficients
            amrex::Real A_cell  = A_default;
            amrex::Real vb_cell = vb_default;
            if (use_fuel_table) {
                const int fuel_code = static_cast<int>(fuel_model(i, j, k));
                if (fuel_code >= 0 && fuel_code < tbl_sz) {
                    A_cell  = tbl[fuel_code].A_coeff;
                    vb_cell = tbl[fuel_code].v_b;
                }
            }
            vb_cell = amrex::max(vb_cell, v_b_min);

            // Flame tilt from wind and slope: tan α = U/v_b + |∇z|
            const amrex::Real dzdx = slopes_arr(i, j, k, 0);
            const amrex::Real dzdy = slopes_arr(i, j, k, 1);
            const amrex::Real tan_alpha = wind_speed / vb_cell
                                        + std::sqrt(dzdx*dzdx + dzdy*dzdy);

            // sin α = tan α / sqrt(1 + tan²α), cos α = 1 / sqrt(1 + tan²α)
            const amrex::Real inv_sec = 1.0 / std::sqrt(1.0 + tan_alpha*tan_alpha);
            const amrex::Real sin_alpha = tan_alpha * inv_sec;
            const amrex::Real cos_alpha = inv_sec;

            ros_arr(i, j, k, 0) = amrex::max(A_cell * (1.0 + sin_alpha - cos_alpha), 0.0);
        });
    }
}

#endif
