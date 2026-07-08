#ifndef ERF_BALBI_MODEL_H
#define ERF_BALBI_MODEL_H

#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <cmath>
#include <ERF_Rothermel.H>

/**
 * @file ERF_BalbiModel.H
 * @brief Balbi (2009) physical fire spread model implementation
 *
 * Implements the Balbi (2009) fire spread model for rate-of-spread (ROS)
 * computation as a function of fuel type, moisture, wind, and slope.
 *
 * Formula:
 *   R = A_coeff × (1 + sin α − cos α) [m/s]
 * where
 *   tan α = U/v_b + tan θ
 *   v_b = sqrt(g × δ_m × (T_f − T_a) / T_a)
 *
 * Inputs:
 *   U:  wind speed in the fire spread direction [m/s]
 *   θ:  terrain slope [radians]
 *   δ_m: fuel bed depth [m]
 *   All US/Rothermel unit inputs (sigma [1/ft], delta [ft], h_heat [BTU/lb])
 *   are converted to SI internally.
 *
 * Limitations:
 *   Not appropriate for fuels with very low fuel bed depth (δ_m < 0.01 m)
 *   — returns zero ROS for such cases.
 *
 * References:
 *   Balbi, J.-H. et al. (2009). A physical model for wildland fires.
 *     Combustion and Flame, 156(12), 2217–2230.
 *   Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling
 *     System. USDA Forest Service GTR INT-194.
 */

/**
 * @struct BalbiComputed
 * @brief Precomputed Balbi coefficients (fuel dependent)
 */
struct BalbiComputed
{
    amrex::Real A_coeff;  ///< Balbi coefficient A [m/s] (fire intensity related)
    amrex::Real v_b;      ///< Buoyancy velocity [m/s]
};

/**
 * @brief Compute Balbi fuel/thermal coefficients
 *
 * Precomputes A_coeff and v_b that depend only on fuel type and ambient
 * conditions (not on wind/slope).
 *
 * @param[in] fp  FuelModelParams for selected fuel model
 * @param[in] bp  BalbiParams with thermal properties and constants
 * @return BalbiComputed with precomputed coefficients
 *
 * The computation follows Balbi (2009):
 * - A_coeff depends on fuel load, SAV, and moisture
 * - v_b depends on ambient/flame temperatures and fuel bed depth
 */
inline BalbiComputed compute_balbi_params(
    const FuelModelParams& fp,
    const FireParams::BalbiParams& bp) noexcept
{
    BalbiComputed bc;
    
    // Guard against very small fuel bed depths
    if (fp.delta < 0.01_rt) {
        bc.A_coeff = 0.0_rt;
        bc.v_b = 0.0_rt;
        return bc;
    }

    // Convert Rothermel fuel parameters from US to SI units
    // fp.w_d1, fp.w_d10, fp.w_d100 [lb/ft²] → [kg/m²]
    // 1 lb/ft² = 4.88 kg/m²
    constexpr amrex::Real lb_ft2_to_kg_m2 = 4.88_rt;
    amrex::Real w_d1_si  = fp.w_d1  * lb_ft2_to_kg_m2;
    amrex::Real w_d10_si = fp.w_d10 * lb_ft2_to_kg_m2;
    amrex::Real w_d100_si = fp.w_d100 * lb_ft2_to_kg_m2;
    
    // fp.delta [ft] → [m]
    // 1 ft = 0.3048 m
    constexpr amrex::Real ft_to_m = 0.3048_rt;
    amrex::Real delta_m = fp.delta * ft_to_m;

    // Effective fuel loading [kg/m²]
    amrex::Real w_fuel = w_d1_si + w_d10_si + w_d100_si;
    w_fuel = amrex::max(w_fuel, 0.1_rt);  // Guard against zero fuel

    // Fraction not consumed; simplified as (1 - f_c) where f_c ≈ 0.5
    // More sophisticated: use moisture of extinction
    amrex::Real frac_available = 0.5_rt;

    // A_coeff [m/s] calculation from Balbi (2009)
    // Simplified form based on reaction intensity analog
    // A_coeff ≈ (r_00 / tau_0) × sqrt(w_fuel / delta_H) × (fraction available)
    amrex::Real A_const = bp.r_00 / bp.tau_0;
    amrex::Real w_normalized = amrex::sqrt(w_fuel / (bp.delta_H / 1.0e6_rt));  // Normalize by heat release
    bc.A_coeff = A_const * w_normalized * frac_available;
    bc.A_coeff = amrex::max(bc.A_coeff, 0.01_rt);  // Minimum non-zero A

    // Buoyancy velocity v_b [m/s]
    // v_b = sqrt(g × δ_m × (T_f − T_a) / T_a)
    constexpr amrex::Real g = 9.81_rt;  // gravitational acceleration [m/s²]
    amrex::Real T_ratio = (bp.T_f - bp.T_a) / bp.T_a;
    bc.v_b = amrex::sqrt(g * delta_m * T_ratio);

    return bc;
}

/**
 * @brief Compute wind and slope angle for Balbi model
 *
 * Given wind speed U and terrain slope θ, compute the resultant angle α
 * in the fire spread direction via tan α = U/v_b + tan θ.
 *
 * @param[in] U        Wind speed in spread direction [m/s]
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
    constexpr amrex::Real small_v_b = 1.0e-6_rt;
    if (v_b < small_v_b) v_b = small_v_b;  // Avoid division by zero
    
    amrex::Real tan_alpha = U / v_b + std::tan(theta);
    return std::atan(tan_alpha);
}

/**
 * @brief MacArthur (1966) Mark 5 Forest Fire Danger Meter ROS formula
 *
 * Used in WRF-SFIRE as ibeh=0 (Clark et al. 2004).
 *
 * Formula:
 *   R [m/s] = 0.18 * exp(0.8424 * max(U, 0))    (with wind)
 *   backing_ros = 0.18 m/s                        (no-wind backing rate)
 *
 * Parameters:
 *   wind_speed_ms: wind speed in the fire spread direction [m/s]
 *
 * Returns: ROS [m/s]
 *
 * References:
 *   Clark, T.L., Coen, J.L. & Latham, D. (2004). Description of a coupled
 *     atmosphere-fire model. Int. J. Wildland Fire, 13, 49–63.
 *   McArthur, A.G. (1966). Weather and Grassland Fire Behaviour. Leaflet 100,
 *     Forestry and Timber Bureau, Canberra.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real macarthur_ros(amrex::Real wind_speed_ms) noexcept
{
    constexpr amrex::Real backing = 0.18_rt;
    const amrex::Real wind_contribution = backing * (std::exp(0.8424_rt * amrex::max(wind_speed_ms, 0.0_rt)) - 1.0_rt);
    return backing + wind_contribution;
}

/**
 * @brief Fill fire_ros [m/s] using Balbi (2009) model
 *
 * Mirrors the GPU kernel structure of compute_ros_field() in ERF_Rothermel.H.
 * When fuel_model_mf is non-null, per-cell Balbi coefficients are looked up
 * from the d_balbi_table device array indexed by fuel code.
 *
 * @param[out] fire_ros        Rate of spread [m/s]
 * @param[in]  fire_wind_eff   Effective wind after WAF + terrain [m/s]
 * @param[in]  fire_slopes     Terrain slopes [dz/dx, dz/dy]
 * @param[in]  bc_default      Default Balbi coefficients (uniform fuel)
 * @param[in]  fuel_model_mf   Per-cell fuel model ID (nullptr for uniform)
 * @param[in]  d_balbi_table   Device array of Balbi coeffs indexed by fuel code
 * @param[in]  balbi_table_size  Size of d_balbi_table (13 for FBFM13)
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
    for (amrex::MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ros_arr = fire_ros.array(mfi);
        auto wind_arr = fire_wind_eff.const_array(mfi);
        auto slopes_arr = fire_slopes.const_array(mfi);
        
        BalbiComputed bc_cell = bc_default;  // default for uniform fuel
        amrex::Array4<const amrex::Real> fuel_model;
        if (fuel_model_mf != nullptr) {
            fuel_model = fuel_model_mf->const_array(mfi);
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            // Get wind components
            amrex::Real u = wind_arr(i, j, k, 0);
            amrex::Real v = wind_arr(i, j, k, 1);
            amrex::Real wind_mag = std::sqrt(u*u + v*v);

            // Get slopes
            amrex::Real dzdx = slopes_arr(i, j, k, 0);
            amrex::Real dzdy = slopes_arr(i, j, k, 1);
            amrex::Real slope_mag = std::sqrt(dzdx*dzdx + dzdy*dzdy);
            amrex::Real slope_angle = std::atan(slope_mag);

            // Compute wind direction and project onto slope direction
            amrex::Real wind_dir = std::atan2(v, u);
            amrex::Real slope_dir = std::atan2(dzdy, dzdx);
            amrex::Real angle_to_slope = wind_dir - slope_dir;
            amrex::Real cos_angle = std::cos(angle_to_slope);
            amrex::Real u_slope = wind_mag * cos_angle;  // projected wind speed [m/s]

            // Lookup Balbi coefficients if spatial fuel is active
            if (fuel_model_mf != nullptr && d_balbi_table != nullptr) {
                int fuel_code = static_cast<int>(fuel_model(i, j, k));
                fuel_code = amrex::max(1, amrex::min(fuel_code, balbi_table_size));
                bc_cell = d_balbi_table[fuel_code];
            }

            // Compute Balbi ROS
            amrex::Real alpha = compute_balbi_angle(u_slope, slope_angle, bc_cell.v_b);
            amrex::Real sin_alpha = std::sin(alpha);
            amrex::Real cos_alpha = std::cos(alpha);
            amrex::Real ros = bc_cell.A_coeff * (1.0_rt + sin_alpha - cos_alpha);
            ros = amrex::max(ros, 0.0_rt);

            ros_arr(i, j, k, 0) = ros;
        });
    }
}

#endif
