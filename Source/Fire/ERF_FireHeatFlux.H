// Sensible heat flux and fuel consumption for the fire layer.
//
// The heat flux Q [W/m²] per cell is:
//   Q = fuel_load_remain [kg/m²] * h_fuel [J/kg] / tau_res [s]
//
// This is the instantaneous power from the exponential burnout rate:
//   dw/dt = -w / tau_res   =>   Q = -h * dw/dt = h * w / tau_res
//
// tau_res is the effective cell residence time:
//   tau_res [s] = dx_fire [m] / ROS [m/s]   (time for front to cross one fire cell)
//
// This is preferred over the SAV-based Rothermel tau (= BURNOUT_ALPHA * rho_p / sigma),
// which describes particle burnout time (~seconds for fine fuels) rather than
// the grid-scale burn duration. The SAV-based value is used only as a floor.
//
// Fuel load is depleted each atmospheric timestep using the exact exponential:
//   w_new = w_old * exp(-dt_atm / tau_res)
//
// References:
//   Albini, F.A. (1976). Estimating wildfire behavior and effects.
//     USDA Forest Service General Technical Report INT-30.
//   Rothermel, R.C. (1983). How to predict the spread and intensity of
//     forest and range fires. USDA Forest Service GTR INT-143.
//   Byram, G.M. (1959). Combustion of forest fuels. In: Forest Fire:
//     Control and Use, pp. 61-89.
//
// Unit conventions:
//   fuel loads        — lb/ft²  (as returned by get_anderson_fuel_params())
//   fuel SAV          — ft⁻¹   (as returned by get_anderson_fuel_params())
//   rho_p             — lb/ft³  (as returned by get_anderson_fuel_params())
//   heat_content      — BTU/lb  (as returned by get_anderson_fuel_params())
//   fire_fuel_load    — kg/m²   (stored on fire grid after init conversion)
//   fire_heat_flux    — W/m²    (output field on fire grid)

#ifndef ERF_FIRE_HEAT_FLUX_H
#define ERF_FIRE_HEAT_FLUX_H

#include <AMReX_MultiFab.H>
#include <ERF_Rothermel.H>
#include <cmath>

// BURNOUT_ALPHA [s·ft/lb from Rothermel (1983) GTR INT-143].
// tau_i [s] = FIRE_BURNOUT_ALPHA * rho_p [lb/ft³] / sigma_i [ft⁻¹]
constexpr amrex::Real FIRE_BURNOUT_ALPHA = 301.0;

// Fixed SAV values [ft⁻¹] for 10-hr and 100-hr dead fuel classes
// (Anderson 1982, FBFM13 table standard values).
constexpr amrex::Real FIRE_SIGMA_D10  = 109.0;
constexpr amrex::Real FIRE_SIGMA_D100 =  30.0;

/**
 * @brief Compute SAV-based flame residence time [s].
 *
 * tau [s] = FIRE_BURNOUT_ALPHA * rho_p [lb/ft³] / sigma [ft⁻¹]
 * Used as a floor on the cell-based tau_res and for fallback when ROS ~ 0.
 * Returns a floor of 1.0 s.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_residence_time_s(amrex::Real sigma_ft, amrex::Real rho_p_lbft3) noexcept
{
    if (sigma_ft <= 0.0 || rho_p_lbft3 <= 0.0) { return 1.0; }
    amrex::Real tau = FIRE_BURNOUT_ALPHA * rho_p_lbft3 / sigma_ft;
    return amrex::max(tau, 1.0);
}

/**
 * @brief Per-class exponential burnout model (Albini 1976 GTR INT-30).
 *
 * Computes the load-weighted mean total combustion fraction:
 *   f_c_i = 1 - exp(-tau_res / tau_i)
 *
 * Used ONLY by the fire diagnostics path (fireline intensity, flame length).
 * NOT used in the heat-flux or fuel-depletion kernels.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_per_class_consumption_frac(
    amrex::Real tau_res,
    amrex::Real w_d1,   amrex::Real sigma_d1,
    amrex::Real w_d10,
    amrex::Real w_d100,
    amrex::Real w_lh,   amrex::Real sigma_lh,
    amrex::Real w_lw,   amrex::Real sigma_lw,
    amrex::Real rho_p) noexcept
{
    amrex::Real tau_d1   = compute_residence_time_s(sigma_d1, rho_p);
    amrex::Real tau_d10  = compute_residence_time_s(FIRE_SIGMA_D10, rho_p);
    amrex::Real tau_d100 = compute_residence_time_s(FIRE_SIGMA_D100, rho_p);
    amrex::Real tau_lh_adj = compute_residence_time_s(sigma_lh, rho_p) * 1.0;
    amrex::Real tau_lw_adj = compute_residence_time_s(sigma_lw, rho_p) * 2.0;

    amrex::Real f_d1   = 1.0 - std::exp(-tau_res / tau_d1);
    amrex::Real f_d10  = 1.0 - std::exp(-tau_res / tau_d10);
    amrex::Real f_d100 = 1.0 - std::exp(-tau_res / tau_d100);
    amrex::Real f_lh   = 1.0 - std::exp(-tau_res / tau_lh_adj);
    amrex::Real f_lw   = 1.0 - std::exp(-tau_res / tau_lw_adj);

    amrex::Real total_w = w_d1 + w_d10 + w_d100 + w_lh + w_lw;
    if (total_w <= 0.0) {
        amrex::Real dead_w = w_d1 + w_d10 + w_d100;
        if (dead_w <= 1.0e-10) { return 0.5; }
        amrex::Real sigma_agg = (w_d1*sigma_d1 + w_d10*FIRE_SIGMA_D10 + w_d100*FIRE_SIGMA_D100) / dead_w;
        amrex::Real tau_agg = compute_residence_time_s(sigma_agg, rho_p);
        return 1.0 - std::exp(-tau_res / tau_agg);
    }

    amrex::Real f_c = (w_d1*f_d1 + w_d10*f_d10 + w_d100*f_d100 + w_lh*f_lh + w_lw*f_lw) / total_w;
    return amrex::max(0.0, amrex::min(f_c, 1.0));
}

/**
 * @brief Compute sensible heat flux at one cell [W/m²].
 *
 * Q = fuel_load_remain [kg/m²] * h_fuel [J/kg] / tau_res_s [s]
 *
 * Instantaneous power from dw/dt = -w/tau_res. Returns 0 for
 * unburned cells (phi >= 0) or when fuel is exhausted.
 *
 * @param[in] phi                  Level-set value (burned if < 0)
 * @param[in] fuel_load_remain     Remaining fuel load [kg/m²]
 * @param[in] h_fuel_J_per_kg      Heat content [J/kg]
 * @param[in] tau_res_s            Cell residence time [s]
 * @return Heat flux [W/m²], non-negative
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_heat_flux_cell(
    amrex::Real phi,
    amrex::Real fuel_load_remain,
    amrex::Real h_fuel_J_per_kg,
    amrex::Real tau_res_s) noexcept
{
    if (phi >= 0.0)              return 0.0;
    if (fuel_load_remain <= 0.0) return 0.0;
    if (tau_res_s <= 0.0)        return 0.0;
    return amrex::max(0.0, fuel_load_remain * h_fuel_J_per_kg / tau_res_s);
}

/**
 * @brief Deplete fuel load using the exact exponential solution.
 *
 * w_new = w_old * exp(-dt_fire_s / tau_res_s)
 *
 * Exact solution to dw/dt = -w/tau_res. Unconditionally stable.
 * Only applied when phi < 0 (burned cell).
 *
 * @param[in] phi              Level-set value (burned if < 0)
 * @param[in] fuel_load_remain Current fuel load [kg/m²]
 * @param[in] tau_res_s        Cell residence time [s]
 * @param[in] dt_fire_s        Atmospheric timestep [s]
 * @return Updated fuel load [kg/m²], minimum 0
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real deplete_fuel_load(
    amrex::Real phi,
    amrex::Real fuel_load_remain,
    amrex::Real tau_res_s,
    amrex::Real dt_fire_s) noexcept
{
    if (phi >= 0.0)                            return fuel_load_remain;
    if (tau_res_s <= 0.0 || dt_fire_s <= 0.0) return fuel_load_remain;
    return amrex::max(0.0, fuel_load_remain * std::exp(-dt_fire_s / tau_res_s));
}

/**
 * @brief Fill fire_heat_flux and deplete fire_fuel_load for one atmospheric timestep.
 *
 * For each burned cell (phi < 0), the per-cell residence time is computed as:
 *   tau_cell = max(dx_fire / ros_cell, tau_sav_floor)
 * where ros_cell is the local ROS from fire_ros (floored at 1e-6 m/s to avoid
 * division by zero). This avoids the global-max-ROS bug where one boundary
 * artifact cell (with spuriously high ROS) would set tau for all burned cells.
 *
 *   Q    = fuel_load * h_fuel / tau_cell          [W/m²]
 *   w_new = fuel_load * exp(-dt_atm / tau_cell)   [kg/m²]
 *
 * dt_fire_s corresponds to the full atmospheric dt, not the CFL substep.
 *
 * Phase 10: When fuel_model_mf != nullptr, reads per-cell fuel code and computes
 * per-cell heat content. When nullptr, uses global fp.heat_content.
 *
 * @param[out] fire_heat_flux     Heat flux [W/m²]
 * @param[in,out] fire_fuel_load  Fuel load [kg/m²], modified in-place
 * @param[in] fire_phi            Level-set field
 * @param[in] fire_ros            Rate-of-spread field [m/s]
 * @param[in] fp                  Fuel model parameters (Rothermel units)
 * @param[in] dx_fire             Fire grid cell size [m]
 * @param[in] tau_sav_floor       SAV-based residence time floor [s]
 * @param[in] dt_fire_s           Atmospheric timestep [s]
 * @param[in] fuel_model_mf       Optional per-cell fuel model codes (nullptr = uniform)
 */
inline void fill_fire_heat_flux(
    amrex::MultiFab& fire_heat_flux,
    amrex::MultiFab& fire_fuel_load,
    const amrex::MultiFab& fire_phi,
    const amrex::MultiFab& fire_ros,
    const FuelModelParams& fp,
    amrex::Real dx_fire,
    amrex::Real tau_sav_floor,
    amrex::Real dt_fire_s,
    const amrex::MultiFab* fuel_model_mf = nullptr)
{
    using namespace amrex;

    // 1 BTU/lb = 2326 J/kg
    // Compute heat content based on whether spatial fuel map is used
    Real h_fuel_J_per_kg;

    if (fuel_model_mf == nullptr) {
        // Uniform fuel: use global parameters
        h_fuel_J_per_kg = fp.heat_content * 2326.0;
    } else {
        h_fuel_J_per_kg = 0.0;  // Will be set per-cell
    }

    for (MFIter mfi(fire_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto const phi_arr  = fire_phi.const_array(mfi);
        auto const ros_arr  = fire_ros.const_array(mfi);
        auto Q_arr          = fire_heat_flux.array(mfi);
        auto fuel_arr       = fire_fuel_load.array(mfi);

        if (fuel_model_mf != nullptr) {
            auto const fuel_code_arr = fuel_model_mf->const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
                int i = iv[0], j = iv[1], k = iv[2];

                Real phi_val = phi_arr(i, j, k);
                Real w_cur   = fuel_arr(i, j, k);

                // Get per-cell fuel code and look up heat content
                int fuel_code = static_cast<int>(fuel_code_arr(i, j, k) + 0.5_rt);
                FuelModelParams fp_cell = get_anderson_fuel_params(fuel_code);
                Real h_fuel_cell = fp_cell.heat_content * 2326.0_rt;

                // Per-cell residence time: floor ROS at 1e-6 m/s to avoid divide-by-zero
                Real ros_cell  = amrex::max(ros_arr(i, j, k), 1.0e-6_rt);
                Real tau_cell  = amrex::max(dx_fire / ros_cell, tau_sav_floor);

                Q_arr(i, j, k)    = compute_heat_flux_cell(phi_val, w_cur, h_fuel_cell, tau_cell);
                fuel_arr(i, j, k) = deplete_fuel_load(phi_val, w_cur, tau_cell, dt_fire_s);
            });
        } else {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
                int i = iv[0], j = iv[1], k = iv[2];

                Real phi_val = phi_arr(i, j, k);
                Real w_cur   = fuel_arr(i, j, k);

                // Per-cell residence time: floor ROS at 1e-6 m/s to avoid divide-by-zero
                Real ros_cell  = amrex::max(ros_arr(i, j, k), 1.0e-6_rt);
                Real tau_cell  = amrex::max(dx_fire / ros_cell, tau_sav_floor);

                Q_arr(i, j, k)    = compute_heat_flux_cell(phi_val, w_cur, h_fuel_J_per_kg, tau_cell);
                fuel_arr(i, j, k) = deplete_fuel_load(phi_val, w_cur, tau_cell, dt_fire_s);
            });
        }
    }
}

#endif