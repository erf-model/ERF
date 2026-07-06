// Sensible heat flux and fuel consumption for the fire layer.
//
// The heat flux Q [W/m²] per cell is:
//   Q = ROS [m/s] * fuel_load_remain [kg/m²] * h_fuel [J/kg] * f_c / tau_res [s]
//
// where f_c is the combustion fraction (Albini 1976, GTR INT-30) and
// tau_res is the flame residence time (Rothermel 1983, GTR INT-143):
//   tau_res [s] = BURNOUT_ALPHA * rho_p [lb/ft³] / sigma [ft⁻¹]
//   BURNOUT_ALPHA = 301.0 s·ft/lb
//
// Fuel load is depleted each fire substep:
//   delta_w = fuel_load * f_c * (dt_fire / tau_res)
//
// This file computes the convective (sensible) heat component only.
//
// References:
//   Albini, F.A. (1976). Estimating wildfire behavior and effects.
//     USDA Forest Service General Technical Report INT-30.
//   Rothermel, R.C. (1983). How to predict the spread and intensity of
//     forest and range fires. USDA Forest Service GTR INT-143.
//   Byram, G.M. (1959). Combustion of forest fuels. In: Forest Fire:
//     Control and Use, pp. 61-89.
//
// Ported from hgopalan/wildfire_levelset src/bulk_fuel_consumption.H.
//
// Unit conventions used in this file:
//   fuel loads        — lb/ft²  (as returned by get_anderson_fuel_params())
//   fuel SAV          — ft⁻¹   (as returned by get_anderson_fuel_params())
//   rho_p             — lb/ft³  (as returned by get_anderson_fuel_params())
//   heat_content      — BTU/lb  (as returned by get_anderson_fuel_params())
//   fire_fuel_load    — kg/m²   (stored on fire grid after init conversion)
//   fire_heat_flux    — W/m²    (output field on fire grid)
//   Unit conversions applied once per call, not per cell.
//
// Assumptions:
//   Only cells where fire_phi < 0 (burned) contribute to heat flux.
//   The fire grid is 2D (C*nx x C*ny x 1); k index is always 0.
//   A uniform tau_res derived from the aggregate SAV of the single fuel model
//   is applied across all cells. Per-cell spatial fuel variation is not
//   included in this phase.
//   FuelModelParams does not carry an aggregate sigma field; it is computed
//   here as the load-weighted mean of dead fuel classes.

#ifndef ERF_FIRE_HEAT_FLUX_H
#define ERF_FIRE_HEAT_FLUX_H

#include <AMReX_MultiFab.H>
#include <ERF_Rothermel.H>
#include <cmath>

// BURNOUT_ALPHA [s·ft/lb from Rothermel (1983) GTR INT-143].
// tau_i [s] = FIRE_BURNOUT_ALPHA * rho_p [lb/ft³] / sigma_i [ft⁻¹]
constexpr amrex::Real FIRE_BURNOUT_ALPHA = 301.0_rt;

// Fixed SAV values [ft⁻¹] for 10-hr and 100-hr dead fuel classes
// (Anderson 1982, FBFM13 table standard values).
constexpr amrex::Real FIRE_SIGMA_D10  = 109.0_rt;
constexpr amrex::Real FIRE_SIGMA_D100 =  30.0_rt;

/**
 * @brief Compute flame residence time [s] from fuel SAV and particle density.
 *
 * tau [s] = FIRE_BURNOUT_ALPHA * rho_p [lb/ft³] / sigma [ft⁻¹]
 * Returns a floor of 1.0 s to avoid division by zero in callers.
 *
 * @param[in] sigma_ft      Surface area to volume ratio [ft⁻¹]
 * @param[in] rho_p_lbft3   Fuel particle density [lb/ft³]
 * @return Residence time [s], minimum 1.0 s
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_residence_time_s(amrex::Real sigma_ft, amrex::Real rho_p_lbft3) noexcept
{
    if (sigma_ft <= 0.0_rt || rho_p_lbft3 <= 0.0_rt) {
        return 1.0_rt;
    }
    amrex::Real tau = FIRE_BURNOUT_ALPHA * rho_p_lbft3 / sigma_ft;
    return amrex::max(tau, 1.0_rt);
}

/**
 * @brief Per-class exponential burnout model (Albini 1976 GTR INT-30).
 *
 * Computes the load-weighted mean combustion fraction across all fuel classes
 * using the exponential burnout function:
 *   f_c_i = 1 - exp(-tau_res / tau_i)  where tau_i = BURNOUT_ALPHA * rho_p / sigma_i
 *
 * The calculation follows Albini (1976) and is ported directly from
 * wildfire_levelset/src/bulk_fuel_consumption.H.
 *
 * For live fuels (herbaceous and woody), calibration multipliers are applied:
 *   tau_lh_adj = tau_lh * 1.0  (1-hr equivalent residence time for herbaceous)
 *   tau_lw_adj = tau_lw * 2.0  (calibration multiplier for woody from wildfire_levelset)
 *
 * Falls back to single-class aggregate when all per-class loads are zero.
 * The aggregate fallback uses load-weighted SAV for dead classes only:
 *   sigma_agg = (w_d1*sigma_d1 + w_d10*FIRE_SIGMA_D10 + w_d100*FIRE_SIGMA_D100)
 *             / max(w_d1+w_d10+w_d100, 1e-10)
 *
 * @param[in] tau_res       Pre-computed aggregate residence time [s]
 * @param[in] w_d1          1-hour dead fuel load [lb/ft²]
 * @param[in] sigma_d1      1-hour dead SAV [ft⁻¹]
 * @param[in] w_d10         10-hour dead fuel load [lb/ft²]
 * @param[in] w_d100        100-hour dead fuel load [lb/ft²]
 * @param[in] w_lh          Live herbaceous fuel load [lb/ft²]
 * @param[in] sigma_lh      Live herbaceous SAV [ft⁻¹]
 * @param[in] w_lw          Live woody fuel load [lb/ft²]
 * @param[in] sigma_lw      Live woody SAV [ft⁻¹]
 * @param[in] rho_p         Fuel particle density [lb/ft³]
 * @return Load-weighted mean combustion fraction [0-1]
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
    // Compute per-class residence times and combustion fractions
    amrex::Real tau_d1  = compute_residence_time_s(sigma_d1, rho_p);
    amrex::Real tau_d10 = compute_residence_time_s(FIRE_SIGMA_D10, rho_p);
    amrex::Real tau_d100= compute_residence_time_s(FIRE_SIGMA_D100, rho_p);

    // Live fuel residence times with calibration multipliers from wildfire_levelset
    amrex::Real tau_lh_adj = compute_residence_time_s(sigma_lh, rho_p) * 1.0_rt;   // herbaceous multiplier = 1.0
    amrex::Real tau_lw_adj = compute_residence_time_s(sigma_lw, rho_p) * 2.0_rt;   // woody multiplier = 2.0

    // Exponential burnout: f_c_i = 1 - exp(-tau_res / tau_i)
    amrex::Real f_d1   = 1.0_rt - std::exp(-tau_res / tau_d1);
    amrex::Real f_d10  = 1.0_rt - std::exp(-tau_res / tau_d10);
    amrex::Real f_d100 = 1.0_rt - std::exp(-tau_res / tau_d100);
    amrex::Real f_lh   = 1.0_rt - std::exp(-tau_res / tau_lh_adj);
    amrex::Real f_lw   = 1.0_rt - std::exp(-tau_res / tau_lw_adj);

    // Compute total load and weighted mean combustion fraction
    amrex::Real total_w = w_d1 + w_d10 + w_d100 + w_lh + w_lw;
    if (total_w <= 0.0_rt) {
        // Fallback: compute aggregate SAV from dead fuels only
        amrex::Real dead_w = w_d1 + w_d10 + w_d100;
        if (dead_w <= 1.0e-10_rt) {
            return 0.5_rt;  // Default combustion fraction when no fuel
        }
        amrex::Real sigma_agg = (w_d1 * sigma_d1 + w_d10 * FIRE_SIGMA_D10 + w_d100 * FIRE_SIGMA_D100) / dead_w;
        amrex::Real tau_agg = compute_residence_time_s(sigma_agg, rho_p);
        return 1.0_rt - std::exp(-tau_res / tau_agg);
    }

    // Load-weighted mean of per-class combustion fractions
    amrex::Real f_c = (w_d1 * f_d1 + w_d10 * f_d10 + w_d100 * f_d100 + w_lh * f_lh + w_lw * f_lw) / total_w;
    return amrex::max(0.0_rt, amrex::min(f_c, 1.0_rt));
}

/**
 * @brief Compute sensible heat flux at one cell [W/m²].
 *
 * Q = ros * fuel_load_remain * h_fuel_J_per_kg * f_c / tau_res_s
 *
 * Returns 0 if phi >= 0 (unburned) or fuel_load_remain <= 0 or ros <= 0.
 *
 * @param[in] phi                  Level-set field value (burned if < 0)
 * @param[in] ros                  Rate of spread [m/s]
 * @param[in] fuel_load_remain     Remaining fuel load [kg/m²]
 * @param[in] h_fuel_J_per_kg      Heat content [J/kg]
 * @param[in] f_c                  Combustion fraction [0-1]
 * @param[in] tau_res_s            Residence time [s]
 * @return Heat flux [W/m²], non-negative
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_heat_flux_cell(
    amrex::Real phi,
    amrex::Real ros,
    amrex::Real fuel_load_remain,
    amrex::Real h_fuel_J_per_kg,
    amrex::Real f_c,
    amrex::Real tau_res_s) noexcept
{
    // Only burned cells contribute to heat flux
    if (phi >= 0.0_rt) return 0.0_rt;

    // Heat flux is zero if fuel, ROS, or combustion is negligible
    if (fuel_load_remain <= 0.0_rt || ros <= 0.0_rt || f_c <= 0.0_rt) {
        return 0.0_rt;
    }

    // Guard against zero or negative tau_res
    if (tau_res_s <= 0.0_rt) return 0.0_rt;

    // Q = ros * fuel_load * h_fuel * f_c / tau_res
    amrex::Real Q = ros * fuel_load_remain * h_fuel_J_per_kg * f_c / tau_res_s;
    return amrex::max(0.0_rt, Q);
}

/**
 * @brief Deplete fuel load for one fire substep using forward Euler.
 *
 * delta_w = fuel_load * f_c * (dt_fire_s / tau_res_s)
 *
 * Returns updated fuel load [kg/m²], floored at 0.
 * Only applied when phi < 0 (burned cell).
 *
 * @param[in] phi              Level-set field value (burned if < 0)
 * @param[in] fuel_load_remain Current fuel load [kg/m²]
 * @param[in] f_c              Combustion fraction [0-1]
 * @param[in] tau_res_s        Residence time [s]
 * @param[in] dt_fire_s        Fire timestep [s]
 * @return Updated fuel load [kg/m²], minimum 0
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real deplete_fuel_load(
    amrex::Real phi,
    amrex::Real fuel_load_remain,
    amrex::Real f_c,
    amrex::Real tau_res_s,
    amrex::Real dt_fire_s) noexcept
{
    // Only burned cells deplete fuel
    if (phi >= 0.0_rt) return fuel_load_remain;

    // No depletion if tau_res or dt are invalid
    if (tau_res_s <= 0.0_rt || dt_fire_s <= 0.0_rt) return fuel_load_remain;

    // Euler: w_new = w_old - w_old * f_c * dt / tau
    amrex::Real delta_w = fuel_load_remain * f_c * (dt_fire_s / tau_res_s);
    amrex::Real w_new = fuel_load_remain - delta_w;

    return amrex::max(0.0_rt, w_new);
}

/**
 * @brief Fill fire_heat_flux and deplete fire_fuel_load for one atmospheric timestep.
 *
 * For each cell where phi < 0 (burned):
 *   1. Compute f_c = compute_per_class_consumption_frac(tau_res_s, ...).
 *   2. h_fuel_J_per_kg = fp.heat_content [BTU/lb] * 2326.0  (1 BTU/lb = 2326 J/kg)
 *   3. Q = compute_heat_flux_cell(phi, ros, fuel_load, h_fuel_J_per_kg, f_c, tau_res_s)
 *   4. fire_heat_flux(i,j,k) = Q
 *   5. fire_fuel_load(i,j,k) = deplete_fuel_load(phi, fuel_load, f_c, tau_res_s, dt_fire_s)
 *
 * For unburned cells (phi >= 0): fire_heat_flux = 0, fuel_load unchanged.
 *
 * dt_fire_s is the fire substep size.
 * fp contains all fuel model parameters in original Rothermel units.
 * tau_res_s is the pre-computed residence time [s].
 * All unit conversions (BTU/lb to J/kg, lb/ft² to kg/m²) are applied
 * once as scalars before the GPU kernel, not per cell.
 *
 * This is a local operation: no MPI communication required.
 *
 * @param[out] fire_heat_flux     Heat flux field [W/m²]
 * @param[in,out] fire_fuel_load  Fuel load field [kg/m²], modified in-place
 * @param[in] fire_phi            Level-set field
 * @param[in] fire_ros            Rate of spread [m/s]
 * @param[in] fp                  Fuel model parameters
 * @param[in] tau_res_s           Residence time [s]
 * @param[in] dt_fire_s           Fire substep duration [s]
 */
inline void fill_fire_heat_flux(
    amrex::MultiFab& fire_heat_flux,
    amrex::MultiFab& fire_fuel_load,
    const amrex::MultiFab& fire_phi,
    const amrex::MultiFab& fire_ros,
    const FuelModelParams& fp,
    amrex::Real tau_res_s,
    amrex::Real dt_fire_s)
{
    using namespace amrex;

    // Convert heat content from BTU/lb to J/kg: 1 BTU/lb = 2326 J/kg
    Real h_fuel_J_per_kg = fp.heat_content * 2326.0_rt;

    // Compute combustion fraction (called once here, reused for all cells)
    Real f_c = compute_per_class_consumption_frac(
        tau_res_s,
        fp.w_d1, fp.sigma_d1,
        fp.w_d10,
        fp.w_d100,
        fp.w_lh, fp.sigma_lh,
        fp.w_lw, fp.sigma_lw,
        fp.rho_p);

    // GPU kernel over all cells
    for (MFIter mfi(fire_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto const phi_arr = fire_phi.const_array(mfi);
        auto const ros_arr = fire_ros.const_array(mfi);
        auto Q_arr = fire_heat_flux.array(mfi);
        auto fuel_arr = fire_fuel_load.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0], j = iv[1], k = iv[2];

            Real phi_val = phi_arr(i, j, k);
            Real ros_val = ros_arr(i, j, k);
            Real w_cur = fuel_arr(i, j, k);

            // Compute heat flux
            Real Q = compute_heat_flux_cell(phi_val, ros_val, w_cur, h_fuel_J_per_kg, f_c, tau_res_s);
            Q_arr(i, j, k) = Q;

            // Deplete fuel load
            Real w_new = deplete_fuel_load(phi_val, w_cur, f_c, tau_res_s, dt_fire_s);
            fuel_arr(i, j, k) = w_new;
        });
    }
}

#endif
