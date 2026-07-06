// Sensible heat flux and fuel consumption for the fire layer.
//
// The heat flux Q [W/m²] per cell is:
//   Q = fuel_load_remain [kg/m²] * h_fuel [J/kg] / tau_res [s]
//
// This is the instantaneous power from the exponential burnout rate:
//   dw/dt = -w / tau_res   =>   Q = -h * dw/dt = h * w / tau_res
//
// tau_res is the flame residence time (Rothermel 1983, GTR INT-143):
//   tau_res [s] = BURNOUT_ALPHA * rho_p [lb/ft³] / sigma [ft⁻¹]
//   BURNOUT_ALPHA = 301.0 s·ft/lb
//
// Fuel load is depleted each fire substep using the exact exponential solution:
//   w_new = w_old * exp(-dt_fire / tau_res)
//
// Note: the combustion fraction f_c = 1 - exp(-tau_res/tau_i) (Albini 1976)
// represents the total fraction consumed over one full residence time and is
// used only for the fire diagnostics path (fireline intensity / flame length),
// NOT in the heat-flux or depletion formulas here. Using f_c as a rate
// multiplier would double-count the asymptotic burnout fraction.
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
constexpr amrex::Real FIRE_BURNOUT_ALPHA = 301.0;

// Fixed SAV values [ft⁻¹] for 10-hr and 100-hr dead fuel classes
// (Anderson 1982, FBFM13 table standard values).
constexpr amrex::Real FIRE_SIGMA_D10  = 109.0;
constexpr amrex::Real FIRE_SIGMA_D100 =  30.0;

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
    if (sigma_ft <= 0.0 || rho_p_lbft3 <= 0.0) {
        return 1.0;
    }
    amrex::Real tau = FIRE_BURNOUT_ALPHA * rho_p_lbft3 / sigma_ft;
    return amrex::max(tau, 1.0);
}

/**
 * @brief Per-class exponential burnout model (Albini 1976 GTR INT-30).
 *
 * Computes the load-weighted mean combustion fraction across all fuel classes
 * using the exponential burnout function:
 *   f_c_i = 1 - exp(-tau_res / tau_i)  where tau_i = BURNOUT_ALPHA * rho_p / sigma_i
 *
 * This function is retained for use by the fire DIAGNOSTICS path (fireline
 * intensity, flame length) which needs the total fraction consumed.
 * It is NOT used in the heat-flux or fuel-depletion kernels.
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
    // Compute per-class residence times
    amrex::Real tau_d1   = compute_residence_time_s(sigma_d1, rho_p);
    amrex::Real tau_d10  = compute_residence_time_s(FIRE_SIGMA_D10, rho_p);
    amrex::Real tau_d100 = compute_residence_time_s(FIRE_SIGMA_D100, rho_p);

    // Live fuel residence times with calibration multipliers from wildfire_levelset
    amrex::Real tau_lh_adj = compute_residence_time_s(sigma_lh, rho_p) * 1.0;  // herbaceous multiplier = 1.0
    amrex::Real tau_lw_adj = compute_residence_time_s(sigma_lw, rho_p) * 2.0;  // woody multiplier = 2.0

    // Exponential burnout: f_c_i = 1 - exp(-tau_res / tau_i)
    amrex::Real f_d1   = 1.0 - std::exp(-tau_res / tau_d1);
    amrex::Real f_d10  = 1.0 - std::exp(-tau_res / tau_d10);
    amrex::Real f_d100 = 1.0 - std::exp(-tau_res / tau_d100);
    amrex::Real f_lh   = 1.0 - std::exp(-tau_res / tau_lh_adj);
    amrex::Real f_lw   = 1.0 - std::exp(-tau_res / tau_lw_adj);

    // Compute total load and weighted mean combustion fraction
    amrex::Real total_w = w_d1 + w_d10 + w_d100 + w_lh + w_lw;
    if (total_w <= 0.0) {
        // Fallback: compute aggregate SAV from dead fuels only
        amrex::Real dead_w = w_d1 + w_d10 + w_d100;
        if (dead_w <= 1.0e-10) {
            return 0.5;  // Default combustion fraction when no fuel data
        }
        amrex::Real sigma_agg = (w_d1 * sigma_d1 + w_d10 * FIRE_SIGMA_D10 + w_d100 * FIRE_SIGMA_D100) / dead_w;
        amrex::Real tau_agg = compute_residence_time_s(sigma_agg, rho_p);
        return 1.0 - std::exp(-tau_res / tau_agg);
    }

    // Load-weighted mean of per-class combustion fractions
    amrex::Real f_c = (w_d1 * f_d1 + w_d10 * f_d10 + w_d100 * f_d100 + w_lh * f_lh + w_lw * f_lw) / total_w;
    return amrex::max(0.0, amrex::min(f_c, 1.0));
}

/**
 * @brief Compute sensible heat flux at one cell [W/m²].
 *
 * Uses the instantaneous power from the exponential burnout rate:
 *   Q = fuel_load_remain [kg/m²] * h_fuel [J/kg] / tau_res_s [s]
 *
 * This follows directly from dw/dt = -w/tau_res, so the instantaneous
 * heat release rate is Q = h * |dw/dt| = h * w / tau_res.
 * f_c is NOT used here; it belongs to the total-consumption diagnostics,
 * not the instantaneous rate formula.
 *
 * Returns 0 if phi >= 0 (unburned) or fuel_load_remain <= 0.
 *
 * @param[in] phi                  Level-set field value (burned if < 0)
 * @param[in] fuel_load_remain     Remaining fuel load [kg/m²]
 * @param[in] h_fuel_J_per_kg      Heat content [J/kg]
 * @param[in] tau_res_s            Residence time [s]
 * @return Heat flux [W/m²], non-negative
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_heat_flux_cell(
    amrex::Real phi,
    amrex::Real fuel_load_remain,
    amrex::Real h_fuel_J_per_kg,
    amrex::Real tau_res_s) noexcept
{
    // Only burned cells contribute to heat flux
    if (phi >= 0.0) return 0.0;

    // Heat flux is zero if fuel or combustion is negligible
    if (fuel_load_remain <= 0.0_rt || f_c <= 0.0_rt) {
        return 0.0_rt;
    }

    // Guard against zero or negative tau_res
    if (tau_res_s <= 0.0) return 0.0;

    // Q = fuel_load * h_fuel * f_c / tau_res
    amrex::Real Q = fuel_load_remain * h_fuel_J_per_kg * f_c / tau_res_s;
    return amrex::max(0.0_rt, Q);
}

/**
 * @brief Deplete fuel load for one timestep using the exact exponential solution.
 *
 * w_new = w_old * exp(-dt_fire_s / tau_res_s)
 *
 * This is the exact solution to dw/dt = -w/tau_res over one timestep.
 * It is unconditionally stable for any dt and avoids the spurious f_c
 * multiplier that the previous forward-Euler formulation carried.
 *
 * Returns updated fuel load [kg/m²], floored at 0.
 * Only applied when phi < 0 (burned cell).
 *
 * @param[in] phi              Level-set field value (burned if < 0)
 * @param[in] fuel_load_remain Current fuel load [kg/m²]
 * @param[in] tau_res_s        Residence time [s]
 * @param[in] dt_fire_s        Fire timestep [s]
 * @return Updated fuel load [kg/m²], minimum 0
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real deplete_fuel_load(
    amrex::Real phi,
    amrex::Real fuel_load_remain,
    amrex::Real tau_res_s,
    amrex::Real dt_fire_s) noexcept
{
    // Only burned cells deplete fuel
    if (phi >= 0.0) return fuel_load_remain;

    // No depletion if tau_res or dt are invalid
    if (tau_res_s <= 0.0 || dt_fire_s <= 0.0) return fuel_load_remain;

    // Exact exponential: w_new = w_old * exp(-dt / tau_res)
    amrex::Real w_new = fuel_load_remain * std::exp(-dt_fire_s / tau_res_s);
    return amrex::max(0.0, w_new);
}

/**
 * @brief Fill fire_heat_flux and deplete fire_fuel_load for one atmospheric timestep.
 *
 * For each cell where phi < 0 (burned):
 *   1. h_fuel_J_per_kg = fp.heat_content [BTU/lb] * 2326.0  (1 BTU/lb = 2326 J/kg)
 *   2. Q = compute_heat_flux_cell(phi, fuel_load, h_fuel_J_per_kg, tau_res_s)
 *        = fuel_load * h_fuel / tau_res   [W/m²]
 *   3. fire_heat_flux(i,j,k) = Q
 *   4. fire_fuel_load(i,j,k) = deplete_fuel_load(phi, fuel_load, tau_res_s, dt_fire_s)
 *        = fuel_load * exp(-dt / tau_res)
 *
 * For unburned cells (phi >= 0): fire_heat_flux = 0, fuel_load unchanged.
 *
 * The combustion fraction f_c (Albini 1976) is NOT used here. It represents
 * the total fraction consumed over one full residence time and is used only
 * by the diagnostics path (ERF_FireDiagnostics.H). Using it as a rate
 * multiplier in depletion double-counts the asymptotic burnout and causes
 * fuel to exhaust ~1.6x faster than the correct exponential rate.
 *
 * dt_fire_s is the atmospheric timestep (passed as `dt` from FireLayer::advance,
 * NOT the CFL-limited fire substep). fill_fire_heat_flux is called once per
 * atmospheric step, so dt_fire_s should be the full atmospheric dt.
 *
 * fp contains all fuel model parameters in original Rothermel units.
 * tau_res_s is the pre-computed residence time [s].
 * All unit conversions (BTU/lb to J/kg) are applied once before the GPU
 * kernel, not per cell.
 *
 * This is a local operation: no MPI communication required.
 *
 * @param[out] fire_heat_flux     Heat flux field [W/m²]
 * @param[in,out] fire_fuel_load  Fuel load field [kg/m²], modified in-place
 * @param[in] fire_phi            Level-set field
 * @param[in] fp                  Fuel model parameters
 * @param[in] tau_res_s           Residence time [s]
 * @param[in] dt_fire_s           Atmospheric timestep duration [s]
 */
inline void fill_fire_heat_flux(
    amrex::MultiFab& fire_heat_flux,
    amrex::MultiFab& fire_fuel_load,
    const amrex::MultiFab& fire_phi,
    const FuelModelParams& fp,
    amrex::Real tau_res_s,
    amrex::Real dt_fire_s)
{
    using namespace amrex;

    // Convert heat content from BTU/lb to J/kg: 1 BTU/lb = 2326 J/kg
    Real h_fuel_J_per_kg = fp.heat_content * 2326.0;

    // GPU kernel over all cells
    for (MFIter mfi(fire_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto const phi_arr  = fire_phi.const_array(mfi);
        auto Q_arr          = fire_heat_flux.array(mfi);
        auto fuel_arr       = fire_fuel_load.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0], j = iv[1], k = iv[2];

            Real phi_val = phi_arr(i, j, k);
            Real w_cur   = fuel_arr(i, j, k);

            // Instantaneous heat flux: Q = w * h / tau_res
            Real Q = compute_heat_flux_cell(phi_val, w_cur, h_fuel_J_per_kg, tau_res_s);
            Q_arr(i, j, k) = Q;

            // Exact exponential depletion: w_new = w * exp(-dt / tau_res)
            Real w_new = deplete_fuel_load(phi_val, w_cur, tau_res_s, dt_fire_s);
            fuel_arr(i, j, k) = w_new;
        });
    }
}

#endif