// Fuel moisture time-lag differential equations and moisture of extinction.
//
// Time-lag model reference:
//   Nelson, R.M. (2000). Prediction of diurnal change in 10-h fuel stick
//     moisture content. Canadian Journal of Forest Research, 30, 1071-1087.
//   Van Wagner, C.E. (1972). Equilibrium moisture contents of some fine
//     forest fuels in eastern Canada. Canadian Forestry Service Report 7.
//
// Moisture of extinction reference:
//   Rothermel, R.C. (1972). A mathematical model for predicting fire spread
//     in wildland fuels. USDA Forest Service Research Paper INT-115.
//
// Ported from hgopalan/wildfire_levelset src/fuel_moisture_timelag_de.H
// and src/moisture_of_extinction.H.
//
// Assumptions:
//   - Temperature is in degrees Celsius.
//   - Relative humidity is in percent (0-100).
//   - Fuel moisture content is a fraction (0.01-0.40).
//   - Precipitation rate is in mm/hr.
//   - Time step is in hours.
//   - The forward Euler scheme is stable for dt_hours << tau (time-lag constant).
//     For a 30-second atmospheric step: dt_hours = 8.3e-3 hr << tau_1hr = 1 hr.

#ifndef ERF_FUEL_MOISTURE_H
#define ERF_FUEL_MOISTURE_H

#include <AMReX_REAL.H>
#include <AMReX_GpuDevice.H>
#include <cmath>

namespace FuelMoistureConst {
    // Time-lag constants [hours] — Van Wagner (1972), Nelson (2000)
    constexpr amrex::Real TAU_1HR   = 1.0;
    constexpr amrex::Real TAU_10HR  = 10.0;
    constexpr amrex::Real TAU_100HR = 100.0;

    // Moisture bounds [fraction]
    constexpr amrex::Real M_MIN = 0.01;
    constexpr amrex::Real M_MAX = 0.40;

    // RH clamp bounds [fraction of 100]
    constexpr amrex::Real RH_MIN = 0.01;
    constexpr amrex::Real RH_MAX = 0.99;

    // Precipitation wetting
    constexpr amrex::Real RAIN_WETTING_COEF = 0.01; // [hr/mm]
    constexpr amrex::Real MIN_RAIN_RATE     = 0.1;  // [mm/hr]

    // Temperature correction
    constexpr amrex::Real TEMP_REFERENCE = 20.0;  // [°C]
    constexpr amrex::Real TEMP_FACTOR_K  = 0.015; // [1/°C]
}

// Equilibrium moisture content (adsorption / wetting path).
// Nelson (2000) 4th-degree polynomial in H = RH/100.
// Valid for H in [0.01, 0.99]. Returns EMC [fraction].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_emc_adsorption(amrex::Real RH) noexcept
{
    // Convert RH% to fraction and clamp
    amrex::Real H = (RH / 100.0_rt);
    H = amrex::min(amrex::max(H, FuelMoistureConst::RH_MIN),
                   FuelMoistureConst::RH_MAX);

    // Nelson (2000) adsorption polynomial
    constexpr amrex::Real a0 = 0.03229_rt;
    constexpr amrex::Real a1 = 0.2810_rt;
    constexpr amrex::Real a2 = 0.4093_rt;
    constexpr amrex::Real a3 = -1.3560_rt;
    constexpr amrex::Real a4 = 1.6596_rt;

    const amrex::Real H2 = H * H;
    const amrex::Real H3 = H2 * H;
    const amrex::Real H4 = H3 * H;

    amrex::Real M_e = a0 + a1*H + a2*H2 + a3*H3 + a4*H4;

    // Clamp to physical range
    M_e = amrex::max(0.0_rt, amrex::min(M_e, 0.35_rt));

    return M_e;
}

// Equilibrium moisture content (desorption / drying path).
// Nelson (2000) 4th-degree polynomial in H = RH/100.
// Returns EMC [fraction].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_emc_desorption(amrex::Real RH) noexcept
{
    // Convert RH% to fraction and clamp
    amrex::Real H = (RH / 100.0_rt);
    H = amrex::min(amrex::max(H, FuelMoistureConst::RH_MIN),
                   FuelMoistureConst::RH_MAX);

    // Nelson (2000) desorption polynomial
    constexpr amrex::Real d0 = 0.05800_rt;
    constexpr amrex::Real d1 = 0.1985_rt;
    constexpr amrex::Real d2 = 0.6250_rt;
    constexpr amrex::Real d3 = -1.1830_rt;
    constexpr amrex::Real d4 = 1.0570_rt;

    const amrex::Real H2 = H * H;
    const amrex::Real H3 = H2 * H;
    const amrex::Real H4 = H3 * H;

    amrex::Real M_e = d0 + d1*H + d2*H2 + d3*H3 + d4*H4;

    // Clamp to physical range
    M_e = amrex::max(0.0_rt, amrex::min(M_e, 0.35_rt));

    return M_e;
}

// Select adsorption or desorption EMC based on current moisture state
// (hysteresis). When M_current < M_e_des the fuel is drying; when
// M_current > M_e_ads it is wetting; when between curves, the average is used.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_emc_with_hysteresis(amrex::Real RH, amrex::Real M_current) noexcept
{
    const amrex::Real M_e_ads = compute_emc_adsorption(RH);
    const amrex::Real M_e_des = compute_emc_desorption(RH);

    // Determine wetting or drying based on current moisture
    if (M_current < M_e_des) {
        // Fuel is drying
        return M_e_des;
    } else if (M_current > M_e_ads) {
        // Fuel is wetting
        return M_e_ads;
    } else {
        // Between curves: use average
        return (M_e_ads + M_e_des) * 0.5_rt;
    }
}

// Temperature correction factor for the effective time-lag constant.
// f_T = exp(-k * (T - T_ref)), clamped to [0.5, 2.0].
// Higher T -> shorter effective time lag -> faster moisture response.
// Reference: Nelson (2000).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_temp_correction(amrex::Real T_C) noexcept
{
    const amrex::Real T_ref = FuelMoistureConst::TEMP_REFERENCE;
    const amrex::Real k     = FuelMoistureConst::TEMP_FACTOR_K;

    amrex::Real f_T = std::exp(-k * (T_C - T_ref));

    // Clamp to reasonable range (0.5× to 2× nominal time lag)
    f_T = amrex::min(amrex::max(f_T, 0.5_rt), 2.0_rt);

    return f_T;
}

// Wetting rate P(t) [1/hr] from precipitation intensity.
// P = k_w * rain_rate_mm_hr  for rain_rate >= MIN_RAIN_RATE
// P = 0                       otherwise
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_precip_wetting_rate(amrex::Real rain_rate_mm_hr) noexcept
{
    if (rain_rate_mm_hr < FuelMoistureConst::MIN_RAIN_RATE) {
        return 0.0_rt;
    }

    const amrex::Real k_w = FuelMoistureConst::RAIN_WETTING_COEF;

    return k_w * rain_rate_mm_hr;
}

// Advance fuel moisture for one size class using the time-lag ODE.
//
// Governing equation (forward Euler):
//   M_new = M + dt_hours * [(M_e - M) / tau_eff + P(t)]
//
// where tau_eff = tau_hours * f_T (temperature-corrected time lag).
//
// Parameters:
//   M_current    — current fuel moisture [fraction]
//   RH           — relative humidity [%]
//   T_C          — temperature [°C]
//   rain_mm_hr   — precipitation rate [mm/hr] (0 if none)
//   dt_hours     — time step [hours]
//   tau_hours    — base time-lag constant [hours]
//
// Returns updated moisture [fraction], clamped to [M_MIN, M_MAX].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real advance_fuel_moisture_one_class(
    amrex::Real M_current,
    amrex::Real RH,
    amrex::Real T_C,
    amrex::Real rain_mm_hr,
    amrex::Real dt_hours,
    amrex::Real tau_hours) noexcept
{
    // Ensure M_current is in valid range
    M_current = amrex::min(amrex::max(M_current, FuelMoistureConst::M_MIN),
                           FuelMoistureConst::M_MAX);

    // Compute equilibrium moisture with hysteresis
    const amrex::Real M_e = compute_emc_with_hysteresis(RH, M_current);

    // Compute temperature correction factor
    const amrex::Real f_T = compute_temp_correction(T_C);

    // Effective time-lag constant [hours]
    const amrex::Real tau_eff = tau_hours * f_T;

    // Compute precipitation wetting rate
    const amrex::Real P = compute_precip_wetting_rate(rain_mm_hr);

    // Time-lag moisture change term
    const amrex::Real dM_dt_timelag = (M_e - M_current) / amrex::max(tau_eff, 0.1_rt);

    // Total moisture change rate
    const amrex::Real dM_dt = dM_dt_timelag + P;

    // Forward Euler update
    amrex::Real M_new = M_current + dt_hours * dM_dt;

    // Clamp to physical range
    M_new = amrex::min(amrex::max(M_new, FuelMoistureConst::M_MIN),
                       FuelMoistureConst::M_MAX);

    return M_new;
}

#endif  // ERF_FUEL_MOISTURE_H
