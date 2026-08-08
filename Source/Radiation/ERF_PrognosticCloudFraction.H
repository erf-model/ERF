#ifndef ERF_PROGNOSTIC_CLOUD_FRACTION_H_
#define ERF_PROGNOSTIC_CLOUD_FRACTION_H_

#include <AMReX_GpuControl.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

/**
 * @file ERF_PrognosticCloudFraction.H
 * @brief Phase 14: Prognostic cloud fraction diagnosis from relative humidity and cloud water.
 *
 * Implements RH/qc-based cloud fraction diagnosis for use in TwoStream radiation.
 * The diagnosed cloud fraction cf(k) can optionally be temporally smoothed.
 *
 * Cloud fraction diagnosis formula:
 *   cf_rh(k) = max(0, min(1, (rh(k) - rh_min) / (rh_max - rh_min)))  [linear ramp]
 *   cf_qc(k) = min(1, qc_scale * qc(k))  [qc contribution]
 *   cf_diag(k) = min(1, cf_rh(k) + cf_qc(k))  [combined with saturation]
 *
 * When smoothing enabled (Phase 14):
 *   cf_smooth(k, t+dt) = alpha * cf_diag(k, t+dt) + (1 - alpha) * cf_smooth(k, t)
 *
 * References:
 * -----------
 * - Phase 14 Radiation Development notes (Source/Radiation/RAD_DEVELOPMENT.md)
 */

/**
 * @brief Compute relative humidity from water vapor mixing ratio.
 *
 * Approximate RH from qv and temperature using saturation mixing ratio.
 * Uses a simple saturation vapor pressure formula (Magnus or similar).
 *
 * @param[in] qv Water vapor mixing ratio [kg/kg]
 * @param[in] T Temperature [K]
 * @param[in] P Pressure [Pa]
 * @return Relative humidity [0, 1]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_relative_humidity(amrex::Real qv, amrex::Real T, amrex::Real P)
{
    // Guard against invalid inputs
    if (T <= 0.0 || P <= 0.0 || qv < 0.0 || !std::isfinite(qv) ||
        !std::isfinite(T) || !std::isfinite(P)) {
        return 0.0;
    }

    // Simple saturation mixing ratio approximation (Magnus formula variant)
    // qsat = 0.622 * e_sat / (P - e_sat)
    // where e_sat = e0 * exp(a * (T - T0) / (T - b))
    // For simplicity, use approximate saturation vapor pressure
    const amrex::Real e0 = 611.2;  // Pa (saturation vapor pressure at 0°C)
    const amrex::Real T0 = 273.15; // K
    const amrex::Real a = 17.62;
    const amrex::Real b = 243.12;  // K

    // Compute saturation vapor pressure (Magnus formula)
    amrex::Real arg = a * (T - T0) / (T - b);
    if (!std::isfinite(arg) || arg > 100.0) {
        // Protect against overflow
        arg = 100.0;
    }
    amrex::Real e_sat = e0 * std::exp(arg);

    // Saturation mixing ratio: qsat = epsilon * e_sat / (P - e_sat)
    const amrex::Real epsilon = 0.622;  // ratio of molecular weights
    amrex::Real qsat = epsilon * e_sat / (P - e_sat);
    if (qsat < 0.0 || !std::isfinite(qsat)) {
        qsat = 1.0e-6;  // fallback
    }

    // RH = qv / qsat, clamped to [0, 1]
    amrex::Real rh = qv / qsat;
    if (rh < 0.0) rh = 0.0;
    if (rh > 1.0) rh = 1.0;

    return rh;
}

/**
 * @brief Diagnose cloud fraction from relative humidity and cloud water.
 *
 * Computes diagnosed cloud fraction from RH and qc:
 *   - RH contribution: linear ramp from 0 at rh_min to 1 at rh_max
 *   - qc contribution: scaled by qc_scale coefficient
 *   - Combined: saturated blend (cf ≤ 1)
 *
 * This is a GPU-safe inline function for device-side kernels.
 *
 * @param[in] rh Relative humidity [0, 1]
 * @param[in] qc Cloud liquid water mixing ratio [kg/kg]
 * @param[in] rh_min Minimum RH threshold [0, 1]
 * @param[in] rh_max Maximum RH threshold [0, 1]
 * @param[in] qc_scale Scaling coefficient for qc contribution
 *
 * @return Diagnosed cloud fraction [0, 1]
 *
 * @note All inputs are assumed to be validated and finite.
 * @note Output is always in [0, 1], finite, and safe for further use.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_cloud_fraction_from_rh_qc(
    amrex::Real rh, amrex::Real qc,
    amrex::Real rh_min, amrex::Real rh_max,
    amrex::Real qc_scale)
{
    // Guard against invalid inputs
    if (!std::isfinite(rh) || !std::isfinite(qc) ||
        !std::isfinite(rh_min) || !std::isfinite(rh_max) ||
        !std::isfinite(qc_scale)) {
        return 0.0;  // Safe fallback on invalid input
    }

    // Clamp RH and qc to valid ranges
    if (rh < 0.0) rh = 0.0;
    if (rh > 1.0) rh = 1.0;
    if (qc < 0.0) qc = 0.0;

    amrex::Real cf = 0.0;

    // RH contribution: linear ramp from 0 at rh_min to 1 at rh_max
    if (rh_max > rh_min) {
        // Linear interpolation
        amrex::Real cf_rh = (rh - rh_min) / (rh_max - rh_min);
        if (cf_rh < 0.0) cf_rh = 0.0;
        if (cf_rh > 1.0) cf_rh = 1.0;
        cf += cf_rh;
    } else if (rh >= rh_min) {
        // rh_max == rh_min: step function
        cf += 1.0;
    }

    // qc contribution: scaled by qc_scale
    if (qc_scale > 0.0) {
        amrex::Real cf_qc = qc_scale * qc;
        if (cf_qc > 1.0) cf_qc = 1.0;
        cf += cf_qc;
    }

    // Saturate combined cloud fraction to [0, 1]
    if (cf < 0.0) cf = 0.0;
    if (cf > 1.0) cf = 1.0;

    // Final sanity check
    if (!std::isfinite(cf)) {
        cf = 0.0;
    }

    return cf;
}

/**
 * @brief Apply exponential moving average (EMA) smoothing to cloud fraction.
 *
 * Temporal smoothing via EMA:
 *   cf_smooth(t+dt) = alpha * cf_new + (1 - alpha) * cf_old
 *
 * alpha = 0.0: no smoothing (cf_old retained)
 * alpha = 1.0: instant update (cf_new applied)
 * alpha ∈ (0, 1): gradual transition
 *
 * This is a GPU-safe inline function for device-side kernels.
 *
 * @param[in] cf_new Newly diagnosed cloud fraction [0, 1]
 * @param[in] cf_old Previous (smoothed) cloud fraction [0, 1]
 * @param[in] alpha Smoothing parameter [0, 1]
 *
 * @return Smoothed cloud fraction [0, 1]
 *
 * @note All inputs are assumed finite; output is always finite and bounded.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real smooth_cloud_fraction_ema(
    amrex::Real cf_new, amrex::Real cf_old,
    amrex::Real alpha)
{
    // Guard against invalid inputs
    if (!std::isfinite(cf_new) || !std::isfinite(cf_old) || !std::isfinite(alpha)) {
        return 0.0;  // Safe fallback
    }

    // Clamp all inputs to [0, 1]
    if (cf_new < 0.0) cf_new = 0.0;
    if (cf_new > 1.0) cf_new = 1.0;
    if (cf_old < 0.0) cf_old = 0.0;
    if (cf_old > 1.0) cf_old = 1.0;
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    // EMA: cf_smooth = alpha * cf_new + (1 - alpha) * cf_old
    amrex::Real cf_smooth = alpha * cf_new + (1.0 - alpha) * cf_old;

    // Final bounds check and finiteness guard
    if (cf_smooth < 0.0) cf_smooth = 0.0;
    if (cf_smooth > 1.0) cf_smooth = 1.0;
    if (!std::isfinite(cf_smooth)) cf_smooth = 0.0;

    return cf_smooth;
}

#endif // ERF_PROGNOSTIC_CLOUD_FRACTION_H_
