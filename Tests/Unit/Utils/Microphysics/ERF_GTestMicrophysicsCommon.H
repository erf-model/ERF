#ifndef ERF_GTEST_MICROPHYSICS_COMMON_H_
#define ERF_GTEST_MICROPHYSICS_COMMON_H_

#include <algorithm>
#include <cmath>
#include <string>

#include <AMReX_Math.H>

#include <ERF_Constants.H>
#include <ERF_HSEUtils.H>
#include <ERF_MicrophysicsUtils.H>

namespace microphysics_test {

#ifdef AMREX_USE_FLOAT
// Test-only relative value tolerance; not a production threshold.
constexpr amrex::Real kValueRelTol = amrex::Real(1.0e-4);
// Test-only derivative tolerance; not a production threshold.
constexpr amrex::Real kDerivativeRelTol = amrex::Real(2.0e-2);
// Test-only empirical-branch agreement tolerance; not a production threshold.
constexpr amrex::Real kEmpiricalAgreementRelTol = amrex::Real(1.0e-2);
// Test-only finite-difference step chosen for stable float sweeps; not a
// production threshold.
constexpr amrex::Real kFiniteDiffStep = amrex::Real(5.0e-2);
#else
// Test-only relative value tolerance; not a production threshold.
constexpr amrex::Real kValueRelTol = amrex::Real(1.0e-12);
// Test-only derivative tolerance; not a production threshold.
constexpr amrex::Real kDerivativeRelTol = amrex::Real(2.0e-5);
// Test-only empirical-branch agreement tolerance; not a production threshold.
constexpr amrex::Real kEmpiricalAgreementRelTol = amrex::Real(2.0e-3);
// Test-only finite-difference step chosen for stable double sweeps; not a
// production threshold.
constexpr amrex::Real kFiniteDiffStep = amrex::Real(1.0e-3);
#endif

// Measured relative value jumps at T = switch +/- 1e-3 K on this branch are
// about 5.12e-5 (lower) and 1.386e-2 (upper). Use about 5x the larger jump to
// characterize the independent-fit mismatch rather than machine precision.
constexpr amrex::Real kSwitchValueRelTol = amrex::Real(7.0e-2);
// Measured relative derivative jumps at T = switch +/- 1e-3 K on this branch
// are about 4.92e-5 (lower) and 2.594e-2 (upper). Use about 5x the larger jump.
constexpr amrex::Real kSwitchDerivativeRelTol = amrex::Real(1.3e-1);

// Analytic reference for Gamma(1/2) = sqrt(pi); not a production threshold.
constexpr amrex::Real kSqrtPi = amrex::Real(1.7724538509055160273);
// Regression reference for erf_esatw(273.15, true) after the empirical Pa to
// mbar conversion inside erf_esatw; not a production threshold.
constexpr amrex::Real kEmpiricalTriplePointMbar = amrex::Real(6.11245202992297);
// Same regression reference as kEmpiricalTriplePointMbar, but in the original
// empirical Pa units before erf_esatw converts them to mbar.
constexpr amrex::Real kEmpiricalTriplePointPa = amrex::Real(611.245202992297);
// MUST MATCH poly_min_dtt = -70.0 in erf_use_positive_esatw_poly:
// 273.16 - 70.0 = 203.16 K.
constexpr amrex::Real kWaterLowerSwitchK = amrex::Real(203.16);
// MUST MATCH the upper dtt < 70.0 branch in erf_use_positive_esatw_poly:
// 273.16 + 70.0 = 343.16 K.
constexpr amrex::Real kWaterUpperSwitchK = amrex::Real(343.16);

inline amrex::Real scaled_tol (const amrex::Real actual,
                               const amrex::Real expected,
                               const amrex::Real factor)
{
    return factor * std::max({amrex::Real(1.0), std::abs(actual), std::abs(expected)});
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real portable_scaled_tol (const amrex::Real actual,
                                 const amrex::Real expected,
                                 const amrex::Real factor)
{
    return factor * amrex::max(amrex::Real(1.0),
                                amrex::max(amrex::Math::abs(actual), amrex::Math::abs(expected)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real normalized_error (const amrex::Real actual,
                              const amrex::Real expected,
                              const amrex::Real factor)
{
    return amrex::Math::abs(actual - expected) / portable_scaled_tol(actual, expected, factor);
}

template <typename Function>
inline amrex::Real central_difference (Function&& func,
                                       const amrex::Real x,
                                       const amrex::Real step = kFiniteDiffStep)
{
    return (func(x + step) - func(x - step)) / (amrex::Real(2.0) * step);
}

inline amrex::Real cc_derivative_formula (const amrex::Real temperature)
{
    // MUST MATCH svp2, svp3, and svpt0 in erf_dtesatw_cc. Despite the `_cc`
    // suffix, the production helper evaluates a Magnus-style closed-form
    // exponential approximation, not a direct Clausius-Clapeyron integration.
    constexpr amrex::Real svp2 = amrex::Real(17.67);
    constexpr amrex::Real svp3 = amrex::Real(29.65);
    constexpr amrex::Real svpt0 = amrex::Real(273.15);
    const amrex::Real esatw = erf_esatw_cc(temperature);

    return esatw * svp2 * (svpt0 - svp3) /
           ((temperature - svp3) * (temperature - svp3));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool uses_normal_qsat_branch (const amrex::Real esat, const amrex::Real pressure)
{
    // MUST MATCH the normal-branch condition (p - esat) >= esat in
    // erf_dtqsat_from_esat in Source/Utils/ERF_MicrophysicsUtils.H.
    return (pressure - esat) >= esat;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real capped_branch_pressure (const amrex::Real esat)
{
    // MUST MATCH the capped branch condition (p - esat) < esat, i.e. p < 2*esat,
    // in erf_dtqsat_from_esat. Choosing 1.5*esat stays strictly inside the
    // capped region whenever esat > 0.
    return amrex::Real(1.5) * esat;
}

inline std::string tp_label (const char* label,
                             const amrex::Real temperature,
                             const amrex::Real pressure)
{
    return std::string(label) + " T=" + std::to_string(temperature) +
           " p=" + std::to_string(pressure);
}

} // namespace microphysics_test

#endif