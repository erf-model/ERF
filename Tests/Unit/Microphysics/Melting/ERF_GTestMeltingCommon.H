#ifndef ERF_GTEST_MELTING_COMMON_H_
#define ERF_GTEST_MELTING_COMMON_H_

#include <cmath>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_SuperDropletPCMassChange.H"

// Unit tests for ERF's ice melt rate
// (Source/Particles/ERF_SuperDropletPCMassChange.H, SDMassChangeUtils_SV::dMdt
// ::meltRate). The rate follows Seifert-Beheng (2006, Eq. 72): heat conducted
// from warm air plus the latent-heat flux of vapour to the 0 C wet surface,
// divided by the latent heat of fusion, with the same capacitance and
// ventilation as deposition. dm_ice/dt is negative when melting. The terminal-
// velocity test covers the melt fall-speed blend; this covers the melt mass rate.
//
// Each check is a non-negative "violation" returned by melting_metric(i): 0 when
// the property holds, > 0 when it fails. Host (Scalar) and device (Kernel) tests
// call the same AMREX_GPU_HOST_DEVICE function.
namespace melting_test {

using SDMassChangeUtils_SV::dMdt;

constexpr amrex::Real kRhoWater = amrex::Real(1000.0);
constexpr amrex::Real kP        = amrex::Real(90000.0);   // ambient pressure [Pa]
constexpr amrex::Real kD        = diffelq;                // vapour diffusivity [m^2 s^-1]

constexpr int kNumMetrics = 5;
constexpr amrex::Real kTol = amrex::Real(1.0e-6);

// dMdt functor built with the water material (per meltRate's contract: L is the
// latent heat of vaporization, K the air thermal conductivity, Rv the vapour gas
// constant; rho is unused by meltRate).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
dMdt<amrex::Real> melt_functor () noexcept
{
    return dMdt<amrex::Real>{ lcond, therco, R_v, kRhoWater };
}

// Saturation pressure over water [Pa].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real e_sat_water (const amrex::Real T) noexcept { return erf_esatw(T) * amrex::Real(100.0); }

// Ice melt rate dm_ice/dt of a spheroid (a, c) falling at vterm in water-saturated
// air at temperature T. Negative when melting.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real melt_rate ( const amrex::Real a, const amrex::Real c,
                        const amrex::Real vterm, const amrex::Real T ) noexcept
{
    const dMdt<amrex::Real> dm = melt_functor();
    const amrex::Real e_inf      = e_sat_water(T);          // water-saturated ambient
    const amrex::Real e_sat_wT0  = e_sat_water(tmelt);      // over water at 0 C
    const amrex::Real rhom       = kP / (R_d * T);
    return dm.meltRate(a, c, vterm, T, e_inf, e_sat_wT0, rhom, kD, lfus);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nonneg (const amrex::Real x) noexcept { return amrex::max(amrex::Real(0.0), x); }

// The i-th property violation (>= 0; 0 when the property holds).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real melting_metric (const int i) noexcept
{
    const amrex::Real a1mm  = amrex::Real(1.0e-3);   // 1 mm graupel sphere
    const amrex::Real v_ref = amrex::Real(3.0);      // representative fall speed [m s^-1]
    const amrex::Real T5    = tmelt + amrex::Real(5.0);
    const amrex::Real T10   = tmelt + amrex::Real(10.0);

    switch (i) {
    case 0: // above 0 C the particle melts: dm_ice/dt < 0
        return nonneg(melt_rate(a1mm, a1mm, v_ref, T5));
    case 1: // melting accelerates with warmth (more negative at 10 C than 5 C)
        return nonneg(melt_rate(a1mm, a1mm, v_ref, T10) - melt_rate(a1mm, a1mm, v_ref, T5));
    case 2: // ventilation by faster fall enhances melt (more negative)
        return nonneg(melt_rate(a1mm, a1mm, amrex::Real(6.0), T5) - melt_rate(a1mm, a1mm, amrex::Real(1.0), T5));
    case 3: // a larger particle melts faster in absolute mass terms (capacitance)
        return nonneg(melt_rate(amrex::Real(3.0e-3), amrex::Real(3.0e-3), v_ref, T5)
                    - melt_rate(a1mm, a1mm, v_ref, T5));
    default: { // no melting at the freezing point in saturated air (rate ~ 0)
        const amrex::Real m0   = std::abs(melt_rate(a1mm, a1mm, v_ref, tmelt));
        const amrex::Real mref = std::abs(melt_rate(a1mm, a1mm, v_ref, T5));
        return nonneg(m0 - amrex::Real(0.01)*mref);
    }
    }
}

inline std::string metric_label (const int i)
{
    switch (i) {
    case 0:  return "melt: dm_ice/dt < 0 above 0 C";
    case 1:  return "melt: accelerates with temperature";
    case 2:  return "melt: ventilation enhances rate";
    case 3:  return "melt: larger particle melts faster";
    default: return "melt: ~0 at the freezing point";
    }
}

} // namespace melting_test

#endif
