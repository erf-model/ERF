#ifndef ERF_GTEST_FREEZING_COMMON_H_
#define ERF_GTEST_FREEZING_COMMON_H_

#include <cmath>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include "ERF_Constants.H"
#include "ERF_SuperDropletPCDefinitions.H"

// Unit tests for ERF's immersion-freezing INAS parameterization
// (Source/Particles/ERF_SuperDropletPCDefinitions.H, INAS_Niemand2012, after
// Niemand et al. 2012). The singular-hypothesis model sets a per-particle
// freezing temperature from its INP surface area A and the active-site density
// n_s(T) = a0 exp(-a1 T + a2): P(T) = 1 - exp(-A n_s(T)), and sample_Tfz inverts
// it. This checks the n_s and P monotonicities, the inverse-CDF round-trip, and
// the homogeneous-freezing and clamping edge cases.
//
// Each check is a non-negative "violation" returned by freezing_metric(i): 0 when
// the property holds, > 0 when it fails. Host (Scalar) and device (Kernel) tests
// call the same AMREX_GPU_HOST_DEVICE function.
namespace freezing_test {

using namespace SDPCDefn;

constexpr int kNumMetrics = 9;
constexpr amrex::Real kTol = amrex::Real(1.0e-6);

// Homogeneous-freezing temperature in Kelvin returned for IN-inactive particles.
constexpr amrex::Real kThomK = amrex::Real(-38.0) + amrex::Real(273.15);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nonneg (const amrex::Real x) noexcept { return amrex::max(amrex::Real(0.0), x); }

// The i-th property violation (>= 0; 0 when the property holds).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real freezing_metric (const int i) noexcept
{
    const INAS_Niemand2012 inas;
    const amrex::Real A = amrex::Real(1.0e-9);   // typical INP surface area [m^2]

    switch (i) {
    case 0: // active-site density rises as temperature falls
        return nonneg(inas.n_s(amrex::Real(-15.0)) - inas.n_s(amrex::Real(-30.0)));
    case 1: { // freezing probability stays in [0,1]
        const amrex::Real P = inas.P_freeze(amrex::Real(-25.0), A);
        return nonneg(-P) + nonneg(P - amrex::Real(1.0));
    }
    case 2: // more INP area raises the freezing probability
        return nonneg( inas.P_freeze(amrex::Real(-25.0), amrex::Real(1.0e-11))
                     - inas.P_freeze(amrex::Real(-25.0), amrex::Real(1.0e-9)) );
    case 3: // colder air raises the freezing probability
        return nonneg( inas.P_freeze(amrex::Real(-15.0), A)
                     - inas.P_freeze(amrex::Real(-30.0), A) );
    case 4: // no INP area -> no heterogeneous freezing probability
        return inas.P_freeze(amrex::Real(-25.0), amrex::Real(0.0));
    case 5: { // inverse-CDF round-trip: P(sample_Tfz(u)) == u in the active range
        const amrex::Real u = amrex::Real(0.5);
        const amrex::Real Tfz_K = inas.sample_Tfz(u, A);
        const amrex::Real P_rt  = inas.P_freeze(Tfz_K - amrex::Real(273.15), A);
        return nonneg(std::abs(P_rt - u) - amrex::Real(1.0e-4));
    }
    case 6: // no INP area -> homogeneous freezing temperature
        return std::abs(inas.sample_Tfz(amrex::Real(0.5), amrex::Real(-1.0)) - kThomK);
    case 7: // IN-inactive particle (u above P_max for tiny area) -> homogeneous
        return std::abs(inas.sample_Tfz(amrex::Real(0.5), amrex::Real(1.0e-20)) - kThomK);
    default: { // sampled freezing temperature is clamped to the warm INAS limit
        const amrex::Real Tmax_K = amrex::Real(-12.0) + amrex::Real(273.15);
        return std::abs(inas.sample_Tfz(amrex::Real(1.0e-6), A) - Tmax_K);
    }
    }
}

inline std::string metric_label (const int i)
{
    switch (i) {
    case 0:  return "INAS: n_s rises as T falls";
    case 1:  return "INAS: P_freeze in [0,1]";
    case 2:  return "INAS: P_freeze rises with INP area";
    case 3:  return "INAS: P_freeze rises as T falls";
    case 4:  return "INAS: P_freeze = 0 with no INP area";
    case 5:  return "INAS: inverse-CDF round-trip P(sample_Tfz)=u";
    case 6:  return "INAS: A<=0 -> homogeneous freezing T";
    case 7:  return "INAS: IN-inactive -> homogeneous freezing T";
    default: return "INAS: sampled T clamped to warm limit";
    }
}

} // namespace freezing_test

#endif
