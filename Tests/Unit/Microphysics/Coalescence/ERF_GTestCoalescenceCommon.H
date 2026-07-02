#ifndef ERF_GTEST_COALESCENCE_COMMON_H_
#define ERF_GTEST_COALESCENCE_COMMON_H_

#include <cmath>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Enum.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include "ERF_Constants.H"
#include "ERF_SuperDropletPCCoalescence.H"

// Unit tests for ERF's collision-coalescence kernels and riming efficiencies
// (Source/Particles/ERF_SuperDropletPCCoalescence.H, struct CollisionKernel):
// the Golovin additive kernel, the gravitational sedimentation kernel, Long's
// (1974) and Hall's (1980) kernels, the Seinfeld-Pandis Brownian coagulation
// kernel, and the Beard-Grover (1974) and Erfani-Mitchell (2017) riming
// efficiencies. The geometric kernels also set the collision rate for riming
// (ice + water) and aggregation (ice + ice).
//
// Each check is a non-negative "violation" returned by coalescence_metric(i): 0
// when the property holds, > 0 when it fails. Host (Scalar) and device (Kernel)
// tests call the same AMREX_GPU_HOST_DEVICE function.
namespace coalescence_test {

using Kernel = CollisionKernel<amrex::Real, AMREX_SPACEDIM>;

constexpr int kNumMetrics = 20;
constexpr amrex::Real kTol = amrex::Real(1.0e-6);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nonneg (const amrex::Real x) noexcept { return amrex::max(amrex::Real(0.0), x); }

// Distance of v outside [lo, hi] (0 when inside).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real band_violation (const amrex::Real v, const amrex::Real lo, const amrex::Real hi) noexcept
{
    return amrex::max(amrex::Real(0.0), lo - v) + amrex::max(amrex::Real(0.0), v - hi);
}

// The i-th property violation (>= 0; 0 when the property holds). Direct calls
// only (no lambdas), so the function is portable to the device path.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real coalescence_metric (const int i) noexcept
{
    const Kernel ck{};

    // Two droplets and a relative fall velocity along x.
    const amrex::Real r1 = amrex::Real(5.0e-5);   // 50 um
    const amrex::Real r2 = amrex::Real(2.0e-5);   // 20 um
    amrex::Real va[AMREX_SPACEDIM] = {};
    amrex::Real vb[AMREX_SPACEDIM] = {};
    amrex::Real vc[AMREX_SPACEDIM] = {};
    va[0] = amrex::Real(1.0);   // particle 1 falls at 1 m/s
    vb[0] = amrex::Real(0.0);   // particle 2 at rest  -> |dv| = 1
    vc[0] = amrex::Real(2.0);   // a larger |dv| = 2 (vs vb)
    const amrex::Real xz[AMREX_SPACEDIM] = {};  // unused positions (absolute_velocity)
    const amrex::Real sumr = r1 + r2;
    const amrex::Real geom = PI*sumr*sumr*amrex::Real(1.0);  // E=1 rate at |dv|=1

    switch (i) {
    // --- Golovin additive kernel (no velocity) ---
    case 0: { // exact additive form b (X1 + X2), b = 1.5e3, X = (4/3)pi r^3
        const amrex::Real ref = amrex::Real(1.5e3)*four_thirds_pi*(r1*r1*r1 + r2*r2*r2);
        return std::abs(ck.golovin(r1, r2) - ref)/ref;
    }
    case 1: // symmetric in the two radii
        return std::abs(ck.golovin(r1, r2) - ck.golovin(r2, r1));
    case 2: // grows with droplet size
        return nonneg(ck.golovin(r1, r2) - ck.golovin(amrex::Real(8.0e-5), r2));

    // --- Gravitational sedimentation kernel ---
    case 3: // vanishes when the fall velocities are equal (no differential settling)
        return std::abs(ck.sedimentation(r1, r2, va, va, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));
    case 4: // symmetric under swapping both radius and velocity
        return std::abs(ck.sedimentation(r1, r2, va, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity) - ck.sedimentation(r2, r1, vb, va, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));
    case 5: // non-negative
        return nonneg(-ck.sedimentation(r1, r2, va, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));
    case 6: // increases with the differential fall speed
        return nonneg(ck.sedimentation(r1, r2, va, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity) - ck.sedimentation(r1, r2, vc, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));

    // --- Long's (1974) kernel ---
    case 7: // vanishes when |dv| = 0
        return std::abs(ck.Longs(r1, r2, va, va, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));
    case 8: // non-negative
        return nonneg(-ck.Longs(r1, r2, va, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));
    case 9: // symmetric under swapping both radius and velocity
        return std::abs(ck.Longs(r1, r2, va, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity) - ck.Longs(r2, r1, vb, va, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));

    // --- Hall's (1980) kernel (efficiency table + bilinear interpolation) ---
    case 10: // vanishes when |dv| = 0
        return std::abs(ck.Halls(r1, r2, va, va, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));
    case 11: // symmetric under swapping both radius and velocity
        return std::abs(ck.Halls(r1, r2, va, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity) - ck.Halls(r2, r1, vb, va, xz, xz, SDKernelRelativeVelocityType::absolute_velocity));
    case 12: // collision efficiency E = K / (pi (r1+r2)^2 |dv|) stays within [0, 1.05]
        return band_violation(ck.Halls(r1, r2, va, vb, xz, xz, SDKernelRelativeVelocityType::absolute_velocity)/geom, amrex::Real(0.0), amrex::Real(1.05));

    // --- Riming collision efficiencies ---
    case 13: // Beard-Grover (1974) efficiency in [0,1]
        return band_violation(ck.BeardGrover1974(amrex::Real(0.5), amrex::Real(100.0), amrex::Real(1.0)),
                              amrex::Real(0.0), amrex::Real(1.0));
    case 14: // Erfani-Mitchell plate: no collection below Re = 1
        return std::abs(ck.ErfaniMitchell2017_Plate(amrex::Real(0.5), amrex::Real(1.0)));
    case 15: // Erfani-Mitchell plate: bounded efficiency for a valid Reynolds number
        return band_violation(ck.ErfaniMitchell2017_Plate(amrex::Real(50.0), amrex::Real(1.0)),
                              amrex::Real(0.0), amrex::Real(1.5));
    case 16: // Erfani-Mitchell column: no collection below Re = 0.2
        return std::abs(ck.ErfaniMitchell2017_Column(amrex::Real(0.1), amrex::Real(1.0)));
    case 17: // Erfani-Mitchell column: bounded efficiency for a valid Reynolds number
        return band_violation(ck.ErfaniMitchell2017_Column(amrex::Real(10.0), amrex::Real(1.0)),
                              amrex::Real(0.0), amrex::Real(1.5));

    // --- Brownian coagulation (Seinfeld & Pandis) ---
    case 18: { // symmetric under swapping both radius and mass
        const amrex::Real m1 = amrex::Real(1.0e-12), m2 = amrex::Real(5.0e-14);
        const amrex::Real P = amrex::Real(9.0e4), T = amrex::Real(283.15);
        return std::abs(ck.Brownian_SeinfeldPandis(amrex::Real(1.0e-7), amrex::Real(2.0e-8), m1, m2, P, T)
                      - ck.Brownian_SeinfeldPandis(amrex::Real(2.0e-8), amrex::Real(1.0e-7), m2, m1, P, T));
    }
    default: { // non-negative
        const amrex::Real m1 = amrex::Real(1.0e-12), m2 = amrex::Real(5.0e-14);
        return nonneg(-ck.Brownian_SeinfeldPandis(amrex::Real(1.0e-7), amrex::Real(2.0e-8), m1, m2,
                                                  amrex::Real(9.0e4), amrex::Real(283.15)));
    }
    }
}

inline std::string metric_label (const int i)
{
    switch (i) {
    case 0:  return "Golovin: exact additive form";
    case 1:  return "Golovin: symmetric in radii";
    case 2:  return "Golovin: grows with size";
    case 3:  return "sedimentation: zero at equal fall speed";
    case 4:  return "sedimentation: symmetric";
    case 5:  return "sedimentation: non-negative";
    case 6:  return "sedimentation: grows with |dv|";
    case 7:  return "Long: zero at |dv| = 0";
    case 8:  return "Long: non-negative";
    case 9:  return "Long: symmetric";
    case 10: return "Hall: zero at |dv| = 0";
    case 11: return "Hall: symmetric";
    case 12: return "Hall: efficiency in [0,1.05]";
    case 13: return "Beard-Grover: efficiency in [0,1]";
    case 14: return "Erfani-Mitchell plate: cutoff below Re=1";
    case 15: return "Erfani-Mitchell plate: bounded efficiency";
    case 16: return "Erfani-Mitchell column: cutoff below Re=0.2";
    case 17: return "Erfani-Mitchell column: bounded efficiency";
    case 18: return "Brownian: symmetric";
    default: return "Brownian: non-negative";
    }
}

} // namespace coalescence_test

#endif
