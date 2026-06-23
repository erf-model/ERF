#ifndef ERF_GTEST_TERMINAL_VELOCITY_COMMON_H_
#define ERF_GTEST_TERMINAL_VELOCITY_COMMON_H_

#include <cmath>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include "ERF_TerminalVelocity.H"

// Unit tests for ERF's hydrometeor terminal-velocity routines
// (Source/MaterialProperties/ERF_TerminalVelocity.H): the Beard (1976)
// cloud-and-rain fall speed with the Gunn-Kinzer 4 mm clamp, the Boehm (1989)
// ice fall speed in the Shima (2020) porous-spheroid form, and the Frick et al.
// (2013) melting-particle blend.
//
// Every check is expressed as a non-negative "violation" returned by
// terminal_velocity_metric(i): 0 when the property holds, > 0 (and larger than
// the tolerance) when it fails. The host (Scalar) and device (Kernel) tests both
// call this same function, so they verify identical physics on both backends.
namespace terminal_velocity_test {

// Material densities (kg m^-3): liquid water and solid ice, as used to build the
// terminal-velocity functor in the SDM advection step.
constexpr amrex::Real kRhoWater = amrex::Real(1000.0);
constexpr amrex::Real kRhoIce   = amrex::Real(916.8);

// Reference ambient states. Rain: sea level, ~20 C. Ice: the rho_a = 1.1 kg m^-3,
// T = 0 C condition under which Boehm (1989, Fig. 3) plots his fall-speed curves,
// so the ice-magnitude bands below are a cross-check against that published figure.
constexpr amrex::Real kTrain    = amrex::Real(293.15);   // [K]
constexpr amrex::Real kPrain    = amrex::Real(101325.0); // [Pa]
constexpr amrex::Real kRhoAirR  = amrex::Real(1.2);      // [kg m^-3]
constexpr amrex::Real kTice     = amrex::Real(273.15);   // [K]
constexpr amrex::Real kRhoAirI  = amrex::Real(1.1);      // [kg m^-3]

constexpr int kNumMetrics = 19;

// Uniform absolute tolerance: each metric is built so that a satisfied property
// returns mathematically zero (equalities, monotonicity gaps, in-band distances),
// so any small slack absorbs only floating-point round-off.
constexpr amrex::Real kTol = amrex::Real(1.0e-6);

// Semi-axes (a equatorial, c polar) of a spheroid of maximum dimension D and
// aspect ratio phi = c/a. For a column (phi > 1) the maximum dimension is the
// polar extent 2c; otherwise it is the equatorial extent 2a.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void spheroid_from_D (const amrex::Real D, const amrex::Real phi,
                      amrex::Real& a, amrex::Real& c) noexcept
{
    if (phi <= amrex::Real(1.0)) { a = amrex::Real(0.5)*D; c = phi*a; }
    else                         { c = amrex::Real(0.5)*D; a = c/phi; }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real spheroid_mass (const amrex::Real a, const amrex::Real c,
                           const amrex::Real rho_app) noexcept
{
    return rho_app * (amrex::Real(4.0)/amrex::Real(3.0)) * PI * a*a*c;
}

// Semi-axes of a spheroid of given mass, aspect ratio, and apparent density:
// m = rho_app (4/3) pi a^2 c with c = phi a  =>  a^3 = 3 m / (4 pi rho_app phi).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void spheroid_from_mass (const amrex::Real m, const amrex::Real phi,
                         const amrex::Real rho_app,
                         amrex::Real& a, amrex::Real& c) noexcept
{
    const amrex::Real a3 = amrex::Real(3.0)*m
                         / (amrex::Real(4.0)*PI*rho_app*phi);
    a = std::cbrt(a3);
    c = phi*a;
}

// Ice (Boehm) fall speed of a spheroid specified by maximum dimension.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ice_vterm_D (const TerminalVelocity<amrex::Real>& tv,
                         const amrex::Real D, const amrex::Real phi,
                         const amrex::Real rho_app) noexcept
{
    amrex::Real a, c;
    spheroid_from_D(D, phi, a, c);
    const amrex::Real m = spheroid_mass(a, c, rho_app);
    return tv.IceBohm(m, a, c, rho_app, kRhoAirI, kTice);
}

// Ice (Boehm) fall speed of a spheroid specified by mass.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ice_vterm_m (const TerminalVelocity<amrex::Real>& tv,
                         const amrex::Real m, const amrex::Real phi,
                         const amrex::Real rho_app) noexcept
{
    amrex::Real a, c;
    spheroid_from_mass(m, phi, rho_app, a, c);
    return tv.IceBohm(m, a, c, rho_app, kRhoAirI, kTice);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real band_violation (const amrex::Real v,
                            const amrex::Real lo, const amrex::Real hi) noexcept
{
    return amrex::max(amrex::Real(0.0), lo - v) + amrex::max(amrex::Real(0.0), v - hi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nonneg (const amrex::Real x) noexcept
{
    return amrex::max(amrex::Real(0.0), x);
}

// Liquid-drop (Beard 1976, with the Gunn-Kinzer clamp) fall speed at the
// reference rain state, as a function of radius (m).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rain_vterm (const TerminalVelocity<amrex::Real>& tv, const amrex::Real R) noexcept
{
    return tv.CloudRainShima(R, kRhoAirR, kPrain, kTrain);
}

// The i-th property violation (>= 0; 0 when the property holds). Direct calls
// only (no lambdas), so the function is portable to the device path.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real terminal_velocity_metric (const int i) noexcept
{
    const TerminalVelocity<amrex::Real> tv{kRhoWater, kRhoIce};

    const amrex::Real v_4mm = rain_vterm(tv, amrex::Real(4.0e-3));

    // Solid-ice sphere and a low-density graupel sphere.
    const amrex::Real v_hail = ice_vterm_D(tv, amrex::Real(1.0e-2), amrex::Real(1.0), kRhoIce);
    const amrex::Real v_grp  = ice_vterm_D(tv, amrex::Real(2.0e-3), amrex::Real(1.0), amrex::Real(400.0));

    switch (i) {
    // --- Beard (1976) liquid drops + Gunn-Kinzer clamp ---
    case 0:  // 4 mm radius clamp: fall speed is held constant beyond 4 mm
        return std::abs(rain_vterm(tv, amrex::Real(5.0e-3)) - v_4mm) / v_4mm;
    case 1:  // ... all the way out to 8 mm
        return std::abs(rain_vterm(tv, amrex::Real(8.0e-3)) - v_4mm) / v_4mm;
    case 2:  // plateau magnitude ~9 m/s at sea level (Gunn-Kinzer)
        return band_violation(v_4mm, amrex::Real(6.0), amrex::Real(12.0));
    case 3:  // monotone increasing with size below the clamp
        return nonneg(rain_vterm(tv, amrex::Real(5.0e-4)) - rain_vterm(tv, amrex::Real(2.0e-3)));
    case 4:  // small cloud droplet (10 um) falls slowly (Stokes regime)
        return nonneg(rain_vterm(tv, amrex::Real(1.0e-5)) - amrex::Real(0.2));
    case 5:  // positive fall speed
        return nonneg(-rain_vterm(tv, amrex::Real(1.0e-3)));

    // --- Boehm (1989) ice, cross-checked against his Fig. 3 (rho_a = 1.1) ---
    case 6:  // monotone increasing with size for a solid-ice sphere
        return nonneg(ice_vterm_D(tv, amrex::Real(1.0e-3), amrex::Real(1.0), kRhoIce)
                    - ice_vterm_D(tv, amrex::Real(5.0e-3), amrex::Real(1.0), kRhoIce));
    case 7:  // hail (10 mm solid sphere) ~10 m/s in Boehm Fig. 3d
        return band_violation(v_hail, amrex::Real(6.0), amrex::Real(15.0));
    case 8:  // graupel (2 mm, 400 kg m^-3 sphere): Boehm 1992 gives ~4 m/s
        return band_violation(v_grp, amrex::Real(3.0), amrex::Real(6.0));
    case 9:  // at fixed mass, a strongly aspherical column falls slower than a sphere
        return nonneg(ice_vterm_m(tv, amrex::Real(1.0e-7), amrex::Real(20.0), kRhoIce)
                    - ice_vterm_m(tv, amrex::Real(1.0e-7), amrex::Real(1.0),  kRhoIce));

    // --- Frick et al. (2013) melting-particle blend ---
    case 10: // Psi(0) = 0 : dry ice keeps its Boehm speed
        return std::abs(tv.MeltBlendWeight(amrex::Real(0.0)));
    case 11: // Psi(1) = 1 : fully melted reaches the raindrop speed
        return std::abs(tv.MeltBlendWeight(amrex::Real(1.0)) - amrex::Real(1.0));
    case 12: // Psi monotone increasing
        return nonneg(tv.MeltBlendWeight(amrex::Real(0.3)) - tv.MeltBlendWeight(amrex::Real(0.7)));
    case 13: // Psi matches Frick Eq. 15 coefficients 0.246 l + 0.754 l^7
    {
        const amrex::Real l = amrex::Real(0.6);
        const amrex::Real ref = amrex::Real(0.246)*l + amrex::Real(0.754)*std::pow(l, amrex::Real(7.0));
        return std::abs(tv.MeltBlendWeight(l) - ref);
    }
    case 14: // Mitra (1990) shape: little speed-up until ~70% melted
        return nonneg(tv.MeltBlendWeight(amrex::Real(0.7)) - amrex::Real(0.30));
    case 15: // ... then a steep rise: the last decile increment is the largest
        return nonneg((tv.MeltBlendWeight(amrex::Real(0.9)) - tv.MeltBlendWeight(amrex::Real(0.8)))
                    - (tv.MeltBlendWeight(amrex::Real(1.0)) - tv.MeltBlendWeight(amrex::Real(0.9))));
    case 16: // blended speed stays between the ice and rain endpoints
    {
        const amrex::Real vi = amrex::Real(2.0), vr = amrex::Real(9.0);
        const amrex::Real vb = tv.MixedPhaseVterm(vi, vr, amrex::Real(0.5));
        return band_violation(vb, vi, vr);
    }
    // --- Stokes regime: the small-ice limit that exercises viscosity directly ---
    case 17: // 40 um ice sphere (Re ~ 0.1): the Boehm speed must track Stokes' law
             // v = (2/9)(rho_p - rho_a) g a^2 / mu. A viscosity unit slip (CGS poise
             // vs SI Pa.s) leaves high-Re ice untouched but makes small ice ~10x slow.
    {
        const amrex::Real a = amrex::Real(4.0e-5), rho_app = amrex::Real(100.0);
        const amrex::Real v = ice_vterm_D(tv, amrex::Real(2.0)*a, amrex::Real(1.0), rho_app);
        const amrex::Real mu = amrex::Real(1.718e-5); // dynamic viscosity at 0 C [Pa s]
        const amrex::Real v_stokes = (amrex::Real(2.0)/amrex::Real(9.0))
                                   * (rho_app - kRhoAirI) * amrex::Real(CONST_GRAV) * a*a / mu;
        return band_violation(v / v_stokes, amrex::Real(0.7), amrex::Real(1.5));
    }

    default: // blend hits the endpoints exactly at ell = 0 and ell = 1
    {
        const amrex::Real vi = amrex::Real(2.0), vr = amrex::Real(9.0);
        return std::abs(tv.MixedPhaseVterm(vi, vr, amrex::Real(0.0)) - vi)
             + std::abs(tv.MixedPhaseVterm(vi, vr, amrex::Real(1.0)) - vr);
    }
    }
}

inline std::string metric_label (const int i)
{
    switch (i) {
    case 0:  return "rain: 4 mm clamp holds at 5 mm";
    case 1:  return "rain: 4 mm clamp holds at 8 mm";
    case 2:  return "rain: plateau magnitude in [6,12] m/s";
    case 3:  return "rain: monotone increasing below clamp";
    case 4:  return "rain: 10 um droplet falls slowly";
    case 5:  return "rain: positive fall speed";
    case 6:  return "ice: solid sphere monotone in size";
    case 7:  return "ice: 10 mm hail in [6,15] m/s (Boehm Fig.3)";
    case 8:  return "ice: 2 mm graupel in [3,6] m/s (Boehm 1992)";
    case 9:  return "ice: column slower than sphere at fixed mass";
    case 10: return "Frick: Psi(0) = 0";
    case 11: return "Frick: Psi(1) = 1";
    case 12: return "Frick: Psi monotone increasing";
    case 13: return "Frick: Psi matches Eq.15 coefficients";
    case 14: return "Mitra: < ~30% speed-up by 70% melt";
    case 15: return "Mitra: steepest rise in final decile";
    case 16: return "blend stays between ice and rain speeds";
    case 17: return "ice: Stokes-regime small sphere matches Stokes' law";
    default: return "blend hits endpoints at ell = 0,1";
    }
}

} // namespace terminal_velocity_test

#endif
