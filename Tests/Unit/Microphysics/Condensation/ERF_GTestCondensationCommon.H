#ifndef ERF_GTEST_CONDENSATION_COMMON_H_
#define ERF_GTEST_CONDENSATION_COMMON_H_

#include <cmath>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_MaterialProperties.H"
#include "ERF_SuperDropletPCMassChange.H"

// Unit tests for ERF's liquid condensation/evaporation
// (Source/Particles/ERF_SuperDropletPCMassChange.H, namespace
// SDMassChangeUtils_LV). The growth equation is integrated in R^2 with the
// Koehler curvature and solute terms and the Bayley et al. (2025) ventilation
// factor, exactly as the production MassChange_LV step builds it
// (dRsqdt -> NewtonSolver -> TI::rk3bs). Each check is replayed through those
// production routines so it verifies the model, not a re-derivation.
//
// Every check is expressed as a non-negative "violation" returned by
// condensation_metric(i): mathematically 0 when the property holds, > 0 when it
// fails. Host (Scalar) and device (Kernel) tests call the same function, so they
// verify identical physics on both backends.
namespace condensation_test {

using SDMassChangeUtils_LV::dRsqdt;
using SDMassChangeUtils_LV::NewtonSolver;
using SDMassChangeUtils_LV::TI;
using SDMassChangeUtils_LV::ventilationFactor;

constexpr amrex::Real kRhoWater = amrex::Real(1000.0);   // liquid water density [kg m^-3]
constexpr amrex::Real kT        = amrex::Real(283.15);   // reference temperature [K] (10 C)
constexpr amrex::Real kP        = amrex::Real(90000.0);  // ambient pressure [Pa] (~1 km)
constexpr amrex::Real kNs       = amrex::Real(1.0e-16);  // solute amount [mol] (~0.1 um (NH4)2SO4 CCN)

// Ventilation-factor defaults (ERF_SuperDropletPCInitializations.cpp).
constexpr amrex::Real kVentA1   = amrex::Real(6.954e7);
constexpr amrex::Real kVentB1   = amrex::Real(1.963);
constexpr amrex::Real kVentA2   = amrex::Real(1.069e3);
constexpr amrex::Real kVentB2   = amrex::Real(0.702);
constexpr amrex::Real kVentFcap = amrex::Real(20.0);

constexpr int kNumMetrics = 14;

// Uniform absolute tolerance. Each metric is built so a satisfied property
// returns mathematically zero (sign tests, monotonicity gaps, band distances),
// leaving only floating-point round-off for this slack to absorb.
constexpr amrex::Real kTol = amrex::Real(1.0e-6);

// Water MaterialPropertiesCore, matching the values ERF assigns the water
// species (ERF_MaterialProperties.cpp): m_is_water = true selects the water
// branches of coeffCurv / coeffVPSolute / coeffMolecularDiffusion.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MaterialPropertiesCore water_material () noexcept
{
    MaterialPropertiesCore m;
    m.m_density          = kRhoWater;
    m.m_mol_weight       = amrex::Real(1.802e-2);
    m.m_lat_vap          = lcond;                  // vaporization
    m.m_lat_fus          = lfus;
    m.m_Rv               = R_v;
    m.m_Tb               = tboil;
    m.m_Nav_by_molweight = MaterialPropertiesCore::s_N_av / amrex::Real(1.802e-2);
    m.m_is_water         = true;
    return m;
}

// Saturation pressure over water [Pa] (erf_esatw returns hPa).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real e_sat_water (const amrex::Real T) noexcept
{
    return erf_esatw(T) * amrex::Real(100.0);
}

// Growth-equation functor with the water constants.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
dRsqdt<amrex::Real> make_drsqdt () noexcept
{
    return dRsqdt<amrex::Real>{ lcond, therco, R_v, kRhoWater };
}

// dR^2/dt from the production routine at radius R, ventilation off (f_v = 1).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real drsqdt_at ( const amrex::Real R, const amrex::Real S,
                        const amrex::Real a_curv, const amrex::Real b_sol,
                        const amrex::Real N_s, const amrex::Real T ) noexcept
{
    const dRsqdt<amrex::Real> ode = make_drsqdt();
    const amrex::Real D   = water_material().coeffMolecularDiffusion(T, kP);
    const amrex::Real e_s = e_sat_water(T);
    return ode.rhs_func(R*R, S, T, e_s, D, a_curv, b_sol, N_s, amrex::Real(1.0));
}

// Koehler equilibrium saturation ratio of a droplet of radius R:
// S_eq(R) = 1 + (a/T)/R - b N_s / R^3. The RHS vanishes at S = S_eq(R).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real koehler_Seq ( const amrex::Real R, const amrex::Real a_curv,
                          const amrex::Real b_sol, const amrex::Real N_s,
                          const amrex::Real T ) noexcept
{
    return amrex::Real(1.0) + (a_curv/T)/R - b_sol*N_s/(R*R*R);
}

// Koehler critical (activation) radius r_c = sqrt(3 b N_s T / a), the maximum of
// S_eq(R).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real koehler_rc ( const amrex::Real a_curv, const amrex::Real b_sol,
                         const amrex::Real N_s, const amrex::Real T ) noexcept
{
    return std::sqrt( amrex::Real(3.0)*b_sol*N_s*T / a_curv );
}

// Integrate one droplet from radius r0 for t_total seconds at fixed S, T through
// the production TI::rk3bs step, returning the final radius (r0 on failure).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real grow_droplet ( const amrex::Real S, const amrex::Real T,
                           const amrex::Real r0, const amrex::Real N_s,
                           const amrex::Real t_total, const bool ventilation ) noexcept
{
    const MaterialPropertiesCore wat = water_material();
    const amrex::Real a_curv = wat.coeffCurv(T);
    const amrex::Real b_sol  = wat.coeffVPSolute();
    const amrex::Real D      = wat.coeffMolecularDiffusion(T, kP);
    const amrex::Real e_s    = e_sat_water(T);

    const dRsqdt<amrex::Real> ode = make_drsqdt();
    const NewtonSolver<dRsqdt<amrex::Real>, amrex::Real> newton{
        ode, amrex::Real(1.0e-6), amrex::Real(1.0e-30), amrex::Real(1.0e-12), 10 };

    TI< dRsqdt<amrex::Real>, NewtonSolver<dRsqdt<amrex::Real>, amrex::Real>, amrex::Real >
        ti{ ode, newton, t_total, amrex::Real(100), S, T, e_s, D, a_curv, b_sol, N_s,
            amrex::Real(1000.0), amrex::Real(1.0e-40), amrex::Real(1.0e-3), amrex::Real(1.0e-6),
            false, false, ventilation, kVentA1, kVentB1, kVentA2, kVentB2, kVentFcap };

    amrex::Real r_sq = r0*r0;
    bool success = false;
    ti.rk3bs(r_sq, success);
    return success ? std::sqrt(r_sq) : r0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nonneg (const amrex::Real x) noexcept { return amrex::max(amrex::Real(0.0), x); }

// The i-th property violation (>= 0; 0 when the property holds). Direct calls
// only (no lambdas), so the function is portable to the device path.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real condensation_metric (const int i) noexcept
{
    const MaterialPropertiesCore wat = water_material();
    const amrex::Real a_curv = wat.coeffCurv(kT);
    const amrex::Real b_sol  = wat.coeffVPSolute();
    const amrex::Real zero   = amrex::Real(0.0);

    switch (i) {
    // --- Diffusional (Maxwellian) growth: curvature and solute off ---
    case 0: { // r^2-law: dR^2/dt is nearly radius-independent for large drops
        const amrex::Real g50  = drsqdt_at(amrex::Real(5.0e-5), amrex::Real(1.01), zero, zero, zero, kT);
        const amrex::Real g100 = drsqdt_at(amrex::Real(1.0e-4), amrex::Real(1.01), zero, zero, zero, kT);
        return nonneg(std::abs(g50 - g100)/std::abs(g100) - amrex::Real(0.02)); // within 2%
    }
    case 1: // supersaturation S > 1 grows the drop (dR^2/dt > 0)
        return nonneg(-drsqdt_at(amrex::Real(2.0e-5), amrex::Real(1.02), zero, zero, zero, kT));
    case 2: // subsaturation S < 1 evaporates the drop (dR^2/dt < 0)
        return nonneg( drsqdt_at(amrex::Real(2.0e-5), amrex::Real(0.98), zero, zero, zero, kT));

    // --- Koehler curvature (Kelvin) and solute (Raoult) terms ---
    case 3: // at S = 1 a small pure droplet evaporates from curvature alone
        return nonneg( drsqdt_at(amrex::Real(5.0e-7), amrex::Real(1.0), a_curv, zero, zero, kT));
    case 4: // at S = 1 solute alone drives growth (vapour-pressure lowering)
        return nonneg(-drsqdt_at(amrex::Real(5.0e-7), amrex::Real(1.0), zero, b_sol, kNs, kT));
    case 5: { // a haze droplet sits in equilibrium at S = S_eq (RHS ~ 0)
        const amrex::Real r_eq = amrex::Real(3.0e-7);
        const amrex::Real S_eq = koehler_Seq(r_eq, a_curv, b_sol, kNs, kT);
        const amrex::Real g_eq  = std::abs(drsqdt_at(r_eq, S_eq, a_curv, b_sol, kNs, kT));
        const amrex::Real g_off = std::abs(drsqdt_at(r_eq, S_eq + amrex::Real(0.001), a_curv, b_sol, kNs, kT));
        return nonneg(g_eq - amrex::Real(0.01)*g_off); // equilibrium RHS < 1% of off-equilibrium
    }
    case 6: { // S_eq is maximal at the critical radius r_c (the activation peak)
        const amrex::Real r_c = koehler_rc(a_curv, b_sol, kNs, kT);
        const amrex::Real Sc  = koehler_Seq(r_c, a_curv, b_sol, kNs, kT);
        return nonneg(koehler_Seq(amrex::Real(0.5)*r_c, a_curv, b_sol, kNs, kT) - Sc)
             + nonneg(koehler_Seq(amrex::Real(2.0)*r_c, a_curv, b_sol, kNs, kT) - Sc);
    }

    // --- Time integration through TI::rk3bs ---
    case 7: // a drop grows under supersaturation
        return nonneg(amrex::Real(5.0e-6) - grow_droplet(amrex::Real(1.02), kT, amrex::Real(5.0e-6), zero, amrex::Real(1.0), false));
    case 8: // a drop shrinks under subsaturation
        return nonneg(grow_droplet(amrex::Real(0.98), kT, amrex::Real(2.0e-5), zero, amrex::Real(1.0), false) - amrex::Real(2.0e-5));
    case 9: { // r^2 grows ~linearly in time: the increment scales with the interval
        const amrex::Real r0 = amrex::Real(2.0e-5);
        const amrex::Real r1 = grow_droplet(amrex::Real(1.02), kT, r0, zero, amrex::Real(0.5), false);
        const amrex::Real r2 = grow_droplet(amrex::Real(1.02), kT, r0, zero, amrex::Real(1.0), false);
        const amrex::Real d1 = r1*r1 - r0*r0;
        const amrex::Real d2 = r2*r2 - r0*r0;
        return nonneg(std::abs(d2 - amrex::Real(2.0)*d1)/d1 - amrex::Real(0.05)); // within 5%
    }

    // --- Ventilation factor (Bayley et al. 2025) ---
    case 10: // f_v >= 1
        return nonneg(amrex::Real(1.0) - ventilationFactor<amrex::Real>(amrex::Real(1.0e-3), kVentA1, kVentB1, kVentA2, kVentB2, kVentFcap));
    case 11: // f_v increases with radius
        return nonneg( ventilationFactor<amrex::Real>(amrex::Real(1.0e-4), kVentA1, kVentB1, kVentA2, kVentB2, kVentFcap)
                     - ventilationFactor<amrex::Real>(amrex::Real(1.0e-3), kVentA1, kVentB1, kVentA2, kVentB2, kVentFcap));
    case 12: // f_v is capped at f_cap
        return nonneg( ventilationFactor<amrex::Real>(amrex::Real(1.0e-1), kVentA1, kVentB1, kVentA2, kVentB2, kVentFcap) - kVentFcap);
    default: // f_v -> 1 as the radius vanishes
        return nonneg( ventilationFactor<amrex::Real>(amrex::Real(1.0e-9), kVentA1, kVentB1, kVentA2, kVentB2, kVentFcap) - amrex::Real(1.001));
    }
}

inline std::string metric_label (const int i)
{
    switch (i) {
    case 0:  return "Maxwell: dR^2/dt radius-independent for large drops";
    case 1:  return "growth: S > 1 increases R";
    case 2:  return "evaporation: S < 1 decreases R";
    case 3:  return "Kelvin: small pure drop evaporates at S = 1";
    case 4:  return "Raoult: solute drives growth at S = 1";
    case 5:  return "Koehler: haze drop in equilibrium at S_eq";
    case 6:  return "Koehler: S_eq peaks at the critical radius";
    case 7:  return "integrator: drop grows under supersaturation";
    case 8:  return "integrator: drop shrinks under subsaturation";
    case 9:  return "integrator: r^2 grows linearly in time";
    case 10: return "ventilation: f_v >= 1";
    case 11: return "ventilation: f_v increases with radius";
    case 12: return "ventilation: f_v capped at f_cap";
    default: return "ventilation: f_v -> 1 as R -> 0";
    }
}

} // namespace condensation_test

#endif
