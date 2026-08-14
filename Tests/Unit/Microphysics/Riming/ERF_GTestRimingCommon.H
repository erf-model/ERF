#ifndef ERF_GTEST_RIMING_COMMON_H_
#define ERF_GTEST_RIMING_COMMON_H_

#include <cmath>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_SuperDropletPCRiming.H"

// Unit tests for ERF's riming helpers (Source/Particles/ERF_SuperDropletPCRiming.H,
// namespace SDRiming, factored out of ERF_SuperDropletPCCoalescence.cpp): the
// dynamic viscosity, the Heymsfield-Pflaum (1985) rime density, and the ice
// surface temperature that sets the rime regime. The rime density is bounded to
// 100-910 kg/m^3 and densifies with impact speed; the ice surface temperature
// rises above ambient under deposition heating.
//
// Each check is a non-negative "violation" returned by riming_metric(i): 0 when
// the property holds, > 0 when it fails. Host (Scalar) and device (Kernel) tests
// call the same AMREX_GPU_HOST_DEVICE function.
namespace riming_test {

using namespace SDRiming;
using SDMassChangeUtils_SV::dMdt;

constexpr amrex::Real kRhoWater = amrex::Real(1000.0);
constexpr amrex::Real kRhoIce   = amrex::Real(916.8);
constexpr amrex::Real kP        = amrex::Real(80000.0);   // ambient pressure [Pa] (~2 km)
constexpr amrex::Real kTice     = amrex::Real(263.15);    // -10 C, for the surface-temperature checks

constexpr int kNumMetrics = 7;
constexpr amrex::Real kTol = amrex::Real(1.0e-6);

// dMdt functor with the ice deposition constants (sublimation latent heat).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
dMdt<amrex::ParticleReal> ice_dmdt () noexcept
{
    return dMdt<amrex::ParticleReal>{ amrex::ParticleReal(lsub), amrex::ParticleReal(therco),
                                      amrex::ParticleReal(R_v), amrex::ParticleReal(kRhoIce) };
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real nonneg (const amrex::Real x) noexcept { return amrex::max(amrex::Real(0.0), x); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real band_violation (const amrex::Real v, const amrex::Real lo, const amrex::Real hi) noexcept
{
    return amrex::max(amrex::Real(0.0), lo - v) + amrex::max(amrex::Real(0.0), v - hi);
}

// Rime density of a 1 mm graupel sphere collecting a 20 um droplet in
// water-saturated air at temperature T, with the droplet and ice falling at
// vz_w and vz_i.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rime_rho (const amrex::Real T, const amrex::Real vz_w, const amrex::Real vz_i) noexcept
{
    using amrex::ParticleReal;
    const dMdt<ParticleReal> dm = ice_dmdt();
    const ParticleReal e_sw = ParticleReal(erf_esatw(T) * amrex::Real(100.0));
    const ParticleReal e_si = ParticleReal(erf_esati(T) * amrex::Real(100.0));
    amrex::Real qsat_i = amrex::Real(0.0); erf_qsati(T, kP*amrex::Real(0.01), qsat_i);  // erf_qsati wants mbar
    const ParticleReal qv   = ParticleReal(qsat_i) * (e_sw/e_si);   // water-saturated mixing ratio
    const ParticleReal rhom = ParticleReal(kP) / ParticleReal(R_d * T);

    return amrex::Real( rimeDensity_HeymsfieldPflaum1985(
        ParticleReal(2.0e-5), ParticleReal(1.0e-3), ParticleReal(1.0e-3),
        ParticleReal(vz_w), ParticleReal(vz_i), ParticleReal(kRhoWater), ParticleReal(diffelq),
        dm, ParticleReal(T), rhom, ParticleReal(kP), qv) );
}

// Ice surface temperature [K] at ice supersaturation ratio s (qv = s * qsat_ice).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ice_surf_T (const amrex::Real T, const amrex::Real s) noexcept
{
    using amrex::ParticleReal;
    const dMdt<ParticleReal> dm = ice_dmdt();
    amrex::Real qsat_i = amrex::Real(0.0); erf_qsati(T, kP*amrex::Real(0.01), qsat_i);  // erf_qsati wants mbar
    const ParticleReal qv = ParticleReal(qsat_i) * ParticleReal(s);
    return amrex::Real( iceSurfaceTemperature(ParticleReal(T), ParticleReal(kP), qv,
                                              ParticleReal(diffelq), dm) );
}

// The i-th property violation (>= 0; 0 when the property holds).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real riming_metric (const int i) noexcept
{
    switch (i) {
    // --- Heymsfield-Pflaum rime density ---
    case 0: // physical density bounds [100, 910] kg/m^3 (cold, slow)
        return band_violation(rime_rho(amrex::Real(258.15), amrex::Real(1.0), amrex::Real(3.0)),
                              amrex::Real(100.0), amrex::Real(910.0));
    case 1: // bounds hold at a warmer, faster-impact state too
        return band_violation(rime_rho(amrex::Real(268.15), amrex::Real(0.5), amrex::Real(4.0)),
                              amrex::Real(100.0), amrex::Real(910.0));
    case 2: // faster impact -> denser rime (low-density regime, cold surface)
        return nonneg(rime_rho(amrex::Real(258.15), amrex::Real(1.0), amrex::Real(2.0))
                    - rime_rho(amrex::Real(258.15), amrex::Real(1.0), amrex::Real(6.0)));

    // --- Ice surface temperature ---
    case 3: // deposition heating: surface warmer than ambient when supersaturated over ice
        return nonneg(kTice - ice_surf_T(kTice, amrex::Real(1.1)));
    case 4: // at exact ice saturation the surface sits at ambient temperature
        return std::abs(ice_surf_T(kTice, amrex::Real(1.0)) - kTice);
    case 5: // stronger supersaturation -> warmer surface
        return nonneg(ice_surf_T(kTice, amrex::Real(1.05)) - ice_surf_T(kTice, amrex::Real(1.2)));

    // --- Dynamic viscosity helper ---
    default: // viscosity is positive and increases with temperature
        return nonneg(viscCoeff(amrex::ParticleReal(250.0)) - viscCoeff(amrex::ParticleReal(280.0)))
             + nonneg(-amrex::Real(viscCoeff(amrex::ParticleReal(250.0))));
    }
}

inline std::string metric_label (const int i)
{
    switch (i) {
    case 0:  return "rime density in [100,910] kg/m^3 (cold)";
    case 1:  return "rime density in [100,910] kg/m^3 (warm)";
    case 2:  return "rime densifies with impact speed";
    case 3:  return "ice surface warmer than ambient when supersaturated";
    case 4:  return "ice surface at ambient when ice-saturated";
    case 5:  return "ice surface warms with supersaturation";
    default: return "viscosity positive and rising with T";
    }
}

} // namespace riming_test

#endif
