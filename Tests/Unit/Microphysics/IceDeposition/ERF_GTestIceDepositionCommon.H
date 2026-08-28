#ifndef ERF_GTEST_ICE_DEPOSITION_COMMON_H_
#define ERF_GTEST_ICE_DEPOSITION_COMMON_H_

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_MaterialProperties.H"
#include "ERF_TerminalVelocity.H"
#include "ERF_SuperDropletPCMassChange.H"

// Unit test for ERF's ice vapour-deposition growth and habit prediction
// (Source/Particles/ERF_SuperDropletPCMassChange.H). It marches a single ice
// crystal through 10 min of growth at constant temperature and water saturation
// -- the Takahashi et al. (1991) chamber condition -- by replaying exactly the
// MassChange_SV sequence (dMdt/TI mass step -> dVolume -> growthRatioStar ->
// dLogRadius) with the ice material constants, then reports mass and the a/c
// semi-axes vs temperature for comparison with TH91 (Welss et al. 2024, Fig.10).
// Doing it as a unit test avoids the full dynamical core, which cannot hold the
// box at a constant supersaturation (kinematic_mode does not isolate the SDM).
namespace ice_deposition_test {

using SDMassChangeUtils_SV::dMdt;
using SDMassChangeUtils_SV::TI;

constexpr amrex::Real kRhoIce  = amrex::Real(916.80);     // solid ice density [kg m^-3]
constexpr amrex::Real kP       = amrex::Real(80000.0);    // ambient pressure [Pa] (~2 km)
constexpr amrex::Real kGrowSec = amrex::Real(600.0);      // 10 min
constexpr amrex::Real kDt      = amrex::Real(0.25);       // micro step [s]

struct CrystalState { amrex::Real mass; amrex::Real a; amrex::Real c; };

// Ice MaterialPropertiesCore, matching the values ERF assigns the ice species
// (ERF_MaterialProperties.cpp): sublimation latent heat, R_v, solid density;
// m_is_water = false so coeffMolecularDiffusion takes the Chapman-Enskog branch.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MaterialPropertiesCore ice_material () noexcept
{
    MaterialPropertiesCore m;
    m.m_density          = kRhoIce;
    m.m_mol_weight       = amrex::Real(1.802e-2);
    m.m_lat_vap          = amrex::Real(2.8342e6);     // sublimation
    m.m_lat_fus          = lfus;
    m.m_Rv               = R_v;
    m.m_Tb               = tboil;
    m.m_Nav_by_molweight = MaterialPropertiesCore::s_N_av / amrex::Real(1.802e-2);
    m.m_is_water         = false;
    return m;
}

// Grow one ice crystal from a 2 um sphere for t_total seconds at constant
// temperature T_K and water saturation. Mirrors MassChange_SV exactly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
CrystalState grow_crystal (const amrex::Real T_K, const amrex::Real t_total = kGrowSec) noexcept
{
    const MaterialPropertiesCore ice = ice_material();
    const dMdt<amrex::Real> dmdt{ amrex::Real(2.8342e6), therco, R_v, kRhoIce };
    const TerminalVelocity<amrex::Real> tv{ amrex::Real(1000.0), kRhoIce };

    // Constant ambient forcing (water saturation): e over ice/water [Pa].
    const amrex::Real e_si  = erf_esati(T_K) * amrex::Real(100.0);   // hPa -> Pa
    const amrex::Real e_sw  = erf_esatw(T_K) * amrex::Real(100.0);
    const amrex::Real S     = e_sw / e_si;          // ice supersaturation ratio at water sat.
    const amrex::Real wi    = e_sw / e_si;          // e_sat ratio water/ice
    const amrex::Real rhom  = kP / (R_d * T_K);     // moist air density
    const amrex::Real D     = ice.coeffMolecularDiffusion(T_K, kP);

    // initial 2 um solid-ice sphere
    amrex::Real a = amrex::Real(2.0e-6), c = amrex::Real(2.0e-6);
    amrex::Real mass = kRhoIce * four_thirds_pi * a * a * c;

    const int nsteps = static_cast<int>(t_total / kDt);
    for (int s = 0; s < nsteps; ++s) {
        const amrex::Real vol_old = four_thirds_pi * a * a * c;
        const amrex::Real rhoi_old = mass / vol_old;
        const amrex::Real vterm = tv.IceBohm(mass, a, c, rhoi_old, rhom, T_K);

        // (1) mass: one forward-Euler step of dm/dt = 4 pi C (S-1) f_v / (Fk+Fd)
        TI<dMdt<amrex::Real>, amrex::Real> ti{ dmdt, kDt, kDt, 100,
                                               a, c, vterm, S, rhom, T_K, e_si, D, false };
        amrex::Real mass_new = mass;
        bool ok = false;
        ti.fe(mass_new, ok);
        mass_new = amrex::max(mass_new, amrex::Real(3.8403e-24));
        const amrex::Real d_mass = mass_new - mass;

        // (2) volume / apparent density via the deposition density
        amrex::Real d_vol = dmdt.dVolume(d_mass, a, c, rhoi_old, S, T_K, e_si, wi, D);
        amrex::Real vol_new = vol_old + d_vol;
        vol_new = amrex::max(vol_new, mass_new / kRhoIce);

        // (3) habit: split the volume change between a and c by the growth ratio
        const amrex::Real gr_star = dmdt.growthRatioStar(d_mass, a, c, vterm, rhom, T_K, D);
        const amrex::Real d_loga = dmdt.dLogRadius(gr_star, std::log(vol_new / vol_old));
        amrex::Real a_new = a * std::exp(d_loga);
        amrex::Real c_new = c * std::exp(gr_star * d_loga);

        // size floors (as in MassChange_SV)
        if (amrex::min(a_new, c_new) <= amrex::Real(1.0e-6)) {
            a_new = c_new = std::cbrt((mass_new / kRhoIce) / four_thirds_pi);
        }
        if (amrex::min(a_new, c_new) <= amrex::Real(1.0e-9)) {
            a_new = c_new = amrex::Real(1.0e-9);
        }
        a = a_new; c = c_new; mass = mass_new;
    }
    return { mass, a, c };
}

// Net ice mass rate dm/dt at saturation ratio S over ice, for a fixed 50 um ice
// sphere -- the deposition/sublimation branch of the production routine
// (dMdt::rhs_func, proportional to S-1). Positive grows the crystal (deposition),
// negative shrinks it (sublimation).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real deposition_mass_rate (const amrex::Real T_K, const amrex::Real S) noexcept
{
    const MaterialPropertiesCore ice = ice_material();
    const dMdt<amrex::Real> dmdt{ amrex::Real(2.8342e6), therco, R_v, kRhoIce };
    const TerminalVelocity<amrex::Real> tv{ amrex::Real(1000.0), kRhoIce };

    const amrex::Real a = amrex::Real(5.0e-5), c = amrex::Real(5.0e-5);
    const amrex::Real mass  = kRhoIce * four_thirds_pi * a * a * c;
    const amrex::Real e_si  = erf_esati(T_K) * amrex::Real(100.0);
    const amrex::Real rhom  = kP / (R_d * T_K);
    const amrex::Real D     = ice.coeffMolecularDiffusion(T_K, kP);
    const amrex::Real vterm = tv.IceBohm(mass, a, c, kRhoIce, rhom, T_K);

    return dmdt.rhs_func(a, c, vterm, S, rhom, T_K, e_si, D);
}

// Convenience reductions.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real aspect_ratio (const CrystalState& x) noexcept { return x.c / x.a; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real mass_ug      (const CrystalState& x) noexcept { return x.mass * amrex::Real(1.0e9); }

// --- TH91 markers (Welss et al. 2024, Fig.10, 10 min) for host-side checks ----
// aspect ratio phi = c/a vs T [C]; the habit signature: plates (phi<1) near -15 C,
// columns (phi>1) near -6 C.
constexpr int kNumTH91Phi = 13;
constexpr amrex::Real kTH91PhiT[kNumTH91Phi]   = {-23,-22,-21,-16,-15,-14,-13,-12,-9,-8,-7,-6,-4};
constexpr amrex::Real kTH91Phi[kNumTH91Phi]    = {1.05,0.82,0.65,0.020,0.010,0.013,0.027,0.076,0.75,0.86,1.34,7.36,0.91};

} // namespace ice_deposition_test

#endif
