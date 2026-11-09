/**
 * \file ERF_MorrisonUtils.H
 * \brief Parameters and pointwise helpers of the Morrison two-moment microphysics.
 *
 * Everything here is a scalar function or a plain struct, so the Morrison advance
 * and any other part of the code (diagnostics, radiation, initialization, ...) can
 * compute Morrison's quantities -- saturation, size distributions, fall speeds,
 * process rates -- exactly the way the scheme does. The formulas follow
 * module_mp_morr_two_moment.F90; the comments name the Fortran variables.
 *
 * Units are those of the Fortran: SI, pressure in Pa, mixing ratios in kg/kg and
 * number concentrations in #/kg.
 */
#ifndef ERF_MORRISON_UTILS_H_
#define ERF_MORRISON_UTILS_H_

#include <algorithm>
#include <cmath>

#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include "ERF_IndexDefines.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_NumericalConstants.H"
#include "ERF_MorrisonGammaFunction.H"
#include "ERF_MorrisonVaporPressure.H"

/**
 * Switches and constants of the Morrison scheme (the module variables of
 * module_mp_morr_two_moment), including the precomputed CONS1..CONS41.
 * Build it with make_morrison_params(); it is cheap to copy into device kernels.
 */
struct MorrisonParams
{
    // Switches
    int iact;    // CCN activation option (1: power-law, 2: lognormal aerosol)
    int inum;    // Droplet number option (0: predict, 1: constant)
    int iliq;    // Liquid-only option (0: include ice, 1: liquid only)
    int inuc;    // Ice nucleation option (0: mid-latitude, 1: arctic)
    int igraup;  // Graupel option (0: include graupel, 1: no graupel)
    int ihail;   // Graupel/hail option (0: graupel, 1: hail)

    amrex::Real ndcnst; // Constant droplet number concentration (cm^-3), used when INUM = 1

    // Physical constants
    amrex::Real pi;     // Pi
    amrex::Real R;      // Gas constant for dry air (J/kg/K)
    amrex::Real Rv;     // Gas constant for water vapor (J/kg/K)
    amrex::Real g;      // Gravitational acceleration (m/s^2)
    amrex::Real ep_2;   // Molecular weight ratio (Rd/Rv)

    // Reference density and species densities (kg/m^3)
    amrex::Real rhosu;  // Standard air density at 850 mb
    amrex::Real rhow;   // Liquid water
    amrex::Real rhoi;   // Bulk density of cloud ice
    amrex::Real rhosn;  // Bulk density of snow
    amrex::Real rhog;   // Bulk density of graupel/hail

    // Fall speed parameters (V = A D^B)
    amrex::Real ai, bi; // Cloud ice
    amrex::Real bc;     // Cloud droplets (A is the temperature-dependent Stokes ACN)
    amrex::Real as, bs; // Snow
    amrex::Real ar, br; // Rain
    amrex::Real ag, bg; // Graupel/hail

    // Microphysical parameters
    amrex::Real aimm;   // Parameter in Bigg immersion freezing
    amrex::Real bimm;   // Parameter in Bigg immersion freezing
    amrex::Real ecr;    // Collection efficiency between droplets/rain and snow/rain
    amrex::Real dcs;    // Threshold size for cloud ice autoconversion (m)
    amrex::Real mi0;    // Initial mass of nucleated ice crystal (kg)
    amrex::Real mg0;    // Mass of embryo graupel (kg)
    amrex::Real f1s;    // Ventilation parameter for snow
    amrex::Real f2s;    // Ventilation parameter for snow
    amrex::Real f1r;    // Ventilation parameter for rain
    amrex::Real f2r;    // Ventilation parameter for rain
    amrex::Real qsmall; // Smallest allowed hydrometeor mixing ratio
    amrex::Real eii;    // Collection efficiency, ice-ice collisions
    amrex::Real eci;    // Collection efficiency, ice-droplet collisions
    amrex::Real cpw;    // Specific heat of liquid water (J/kg/K)
    amrex::Real rin;    // Radius of contact nuclei (m)
    amrex::Real mmult;  // Mass of splintered ice particle (kg)

    // Size distribution parameters (M = C D^D)
    amrex::Real ci, di; // Cloud ice
    amrex::Real cs, ds; // Snow
    amrex::Real cg, dg; // Graupel

    // Limits of the slope parameter lambda (1/m)
    amrex::Real lammaxi, lammini; // Cloud ice
    amrex::Real lammaxr, lamminr; // Rain
    amrex::Real lammaxs, lammins; // Snow
    amrex::Real lammaxg, lamming; // Graupel

    // Aerosol concentrations (IACT = 2), used to bound predicted droplet number (m^-3)
    amrex::Real nanew1; // Mode 1
    amrex::Real nanew2; // Mode 2

    // Precomputed constants
    amrex::Real cons1, cons2, cons3, cons4, cons5;
    amrex::Real cons6, cons7, cons8, cons9, cons10;
    amrex::Real cons11, cons12, cons13, cons14, cons15;
    amrex::Real cons16, cons17, cons18, cons19, cons20;
    amrex::Real cons21, cons22, cons23, cons24, cons25;
    amrex::Real cons26, cons27, cons28, cons29;
    amrex::Real cons31, cons32, cons34, cons35;
    amrex::Real cons36, cons37, cons38, cons39, cons40;
    amrex::Real cons41;
};

/**
 * Set the Morrison switches and constants (MORR_TWO_MOMENT_INIT).
 *
 * @param[in] no_ice  true for Morrison_NoIce (liquid only, no graupel)
 * @param[in] ndcnst  constant droplet number concentration (cm^-3)
 */
inline MorrisonParams
make_morrison_params (const bool no_ice, const amrex::Real ndcnst)
{
    using amrex::Real;
    MorrisonParams mp;

    mp.iact   = 2;              // Lognormal aerosol activation
    mp.inum   = 1;              // Constant droplet number concentration
    mp.iliq   = no_ice ? 1 : 0; // Liquid only for Morrison_NoIce
    mp.inuc   = 0;              // Mid-latitude ice nucleation (Cooper)
    mp.igraup = no_ice ? 1 : 0; // No graupel for Morrison_NoIce
    mp.ihail  = 0;              // Graupel (0) rather than hail (1)

    mp.ndcnst = ndcnst;

    mp.pi   = Real(3.1415926535897932384626434);
    mp.R    = Real(287.0);
    mp.Rv   = Real(461.6);
    mp.g    = Real(9.81);
    mp.ep_2 = Real(287.0) / mp.Rv;

    mp.rhosu = Real(85000.0)/(Real(287.15)*Real(273.15));
    mp.rhow  = Real(997.0);
    mp.rhoi  = Real(500.0);
    mp.rhosn = Real(100.0);
    mp.rhog  = (mp.ihail == 0) ? Real(400.0) : Real(900.0);

    mp.ai = Real(700.0);
    mp.bi = one;
    mp.bc = Real(2);
    mp.as = Real(11.72);
    mp.bs = Real(0.41);
    mp.ar = Real(841.99667);
    mp.br = Real(0.8);
    if (mp.ihail == 0) {
        mp.ag = Real(19.3);
        mp.bg = Real(0.37);
    } else {
        // Hail (Matsun and Huggins 1980)
        mp.ag = Real(114.5);
        mp.bg = myhalf;
    }

    mp.aimm   = Real(0.66);
    mp.bimm   = Real(100.0);
    mp.ecr    = one;
    mp.dcs    = Real(125.0E-6);
    mp.mi0    = Real(4.0)/three*mp.pi*mp.rhoi*amrex::Math::powi<3>(Real(10.0E-6));
    mp.mg0    = Real(1.6E-10);
    mp.f1s    = Real(0.86);
    mp.f2s    = Real(0.28);
    mp.f1r    = Real(0.78);
    mp.f2r    = Real(0.308);
    mp.qsmall = Real(1.0E-14);
    mp.eii    = Real(0.1);
    mp.eci    = Real(0.7);
    mp.cpw    = Real(4187.0);

    mp.ci = mp.rhoi * mp.pi / Real(6.0);
    mp.di = three;
    mp.cs = mp.rhosn * mp.pi / Real(6.0);
    mp.ds = three;
    mp.cg = mp.rhog * mp.pi / Real(6.0);
    mp.dg = three;

    mp.rin   = Real(0.1E-6);
    mp.mmult = Real(4.0)/three*mp.pi*mp.rhoi*amrex::Math::powi<3>(Real(5.0E-6));

    mp.lammaxi = one/Real(1.0E-6);
    mp.lammini = one/(Real(2)*mp.dcs + Real(100.0E-6));
    mp.lammaxr = one/Real(20.0E-6);
    mp.lamminr = one/Real(2800.0E-6);
    mp.lammaxs = one/Real(10.0E-6);
    mp.lammins = one/Real(2000.0E-6);
    mp.lammaxg = one/Real(20.0E-6);
    mp.lamming = one/Real(2000.0E-6);

    // Aerosol size distribution for MPACE (Morrison et al. 2007, JGR). The other
    // activation parameters of the Fortran are not needed while INUM = 1.
    mp.nanew1 = Real(72.2E6);
    mp.nanew2 = Real(1.8E6);

    mp.cons1  = gamma_function(one + mp.ds) * mp.cs;
    mp.cons2  = gamma_function(one + mp.dg) * mp.cg;
    mp.cons3  = gamma_function(Real(4.0) + mp.bs) / Real(6.0);
    mp.cons4  = gamma_function(Real(4.0) + mp.br) / Real(6.0);
    mp.cons5  = gamma_function(one + mp.bs);
    mp.cons6  = gamma_function(one + mp.br);
    mp.cons7  = gamma_function(Real(4.0) + mp.bg) / Real(6.0);
    mp.cons8  = gamma_function(one + mp.bg);
    mp.cons9  = gamma_function(Real(5.0)/Real(2) + mp.br/Real(2));
    mp.cons10 = gamma_function(Real(5.0)/Real(2) + mp.bs/Real(2));
    mp.cons11 = gamma_function(Real(5.0)/Real(2) + mp.bg/Real(2));
    mp.cons12 = gamma_function(one + mp.di) * mp.ci;
    mp.cons13 = gamma_function(mp.bs + three) * mp.pi / Real(4.0) * mp.eci;
    mp.cons14 = gamma_function(mp.bg + three) * mp.pi / Real(4.0) * mp.eci;
    mp.cons15 = -Real(1108.0) * mp.eii * std::pow(mp.pi, (one-mp.bs)/three) *
        std::pow(mp.rhosn, (-Real(2)-mp.bs)/three) / (Real(4.0)*Real(720.0));
    mp.cons16 = gamma_function(mp.bi + three) * mp.pi / Real(4.0) * mp.eci;
    mp.cons17 = Real(4.0) * Real(2) * three * mp.rhosu * mp.pi * mp.eci * mp.eci *
        gamma_function(Real(2)*mp.bs + Real(2)) / (Real(8.0)*(mp.rhog-mp.rhosn));
    mp.cons18 = mp.rhosn * mp.rhosn;
    mp.cons19 = mp.rhow * mp.rhow;
    mp.cons20 = Real(20.0) * mp.pi * mp.pi * mp.rhow * mp.bimm;
    mp.cons21 = Real(4.0) / (mp.dcs * mp.rhoi);
    mp.cons22 = mp.pi * mp.rhoi * amrex::Math::powi<3>(mp.dcs) / Real(6.0);
    mp.cons23 = mp.pi / Real(4.0) * mp.eii * gamma_function(mp.bs + three);
    mp.cons24 = mp.pi / Real(4.0) * mp.ecr * gamma_function(mp.br + three);
    mp.cons25 = mp.pi * mp.pi / Real(24.0) * mp.rhow * mp.ecr * gamma_function(mp.br + Real(6.0));
    mp.cons26 = mp.pi / Real(6.0) * mp.rhow;
    mp.cons27 = gamma_function(one + mp.bi);
    mp.cons28 = gamma_function(Real(4.0) + mp.bi) / Real(6.0);
    mp.cons29 = Real(4.0)/three * mp.pi * mp.rhow * amrex::Math::powi<3>(Real(25.0E-6));
    mp.cons31 = mp.pi * mp.pi * mp.ecr * mp.rhosn;
    mp.cons32 = mp.pi / Real(2) * mp.ecr;
    mp.cons34 = Real(5.0)/Real(2) + mp.br/Real(2);
    mp.cons35 = Real(5.0)/Real(2) + mp.bs/Real(2);
    mp.cons36 = Real(5.0)/Real(2) + mp.bg/Real(2);
    mp.cons37 = Real(4.0) * mp.pi * Real(1.38E-23) / (Real(6.0) * mp.pi * mp.rin);
    mp.cons38 = mp.pi * mp.pi / three * mp.rhow;
    mp.cons39 = mp.pi * mp.pi / Real(36.0) * mp.rhow * mp.bimm;
    mp.cons40 = mp.pi / Real(6.0) * mp.bimm;
    mp.cons41 = mp.pi * mp.pi * mp.ecr * mp.rhow;

    return mp;
}

// ---------------------------------------------------------------------------
// State conversion
// ---------------------------------------------------------------------------

/**
 * Fill one cell of Morrison's microphysics variables from the conserved state:
 * density, potential temperature, non-negative mixing ratios, temperature and
 * pressure (in Pa, Morrison's native unit).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_copy_state_to_micro_cell (
    const amrex::Array4<const amrex::Real>& states,
    const amrex::Array4<const amrex::Real>& base,
    const amrex::Array4<amrex::Real>& rho,
    const amrex::Array4<amrex::Real>& theta,
    const amrex::Array4<amrex::Real>& qv,
    const amrex::Array4<amrex::Real>& qc,
    const amrex::Array4<amrex::Real>& qi,
    const amrex::Array4<amrex::Real>& qn,
    const amrex::Array4<amrex::Real>& qt,
    const amrex::Array4<amrex::Real>& qpr,
    const amrex::Array4<amrex::Real>& qps,
    const amrex::Array4<amrex::Real>& qpg,
    const amrex::Array4<amrex::Real>& qp,
    const amrex::Array4<amrex::Real>& tabs,
    const amrex::Array4<amrex::Real>& pres,
    const amrex::Real rdOcp,
    const bool use_anelastic_reference_pressure,
    const int i, const int j, const int k) noexcept
{
    const amrex::Real rho_value = states(i,j,k,Rho_comp);
    const amrex::Real rho_theta = states(i,j,k,RhoTheta_comp);
    rho(i,j,k) = rho_value;
    theta(i,j,k) = rho_theta / rho_value;
    qv(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ1_comp) / rho_value);
    qc(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ2_comp) / rho_value);
    qi(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ3_comp) / rho_value);
    qn(i,j,k) = qc(i,j,k) + qi(i,j,k);
    qt(i,j,k) = qv(i,j,k) + qn(i,j,k);
    qpr(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ4_comp) / rho_value);
    qps(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ5_comp) / rho_value);
    qpg(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ6_comp) / rho_value);
    qp(i,j,k) = qpr(i,j,k) + qps(i,j,k) + qpg(i,j,k);

    const amrex::Real p0 = use_anelastic_reference_pressure
        ? base(i,j,k,BaseState::p0_comp) : amrex::Real(0.0);
    const MicrophysicsThermoState thermo = diagnose_microphysics_thermo_state(
        rho_value, rho_theta, qv(i,j,k), rdOcp,
        use_anelastic_reference_pressure, p0);
    tabs(i,j,k) = thermo.temperature;
    // Morrison's native pressure unit is Pa.
    pres(i,j,k) = thermo.pressure_pa;
}

/**
 * Droplet number mixing ratio (#/kg) of the constant droplet concentration
 * NDCNST (#/cm^3) in air of density rho (kg/m^3).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_ndcnst_to_number_mixing_ratio (const amrex::Real ndcnst,
                                                    const amrex::Real rho) noexcept
{
    return ndcnst * amrex::Real(1.0e6) / rho;
}

// ---------------------------------------------------------------------------
// Thermodynamics and air properties
// ---------------------------------------------------------------------------

//! Latent heat of vaporization XXLV (J/kg)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_latent_heat_vaporization (const amrex::Real t) noexcept
{
    return amrex::Real(3.1484E6) - amrex::Real(2370.0) * t;
}

//! Latent heat of sublimation XXLS (J/kg)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_latent_heat_sublimation (const amrex::Real t) noexcept
{
    return amrex::Real(3.15E6) - amrex::Real(2370.0) * t + amrex::Real(0.3337E6);
}

//! Specific heat at constant pressure of moist air CPM (J/kg/K)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_cpm (const amrex::Real qv) noexcept
{
    return amrex::Real(1004.5) * (one + amrex::Real(0.887) * qv);
}

//! Saturation vapor pressure over liquid (Pa), capped at 0.99 p for low pressure
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_svp_liquid (const amrex::Real t, const amrex::Real pres) noexcept
{
    return std::min(amrex::Real(0.99) * pres, calc_saturation_vapor_pressure(t, 0));
}

//! Saturation vapor pressure over ice (Pa), capped at 0.99 p for low pressure
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_svp_ice (const amrex::Real t, const amrex::Real pres) noexcept
{
    return std::min(amrex::Real(0.99) * pres, calc_saturation_vapor_pressure(t, 1));
}

//! Mixing ratio (kg/kg) of vapor at partial pressure e (Pa)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_mixing_ratio (const amrex::Real ep_2, const amrex::Real e,
                                   const amrex::Real pres) noexcept
{
    return ep_2 * e / (pres - e);
}

//! Saturation mixing ratios and saturation ratios over liquid and ice
struct MorrisonSaturation
{
    amrex::Real qvs;    // QVS: saturation mixing ratio over liquid
    amrex::Real qvi;    // QVI: saturation mixing ratio over ice
    amrex::Real qvqvs;  // QVQVS: saturation ratio over liquid
    amrex::Real qvqvsi; // QVQVSI: saturation ratio over ice
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonSaturation morrison_saturation (const amrex::Real t, const amrex::Real pres,
                                        const amrex::Real qv, const amrex::Real ep_2) noexcept
{
    const amrex::Real evs = morrison_svp_liquid(t, pres);
    amrex::Real eis = morrison_svp_ice(t, pres);
    // Make sure ice saturation doesn't exceed water saturation near freezing
    if (eis > evs) {
        eis = evs;
    }
    MorrisonSaturation s;
    s.qvs = morrison_mixing_ratio(ep_2, evs, pres);
    s.qvi = morrison_mixing_ratio(ep_2, eis, pres);
    s.qvqvs  = qv / s.qvs;
    s.qvqvsi = qv / s.qvi;
    return s;
}

//! Dynamic viscosity of air MU (kg/m/s)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_air_viscosity (const amrex::Real t) noexcept
{
    return amrex::Real(1.496e-6) * std::pow(t, amrex::Real(1.5)) / (t + amrex::Real(120.0));
}

//! Thermal conductivity of air KAP (J/m/s/K)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_thermal_conductivity (const amrex::Real mu) noexcept
{
    return amrex::Real(1.414e3) * mu;
}

//! Diffusivity of water vapor in air DV (m^2/s)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_vapor_diffusivity (const amrex::Real t, const amrex::Real pres) noexcept
{
    return amrex::Real(8.794e-5) * std::pow(t, amrex::Real(1.81)) / pres;
}

/**
 * Air-density correction of fall speeds, (RHOSU/RHO)^0.54. Also scales the
 * upper limits on the rain, snow and graupel fall speeds.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_fallspeed_density_factor (const amrex::Real rhosu,
                                               const amrex::Real rho) noexcept
{
    return std::pow(rhosu / rho, amrex::Real(0.54));
}

//! Air-density corrected fall speed coefficients A of V = A D^B
struct MorrisonFallSpeedCoefficients
{
    amrex::Real ain; // Cloud ice (Ikawa and Saito 1991 correction)
    amrex::Real arn; // Rain
    amrex::Real asn; // Snow
    amrex::Real acn; // Cloud droplets (temperature-dependent Stokes)
    amrex::Real agn; // Graupel
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonFallSpeedCoefficients
morrison_fallspeed_coefficients (const amrex::Real rho, const amrex::Real mu,
                                 const MorrisonParams& mp) noexcept
{
    const amrex::Real dum = morrison_fallspeed_density_factor(mp.rhosu, rho);
    MorrisonFallSpeedCoefficients c;
    c.ain = std::pow(mp.rhosu / rho, amrex::Real(0.35)) * mp.ai;
    c.arn = dum * mp.ar;
    c.asn = dum * mp.as;
    c.acn = mp.g * mp.rhow / (amrex::Real(18.0) * mu);
    c.agn = dum * mp.ag;
    return c;
}

// ---------------------------------------------------------------------------
// Size distributions
// ---------------------------------------------------------------------------

/**
 * Slope LAM and intercept N0 of the exponential size distribution of rain, snow,
 * cloud ice or graupel, with mass-diameter coefficient c (PI*RHOW, CONS1, CONS12
 * or CONS2) and exponent 1/d. When LAM falls outside [lammin, lammax] it is
 * limited and the number concentration n is adjusted to match.
 * Call only for q >= QSMALL.
 *
 * TODO(MorrisonUtils): the StdPow4 switch only keeps answers bit-for-bit while the
 * call sites are converted; it is removed in the next commit.
 */
template <bool StdPow4 = false>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_exponential_psd (const amrex::Real q, amrex::Real& n,
                               const amrex::Real c, const amrex::Real inv_d,
                               const amrex::Real lammin, const amrex::Real lammax,
                               amrex::Real& lam, amrex::Real& n0) noexcept
{
    auto pow4 = [] (amrex::Real x) {
        if constexpr (StdPow4) {
            return std::pow(x, amrex::Real(4.0));
        } else {
            return amrex::Math::powi<4>(x);
        }
    };
    lam = std::pow(c * n / q, inv_d);
    n0 = n * lam;
    if (lam < lammin) {
        lam = lammin;
        n0 = pow4(lam) * q / c;
        n = n0 / lam;
    } else if (lam > lammax) {
        lam = lammax;
        n0 = pow4(lam) * q / c;
        n = n0 / lam;
    }
}

/**
 * Slope of the exponential size distribution limited to [lammin, lammax],
 * without adjusting the number concentration (used for sedimentation).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_limited_slope (const amrex::Real q, const amrex::Real n,
                                    const amrex::Real c, const amrex::Real inv_d,
                                    const amrex::Real lammin, const amrex::Real lammax) noexcept
{
    amrex::Real lam = std::pow(c * n / q, inv_d);
    lam = amrex::max(lam, lammin);
    lam = amrex::min(lam, lammax);
    return lam;
}

//! Shape parameter PGAM of the cloud droplet gamma distribution (Martin et al. 1994)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_cloud_droplet_pgam (const amrex::Real nc, const amrex::Real pres,
                                         const amrex::Real t) noexcept
{
    const amrex::Real rho_air = pres / (amrex::Real(287.15) * t);
    amrex::Real pgam = amrex::Real(0.0005714) * (nc / amrex::Real(1.0e6) * rho_air) + amrex::Real(0.2714);
    pgam = one / (pgam * pgam) - one;
    pgam = amrex::max(pgam, amrex::Real(2));
    pgam = amrex::min(pgam, amrex::Real(10.0));
    return pgam;
}

/**
 * Slope LAMC of the cloud droplet distribution before limiting, given
 * g1 = GAMMA(PGAM+1) and g4 = GAMMA(PGAM+4).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_cloud_droplet_lambda (const amrex::Real qc, const amrex::Real nc,
                                           const amrex::Real cons26,
                                           const amrex::Real g1, const amrex::Real g4) noexcept
{
    return std::pow((cons26 * nc * g4) / (qc * g1), one/three);
}

//! Smallest LAMC, for a 60 micron diameter
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_cloud_droplet_lammin (const amrex::Real pgam) noexcept
{
    return (pgam + one) / amrex::Real(60.0e-6);
}

//! Largest LAMC, for a 1 micron diameter
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_cloud_droplet_lammax (const amrex::Real pgam) noexcept
{
    return (pgam + one) / amrex::Real(1.0e-6);
}

/**
 * Shape PGAM, slope LAMC and intercept CDIST1 of the cloud droplet gamma
 * distribution. When LAMC is limited the droplet number nc is adjusted to match.
 * Call only for qc >= QSMALL.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_cloud_droplet_psd (const amrex::Real qc, amrex::Real& nc,
                                 const amrex::Real pres, const amrex::Real t,
                                 const amrex::Real cons26,
                                 amrex::Real& pgam, amrex::Real& lamc, amrex::Real& cdist1) noexcept
{
    pgam = morrison_cloud_droplet_pgam(nc, pres, t);

    const amrex::Real g1 = gamma_function(pgam + one);
    const amrex::Real g4 = gamma_function(pgam + amrex::Real(4.0));

    lamc = morrison_cloud_droplet_lambda(qc, nc, cons26, g1, g4);

    const amrex::Real lammin = morrison_cloud_droplet_lammin(pgam);
    const amrex::Real lammax = morrison_cloud_droplet_lammax(pgam);
    if (lamc < lammin) {
        lamc = lammin;
        nc = std::exp(three*std::log(lamc) + std::log(qc) + std::log(g1) - std::log(g4)) / cons26;
    } else if (lamc > lammax) {
        lamc = lammax;
        nc = std::exp(three*std::log(lamc) + std::log(qc) + std::log(g1) - std::log(g4)) / cons26;
    }

    cdist1 = nc / g1;
}

//! Effective radius (microns) of an exponential distribution with slope lam
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_effective_radius (const amrex::Real lam) noexcept
{
    return three / lam / amrex::Real(2) * amrex::Real(1.0e6);
}

//! Effective radius (microns) of the cloud droplet gamma distribution
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_cloud_droplet_effective_radius (const amrex::Real pgam,
                                                     const amrex::Real lamc) noexcept
{
    return gamma_function(pgam + amrex::Real(4.0)) / gamma_function(pgam + three) / lamc
        / amrex::Real(2) * amrex::Real(1.0e6);
}

// ---------------------------------------------------------------------------
// Fall speeds
// ---------------------------------------------------------------------------

//! Mass- or number-weighted fall speed a*cons/lam^b of an exponential distribution
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_fallspeed (const amrex::Real a, const amrex::Real cons,
                                const amrex::Real lam, const amrex::Real b) noexcept
{
    return a * cons / std::pow(lam, b);
}

//! Mass- and number-weighted fall speeds (m/s)
struct MorrisonFallSpeed
{
    amrex::Real mass;   // UM*
    amrex::Real number; // UN*
};

//! Rain fall speeds, limited to 9.1 m/s times the density factor
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonFallSpeed morrison_rain_fallspeed (const amrex::Real arn, const amrex::Real lamr,
                                           const amrex::Real rhofac, const MorrisonParams& mp) noexcept
{
    return {std::min(morrison_fallspeed(arn, mp.cons4, lamr, mp.br), amrex::Real(9.1) * rhofac),
            std::min(morrison_fallspeed(arn, mp.cons6, lamr, mp.br), amrex::Real(9.1) * rhofac)};
}

//! Snow fall speeds, limited to 1.2 m/s times the density factor
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonFallSpeed morrison_snow_fallspeed (const amrex::Real asn, const amrex::Real lams,
                                           const amrex::Real rhofac, const MorrisonParams& mp) noexcept
{
    return {std::min(morrison_fallspeed(asn, mp.cons3, lams, mp.bs), amrex::Real(1.2) * rhofac),
            std::min(morrison_fallspeed(asn, mp.cons5, lams, mp.bs), amrex::Real(1.2) * rhofac)};
}

//! Graupel fall speeds, limited to 20 m/s times the density factor
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonFallSpeed morrison_graupel_fallspeed (const amrex::Real agn, const amrex::Real lamg,
                                              const amrex::Real rhofac, const MorrisonParams& mp) noexcept
{
    return {std::min(morrison_fallspeed(agn, mp.cons7, lamg, mp.bg), amrex::Real(20.0) * rhofac),
            std::min(morrison_fallspeed(agn, mp.cons8, lamg, mp.bg), amrex::Real(20.0) * rhofac)};
}

//! Cloud ice fall speeds, limited to 1.2 m/s times (RHOSU/RHO)^0.35
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonFallSpeed morrison_ice_fallspeed (const amrex::Real ain, const amrex::Real lami,
                                          const amrex::Real rho, const MorrisonParams& mp) noexcept
{
    const amrex::Real umax = amrex::Real(1.2) * std::pow(mp.rhosu / rho, amrex::Real(0.35));
    return {std::min(morrison_fallspeed(ain, mp.cons28, lami, mp.bi), umax),
            std::min(morrison_fallspeed(ain, mp.cons27, lami, mp.bi), umax)};
}

//! Cloud droplet fall speeds of the gamma distribution (not limited)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonFallSpeed morrison_cloud_droplet_fallspeed (const amrex::Real acn, const amrex::Real lamc,
                                                    const amrex::Real pgam, const MorrisonParams& mp) noexcept
{
    return {acn * gamma_function(amrex::Real(4.) + mp.bc + pgam) /
                (std::pow(lamc, mp.bc) * gamma_function(pgam + amrex::Real(4.))),
            acn * gamma_function(one + mp.bc + pgam) /
                (std::pow(lamc, mp.bc) * gamma_function(pgam + one))};
}

// ---------------------------------------------------------------------------
// Process rates (per second)
// ---------------------------------------------------------------------------

/**
 * Ventilation term of vapor diffusion and melting (Rutledge and Hobbs),
 * F1/LAM^2 + F2 (A RHO/MU)^1/2 SC^1/3 CONS/LAM^EXP, for fall speed coefficient a
 * and the ventilation constants of the species (CONS9/CONS34 for rain,
 * CONS10/CONS35 for snow, CONS11/CONS36 for graupel).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_ventilation (const amrex::Real f1, const amrex::Real f2,
                                  const amrex::Real a, const amrex::Real rho,
                                  const amrex::Real mu, const amrex::Real sc,
                                  const amrex::Real cons, const amrex::Real lam,
                                  const amrex::Real cons_exp) noexcept
{
    return f1/(lam*lam) + f2 * std::sqrt(a*rho/mu) * std::pow(sc, one/three) * cons / std::pow(lam, cons_exp);
}

//! Inverse phase relaxation time EPS* of a ventilated species (see M2005)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_phase_relaxation (const amrex::Real pi, const amrex::Real n0,
                                       const amrex::Real rho, const amrex::Real dv,
                                       const amrex::Real ventilation) noexcept
{
    return amrex::Real(2) * pi * n0 * rho * dv * ventilation;
}

//! Rain evaporation PRE (Rutledge and Hobbs 1983); no condensation onto rain
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_rain_evaporation (const amrex::Real qv, const amrex::Real qvs,
                                       const amrex::Real epsr, const amrex::Real ab) noexcept
{
    if (qv < qvs) {
        return std::min(epsr * (qv - qvs) / ab, amrex::Real(0));
    }
    return amrex::Real(0);
}

//! Autoconversion of cloud water to rain (Khairoutdinov and Kogan 2000)
struct MorrisonAutoconversion
{
    amrex::Real prc;   // PRC: cloud water mass
    amrex::Real nprc;  // NPRC: change in droplet number
    amrex::Real nprc1; // NPRC1: change in rain number
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonAutoconversion morrison_cloud_autoconversion (const amrex::Real qc, const amrex::Real nc,
                                                      const amrex::Real rho, const amrex::Real cons29,
                                                      const amrex::Real dt) noexcept
{
    MorrisonAutoconversion r{amrex::Real(0), amrex::Real(0), amrex::Real(0)};
    // Use a minimum value of 1e-6 to prevent floating point error
    if (qc >= amrex::Real(1.0e-6)) {
        r.prc = amrex::Real(1350.0) * std::pow(qc, amrex::Real(2.47)) *
            std::pow((nc/amrex::Real(1.0e6)*rho), -amrex::Real(1.79));
        r.nprc1 = r.prc / cons29;
        r.nprc  = r.prc / (qc / nc);
        // hm bug fix 3/20/12
        r.nprc  = std::min(r.nprc, nc / dt);
        r.nprc1 = std::min(r.nprc1, r.nprc);
    }
    return r;
}

//! Mass and number rates of a process
struct MorrisonRates
{
    amrex::Real mass;   // P*
    amrex::Real number; // N*
};

//! Accretion of cloud water by rain, PRA and NPRA (Khairoutdinov and Kogan 2000)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonRates morrison_cloud_accretion_by_rain (const amrex::Real qc, const amrex::Real qr,
                                                     const amrex::Real nc) noexcept
{
    MorrisonRates r{amrex::Real(0), amrex::Real(0)};
    if (qr >= amrex::Real(1.0e-8) && qc >= amrex::Real(1.0e-8)) {
        const amrex::Real dum = qc * qr;
        r.mass = amrex::Real(67.0) * std::pow(dum, amrex::Real(1.15));
        r.number = r.mass / (qc / nc);
    }
    return r;
}

//! Self-collection and breakup of rain drops NRAGG (Beheng 1994)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_rain_self_collection (const amrex::Real qr, const amrex::Real nr,
                                           const amrex::Real lamr, const amrex::Real rho) noexcept
{
    if (qr >= amrex::Real(1.0e-8)) {
        // include breakup add 10/09/09
        const amrex::Real dum1 = amrex::Real(300.0e-6);
        const amrex::Real dum = (one/lamr < dum1)
            ? one : amrex::Real(2) - std::exp(amrex::Real(2300.0) * (one/lamr - dum1));
        return -amrex::Real(5.78) * dum * nr * qr * rho;
    }
    return amrex::Real(0);
}

/**
 * Mass of rain collected by snow or graupel, or of snow collected by rain
 * (Ikawa and Saito 1991): PRACS and PRACG with (lam1, lam2) = (LAMR, LAMx) and
 * CONS41, PSACR with (LAMS, LAMR) and CONS31. umr and umx are the mass-weighted
 * fall speeds of rain and of the other species.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_rain_collection_mass (const amrex::Real cons,
                                           const amrex::Real umr, const amrex::Real umx,
                                           const amrex::Real rho,
                                           const amrex::Real n0r, const amrex::Real n0x,
                                           const amrex::Real lam1, const amrex::Real lam2) noexcept
{
    using amrex::Math::powi;
    return cons * (std::sqrt(powi<2>(amrex::Real(1.2)*umr - amrex::Real(0.95)*umx) +
                             amrex::Real(0.08)*umx*umr) * rho * n0r * n0x / powi<3>(lam1) *
                   (amrex::Real(5.0)/(powi<3>(lam1) * lam2) +
                    amrex::Real(2)/(powi<2>(lam1) * powi<2>(lam2)) +
                    myhalf/(lam1 * powi<3>(lam2))));
}

/**
 * Number of rain drops collected by snow or graupel, NPRACS and NPRACG; unr and
 * unx are the number-weighted fall speeds of rain and of the other species.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_rain_collection_number (const amrex::Real cons32,
                                             const amrex::Real unr, const amrex::Real unx,
                                             const amrex::Real rho,
                                             const amrex::Real n0r, const amrex::Real n0x,
                                             const amrex::Real lamr, const amrex::Real lamx) noexcept
{
    using amrex::Math::powi;
    return cons32 * rho * std::sqrt(amrex::Real(1.7) * powi<2>(unr - unx) +
                                    amrex::Real(0.3) * unr * unx) * n0r * n0x *
        (one / (powi<3>(lamr) * lamx) +
         one / (powi<2>(lamr) * powi<2>(lamx)) +
         one / (lamr * powi<3>(lamx)));
}

/**
 * Continuous collection with a gravitational kernel by an exponential
 * distribution (fall speed a D^b, intercept n0, slope lam) of a species of
 * mixing ratio or number x: PSACWS, PSACWG, PSACWI, PRAI and their numbers.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_continuous_collection (const amrex::Real cons, const amrex::Real a,
                                            const amrex::Real x, const amrex::Real rho,
                                            const amrex::Real n0, const amrex::Real lam,
                                            const amrex::Real b) noexcept
{
    return cons * a * x * rho * n0 / std::pow(lam, (b + three));
}

/**
 * Melting of snow or graupel above freezing, PSMLT or PGMLT (Rutledge and Hobbs
 * 1984), accelerated by the collection of rain prac (PRACS or PRACG).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_melting (const amrex::Real n0, const amrex::Real kap,
                              const amrex::Real t, const amrex::Real xlf,
                              const amrex::Real ventilation, const amrex::Real prac,
                              const MorrisonParams& mp) noexcept
{
    const amrex::Real dum = -mp.cpw/xlf * (t - amrex::Real(273.15)) * prac;
    return amrex::Real(2) * mp.pi * n0 * kap * (amrex::Real(273.15) - t) / xlf * ventilation + dum;
}

/**
 * Snow or graupel melting above freezing, PSMLT/PGMLT, and the evaporation of the
 * melt in water subsaturation, EVPMS/EVPMG. q, n0, lam and a are those of the
 * species, cons and cons_exp its ventilation constants, prac the rain it collects.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_melting_and_evaporation (const amrex::Real q, const amrex::Real n0,
                                       const amrex::Real lam, const amrex::Real a,
                                       const amrex::Real cons, const amrex::Real cons_exp,
                                       const amrex::Real prac,
                                       const amrex::Real t, const amrex::Real qv,
                                       const amrex::Real rho, const amrex::Real mu,
                                       const amrex::Real xlf, const amrex::Real kap,
                                       const amrex::Real dv, const amrex::Real sc,
                                       const amrex::Real ab,
                                       const MorrisonSaturation& sat,
                                       const MorrisonParams& mp,
                                       amrex::Real& pmlt, amrex::Real& evpm) noexcept
{
    if (q >= amrex::Real(1.0e-8)) {
        const amrex::Real vent = morrison_ventilation(mp.f1s, mp.f2s, a, rho, mu, sc, cons, lam, cons_exp);
        pmlt = morrison_melting(n0, kap, t, xlf, vent, prac, mp);

        // In water subsaturation, snow/graupel melts and evaporates
        if (sat.qvqvs < one) {
            const amrex::Real eps = morrison_phase_relaxation(mp.pi, n0, rho, dv, vent);
            evpm = (qv - sat.qvs) * eps / ab;
            evpm = std::max(evpm, pmlt);
            pmlt = pmlt - evpm;
        }
    }
}

/**
 * Number tendency that goes with the loss rate of mass, assuming a 1:1 ratio of
 * relative mass and number losses (NSUBR, NSUBI, NSMLTS, ...).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_number_loss (const amrex::Real mass_rate, const amrex::Real q,
                                  const amrex::Real n, const amrex::Real dt) noexcept
{
    amrex::Real dum = mass_rate * dt / q;
    dum = std::max(-one, dum);
    return dum * n / dt;
}

//! Contact (Meyers et al. 1992) and immersion (Bigg 1953) freezing of cloud droplets
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonRates morrison_cloud_droplet_freezing (const amrex::Real qc, const amrex::Real nc,
                                                    const amrex::Real t, const amrex::Real pres,
                                                    const amrex::Real mu, const amrex::Real cdist1,
                                                    const amrex::Real pgam, const amrex::Real lamc,
                                                    const MorrisonParams& mp, const amrex::Real dt) noexcept
{
    using amrex::Real;
    MorrisonRates r{Real(0), Real(0)}; // MNUCCC, NNUCCC
    // Only allowed below -4C
    if (qc >= mp.qsmall && t < Real(269.15)) {
        // Number of contact nuclei (m^-3) from Meyers et al., 1992; factor of 1000 converts L^-1 to m^-3
        const Real nacnt = std::exp(-Real(2.80) + Real(0.262) * (Real(273.15) - t)) * Real(1000.0);

        // Mean free path
        const Real dum = Real(7.37) * t / (Real(288.0) * Real(10.0) * pres) / Real(100.0);

        // Effective diffusivity of contact nuclei, based on Brownian diffusion
        const Real dap = mp.cons37 * t * (one + dum / mp.rin) / mu;

        // Contact freezing
        r.mass = mp.cons38 * dap * nacnt * std::exp(std::log(cdist1) +
                     std::log(gamma_function(pgam + Real(5.0))) - Real(4.0) * std::log(lamc));
        r.number = Real(2) * mp.pi * dap * nacnt * cdist1 * gamma_function(pgam + Real(2)) / lamc;

        // Immersion freezing (Bigg 1953); hm 7/15/13 fix for consistency w/ original formula
        r.mass = r.mass + mp.cons39 *
            std::exp(std::log(cdist1) + std::log(gamma_function(Real(7.0) + pgam)) - Real(6.0) * std::log(lamc)) *
            (std::exp(mp.aimm * (Real(273.15) - t)) - one);
        r.number = r.number +
            mp.cons40 * std::exp(std::log(cdist1) + std::log(gamma_function(pgam + Real(4.0))) - three * std::log(lamc)) *
            (std::exp(mp.aimm * (Real(273.15) - t)) - one);

        // Prevent divergence between number and mixing ratio, since strict
        // conservation is not checked for number
        r.number = std::min(r.number, nc / dt);
    }
    return r;
}

//! Immersion freezing of rain (Bigg 1953), MNUCCR and NNUCCR; allowed below -4C
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonRates morrison_rain_freezing (const amrex::Real qr, const amrex::Real nr,
                                           const amrex::Real lamr, const amrex::Real t,
                                           const MorrisonParams& mp, const amrex::Real dt) noexcept
{
    using amrex::Real;
    using amrex::Math::powi;
    MorrisonRates r{Real(0), Real(0)};
    if (t < Real(269.15) && qr >= mp.qsmall) {
        // hm fix 7/15/13 for consistency w/ original formula
        r.mass = mp.cons20 * nr * (std::exp(mp.aimm * (Real(273.15) - t)) - one) /
            powi<3>(lamr) / powi<3>(lamr);
        r.number = mp.pi * nr * mp.bimm * (std::exp(mp.aimm * (Real(273.15) - t)) - one) /
            powi<3>(lamr);
        // Prevent divergence between mixing ratio and number
        r.number = std::min(r.number, nr / dt);
    }
    return r;
}

//! Snow aggregation NSAGG (Passarelli 1978, used by Reisner 1998), hard-wired for BS = 0.4
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_snow_aggregation (const amrex::Real qni, const amrex::Real ns,
                                       const amrex::Real asn, const amrex::Real rho,
                                       const MorrisonParams& mp) noexcept
{
    if (qni >= amrex::Real(1.0e-8)) {
        return mp.cons15 * asn * std::pow(rho, ((amrex::Real(2) + mp.bs) / three)) *
            std::pow(qni, ((amrex::Real(2) + mp.bs) / three)) *
            std::pow((ns * rho), ((amrex::Real(4.0) - mp.bs) / three)) / rho;
    }
    return amrex::Real(0);
}

/**
 * Autoconversion of cloud ice to snow, PRCI and NPRCI (Harrington et al. 1995,
 * modified); only while ice grows, i.e. in ice supersaturation.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonRates morrison_ice_autoconversion (const amrex::Real qi, const amrex::Real ni,
                                                const amrex::Real n0i, const amrex::Real lami,
                                                const amrex::Real qv, const amrex::Real rho,
                                                const amrex::Real dv, const amrex::Real abi,
                                                const MorrisonSaturation& sat,
                                                const MorrisonParams& mp, const amrex::Real dt) noexcept
{
    MorrisonRates r{amrex::Real(0), amrex::Real(0)};
    if (qi >= amrex::Real(1.0e-8) && sat.qvqvsi >= one) {
        r.number = mp.cons21 * (qv - sat.qvi) * rho * n0i * std::exp(-lami * mp.dcs) * dv / abi;
        r.mass = mp.cons22 * r.number;
        r.number = std::min(r.number, ni / dt);
    }
    return r;
}

//! Rain-cloud ice collisions (Reisner et al. 1998), producing snow or graupel
struct MorrisonRainIceCollision
{
    amrex::Real niacr; // NIACR: ice number
    amrex::Real piacr; // PIACR: rain mass
    amrex::Real praci; // PRACI: ice mass
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonRainIceCollision morrison_rain_ice_collision (const amrex::Real qi, const amrex::Real ni,
                                                      const amrex::Real nr, const amrex::Real n0r,
                                                      const amrex::Real lamr, const amrex::Real arn,
                                                      const amrex::Real rho,
                                                      const MorrisonParams& mp, const amrex::Real dt) noexcept
{
    MorrisonRainIceCollision r;
    r.niacr = mp.cons24 * ni * n0r * arn / std::pow(lamr, (mp.br + three)) * rho;
    r.piacr = mp.cons25 * ni * n0r * arn / std::pow(lamr, (mp.br + three)) /
        amrex::Math::powi<3>(lamr) * rho;
    r.praci = mp.cons24 * qi * n0r * arn / std::pow(lamr, (mp.br + three)) * rho;
    r.niacr = std::min(r.niacr, nr / dt);
    r.niacr = std::min(r.niacr, ni / dt);
    return r;
}

/**
 * Nucleation of cloud ice from homogeneous and heterogeneous freezing on aerosol,
 * MNUCCD and NNUCCD (Cooper curve for INUC = 0, constant for INUC = 1).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonRates morrison_ice_nucleation (const amrex::Real t, const amrex::Real rho,
                                            const amrex::Real ni, const amrex::Real ns,
                                            const amrex::Real ng,
                                            const MorrisonSaturation& sat,
                                            const MorrisonParams& mp, const amrex::Real dt) noexcept
{
    using amrex::Real;
    MorrisonRates r{Real(0), Real(0)};
    if (mp.inuc == 0) {
        // Add threshold according to Greg Thompson
        if ((sat.qvqvs >= Real(0.999) && t <= Real(265.15)) || sat.qvqvsi >= Real(1.08)) {
            // hm, modify dec. 5, 2006, replace with cooper curve
            Real kc2 = Real(0.005) * std::exp(Real(0.304) * (Real(273.15) - t)) * Real(1000.0); // L^-1 to m^-3
            kc2 = std::min(kc2, Real(500.0e3));     // Limit to 500 L^-1
            kc2 = std::max(kc2 / rho, Real(0));    // Convert to kg^-1
            if (kc2 > ni + ns + ng) {
                r.number = (kc2 - ni - ns - ng) / dt;
                r.mass = r.number * mp.mi0;
            }
        }
    } else if (mp.inuc == 1) {
        if (t < Real(273.15) && sat.qvqvsi > one) {
            const Real kc2 = Real(0.16) * Real(1000.0) / rho; // L^-1 to kg^-1
            if (kc2 > ni + ns + ng) {
                r.number = (kc2 - ni - ns - ng) / dt;
                r.mass = r.number * mp.mi0;
            }
        }
    }
    return r;
}

//! Fraction of rime that produces Hallet-Mossop splinters, peaking at -5C
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_rime_splinter_factor (const amrex::Real t) noexcept
{
    using amrex::Real;
    Real fmult = Real(0);
    if (t > Real(270.16)) {
        fmult = Real(0);
    } else if (t <= Real(270.16) && t > Real(268.16)) {
        fmult = (Real(270.16) - t) / Real(2);
    } else if (t >= Real(265.16) && t <= Real(268.16)) {
        fmult = (t - Real(265.16)) / three;
    } else if (t < Real(265.16)) {
        fmult = Real(0);
    }
    return fmult;
}

/**
 * Rime-splintering (Hallet and Mossop 1974) on snow or graupel of mixing ratio
 * qx, which collects cloud water at psacw (PSACWS or PSACWG) and rain at prac
 * (PRACS or PRACG). The number of splinters follows the mass of rimed water;
 * the splinter mass is taken out of the riming rates, which it cannot exceed.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_rime_splintering (const amrex::Real qx, const amrex::Real qc,
                                const amrex::Real qr, const amrex::Real t,
                                const amrex::Real mmult,
                                amrex::Real& psacw, amrex::Real& prac,
                                amrex::Real& nmult_c, amrex::Real& qmult_c,
                                amrex::Real& nmult_r, amrex::Real& qmult_r) noexcept
{
    using amrex::Real;
    // Thresholds limit rime-splintering in stratiform clouds; they correspond with
    // the graupel thresholds in Rutledge and Hobbs 1984
    if (qx >= Real(0.1e-3)) {
        if (qc >= Real(0.5e-3) || qr >= Real(0.1e-3)) {
            if (psacw > Real(0) || prac > Real(0)) {
                if (t < Real(270.16) && t > Real(265.16)) {
                    const Real fmult = morrison_rime_splinter_factor(t);

                    // 1000 is to convert from kg to g
                    // Splintering from droplets accreted onto snow/graupel
                    if (psacw > Real(0)) {
                        nmult_c = Real(35.0e4) * psacw * fmult * Real(1000.0);
                        qmult_c = nmult_c * mmult;
                        qmult_c = std::min(qmult_c, psacw);
                        psacw = psacw - qmult_c;
                    }

                    // Riming and splintering from accreted raindrops
                    if (prac > Real(0)) {
                        nmult_r = Real(35.0e4) * prac * fmult * Real(1000.0);
                        qmult_r = nmult_r * mmult;
                        qmult_r = std::min(qmult_r, prac);
                        prac = prac - qmult_r;
                    }
                }
            }
        }
    }
}

/**
 * Condensation rate PCC of the saturation adjustment that condenses the vapor
 * above water saturation after the tendencies so far have been applied over dt.
 * It cannot evaporate more cloud water than there is.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real morrison_condensation_rate (const amrex::Real t, const amrex::Real qv,
                                        const amrex::Real qc,
                                        const amrex::Real t_ten, const amrex::Real qv_ten,
                                        const amrex::Real qc_ten,
                                        const amrex::Real pres, const amrex::Real xxlv,
                                        const amrex::Real cpm,
                                        const MorrisonParams& mp, const amrex::Real dt) noexcept
{
    const amrex::Real dumt  = t + dt * t_ten;
    const amrex::Real dumqv = qv + dt * qv_ten;

    // hm, add fix for low pressure, 5/12/10
    const amrex::Real dumqss = morrison_mixing_ratio(mp.ep_2, morrison_svp_liquid(dumt, pres), pres);
    amrex::Real dumqc = qc + dt * qc_ten;
    dumqc = std::max(dumqc, amrex::Real(0));

    const amrex::Real dums = dumqv - dumqss;
    amrex::Real pcc = dums / (one + amrex::Math::powi<2>(xxlv) * dumqss /
                              (cpm * mp.Rv * amrex::Math::powi<2>(dumt))) / dt;
    if (pcc * dt + dumqc < amrex::Real(0)) {
        pcc = -dumqc / dt;
    }
    return pcc;
}

// ---------------------------------------------------------------------------
// Instantaneous adjustments
// ---------------------------------------------------------------------------

/**
 * Move all of species src (mass and number) into species dst. The temperature
 * changes by q_src*latent/cpm: pass latent = XLF for freezing, -XLF for melting.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_convert_all (amrex::Real& q_src, amrex::Real& n_src,
                           amrex::Real& q_dst, amrex::Real& n_dst,
                           amrex::Real& t, const amrex::Real latent,
                           const amrex::Real cpm) noexcept
{
    q_dst = q_dst + q_src;
    n_dst = n_dst + n_src;
    t = t + q_src * latent / cpm;
    q_src = amrex::Real(0);
    n_src = amrex::Real(0);
}

//! Evaporate (or sublimate) a mixing ratio q below 1e-8 into vapor
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_evaporate_small (amrex::Real& q, amrex::Real& qv, amrex::Real& t,
                               const amrex::Real latent, const amrex::Real cpm) noexcept
{
    if (q < amrex::Real(1.0e-8)) {
        qv += q;
        t -= q * latent / cpm;
        q = amrex::Real(0);
    }
}

/**
 * At subsaturation, remove small amounts of cloud and precipitation: liquid below
 * 90% of water saturation, ice below 90% of ice saturation.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_evaporate_small_at_subsaturation (const MorrisonSaturation& sat,
                                                const amrex::Real xxlv, const amrex::Real xxls,
                                                const amrex::Real cpm,
                                                amrex::Real& qv, amrex::Real& t,
                                                amrex::Real& qc, amrex::Real& qr,
                                                amrex::Real& qi, amrex::Real& qni,
                                                amrex::Real& qg) noexcept
{
    // hm modify 7/0/09 change limit to 1.e-8
    if (sat.qvqvs < amrex::Real(0.9)) {
        morrison_evaporate_small(qr, qv, t, xxlv, cpm);
        morrison_evaporate_small(qc, qv, t, xxlv, cpm);
    }
    if (sat.qvqvsi < amrex::Real(0.9)) {
        morrison_evaporate_small(qi,  qv, t, xxls, cpm);
        morrison_evaporate_small(qni, qv, t, xxls, cpm);
        morrison_evaporate_small(qg,  qv, t, xxls, cpm);
    }
}

//! Below QSMALL, set the mixing ratio, number and effective radius to zero
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_zero_below_qsmall (amrex::Real& q, amrex::Real& n, amrex::Real& eff,
                                 const amrex::Real qsmall) noexcept
{
    if (q < qsmall) {
        q   = amrex::Real(0);
        n   = amrex::Real(0);
        eff = amrex::Real(0);
    }
}

#endif
