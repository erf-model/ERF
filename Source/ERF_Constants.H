#ifndef ERF_CONSTANTS_H_
#define ERF_CONSTANTS_H_

#include <cmath>
#include <limits>
#include <AMReX_REAL.H>

// Dimensionless numeric literals (zero, one, myhalf, ...) and the pure
// mathematical constants (PI, PIoTwo, four_thirds_pi) live here.  This header
// is included so that every existing consumer of ERF_Constants.H keeps seeing
// them; code that needs only numeric literals should include it directly.
#include "ERF_NumericalConstants.H"

// Safe overly large value to fill data so that it is not uninitialized
//      but we can easily identify if the value is used in a computation
// Would like to use the following but std::sqrt will not be allowable
//       in constexpr until C++26
// const amrex::Real bogus_large_value = std::sqrt(std::numeric_limits<amrex::Real>::max());
#ifdef AMREX_USE_FLOAT
constexpr amrex::Real bogus_large_value = amrex::Real(1.e18);
#else
constexpr amrex::Real bogus_large_value = amrex::Real(1.e150);
#endif

// lsm_undefined: sentinel written into the LSM MultiFabs for cells
// the land model did NOT process (e.g. sea-ice / open-water points, where Noah-MP returns the
// -9999 fill value instead of a computed flux). The surface layer detects this
// and falls back to the Monin-Obukhov (MOST) flux for those cells -- analogous
// to how WRF computes the atmospheric surface flux from a skin temperature via
// the surface-layer scheme rather than taking an undefined LSM flux.
constexpr amrex::Real lsm_undefined = bogus_large_value;

// Physical ceiling [mm] on interval precip handed to Noah-MP over one land call.
// 500 mm is far above any real rate, so a larger value signals noisy IC/BC or a
// corrupt accumulator snapshot; the coupling clamps to this and warns (ERF_NOAHMP.cpp).
constexpr amrex::Real lsm_max_precip_interval = amrex::Real(500.0);

// Physical Constants
constexpr amrex::Real R_d          = amrex::Real(287.0);    // dry air constant for dry air [J/(kg-K)]
constexpr amrex::Real R_v          = amrex::Real(461.505);  // water vapor constant for water vapor [J/(kg-K)]
constexpr amrex::Real Cp_d         = amrex::Real(1004.5);   // We have set this so that with qv=0 we get identically gamma = amrex::Real(1.4)
constexpr amrex::Real Cp_v         = amrex::Real(1859.0);
constexpr amrex::Real Cp_l         = amrex::Real(4200.0);

constexpr amrex::Real epsv         = R_v / R_d - one;
constexpr amrex::Real RdoCp        = R_d / Cp_d;
constexpr amrex::Real CpoRd        = Cp_d / R_d;
constexpr amrex::Real RvoRd        = R_v / R_d;
constexpr amrex::Real RdoRv        = R_d / R_v;

// epsilo / omeps are the CAM spellings of RdoRv and its complement; they appear
// in the saturation-vapor-pressure formulas in Utils/ERF_SatMethods.H.
constexpr amrex::Real epsilo       = RdoRv;
constexpr amrex::Real omeps        = one - epsilo;

constexpr amrex::Real L_v          = amrex::Real(2.5e6);    // latent heat of vaporization (J / kg)

constexpr amrex::Real p_0          = amrex::Real(1.0e5);    // reference surface pressure [Pa]
constexpr amrex::Real Gamma        = amrex::Real(1.4);      // c_p / c_v [-]
constexpr amrex::Real KAPPA        = amrex::Real(0.41);     // von Karman constant
constexpr amrex::Real CONST_GRAV   = amrex::Real(9.81);

// Derived Constants
constexpr amrex::Real ip_0         = one/p_0;
constexpr amrex::Real iR_d         = one/R_d;
constexpr amrex::Real iGamma       = one/Gamma;

// Molecular / kinetic-theory constants.  These duplicate, in a different unit
// convention, the s_kb / s_N_av / s_mol_weight_air statics on
// MaterialPropertiesCore; consolidating them is left to a separate change.
// Consumers: Particles/ERF_SuperDropletPCCoalescence.H and
// PhysicsInterfaces/Radiation/ERF_Radiation.cpp.
constexpr amrex::Real boltz    = amrex::Real(1.38065e-23);
constexpr amrex::Real avogadro = amrex::Real(6.02214e26);
constexpr amrex::Real mwdair   = amrex::Real(28.966);

// Constants used solely by the moisture and cloud-physics code -- hydrometeor
// densities, phase-partition temperature thresholds, fall speeds,
// autoconversion and collection coefficients, size-distribution intercepts and
// latent heats -- now live in Source/Microphysics/ERF_MicrophysicsConstants.H.

static constexpr int  ORB_UNDEF_INT  = 2000000000;  // undefined int
#endif
