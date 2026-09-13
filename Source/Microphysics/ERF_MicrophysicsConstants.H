/*! @file ERF_MicrophysicsConstants.H
 *  \brief Physical constants and tuning parameters used only by the moisture
 *         and cloud-physics code.
 *
 *  These constants were previously carried in Source/ERF_Constants.H, which is
 *  included (directly or transitively) by essentially every translation unit in
 *  ERF.  None of them is referenced by the dycore, advection, diffusion, I/O or
 *  boundary-condition code, so they live here instead: the microphysics schemes,
 *  the superdroplet particle code, the material-property database, and the
 *  moisture utilities are the only consumers.
 *
 *  Source/Microphysics is unconditionally on the include path for both the GNU
 *  Make build (Exec/Make.ERF) and the CMake build (CMake/BuildERFExe.cmake), so
 *  any file in ERF may include this header.
 *
 *  Note on naming: several latent heats appear here under names inherited from
 *  the schemes they were ported from.  `lcond` (2.501e6) comes with the SAM
 *  lineage, while `lat_vap` (2.5e6) comes with the CAM/WRF lineage used by
 *  WSM6/WDM6 and is numerically identical to `L_v` in ERF_Constants.H.
 *  Likewise `lat_ice` and `lfus` are the same number under two names.  They are
 *  kept distinct for now so that this reorganization does not change any
 *  answers; consolidating them is a separate, answer-affecting change.
 */

#ifndef ERF_MICROPHYSICS_CONSTANTS_H_
#define ERF_MICROPHYSICS_CONSTANTS_H_

#include <AMReX_REAL.H>

// This header deliberately depends only on the dimensionless numeric literals,
// not on ERF_Constants.H.  That keeps the layering one-directional: numeric ->
// physical -> scheme-specific.  A file that needs both the thermodynamic
// constants and these should include ERF_Constants.H itself.
#include "ERF_NumericalConstants.H"  // one, fourth

//
// Densities of the condensed phases
//
constexpr amrex::Real rhor         = amrex::Real(1000.); // Density of water, kg/m3
constexpr amrex::Real rhos         = amrex::Real(100.);  // Density of snow, kg/m3
constexpr amrex::Real rhog         = amrex::Real(400.);  // Density of graupel, kg/m3
constexpr amrex::Real rhoh2o       = amrex::Real(1.000e3);  // density of liquid water

//
// Temperature thresholds that partition condensate among the hydrometeor
// species, and the inverse widths of the corresponding transition ranges
//
constexpr amrex::Real tbgmin       = amrex::Real(253.16);    // Minimum temperature for cloud water., K
constexpr amrex::Real tbgmax       = amrex::Real(273.16);    // Maximum temperature for cloud ice, K
constexpr amrex::Real tprmin       = amrex::Real(268.16);    // Minimum temperature for rain, K
constexpr amrex::Real tprmax       = amrex::Real(283.16);    // Maximum temperature for snow+graupel, K
constexpr amrex::Real tgrmin       = amrex::Real(223.16);    // Minimum temperature for graupel, K
constexpr amrex::Real tgrmax       = amrex::Real(283.16);    // Maximum temperature for graupel, K

constexpr amrex::Real a_bg = one/(tbgmax-tbgmin);
constexpr amrex::Real a_pr = one/(tprmax-tprmin);
constexpr amrex::Real a_gr = one/(tgrmax-tgrmin);

//
// Terminal fall speed coefficients and exponents
//
constexpr amrex::Real a_rain       = amrex::Real(842.);   // Coeff.for rain term vel
constexpr amrex::Real b_rain       = amrex::Real(0.8);    // Fall speed exponent for rain
constexpr amrex::Real a_snow       = amrex::Real(4.84);   // Coeff.for snow term vel
constexpr amrex::Real b_snow       = fourth;   // Fall speed exponent for snow
constexpr amrex::Real a_grau       = amrex::Real(94.5);   // Lin (1983) (rhog=400)
constexpr amrex::Real b_grau       = amrex::Real(0.5);    // Fall speed exponent for graupel

constexpr amrex::Real crain = b_rain / amrex::Real(4.0);
constexpr amrex::Real csnow = b_snow / amrex::Real(4.0);
constexpr amrex::Real cgrau = b_grau / amrex::Real(4.0);

//
// Autoconversion
//
constexpr amrex::Real qcw0         = amrex::Real(1.e-3);  // Threshold for water autoconversion, g/g
constexpr amrex::Real qci0         = amrex::Real(1.e-4);  // Threshold for ice autoconversion, g/g
constexpr amrex::Real alphaelq     = amrex::Real(1.e-3);  // autoconversion of cloud water rate coef
constexpr amrex::Real betaelq      = amrex::Real(1.e-3);  // autoconversion of cloud ice rate coef

//
// Collection efficiencies
//
constexpr amrex::Real erccoef      = one;    // Rain/Cloud water collection efficiency
constexpr amrex::Real esccoef      = one;    // Snow/Cloud water collection efficiency
constexpr amrex::Real esicoef      = amrex::Real(0.1);    // Snow/cloud ice collection efficiency
constexpr amrex::Real egccoef      = one;    // Graupel/Cloud water collection efficiency
constexpr amrex::Real egicoef      = amrex::Real(0.1);    // Graupel/Cloud ice collection efficiency

//
// Marshall-Palmer size distribution intercepts, and the precipitate threshold
//
constexpr amrex::Real nzeror       = amrex::Real(8.e6);   // Intercept coeff. for rain
constexpr amrex::Real nzeros       = amrex::Real(3.e6);   // Intersept coeff. for snow
constexpr amrex::Real nzerog       = amrex::Real(4.e6);   // Intersept coeff. for graupel
constexpr amrex::Real qp_threshold = amrex::Real(1.e-8);  // minimal rain/snow water content

//
// Molecular transport coefficients entering the diffusional growth and
// ventilation terms
//
constexpr amrex::Real diffelq = amrex::Real(2.21e-05);     // Diffusivity of water vapor, m2/s
constexpr amrex::Real therco  = amrex::Real(2.40e-02);     // Thermal conductivity of air, J/m/s/K
constexpr amrex::Real muelq   = amrex::Real(1.717e-05);    // Dynamic viscosity of air

//
// Latent heats.  See the note in the file header on the duplicate names.
//
constexpr amrex::Real lcond    = amrex::Real(2.501e6);
constexpr amrex::Real lfus     = amrex::Real(3.337e5);
constexpr amrex::Real lsub     = lcond+lfus;

constexpr amrex::Real lat_vap  = amrex::Real(2.5e6);     // Latent heat of vaporization (J/kg)
constexpr amrex::Real lat_ice  = amrex::Real(3.337e5);   // latent heat of fusion (J/kg)

//
// Water phase-change reference temperatures (CAM lineage, used by
// ERF_SatMethods.H).  The companion parameters epsilo and omeps are plain
// R_d/R_v ratios rather than microphysics tuning values, so they live in
// ERF_Constants.H alongside RdoRv.
//
constexpr amrex::Real tmelt    = amrex::Real(273.15);  // melting temp.
constexpr amrex::Real h2otrip  = amrex::Real(273.15);  // Triple point temperature of water (K)
constexpr amrex::Real tboil    = amrex::Real(373.16);  // Boiling point of water at 1 atm (K)
constexpr amrex::Real ttrice   = amrex::Real(20.00);   // transition range from es over H2O to es over ice

#endif
