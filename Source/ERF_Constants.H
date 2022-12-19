#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <AMReX_REAL.H>

constexpr amrex::Real PI           = 3.14159265358979323846264338327950288;
constexpr amrex::Real PIoTwo       = PI/2.0;

// Physical Constants
constexpr amrex::Real R_d          = 287.0; // gas constant for dry air [J/(kg-K)]
constexpr amrex::Real R_v          = 461.505; // gas constant for water vapor [J/(kg-K)]
constexpr amrex::Real c_p          = 1004.5; // specific heat at constant pressure for dry air [J/(kg-K)]
constexpr amrex::Real cp_v         = 1859.0; // specific heat at constant pressure for water vapor [J/(kg-K)]
constexpr amrex::Real p_0          = 1.0e5; // reference surface pressure [Pa]
constexpr amrex::Real Gamma        = 1.4;   // C_p/C_v [-]
constexpr amrex::Real KAPPA        = 0.41;  // von Karman constant
constexpr amrex::Real CONST_GRAV   = 9.81;

// Derived Constants
constexpr amrex::Real ip_0         = 1./p_0;
constexpr amrex::Real iR_d         = 1./R_d;
constexpr amrex::Real iGamma       = 1./Gamma;
constexpr amrex::Real rdOcp        = R_d/c_p;
constexpr amrex::Real inv_cp       = 1./c_p;
constexpr amrex::Real cv_v         = R_v - cp_v;

constexpr amrex::Real rhor         = 1000.; // Density of water, kg/m3
constexpr amrex::Real rhos         = 100.;  // Density of snow, kg/m3
constexpr amrex::Real rhog         = 400.;  // Density of graupel, kg/m3
constexpr amrex::Real tbgmin       = 253.16;    // Minimum temperature for cloud water., K
constexpr amrex::Real tbgmax       = 273.16;    // Maximum temperature for cloud ice, K
constexpr amrex::Real tprmin       = 268.16;    // Minimum temperature for rain, K
constexpr amrex::Real tprmax       = 283.16;    // Maximum temperature for snow+graupel, K
constexpr amrex::Real tgrmin       = 223.16;    // Minimum temperature for snow, K
constexpr amrex::Real tgrmax       = 283.16;    // Maximum temperature for graupel, K

constexpr amrex::Real a_rain       = 842.;   // Coeff.for rain term vel
constexpr amrex::Real b_rain       = 0.8;    // Fall speed exponent for rain
constexpr amrex::Real a_snow       = 4.84;   // Coeff.for snow term vel
constexpr amrex::Real b_snow       = 0.25;   // Fall speed exponent for snow
constexpr amrex::Real a_grau       = 94.5;   // Lin (1983) (rhog=400)
constexpr amrex::Real b_grau       = 0.5;    // Fall speed exponent for graupel

// Autoconversion
constexpr amrex::Real qcw0         = 1.e-3;  // Threshold for water autoconversion, g/g
constexpr amrex::Real qci0         = 1.e-4;  // Threshold for ice autoconversion, g/g
constexpr amrex::Real alphaelq     = 1.e-3;  // autoconversion of cloud water rate coef
constexpr amrex::Real betaelq      = 1.e-3;  // autoconversion of cloud ice rate coef

constexpr amrex::Real erccoef      = 1.0;    // Rain/Cloud water collection efficiency
constexpr amrex::Real esccoef      = 1.0;    // Snow/Cloud water collection efficiency
constexpr amrex::Real esicoef      = 0.1;    // Snow/cloud ice collection efficiency
constexpr amrex::Real egccoef      = 1.0;    // Graupel/Cloud water collection efficiency
constexpr amrex::Real egicoef      = 0.1;    // Graupel/Cloud ice collection efficiency

constexpr amrex::Real nzeror       = 8.e6;   // Intercept coeff. for rain
constexpr amrex::Real nzeros       = 3.e6;   // Intersept coeff. for snow
constexpr amrex::Real nzerog       = 4.e6;   // Intersept coeff. for graupel
constexpr amrex::Real qp_threshold = 1.e-8;  // minimal rain/snow water content

constexpr amrex::Real boltz  = 1.38065e-23;
constexpr amrex::Real avogad = 6.02214e26;
constexpr amrex::Real mwdair = 28.966;
constexpr amrex::Real mwwv   = 18.016;
constexpr amrex::Real lcond  = 2.501e6;
constexpr amrex::Real lfus   = 2.11727e3;
constexpr amrex::Real lsub   = lcond+lfus;
constexpr amrex::Real rgas   = boltz*avogad/mwdair;
constexpr amrex::Real rv     = rgas/mwwv;

constexpr amrex::Real diffelq = 2.21e-05;     // Diffusivity of water vapor, m2/s
constexpr amrex::Real therco  = 2.40e-02;     // Thermal conductivity of air, J/m/s/K
constexpr amrex::Real muelq   = 1.717e-05;    // Dynamic viscosity of air

constexpr amrex::Real fac_cond = lcond/c_p;
constexpr amrex::Real fac_fus  = lfus/c_p;
constexpr amrex::Real fac_sub  = lsub/c_p;

constexpr amrex::Real a_bg = 1.0/(tbgmax-tbgmin);
constexpr amrex::Real a_pr = 1.0/(tprmax-tprmin);
constexpr amrex::Real a_gr = 1.0/(tgrmax-tgrmin);

constexpr amrex::Real crain = b_rain / 4.0;
constexpr amrex::Real csnow = b_snow / 4.0;
constexpr amrex::Real cgrau = b_grau / 4.0;

#endif
