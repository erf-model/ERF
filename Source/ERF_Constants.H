#ifndef ERF_CONSTANTS_H_
#define ERF_CONSTANTS_H_

#include <AMReX_REAL.H>

constexpr amrex::Real zero         = amrex::Real(0.0);
constexpr amrex::Real one          = amrex::Real(1.0);
constexpr amrex::Real two          = amrex::Real(2.0);
constexpr amrex::Real three        = amrex::Real(3.0);
constexpr amrex::Real four         = amrex::Real(4.0);
constexpr amrex::Real myhalf         = amrex::Real(0.5);
constexpr amrex::Real fourth       = amrex::Real(0.25);
constexpr amrex::Real third        = amrex::Real(1.0/3.0);
constexpr amrex::Real eightthirds  = amrex::Real(8.0/3.0);

constexpr amrex::Real PI           = amrex::Real(3.14159265358979323846264338327950288);
constexpr amrex::Real PIoTwo       = PI/two;

// Physical Constants
constexpr amrex::Real R_d          = amrex::Real(287.0);    // dry air constant for dry air [J/(kg-K)]
constexpr amrex::Real R_v          = amrex::Real(461.505);  // water vapor constant for water vapor [J/(kg-K)]
constexpr amrex::Real Cp_d         = amrex::Real(1004.5);   // We have set this so that with qv=0 we get identically gamma = amrex::Real(1.4)
constexpr amrex::Real Cp_v         = amrex::Real(1859.0);
constexpr amrex::Real Cp_l         = amrex::Real(4200.0);

constexpr amrex::Real L_v          = amrex::Real(2.5e6);    // latent heat of vaporization (J / kg)

constexpr amrex::Real p_0          = amrex::Real(1.0e5);    // reference surface pressure [Pa]
constexpr amrex::Real Gamma        = amrex::Real(1.4);      // c_p / c_v [-]
constexpr amrex::Real KAPPA        = amrex::Real(0.41);     // von Karman constant
constexpr amrex::Real CONST_GRAV   = amrex::Real(9.81);

// Derived Constants
constexpr amrex::Real ip_0         = one/p_0;
constexpr amrex::Real iR_d         = one/R_d;
constexpr amrex::Real iGamma       = one/Gamma;

constexpr amrex::Real rhor         = amrex::Real(1000.); // Density of water, kg/m3
constexpr amrex::Real rhos         = amrex::Real(100.);  // Density of snow, kg/m3
constexpr amrex::Real rhog         = amrex::Real(400.);  // Density of graupel, kg/m3
constexpr amrex::Real tbgmin       = amrex::Real(253.16);    // Minimum temperature for cloud water., K
constexpr amrex::Real tbgmax       = amrex::Real(273.16);    // Maximum temperature for cloud ice, K
constexpr amrex::Real tprmin       = amrex::Real(268.16);    // Minimum temperature for rain, K
constexpr amrex::Real tprmax       = amrex::Real(283.16);    // Maximum temperature for snow+graupel, K
constexpr amrex::Real tgrmin       = amrex::Real(223.16);    // Minimum temperature for graupel, K
constexpr amrex::Real tgrmax       = amrex::Real(283.16);    // Maximum temperature for graupel, K

constexpr amrex::Real a_rain       = amrex::Real(842.);   // Coeff.for rain term vel
constexpr amrex::Real b_rain       = amrex::Real(0.8);    // Fall speed exponent for rain
constexpr amrex::Real a_snow       = amrex::Real(4.84);   // Coeff.for snow term vel
constexpr amrex::Real b_snow       = fourth;   // Fall speed exponent for snow
constexpr amrex::Real a_grau       = amrex::Real(94.5);   // Lin (1983) (rhog=400)
constexpr amrex::Real b_grau       = amrex::Real(0.5);    // Fall speed exponent for graupel

// Autoconversion
constexpr amrex::Real qcw0         = amrex::Real(1.e-3);  // Threshold for water autoconversion, g/g
constexpr amrex::Real qci0         = amrex::Real(1.e-4);  // Threshold for ice autoconversion, g/g
constexpr amrex::Real alphaelq     = amrex::Real(1.e-3);  // autoconversion of cloud water rate coef
constexpr amrex::Real betaelq      = amrex::Real(1.e-3);  // autoconversion of cloud ice rate coef

constexpr amrex::Real erccoef      = one;    // Rain/Cloud water collection efficiency
constexpr amrex::Real esccoef      = one;    // Snow/Cloud water collection efficiency
constexpr amrex::Real esicoef      = amrex::Real(0.1);    // Snow/cloud ice collection efficiency
constexpr amrex::Real egccoef      = one;    // Graupel/Cloud water collection efficiency
constexpr amrex::Real egicoef      = amrex::Real(0.1);    // Graupel/Cloud ice collection efficiency

constexpr amrex::Real nzeror       = amrex::Real(8.e6);   // Intercept coeff. for rain
constexpr amrex::Real nzeros       = amrex::Real(3.e6);   // Intersept coeff. for snow
constexpr amrex::Real nzerog       = amrex::Real(4.e6);   // Intersept coeff. for graupel
constexpr amrex::Real qp_threshold = amrex::Real(1.e-8);  // minimal rain/snow water content

constexpr amrex::Real boltz    = amrex::Real(1.38065e-23);
constexpr amrex::Real avogadro = amrex::Real(6.02214e26);
constexpr amrex::Real mwdair   = amrex::Real(28.966);
constexpr amrex::Real mwwv     = amrex::Real(18.016);
constexpr amrex::Real lcond    = amrex::Real(2.501e6);
constexpr amrex::Real lfus     = amrex::Real(3.337e5);
constexpr amrex::Real lsub     = lcond+lfus;
constexpr amrex::Real rair     = boltz*avogadro/mwdair;
constexpr amrex::Real rh20     = rair/mwwv;
constexpr amrex::Real rga      = one/CONST_GRAV;

constexpr amrex::Real diffelq = amrex::Real(2.21e-05);     // Diffusivity of water vapor, m2/s
constexpr amrex::Real therco  = amrex::Real(2.40e-02);     // Thermal conductivity of air, J/m/s/K
constexpr amrex::Real muelq   = amrex::Real(1.717e-05);    // Dynamic viscosity of air

constexpr amrex::Real a_bg = one/(tbgmax-tbgmin);
constexpr amrex::Real a_pr = one/(tprmax-tprmin);
constexpr amrex::Real a_gr = one/(tgrmax-tgrmin);

constexpr amrex::Real crain = b_rain / amrex::Real(4.0);
constexpr amrex::Real csnow = b_snow / amrex::Real(4.0);
constexpr amrex::Real cgrau = b_grau / amrex::Real(4.0);

constexpr amrex::Real lat_vap  = amrex::Real(2.5e6);     // Latent heat of vaporization (J/kg)
constexpr amrex::Real lat_ice  = amrex::Real(3.337e5);   // latent heat of fusion (J/kg)
constexpr amrex::Real Rd_on_Rv = R_d/R_v;
constexpr amrex::Real tmelt    = amrex::Real(273.15);  // melting temp.
constexpr amrex::Real h2otrip  = amrex::Real(273.15);  // Triple point temperature of water (K)
constexpr amrex::Real tboil    = amrex::Real(373.16);  // Boiling point of water at 1 atm (K)
constexpr amrex::Real ttrice   = amrex::Real(20.00);   // transition range from es over H2O to es over ice
constexpr amrex::Real epsilo   = Rd_on_Rv;
constexpr amrex::Real omeps    = one - epsilo;
constexpr amrex::Real rhoh2o   = amrex::Real(1.000e3);  // density of liquid water

static constexpr int  ORB_UNDEF_INT  = 2000000000;  // undefined int
#endif
