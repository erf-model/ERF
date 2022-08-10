#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <AMReX_REAL.H>

constexpr amrex::Real PI = 3.14159265358979323846264338327950288;

//TODO: Change these types of macros to 'const'
#define R_d      287.0   // gas constant for dry air [J/(kg-K)]
#define R_v      461.6   // gas constant for dry air [J/(kg-K)]
#define c_p     1004.5   // specific heat at constant pressure for dry air [J/(kg-K)]
#define p_0        1.0e5 // reference surface pressure [Pa]
#define Gamma      1.4   // C_p/C_v [-]
#define KAPPA      0.41  // von Karman constant
#define CONST_GRAV 9.81

// Derived Constants
#define ip_0    1./p_0
#define iR_d    1./R_d
#define iGamma  1./Gamma
#define rdOcp   R_d/c_p

#endif
