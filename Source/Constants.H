#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#define PI 3.14159265358979323846264338327950288

#define SMALL 1.0e-8
#define SMALL_DENS 1.0e-200
#define SMALL_PRES 1.0e-200
#define SMALL_TEMP 1.0e-200
#define SMALL_E 1.0e-200
#define SMALLU 1.0e-12

//TODO: Change these types of macros to 'const'
#define R_d     287.0 // gas constant for dry air [J/(kg-K)]
#define c_p    1004.5 // specific heat at constant pressure for dry air [J/(kg-K)]
#define p_0    1.0e5 // reference surface pressure [Pa]
#define Gamma  1.4   // C_p/C_v [-]

#define CONST_GRAV 9.81

#define Outflow 2

#endif
