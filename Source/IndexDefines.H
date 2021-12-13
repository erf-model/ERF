#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// cell-centered state variables
//TODO: Change these types of macros to 'const'
#define Rho_comp  0
#define RhoTheta_comp  1
#define RhoKE_comp     2
#define RhoScalar_comp 3
#define NVAR           4

#define NUM_ADV 1

#endif
