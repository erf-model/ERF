#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// cell-centered state variables

const int Rho_comp = 0;
const int RhoTheta_comp = 1;
const int RhoScalar_comp = 2;
const int NVAR = 3;

//TODO: Change these types of macros to 'const'
#define SlipWall 4
#define NoSlipWall 5
#define SimSlipWall 6

const int NUM_ADV = 1;

// This value was hardcoded to 4 in get_method_params in the old Fortran code
const int NUM_GROW = 4;

#endif
