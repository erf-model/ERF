#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// cell-centered state variables
#define Density_comp  0
#define RhoTheta_comp 1
#define Scalar_comp   2
#define NVAR          3

#define SlipWall 4
#define NoSlipWall 5

#define NUM_ADV 1

// This value was hardcoded to 4 in get_method_params in the old Fortran code
#define NUM_GROW 4

#endif
