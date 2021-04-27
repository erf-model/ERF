#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// cell-centered state variables
#define Density_comp  0
#define RhoTheta_comp 1
#define Scalar_comp   2
#define NVAR          3

#define URHO 0
#define UMX 1
#define UMY 2
#define UMZ 3
#define UEDEN 4
#define UEINT 5
#define UTEMP 6
#define UFA 7

#define QRHO 0
#define QU 1
#define QV 2
#define QW 3
#define QGAME 4
#define QPRES 5
#define QREINT 6
#define QTEMP 7
#define QFA 8

#define dComp_rhoD 0
#define dComp_mu dComp_rhoD + 1
#define dComp_xi dComp_mu + 1
#define dComp_lambda dComp_xi + 1

#define comp_Cs2 0
#define comp_CI comp_Cs2 + 1
#define comp_PrT comp_CI + 1
#define comp_Cs2ovPrT comp_PrT
#define nCompC comp_PrT + 1

#define SlipWall 4
#define NoSlipWall 5

#define QTHERM 8
#define NUM_ADV 1

#define QVAR QTHERM + NUM_ADV

#define NQ QVAR
#define NGDNV 6

#define NPASSIVE NUM_ADV 

// This value was hardcoded to 4 in get_method_params in the old Fortran code
#define NUM_GROW 4

namespace indxmap {

extern AMREX_GPU_DEVICE_MANAGED int qpassMap[NPASSIVE];
extern AMREX_GPU_DEVICE_MANAGED int upassMap[NPASSIVE];

void init();

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
qpass_map(const int i)
{
  return qpassMap[i];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
upass_map(const int i)
{
  return upassMap[i];
}
} // namespace indxmap

#endif
