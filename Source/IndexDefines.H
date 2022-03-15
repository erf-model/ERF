#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// Cell-centered state variables
#define Rho_comp  0
#define RhoTheta_comp  1
#define RhoKE_comp     2 // for Deardorff LES Model
#define RhoQKE_comp    3 // for MYNN PBL Model
#define RhoScalar_comp 4
#define RhoQv_comp     5
#define RhoQc_comp     6
#define NVAR           7

// Cell-centered primitive variables
#define PrimTheta_comp   RhoTheta_comp -1
#define PrimKE_comp      RhoKE_comp    -1
#define PrimQKE_comp     RhoQKE_comp   -1
#define PrimScalar_comp  RhoScalar_comp-1
#define PrimQv_comp      RhoQv_comp-1
#define PrimQc_comp      RhoQc_comp-1
#define NUM_PRIM         NVAR-1

namespace BCVars {
    enum {
        cons_bc = 0,
        Rho_bc_comp = 0,
        RhoTheta_bc_comp,
        RhoKE_bc_comp,
        RhoQKE_bc_comp,
        RhoScalar_bc_comp,
        RhoQv_bc_comp,
        RhoQc_bc_comp,
        xvel_bc = NVAR,
        yvel_bc,
        zvel_bc,
        NumTypes
    };
}

namespace Vars {
    enum {
        cons = 0,
        xvel,
        yvel,
        zvel,
        NumTypes,
        xmom,
        ymom,
        zmom
    };
}

namespace Cons {
    enum {
        Rho = 0,
        RhoTheta,
        RhoKE,
        RhoQKE,
        RhoScalar,
        RhoQv,
        RhoQc,
        NumVars
    };
}

namespace Prim {
    enum {
        Theta = 0,
        KE,
        QKE,
        Scalar,
        Qv,
        Qc,
        NumVars
    };
}

// We separate out horizontal and vertical turbulent diffusivities
// These are the same for LES, but different for PBL models
namespace EddyDiff {
    enum {
        Mom_h = 0,
        Mom_v,
        Theta_h,
        Theta_v,
	KE_h,
	KE_v,
        QKE_h,
        QKE_v,
        Scalar_h,
        Scalar_v,
	NumDiffs
    };
}

enum struct BC {
    symmetry, inflow, outflow, no_slip_wall, slip_wall, periodic, MOST, undefined
};

// NOTE: the first of these must match up with the BCType enum
//     in amrex/Src/Base/AMReX_BC_TYPES.H.  We had extras at
//     the end to use locally
namespace ERFBCType {
enum mathematicalBndryTypes : int { bogus        = -666,
    reflect_odd  = -1,
    int_dir      =  0,
    reflect_even =  1,
    foextrap     =  2,
    ext_dir      =  3,
    hoextrap     =  4,
    hoextrapcc   =  5,
    MOST         =  101,
    ext_dir_ingested =  102
};
}

#endif
