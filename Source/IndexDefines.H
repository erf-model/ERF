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

// cell-centered primitive variables
#define PrimTheta_comp   0
#define PrimKE_comp      1
#define PrimScalar_comp  2
#define NUM_PRIM         2 + NUM_ADV // (assumes NUM_ADV >= 1)

namespace BCVars {
    enum {
        cons_bc = 0,
        Rho_bc_comp = 0,
        RhoTheta_bc_comp,
        RhoKE_bc_comp,
        RhoScalar_bc_comp,
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
        NumTypes
    };
}

namespace Cons {
    enum {
        Rho = 0,
        RhoTheta,
        RhoKE,
        RhoScalar,
        NumVars
    };
}

namespace Prim {
    enum {
        Theta = 0,
        KE,
        Scalar,
        NumVars
    };
}

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
    MOST         =  101
};
}

#endif
