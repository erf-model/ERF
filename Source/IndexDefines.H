#ifndef INDEX_DEFINES_H_
#define INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

/**
 * Definition of indexing parameters
*/

// Cell-centered state variables
#define Rho_comp       0
#define RhoTheta_comp  1
#define RhoKE_comp     2 // for Deardorff LES Model
#define RhoQKE_comp    3 // for MYNN or YSU PBL Model
#define RhoScalar_comp 4
#define RhoQ1_comp     5
#define RhoQ2_comp     6
#define RhoQ3_comp     7
#define RhoQ4_comp     8
#define RhoQ5_comp     9
#define RhoQ6_comp     10

// This is the number of components if using moisture
// We use this to allocate the 1d arrays of boundary condition types,
//    but not to allocate actual solution data
#define NVAR_max 11

// This defines the maximum number of moisture vars
#define NMOIST_max 6

// Cell-centered primitive variables
#define PrimTheta_comp   (RhoTheta_comp -1)
#define PrimKE_comp      (RhoKE_comp    -1)
#define PrimQKE_comp     (RhoQKE_comp   -1)
#define PrimScalar_comp  (RhoScalar_comp-1)
#define PrimQ1_comp      (RhoQ1_comp-1)
#define PrimQ2_comp      (RhoQ2_comp-1)
#define PrimQ3_comp      (RhoQ3_comp-1)
#define PrimQ4_comp      (RhoQ4_comp-1)
#define PrimQ5_comp      (RhoQ5_comp-1)
#define PrimQ6_comp      (RhoQ6_comp-1)

// NOTE: we still use this indexing even if no moisture
namespace BCVars {
    enum {
        cons_bc = 0,
        Rho_bc_comp = 0,
        RhoTheta_bc_comp,
        RhoKE_bc_comp,
        RhoQKE_bc_comp,
        RhoScalar_bc_comp,
        RhoQ1_bc_comp,
        RhoQ2_bc_comp,
        xvel_bc = NVAR_max,
        yvel_bc = NVAR_max+1,
        zvel_bc = NVAR_max+2,
        NumTypes
    };
}

namespace RealBdyVars {
    enum {
        U  = 0,
        V  = 1,
        R  = 2,
        T  = 3,
        QV,
        NumTypes
    };
}

namespace WRFBdyVars {
    enum {
        U  = 0,
        V  = 1,
        R  = 2,
        T  = 3,
        QV    , // water vapor
        MU    , // bdy perturbation dry air mass in column (we will get mub from the initial data)
        PC    , // p_s - p_top = dry hydrostatic pressure diference between the surface and the model top
        NumTypes
    };
}

namespace MetGridBdyVars {
    enum {
        U = 0,
        V = 1,
        R = 2,
        T = 3,
        QV,
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

namespace IntVars {
    enum {
        cons = 0,
        xmom,
        ymom,
        zmom,
        NumTypes
    };
}

// We separate out horizontal and vertical turbulent diffusivities
// These are the same for LES, but different for PBL models
namespace EddyDiff {
    enum {
        Mom_h = 0,
        Theta_h,
        KE_h,
        QKE_h,
        Scalar_h,
        Q_h,
        Mom_v,
        Theta_v,
        KE_v,
        QKE_v,
        Scalar_v,
        Q_v,
        PBL_lengthscale,
        NumDiffs
    };
}

enum struct ERF_BC {
    symmetry, inflow, outflow, ho_outflow, open, no_slip_wall, slip_wall, periodic, MOST, undefined
};

// NOTE: the first of these must match up with the BCType enum
//     in amrex/Src/Base/AMReX_BC_TYPES.H.  We had extras at
//     the end to use locally
namespace ERFBCType {
enum mathematicalBndryTypes : int {
    bogus        = -666,
    reflect_odd  = -1,
    int_dir      =  0,
    reflect_even =  1,
    foextrap     =  2,
    ext_dir      =  3,
    MOST         =  101,
    ext_dir_ingested =  102,
    neumann          =  103,
    neumann_int      =  104,
    open             =  105,
    hoextrapcc         =  106
};
}

enum struct AdvType : int {
    Centered_2nd = 101,
    Upwind_3rd   = 102,
    Centered_4th = 103,
    Upwind_5th   = 104,
    Centered_6th = 105,
    Weno_3        = 106,
    Weno_3Z       = 107,
    Weno_5        = 108,
    Weno_5Z       = 109,
    Weno_3MZQ     = 110,
    Unknown       = 111
};
#endif
