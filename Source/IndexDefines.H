#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// Cell-centered state variables
#define Rho_comp       0
#define RhoTheta_comp  1
#define RhoKE_comp     2 // for Deardorff LES Model
#define RhoQKE_comp    3 // for MYNN PBL Model
#define RhoScalar_comp 4

#ifdef ERF_USE_MOISTURE
#define RhoQv_comp     5
#define RhoQc_comp     6
#define NVAR           7
#else
#define NVAR           5
#endif

// Cell-centered primitive variables
#define PrimTheta_comp   RhoTheta_comp -1
#define PrimKE_comp      RhoKE_comp    -1
#define PrimQKE_comp     RhoQKE_comp   -1
#define PrimScalar_comp  RhoScalar_comp-1

#ifdef ERF_USE_MOISTURE
#define PrimQv_comp      RhoQv_comp-1
#define PrimQc_comp      RhoQc_comp-1
#endif

#define NUM_PRIM         NVAR-1

namespace BCVars {
    enum {
        cons_bc = 0,
        Rho_bc_comp = 0,
        RhoTheta_bc_comp,
        RhoKE_bc_comp,
        RhoQKE_bc_comp,
        RhoScalar_bc_comp,
#ifdef ERF_USE_MOISTURE
        RhoQv_bc_comp,
        RhoQc_bc_comp,
#endif
        xvel_bc = NVAR,
        yvel_bc,
        zvel_bc,
        NumTypes
    };
}

namespace WRFBdyVars {
    enum {
        U  = 0,
        V  = 1,
        T  = 2,
        QV    , // water vapor
        MU    , // bdy perturbation dry air mass in column (we will get mub from the initial data)
        PC    , // p_s - p_top = dry hydrostatic pressure diference between the surface and the model top
        NumTypes
    };
}

namespace WRFBdyTypes {
    enum {
        x_lo,
        x_hi,
        y_lo,
        y_hi
    };
}
namespace RHSVar {
    enum {
        slow = 0,
        fast,
        all
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
#ifdef ERF_USE_MOISTURE
        RhoQv,
        RhoQc,
#endif
        NumVars
    };
}

namespace Prim {
    enum {
        Theta = 0,
        KE,
        QKE,
        Scalar,
#ifdef ERF_USE_MOISTURE
        Qv,
        Qc,
#endif
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
        PBL_lengthscale,
        NumDiffs
    };
}

enum struct ERF_BC {
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
    MOST         =  101,
    ext_dir_ingested =  102
};
}

namespace TerrainMet {
  enum {
        h_xi=0,
        h_eta,
        h_zeta,
        h_xi_eta,
        h_xi_zeta,
        h_eta_zeta,
        all
  };
}
#endif
