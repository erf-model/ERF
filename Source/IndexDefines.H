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
#if defined(ERF_USE_MOISTURE)
  #define RhoQt_comp   5
  #define RhoQp_comp   6
  #define NVAR         7
#elif defined(ERF_USE_WARM_NO_PRECIP)
  #define RhoQv_comp   5
  #define RhoQc_comp   6
  #define NVAR         7
#else
  #define NVAR         5
#endif

// Cell-centered primitive variables
#define PrimTheta_comp   RhoTheta_comp -1
#define PrimKE_comp      RhoKE_comp    -1
#define PrimQKE_comp     RhoQKE_comp   -1
#define PrimScalar_comp  RhoScalar_comp-1
#if defined(ERF_USE_MOISTURE)
  #define PrimQt_comp    RhoQt_comp-1
  #define PrimQp_comp    RhoQp_comp-1
#elif defined(ERF_USE_WARM_NO_PRECIP)
  #define PrimQv_comp    RhoQv_comp-1
  #define PrimQc_comp    RhoQc_comp-1
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
#if defined(ERF_USE_MOISTURE)
        RhoQt_bc_comp,
        RhoQp_bc_comp,
#elif defined(ERF_USE_WARM_NO_PRECIP)
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
        R  = 2,
        T  = 3,
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
#if defined(ERF_USE_MOISTURE)
        RhoQt,
        RhoQp,
#elif defined(ERF_USE_WARM_NO_PRECIP)
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
#if defined(ERF_USE_MOISTURE)
        Qt,
        Qp,
#elif defined(ERF_USE_WARM_NO_PRECIP)
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
        Theta_h,
        Scalar_h,
        KE_h,
        QKE_h,
#if defined(ERF_USE_MOISTURE)
        Qt_h,
        Qp_h,
#elif defined(ERF_USE_WARM_NO_PRECIP)
        Qv_h,
        Qc_h,
#endif
        Mom_v,
        Theta_v,
        Scalar_v,
        KE_v,
        QKE_v,
#if defined(ERF_USE_MOISTURE)
        Qt_v,
        Qp_v,
#elif defined(ERF_USE_WARM_NO_PRECIP)
        Qv_v,
        Qc_v,
#endif
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
    ext_dir_ingested =  102,
    neumann          =  103,
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

namespace IntVar {
    enum {
        cons = 0,
        xmom,
        ymom,
        zmom,
        NumVars
    };
}

namespace MicVar {
   enum {
      // independent variables
      qt = 0,
      qp,
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      rho,   // density
      pres,  // pressure
      // derived variables
      qr,   // rain water
      qv,   // water vapor
      qn,   // cloud condensate (liquid+ice), initial to zero
      qci,  // cloud ice
      qcl,  // cloud water
      qpl,  // precip rain
      qpi,  // precip ice
      // temporary variable
      omega,
      NumVars
  };
}
#endif
