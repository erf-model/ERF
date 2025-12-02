#ifndef ERF_DAMPING_STRUCT_H_
#define ERF_DAMPING_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

enum struct RayleighDampingType {
    None, SlowExplicit, FastExplicit, FastImplicit
};

/**
 * Container holding damping-related choices
 */

struct DampingChoice {
  public:
    void init_params(std::string pp_prefix)
    {
        amrex::ParmParse pp(pp_prefix);

        static std::string rayleigh_type_string = "SlowExplicit";
        pp.query("rayleigh_damping_type",rayleigh_type_string);

        if (!rayleigh_type_string.compare("SlowExplicit")) {
            rayleigh_damping_type = RayleighDampingType::SlowExplicit;
        } else if (!rayleigh_type_string.compare("FastExplicit")) {
            rayleigh_damping_type = RayleighDampingType::FastExplicit;
        } else if (!rayleigh_type_string.compare("FastImplicit")) {
            rayleigh_damping_type = RayleighDampingType::FastImplicit;
        } else {
            amrex::Error("Don't know this rayleigh_damping_type");
        }

        // Include Rayleigh damping (separate flags for each variable)
        pp.query("rayleigh_damp_U", rayleigh_damp_U);
        pp.query("rayleigh_damp_V", rayleigh_damp_V);
        pp.query("rayleigh_damp_W", rayleigh_damp_W);
        pp.query("rayleigh_damp_T", rayleigh_damp_T);
        pp.query("rayleigh_dampcoef", rayleigh_dampcoef);
        pp.query("rayleigh_zdamp", rayleigh_zdamp);

        if (pp.countval("rayleigh_damp_substep") > 0) {
            amrex::Abort("The input rayleigh_damp_substep is deprecated.  Set rayleigh_damping_type instead.");
        }

        // Include vertical-velocity damping to improve robustness
        pp.query("w_damping", w_damping);
        pp.query("w_damping_cfl", w_damping_cfl);
        pp.query("w_damping_const", w_damping_const);
        pp.query("w_damping_coeff", w_damping_coeff);

        if (w_damping && (w_damping_const < 0 && w_damping_coeff < 0)) {
            amrex::Abort("Need to specify vertical damping coefficient w_damping_const (like WRF) or w_damping_coeff (~ dz/dt**2)");
        }
    }

    void display()
    {
        amrex::Print() << "Rayleigh damping            :";
        if (rayleigh_damp_U || rayleigh_damp_V || rayleigh_damp_W || rayleigh_damp_T) {
            if (rayleigh_damp_U) amrex::Print() << " U";
            if (rayleigh_damp_V) amrex::Print() << " V";
            if (rayleigh_damp_W) amrex::Print() << " W";
            if (rayleigh_damp_T) amrex::Print() << " T";
            amrex::Print() << " (coef=" << rayleigh_dampcoef << " 1/s,"
                           << " depth=" << rayleigh_zdamp << " m)"
                           << std::endl;
        } else {
            amrex::Print() << " None" << std::endl;
        }

        amrex::Print() << "w damping                   : " << w_damping;
        if (w_damping_const > 0) {
            amrex::Print() << " (coeff = " << w_damping_const << ")" << std::endl;
        } else {
            amrex::Print() << " (coeff = " << w_damping_coeff << " * dz/dt**2)" << std::endl;
        }
    }

    bool        rayleigh_damp_U        = false;
    bool        rayleigh_damp_V        = false;
    bool        rayleigh_damp_W        = false;
    bool        rayleigh_damp_T        = false;
    amrex::Real rayleigh_dampcoef      = 0.2; // inverse time scale [1/s]
    amrex::Real rayleigh_zdamp         = 500.0; // damping layer depth [m]
    amrex::Real rayleigh_ztop;

    bool        w_damping              = false;
    amrex::Real w_damping_cfl          = 1.0; // from WRF manual -- WRF source code has default
                                              //   w_crit_cfl=1.2
    amrex::Real w_damping_const        = -1;  // damping coefficient -- to used a constant value
                                              //   (WRF uses 0.3 m/s2), set to a positive value;
                                              //   otherwise, the damping coefficient will be
    amrex::Real w_damping_coeff        = -1;  //   grid-dependent: coeff * dz/dt**2

    // Molecular transport model
    RayleighDampingType rayleigh_damping_type = RayleighDampingType::SlowExplicit;
};
#endif
