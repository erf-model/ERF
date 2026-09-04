#ifndef ERF_IBSEB_PARAMS_H
#define ERF_IBSEB_PARAMS_H

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX.H>
#include <string>

/**
 * @file ERF_IBSEBParams.H
 * @brief Inputs of the immersed-boundary surface energy balance (erf.ibseb.*)
 *
 * The balance runs on the faces of resolved buildings (immersed forcing).
 * Everything is off unless erf.ibseb.enable is set; see
 * Source/ImmersedBoundarySEB/IBSEB_DEVELOPMENT.md for the phases.
 */
struct IBSEBParams
{
    bool        enable        = false;                 ///< Master switch
    int         n_slab_layers = 4;                     ///< Conduction layers per face
    amrex::Real T_skin_init   = 300.0;                 ///< Initial skin temperature [K]
    amrex::Real T_interior    = 293.0;                 ///< Interior (deep) temperature [K]
    std::string csv_file      = "ibseb_buildings.csv"; ///< Per-building report
    int         csv_int       = 100;                   ///< Steps between report rows; <= 0 disables
    bool        debug         = false;                 ///< Extra prints

    void init_params ()
    {
        amrex::ParmParse pp("erf.ibseb");
        pp.query("enable",        enable);
        pp.query("n_slab_layers", n_slab_layers);
        pp.query("T_skin_init",   T_skin_init);
        pp.query("T_interior",    T_interior);
        pp.query("csv_file",      csv_file);
        pp.query("csv_int",       csv_int);
        pp.query("debug",         debug);
        if (enable) {
            if (n_slab_layers < 1) { amrex::Abort("erf.ibseb.n_slab_layers must be >= 1"); }
            if (T_skin_init <= 0.0 || T_interior <= 0.0) {
                amrex::Abort("erf.ibseb.T_skin_init and T_interior must be positive [K]");
            }
        }
    }
};

#endif
