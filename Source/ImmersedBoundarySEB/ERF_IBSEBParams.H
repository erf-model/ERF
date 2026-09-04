#ifndef ERF_IBSEB_PARAMS_H
#define ERF_IBSEB_PARAMS_H

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX.H>
#include <string>

/**
 * \file ERF_IBSEBParams.H
 * \brief Inputs of the immersed-boundary surface energy balance (``erf.ibseb.*``).
 *
 * The balance runs on the faces of buildings represented by immersed forcing
 * (``erf.buildings_type = ImmersedForcing``). Every input lives under the
 * ``erf.ibseb`` prefix and everything is off unless ``erf.ibseb.enable`` is set,
 * in which case ERF::init_ibseb() builds one IBFaceSet per level.
 *
 * The parameters are parsed once, in ERF::ReadParameters(), and copied into
 * each IBFaceSet, so a face set never reads ParmParse itself. New inputs of
 * later phases go here, with their defaults, so that the whole input surface
 * of the balance is in one place. The phases and their inputs are listed in
 * Source/ImmersedBoundarySEB/IBSEB_DEVELOPMENT.md.
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

    /**
     * Read the ``erf.ibseb.*`` inputs and validate them.
     *
     * Only queried, never required: a deck without any ``erf.ibseb`` key gets
     * the defaults above with the balance off. Validation is skipped when the
     * balance is off so an unrelated deck cannot abort on an unused key.
     */
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
