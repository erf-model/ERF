/**
 * @file ERF_DustPlotfileCatalog.H
 * @brief Variable name list for dust 2D plotfile output.
 *
 * Pattern follows ERF_FirePlotfileCatalog.H in Source/Fire/:
 *   https://github.com/hgopalan/ERF/blob/ERF-Fire/Source/Fire/ERF_FirePlotfileCatalog.H
 */
#ifdef ERF_USE_DUST
#include <AMReX_Vector.H>
#include <string>

inline amrex::Vector<std::string> dust_plotfile_var_names()
{
    return {
        "dust_emission_flux",  ///< Total emission flux [kg/m^2/s]. Marticorena & Bergametti (1995).
        "dust_ustar_in",       ///< Friction velocity from SurfaceLayer [m/s].
        "dust_ustar_t",        ///< Threshold friction velocity [m/s]. Shao & Lu (2000).
        "dust_deposition_rate",///< Accumulated dry deposition [kg/m^2]. Zhang et al. (2001).
        "dust_conc_sfc",       ///< Near-surface dust concentration [kg/m^3].
        "dust_surf_moist",     ///< Surface moisture flux [kg/m^2/s]. Null when no moisture scheme.
        "dust_suppression",    ///< Suppression agent coverage [0-1].
        "dust_retreat_flag"    ///< Re-treatment flag [0/1].
    };
}

inline int dust_plotfile_ncomp() { return 8; }
#endif
