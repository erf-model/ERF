#ifndef ERF_NOAHMP_LEVEL_DRIVER_H_
#define ERF_NOAHMP_LEVEL_DRIVER_H_

/**
 * @brief Whether Noah-MP at AMR level lev runs the land model on a setup file of its own.
 *
 * Level 0 always does. A finer level does if it has a land setup file of its own; it
 * otherwise takes its land state from level 0 (NOAHMP::interp_from_lev0). Which file counts:
 *
 * - a run initialized from WRF or metgrid files: erf.nc_init_file_<lev>, which also
 *   supplies the level's atmosphere. ERF_SETUP_FILE_0<lev+1> in namelist.erf alone does
 *   not, so such runs behave as they always have;
 * - any other run (input_sounding, a custom problem), which has no such file:
 *   ERF_SETUP_FILE_0<lev+1> in namelist.erf.
 *
 * @param lev                    AMR level
 * @param has_init_file          erf.nc_init_file_<lev> is given
 * @param idealized_init         the run is not initialized from WRF or metgrid files
 * @param namelist_names_setup   namelist.erf names ERF_SETUP_FILE_0<lev+1>
 */
[[nodiscard]] inline bool
noahmp_level_runs_driver (int lev, bool has_init_file, bool idealized_init,
                          bool namelist_names_setup) noexcept
{
    if (lev == 0 || has_init_file) { return true; }
    return idealized_init && namelist_names_setup;
}

#endif
