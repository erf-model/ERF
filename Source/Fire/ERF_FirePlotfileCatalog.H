#ifndef ERF_FIRE_PLOTFILE_CATALOG_H
#define ERF_FIRE_PLOTFILE_CATALOG_H

#include <AMReX_Vector.H>
#include <string>

/**
 * @file ERF_FirePlotfileCatalog.H
 * @brief Header-only catalog of fire plotfile variable names
 *
 * Defines the canonical output component names for fire plotfiles.
 * All variables are on the native fire grid (C x finer than atmospheric grid).
 *
 * Output order:
 *   1. fire_phi           (1 component  — level-set)
 *   2. fire_ros           (1 component  — rate of spread [m/s])
 *   3. fire_wind_eff_u    (component 0 of fire_wind_eff)
 *   4. fire_wind_eff_v    (component 1 of fire_wind_eff)
 *   5. fire_wind_ref_u    (component 0 of fire_wind_ref)
 *   6. fire_wind_ref_v    (component 1 of fire_wind_ref)
 *   7. fire_slope_x       (component 0 of fire_slopes)
 *   8. fire_slope_y       (component 1 of fire_slopes)
 *   9. fire_fuel_mc_1hr   (component 0 of fire_fuel_mc)
 *  10. fire_fuel_mc_10hr  (component 1 of fire_fuel_mc)
 *  11. fire_fuel_mc_100hr (component 2 of fire_fuel_mc)
 */

/**
 * @brief Get the canonical fire plotfile variable names.
 * @return Vector of 11 variable name strings in canonical write order.
 */
inline amrex::Vector<std::string> fire_plotfile_var_names()
{
    return amrex::Vector<std::string>{
        "fire_phi",
        "fire_ros",
        "fire_wind_eff_u",
        "fire_wind_eff_v",
        "fire_wind_ref_u",
        "fire_wind_ref_v",
        "fire_slope_x",
        "fire_slope_y",
        "fire_fuel_mc_1hr",
        "fire_fuel_mc_10hr",
        "fire_fuel_mc_100hr"
    };
}

/**
 * @brief Get the total number of fire plotfile components.
 * @return 11 (the number of scalar fire output variables).
 */
inline int fire_plotfile_ncomp()
{
    return 11;
}

#endif
