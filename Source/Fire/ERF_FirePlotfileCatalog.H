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
 * Base output order (19 components, indices 0–18):
 *   0. fire_phi           (level-set)
 *   1. fire_ros           (rate of spread [m/s])
 *   2. fire_wind_eff_u    (component 0 of fire_wind_eff)
 *   3. fire_wind_eff_v    (component 1 of fire_wind_eff)
 *   4. fire_wind_ref_u    (component 0 of fire_wind_ref)
 *   5. fire_wind_ref_v    (component 1 of fire_wind_ref)
 *   6. fire_extract_z      (wind extraction height [m])
 *   7. fire_slope_x       (component 0 of fire_slopes)
 *   8. fire_slope_y       (component 1 of fire_slopes)
 *   9. fire_fuel_mc_1hr   (component 0 of fire_fuel_mc)
 *  10. fire_fuel_mc_10hr  (component 1 of fire_fuel_mc)
 *  11. fire_fuel_mc_100hr (component 2 of fire_fuel_mc)
 *  12. fire_surface_temp_K (surface temperature [K])
 *  13. fire_surface_rh    (surface relative humidity [0-1])
 *  14. fire_heat_flux     (heat flux [W/m²])
 *  15. fire_fuel_load     (fuel load [kg/m²])
 *  16. fire_fireline_intensity (Byram fireline intensity [kW/m])
 *  17. fire_flame_length  (Thomas flame length [m])
 *  18. fire_arrival_time  (fire arrival time [s])
 *
 * Phase 8 spotting diagnostics append 4 components when spotting.enable = true:
 *   fire_lofting_height, fire_spot_brand_count, fire_spot_max_dist, fire_spot_active
 *
 * Phase 9 crown-fire diagnostics append 4 or 5 components when crown fire output
 * is enabled:
 *   fire_crown_active, fire_crown_load, fire_crown_fraction_burned,
 *   [optional] fire_flame_tilt, fire_flame_temp
 */

/**
 * @brief Get the canonical fire plotfile variable names.
 * @param[in] include_spotting   Append 4 Albini spotting components when true.
 * @param[in] include_crown      Append crown-fire diagnostics when true.
 * @param[in] include_flame_tilt Append fire_flame_tilt within the crown block when true.
 * @return Vector of variable name strings in canonical write order.
 */
inline amrex::Vector<std::string> fire_plotfile_var_names(
    bool include_spotting,
    bool include_crown,
    bool include_flame_tilt)
{
    amrex::Vector<std::string> names{
        "fire_phi",
        "fire_ros",
        "fire_wind_eff_u",
        "fire_wind_eff_v",
        "fire_wind_ref_u",
        "fire_wind_ref_v",
        "fire_extract_z",
        "fire_slope_x",
        "fire_slope_y",
        "fire_fuel_mc_1hr",
        "fire_fuel_mc_10hr",
        "fire_fuel_mc_100hr",
        "fire_surface_temp_K",
        "fire_surface_rh",
        "fire_heat_flux",
        "fire_fuel_load",
        "fire_fireline_intensity",
        "fire_flame_length",
        "fire_arrival_time"
    };
    if (include_spotting) {
        names.push_back("fire_lofting_height");
        names.push_back("fire_spot_brand_count");
        names.push_back("fire_spot_max_dist");
        names.push_back("fire_spot_active");
    }
    if (include_crown) {
        names.push_back("fire_crown_active");
        names.push_back("fire_crown_load");
        names.push_back("fire_crown_fraction_burned");
        if (include_flame_tilt) {
            names.push_back("fire_flame_tilt");
        }
        names.push_back("fire_flame_temp");
    }
    return names;
}


inline amrex::Vector<std::string> fire_plotfile_var_names(
    bool include_spotting = false,
    bool include_crown = false)
{
    return fire_plotfile_var_names(include_spotting, include_crown, false);
}

/**
 * @brief Get the total number of fire plotfile components.
 * @param[in] include_spotting   Count 4 additional Albini spotting components when true.
 * @param[in] include_crown      Count crown-fire diagnostics when true.
 * @param[in] include_flame_tilt Count fire_flame_tilt within the crown block when true.
 * @return Number of scalar fire output variables.
 */
inline int fire_plotfile_ncomp(
    bool include_spotting,
    bool include_crown,
    bool include_flame_tilt)
{
    return 19
        + (include_spotting ? 4 : 0)
        + (include_crown ? (include_flame_tilt ? 5 : 4) : 0);
}


inline int fire_plotfile_ncomp(
    bool include_spotting = false,
    bool include_crown = false)
{
    return fire_plotfile_ncomp(include_spotting, include_crown, false);
}

#endif
