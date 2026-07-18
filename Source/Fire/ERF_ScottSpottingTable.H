#ifndef ERF_SCOTT_SPOTTING_TABLE_H
#define ERF_SCOTT_SPOTTING_TABLE_H

#include <AMReX_REAL.H>
#include <string>

/**
 * @file ERF_ScottSpottingTable.H
 * @brief Maximum firebrand spotting distance table for Anderson FBFM13 fuel models.
 *
 * Provides per-fuel-model maximum spotting distances derived from:
 *   Scott, J.H. (2006). "Comparison of crown fire behavior characteristics between
 *   prescribed fires and wildfires in ponderosa pine forests." USDA Forest Service
 *   Proceedings RMRS-P-41.
 *   Albini, F.A. (1983). "Potential spotting distance from wind-driven surface fires."
 *   USDA Forest Service Research Paper INT-309.
 *
 * Two fuel system tables are available: "13" (Anderson FBFM13) and "40" (Scott-Burgan).
 * All distances are in metres.
 *
 * Ported from: wildfire_levelset/src/scott_spotting_table.H (no logic changes).
 */

/**
 * @brief Maximum spotting distance [m] for Anderson FBFM13 fuel models (fuel_system="13").
 *
 * Indexed 0–12 corresponding to FM1–FM13. Values represent empirical
 * maximum spotting distances for each fuel model under extreme fire conditions.
 *
 * @param[in] model_id  Anderson fuel model ID (1–13)
 * @return Maximum spotting distance [m]; returns 500 m for out-of-range IDs.
 */
inline amrex::Real
scott_max_spot_dist_13(int model_id)
{
    // Maximum spotting distances [m] for FM1–FM13 (index 0 = FM1).
    // Source: Scott (2006) and Albini (1983) INT-309.
    static const amrex::Real table13[13] = {
        200.0,   // FM1  Short Grass
        400.0,   // FM2  Timber Grass and Understory
        300.0,   // FM3  Tall Grass
        500.0,   // FM4  Chaparral
        400.0,   // FM5  Brush
        300.0,   // FM6  Dormant Brush, Hardwood Slash
        400.0,   // FM7  Southern Rough
        200.0,   // FM8  Closed Timber Litter
        300.0,   // FM9  Hardwood Litter
        600.0,   // FM10 Timber Litter and Understory
        700.0,   // FM11 Light Logging Slash
        800.0,   // FM12 Medium Logging Slash
        1000.0   // FM13 Heavy Logging Slash
    };

    if (model_id < 1 || model_id > 13) { return 500.0; }
    return table13[model_id - 1];
}

/**
 * @brief Maximum spotting distance [m] for Scott-Burgan 40 fuel models (fuel_system="40").
 *
 * Indexed 0–39 corresponding to GR1–TU5. Returns 500 m for unrecognised IDs.
 *
 * @param[in] model_id  Scott-Burgan fuel model ID (101–204)
 * @return Maximum spotting distance [m]; returns 500 m for out-of-range IDs.
 */
inline amrex::Real
scott_max_spot_dist_40(int model_id)
{
    // Simplified Scott-Burgan 40 table: grass, shrub, timber-understory, timber-litter groups.
    // GR: 101-109, GS: 121-124, SH: 141-149, TU: 161-165, TL: 181-189, SB: 201-204
    if (model_id >= 101 && model_id <= 109) { return 250.0 + (model_id - 101) * 50.0; }  // GR group
    if (model_id >= 121 && model_id <= 124) { return 400.0 + (model_id - 121) * 50.0; }  // GS group
    if (model_id >= 141 && model_id <= 149) { return 300.0 + (model_id - 141) * 50.0; }  // SH group
    if (model_id >= 161 && model_id <= 165) { return 400.0 + (model_id - 161) * 50.0; }  // TU group
    if (model_id >= 181 && model_id <= 189) { return 500.0 + (model_id - 181) * 50.0; }  // TL group
    if (model_id >= 201 && model_id <= 204) { return 800.0 + (model_id - 201) * 50.0; }  // SB group
    return 500.0;
}

/**
 * @brief Get maximum spotting distance [m] for a given fuel model and fuel system.
 *
 * @param[in] model_id    Fuel model ID (1–13 for FBFM13, 101–204 for Scott-Burgan 40).
 * @param[in] fuel_system "13" for Anderson FBFM13, "40" for Scott-Burgan 40.
 * @return Maximum spotting distance [m].
 */
inline amrex::Real
scott_max_spot_dist(int model_id, const std::string& fuel_system)
{
    if (fuel_system == "40") { return scott_max_spot_dist_40(model_id); }
    return scott_max_spot_dist_13(model_id);
}

#endif
