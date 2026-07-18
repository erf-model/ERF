#ifndef ERF_FUEL_WIND_HEIGHT_H
#define ERF_FUEL_WIND_HEIGHT_H

#include <vector>
#include <AMReX_Real.H>

/**
 * @file ERF_FuelWindHeight.H
 * @brief Per-fuel wind height (fcwh) and roughness length (fcz0) tables
 *
 * Per-fuel wind height (fcwh) and roughness length (fcz0) for the Anderson
 * FBFM13 fuel models.
 *
 * fcwh[k]: height [m] above terrain at which wind is extracted for fuel model k
 *          (equivalent to WRF-SFIRE fcwh per fuel category).
 *          Index 0 unused; valid range 1-13.
 * fcz0[k]: surface roughness length [m] for log-profile WAF interpolation.
 *
 * Default values follow WRF-SFIRE namelist.fire defaults:
 *   fcwh: all 6.096 m (20 ft) — consistent with BEHAVE midflame height default
 *   fcz0: from WRF-SFIRE data statement:
 *         0.0396, 0.0396, 0.100, 0.2378, 0.0793, 0.0991, 0.0991,
 *         0.0079, 0.0079, 0.0396, 0.0396, 0.0911, 0.1188 for FM1-FM13
 *
 * References:
 *   Andrews, P.L. (1986). BEHAVE: Fire Behavior Prediction and Fuel Modeling
 *     System. USDA Forest Service GTR INT-194.
 *   WRF-SFIRE module_fr_sfire_phys.F, fcwh and fcz0 data statements.
 */

/**
 * @brief Build per-fuel wind height (fcwh) table indexed 0..13
 *
 * @param[in] global_z_ref    Global fallback wind reference height [m]
 * @param[in] use_per_fuel    When true, use WRF-SFIRE fcwh defaults per fuel;
 *                            when false, all entries equal global_z_ref
 * @return host-side vector of size 14; index 0 unused, 1-13 valid
 *
 * Returns vector indexed [0..13]. Index 0 is unused (kept for array alignment).
 * When use_per_fuel=false, all entries 1-13 are set to global_z_ref
 * (backward-compatible, zero cost).
 * When use_per_fuel=true, entries 1-13 are set to WRF-SFIRE defaults (all 6.096 m).
 */
inline std::vector<amrex::Real> build_fcwh_table(
    amrex::Real global_z_ref,
    bool use_per_fuel = false)
{
    std::vector<amrex::Real> fcwh(14, 0.0);
    
    if (use_per_fuel) {
        // WRF-SFIRE defaults: all fuel models use 6.096 m (20 ft)
        for (int i = 1; i <= 13; ++i) {
            fcwh[i] = 6.096;
        }
    } else {
        // Uniform: all entries use global_z_ref
        for (int i = 1; i <= 13; ++i) {
            fcwh[i] = global_z_ref;
        }
    }
    
    return fcwh;
}

/**
 * @brief Build per-fuel roughness length (fcz0) table indexed 0..13
 *
 * @return host-side vector of size 14; index 0 unused, 1-13 valid
 *
 * Returns vector indexed [0..13]. Index 0 is unused (kept for array alignment).
 * Entries 1-13 are WRF-SFIRE defaults [m].
 */
inline std::vector<amrex::Real> build_fcz0_table()
{
    std::vector<amrex::Real> fcz0(14, 0.0);
    
    // WRF-SFIRE data statement: fcz0(k) for k=1..13
    fcz0[1]  = 0.0396;   // FM1
    fcz0[2]  = 0.0396;   // FM2
    fcz0[3]  = 0.100;    // FM3
    fcz0[4]  = 0.2378;   // FM4
    fcz0[5]  = 0.0793;   // FM5
    fcz0[6]  = 0.0991;   // FM6
    fcz0[7]  = 0.0991;   // FM7
    fcz0[8]  = 0.0079;   // FM8
    fcz0[9]  = 0.0079;   // FM9
    fcz0[10] = 0.0396;   // FM10
    fcz0[11] = 0.0396;   // FM11
    fcz0[12] = 0.0911;   // FM12
    fcz0[13] = 0.1188;   // FM13
    
    return fcz0;
}

#endif
