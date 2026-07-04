#ifndef ERF_FIRE_PARAMS_H
#define ERF_FIRE_PARAMS_H

#include <AMReX_ParmParse.H>

/**
 * @file ERF_FireParams.H
 * @brief Fire simulation parameters container
 *
 * Loads fire configuration from ParmParse "erf.fire" prefix at runtime.
 * Stores all tunable parameters for the wildfire spread module including
 * fuel model selection, moisture conditions, wind reference height, and
 * wind adjustment factors.
 */

struct FireParams
{
    bool              enable              = false;       ///< Enable fire module
    bool              fire_debug          = false;       ///< Enable debug output for fire calculations
    int               grid_ratio          = 5;           ///< Fire grid refinement factor (C)
    int               fuel_model_id       = 1;           ///< Anderson FBFM13 model (1-13)
    amrex::Real       moisture_1hr        = 0.08;        ///< 1-hour fuel moisture [fraction]
    amrex::Real       moisture_10hr       = 0.08;        ///< 10-hour fuel moisture [fraction]
    amrex::Real       moisture_100hr      = 0.10;        ///< 100-hour fuel moisture [fraction]
    amrex::Real       ignition_x          = 0.0;         ///< Ignition center x [m]
    amrex::Real       ignition_y          = 0.0;         ///< Ignition center y [m]
    amrex::Real       ignition_r          = 20.0;        ///< Ignition radius [m]
    amrex::Real       wind_ref_ht         = 6.1;         ///< MOST wind reference height [m]
    bool              use_waf             = true;        ///< Apply Wind Adjustment Factor
    std::string       waf_formula         = "andrews";   ///< "andrews" or "behaviorplus"
    bool              use_wind_limit      = true;        ///< Apply MEWS cap to wind
    bool              use_terrain_wind    = true;        ///< Apply terrain corrections
    amrex::Real       k_ridge             = 1.5;         ///< Ridge speed-up factor
    amrex::Real       k_shelter           = 0.6;         ///< Sheltered wind factor
    amrex::Real       k_valley            = 0.8;         ///< Valley channeling factor
    amrex::Real       k_deflect           = 0.3;         ///< Wind deflection factor
    
    // FARSITE parameters
    amrex::Real       farsite_phi_threshold   = 0.0;  ///< front detection level in phi
    int               farsite_use_anderson_lw = 1;       ///< derive coeff_a/b/c from Anderson L/W
    amrex::Real       farsite_coeff_a         = 0.5;  ///< Richards head coefficient
    amrex::Real       farsite_coeff_b         = 0.25; ///< Richards flank coefficient
    amrex::Real       farsite_coeff_c         = 0.1;  ///< Richards backing coefficient
    amrex::Real       farsite_gaussian_sigma  = -1.0; ///< phi stamping radius [m]
    amrex::Real       farsite_cfl_fire        = 0.5;  ///< fire CFL number

    /**
     * @brief Construct FireParams from ParmParse
     *
     * Reads from "erf.fire.*" parameters. Unspecified values use defaults above.
     */
    FireParams()
    {
       amrex::ParmParse pp("erf.fire");
       pp.query("enable", enable);
       pp.query("fire_debug", fire_debug);
       pp.query("grid_ratio", grid_ratio);
       pp.query("fuel_model_id", fuel_model_id);
       pp.query("moisture_1hr", moisture_1hr);
       pp.query("moisture_10hr", moisture_10hr);
       pp.query("moisture_100hr", moisture_100hr);
       pp.query("ignition_x", ignition_x);
       pp.query("ignition_y", ignition_y);
       pp.query("ignition_r", ignition_r);
       pp.query("wind_ref_ht", wind_ref_ht);
       pp.query("use_waf", use_waf);
       pp.query("waf_formula", waf_formula);
       pp.query("use_wind_limit", use_wind_limit);
       pp.query("use_terrain_wind", use_terrain_wind);
       pp.query("k_ridge", k_ridge);
       pp.query("k_shelter", k_shelter);
       pp.query("k_valley", k_valley);
       pp.query("k_deflect", k_deflect);
        
       // FARSITE parameters
       pp.query("farsite.phi_threshold", farsite_phi_threshold);
       pp.query("farsite.use_anderson_lw", farsite_use_anderson_lw);
       pp.query("farsite.coeff_a", farsite_coeff_a);
       pp.query("farsite.coeff_b", farsite_coeff_b);
       pp.query("farsite.coeff_c", farsite_coeff_c);
       pp.query("farsite.gaussian_sigma", farsite_gaussian_sigma);
       pp.query("farsite.cfl_fire", farsite_cfl_fire);
    }
};

#endif
