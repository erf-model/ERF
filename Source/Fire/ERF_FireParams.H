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
    
    // Fuel moisture control (Phase 4)
    bool              moisture_dynamic    = true;        ///< Advance moisture each step via ODEs
    amrex::Real       precip_rate_mm_hr   = 0.0;         ///< Uniform precipitation rate [mm/hr]
    amrex::Real       moisture_live       = 0.60;        ///< Live fuel moisture [fraction]
    
    // Moisture of extinction (Phase 4)
    bool              use_dynamic_mext    = true;        ///< Compute M_x from fuel SAV
    
    // Heat flux parameters (Phase 5)
    amrex::Real       tau_residence_s     = 0.0;         ///< Residence time [s]; 0 = derive from fuel SAV
    amrex::Real       f_c_min             = 0.2;         ///< Minimum combustion fraction
    amrex::Real       f_c_max             = 0.9;         ///< Maximum combustion fraction
    
    // FARSITE parameters
    amrex::Real       farsite_phi_threshold   = 0.1;  ///< front detection level in phi
    int               farsite_use_anderson_lw = 1;       ///< derive coeff_a/b/c from Anderson L/W
    amrex::Real       farsite_coeff_a         = 0.5;  ///< Richards head coefficient
    amrex::Real       farsite_coeff_b         = 0.25; ///< Richards flank coefficient
    amrex::Real       farsite_coeff_c         = 0.1;  ///< Richards backing coefficient
    amrex::Real       farsite_gaussian_sigma  = -1.0; ///< phi stamping radius [m]
    amrex::Real       farsite_cfl_fire        = 0.5;  ///< fire CFL number

    // Terrain file for fine-grid terrain data
    std::string       terrain_file_name       = "";   ///< Path to terrain file for fine-grid interpolation
    
    // CSV output for fire statistics (burned area, perimeter, etc.)
    bool              write_fire_stats_csv    = true; ///< Write burned area and perimeter CSV
    std::string       fire_stats_csv_file     = "fire_stats.csv"; ///< Fire statistics output filename
    
    // Phase 6: fire-atmosphere coupling parameters
    amrex::Real       heat_flux_alfg          = 45.0;  ///< e-folding height for sensible heat distribution [m]
    amrex::Real       fire_atm_feedback       = 1.0;   ///< multiplier on fire flux (0 = one-way, 1 = full)
    // Coupling mode between fire and atmosphere.
    //   "passive"      — fire heat computed but not injected into atmosphere
    //   "lagged"       — lagged coupling: fire heat from step n enters atmosphere
    //                    at step n+1 via RK source term (default, Phase 6 behaviour)
    //   "synchronous"  — synchronous coupling: fire advances using post-dycore wind
    //                    (vars_new), fire heat still enters atmosphere at step n+1
    std::string       coupling_type           = "lagged";
    bool              inject_latent           = true;  ///< inject latent heat alongside sensible (WRF-Fire parity)

    // Phase 8: Albini (1983) ember spotting parameters
    struct SpottingParams {
        bool              enable              = false;   ///< Master on/off switch
        amrex::Real       I_B_min             = 100.0;   ///< Minimum Byram intensity to launch [kW/m]
        amrex::Real       P_base              = 0.01;    ///< Base launch probability per front cell per step
        amrex::Real       terminal_velocity   = 0.5;     ///< Firebrand terminal fall velocity [m/s]
        int               n_traj_steps        = 20;      ///< Trajectory integration sub-steps
        amrex::Real       spot_radius         = 10.0;    ///< Ignition disk radius at landing [m]
        amrex::Real       P_catch             = 1.0;     ///< Catching probability at landing [0–1]
        amrex::Real       reentry_fuel_thresh = 0.05;    ///< Minimum residual fuel fraction for re-entry
        int               random_seed         = 0;       ///< 0 = time-based seed, >0 = fixed seed
        int               spotting_interval   = 1;       ///< Apply spotting every N fire subcycles
        std::string       fuel_system         = "13";    ///< "13" (FBFM13) or "40" (Scott-Burgan 40)
    } spotting;

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
        
       // Fuel moisture control (Phase 4)
       pp.query("moisture_dynamic", moisture_dynamic);
       pp.query("precip_rate_mm_hr", precip_rate_mm_hr);
       pp.query("moisture_live", moisture_live);
        
       // Moisture of extinction (Phase 4)
       pp.query("use_dynamic_mext", use_dynamic_mext);
          
       // Heat flux parameters (Phase 5)
       pp.query("tau_residence_s", tau_residence_s);
       pp.query("f_c_min", f_c_min);
       pp.query("f_c_max", f_c_max);

       // FARSITE parameters
       pp.query("farsite.phi_threshold", farsite_phi_threshold);
       pp.query("farsite.use_anderson_lw", farsite_use_anderson_lw);
       pp.query("farsite.coeff_a", farsite_coeff_a);
       pp.query("farsite.coeff_b", farsite_coeff_b);
       pp.query("farsite.coeff_c", farsite_coeff_c);
       pp.query("farsite.gaussian_sigma", farsite_gaussian_sigma);
       pp.query("farsite.cfl_fire", farsite_cfl_fire);

       // Terrain file
       amrex::ParmParse pp_erf("erf");
       pp_erf.query("terrain_file_name", terrain_file_name);
        
       // CSV output parameters
       pp.query("write_fire_stats_csv", write_fire_stats_csv);
       pp.query("fire_stats_csv_file", fire_stats_csv_file);
        
       // Phase 6: fire-atmosphere coupling parameters
       pp.query("heat_flux_alfg", heat_flux_alfg);
       pp.query("fire_atm_feedback", fire_atm_feedback);
       pp.query("coupling_type", coupling_type);
       pp.query("inject_latent", inject_latent);

       // Phase 8: Albini (1983) ember spotting parameters
       auto& sp = spotting;
       pp.query("spotting.enable",              sp.enable);
       pp.query("spotting.I_B_min",             sp.I_B_min);
       pp.query("spotting.P_base",              sp.P_base);
       pp.query("spotting.terminal_velocity",   sp.terminal_velocity);
       pp.query("spotting.n_traj_steps",        sp.n_traj_steps);
       pp.query("spotting.spot_radius",         sp.spot_radius);
       pp.query("spotting.P_catch",             sp.P_catch);
       pp.query("spotting.reentry_fuel_thresh", sp.reentry_fuel_thresh);
       pp.query("spotting.random_seed",         sp.random_seed);
       pp.query("spotting.spotting_interval",   sp.spotting_interval);
       pp.query("spotting.fuel_system",         sp.fuel_system);
    }

    // Convenience accessors for coupling mode
    bool is_passive()     const { return coupling_type == "passive"; }
    bool is_lagged()      const { return coupling_type == "lagged"; }
    bool is_synchronous() const { return coupling_type == "synchronous"; }
    bool injects_flux()   const { return is_lagged() || is_synchronous(); }
};

#endif
