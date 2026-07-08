#ifndef ERF_FIRE_PARAMS_H
#define ERF_FIRE_PARAMS_H

#include <AMReX_ParmParse.H>
#include <ERF_FireBreak.H>
#include <ERF_FireAcceleration.H>
#include <string>
#include <vector>

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
    // Phase 13A: Per-fuel wind height
    bool              use_per_fuel_wind_ht = false; ///< Use per-fuel fcwh table for wind height
    amrex::Real       waf_fcz0_scale       = 1.0;   ///< Scale factor on fcz0 roughness for WAF [dimensionless]
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

    // Phase 9: Crown fire initiation parameters.
    // Van Wagner (1977) surface-to-crown criterion + Cruz et al. (2005) active crown ROS.
    // Applicable to conifer forest stands. Not appropriate for grass or shrub fuels.
    struct CrownFireParams {
        bool          enable          = false;   ///< Enable crown fire initiation
        amrex::Real   canopy_base_ht  = 5.0;    ///< Canopy base height C_BH [m]
        amrex::Real   canopy_bulk_den = 0.2;    ///< Canopy bulk density C_BD [kg/m³]
        amrex::Real   canopy_depth    = 10.0;   ///< Canopy depth [m]
        amrex::Real   foliar_moisture = 1.0;    ///< Foliar moisture M_fold [fraction]
        amrex::Real   M_c             = 0.30;   ///< Critical foliar moisture [fraction] (Van Wagner 1977)
        amrex::Real   h_crown_BTU_lb  = 8000.0; ///< Canopy heat content [BTU/lb]
        amrex::Real   wind_10m_kmh    = 0.0;    ///< 10-m wind speed [km/h]; 0 = derive from fire_wind_eff
        std::string   ros_model       = "cruz"; ///< Crown ROS model: "cruz", "rothermel1991", "van_wagner_proxy"
        bool          use_passive_blend = false; ///< Use Van Wagner passive-to-active ROS blend
    } crown;

    // Phase 9: Flame temperature diagnostic
    std::string       flame_temp_method   = "byram_radiant"; ///< "byram_radiant", "mcalpine_heat", "nelson_emissivity"
    amrex::Real       flame_temp_T_amb    = 300.0;           ///< Ambient temperature for flame temp formula [K]
    
    // Phase 9: Flame tilt angle diagnostic
    bool              compute_flame_tilt  = false;           ///< Compute flame tilt angle diagnostic [degrees]
    amrex::Real       flame_tilt_T_amb    = 300.0;           ///< Ambient temperature for flame tilt [K]
    amrex::Real       flame_tilt_rho_air  = 1.2;             ///< Air density for flame tilt [kg/m³]

    // Phase 10: Spatial fuel map parameters.
    struct FuelMapParams {
        std::string   fuel_map_file       = "";        ///< Path to fuel map file; empty = uniform
        std::string   fuel_map_format     = "ascii";   ///< "ascii" or "lcp"
        amrex::Real   blending_fraction   = 0.0;       ///< ROS blending at fuel boundaries [0-1]
    } fuel_map;

    // Phase 11: Multi-ignition schedule parameters.
    // When ignition_schedule_file is non-empty, a CSV schedule of timed ignition
    // events is loaded and applied during advance(). The primary ignition from
    // ignition_x/y/r is always applied at t=0 regardless of this setting.
    //
    // CSV format: time_s  cx  cy  radius  [source_type]  [priority]
    //   source_type: "sphere" (default) or "polyline"
    //   priority: integer; higher priority events fire first within the same timestep
    //   Comment lines start with # or !
    //
    // Polygon ignition (initial perimeter from a vertex file):
    //   polygon_file: path to CSV of polygon vertices (x y per line)
    //   polygon_type: "polygon" (closed) or "polyline" (open line fire)
    //   polyline_width: half-width of line ignition zone [m]
    struct IgnitionParams {
        std::string   ignition_schedule_file  = "";         ///< Path to ignition schedule CSV; empty = disabled
        std::string   polygon_file            = "";         ///< Path to polygon vertex CSV; empty = disabled
        std::string   polygon_type            = "polygon";  ///< "polygon" or "polyline"
        amrex::Real   polyline_width          = 10.0;       ///< Half-width of polyline ignition zone [m]
    } ignition;

    // Phase 12: Fire acceleration at startup.
    // Small fires do not immediately achieve their quasi-steady-state (QSS) ROS.
    // Two models are available, selected by use_temporal:
    //
    //   Size-based (use_temporal = false, default):
    //     alpha = 1 - exp(-r_fire / L_acc)
    //     R_effective = alpha * R_Rothermel
    //   where r_fire is the effective radius from burned area and L_acc is the
    //   acceleration length scale [m].
    //   Reference: Rothermel (1983) GTR INT-143.
    //
    //   Temporal (use_temporal = true):
    //     R(t) = R_E * (1 - exp(-A * t))
    //   where t is elapsed time since ignition [s] and A [1/min] is the
    //   acceleration constant, selected based on perimeter length relative to
    //   perim_limit (point vs line fire).
    //   Reference: McAlpine & Wakimoto (1991) Forest Science 37(5).
    //
    // Disabled by default (enable = false). When disabled, no MultiFab is
    // allocated and the call in advance() returns immediately.
    AccelerationParams accel;

    // Phase 13B: Rate-of-spread model selection and parameters
    
    // Phase 13B: Rate-of-spread model selection.
    // "rothermel"    — Rothermel (1972) single-class model (default, Phases 1-12)
    // "macarthur"    — MacArthur (1966) Australian formula R = 0.18*exp(0.8424*U)
    // "balbi"        — Balbi (2009) physical model R = A*(1+sin α − cos α)
    // "cheney_gould" — Cheney & Gould (1998) empirical grassland model
    std::string       ros_model = "rothermel";  ///< ROS model selector

    // Phase 13B: Balbi (2009) physical fire spread model parameters.
    // Only used when ros_model = "balbi".
    // Defaults from Balbi (2009) Table 1.
    struct BalbiParams {
        amrex::Real T_a      = 300.0;    ///< Ambient temperature [K]
        amrex::Real T_f      = 1000.0;   ///< Mean flame temperature [K]
        amrex::Real T_i      = 600.0;    ///< Ignition temperature [K]
        amrex::Real delta_H  = 2.26e6;   ///< Latent heat of vaporisation [J/kg]
        amrex::Real C_pf     = 1800.0;   ///< Specific heat of dry fuel [J/(kg·K)]
        amrex::Real r_00     = 2.5e-4;   ///< Radiation length scale [m]
        amrex::Real tau_0    = 75591.0;  ///< Residence time coefficient [s/m]
    } balbi;

    // Phase 13B: Cheney-Gould (1998) grassland model parameters.
    // Only used when ros_model = "cheney_gould".
    struct CheneyGouldParams {
        amrex::Real moisture = 10.0;  ///< Dead fine fuel moisture [%]
        amrex::Real curing   = 1.0;   ///< Degree of curing [0-1]; 1 = fully cured
    } cheney_gould;

    // Phase 15: BEHAVE multi-class Rothermel parameters.
    // Only used when ros_model = "behave".
    // Controls the dynamic live-to-dead load transfer threshold.
    // All other fuel properties come from FuelModelParams via get_anderson_fuel_params().
    struct BehaveParams {
        amrex::Real dynamic_transfer_lo  = 0.30;  ///< Live herba moisture [fraction] below which full transfer occurs
        amrex::Real dynamic_transfer_hi  = 1.20;  ///< Live herba moisture [fraction] above which no transfer occurs
    } behave;

    // Phase 10: Firebreak barriers (up to 10, indexed by N=0..9)
    int n_firebreaks = 0;
    std::vector<FireBreak> firebreaks;

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
       // Phase 13A: Per-fuel wind height
       pp.query("use_per_fuel_wind_ht", use_per_fuel_wind_ht);
       pp.query("waf_fcz0_scale", waf_fcz0_scale);
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

       // Phase 9: Crown fire initiation parameters
       auto& cr = crown;
       pp.query("crown.enable",          cr.enable);
       pp.query("crown.canopy_base_ht",  cr.canopy_base_ht);
       pp.query("crown.canopy_bulk_den", cr.canopy_bulk_den);
       pp.query("crown.canopy_depth",    cr.canopy_depth);
       pp.query("crown.foliar_moisture", cr.foliar_moisture);
       pp.query("crown.M_c",             cr.M_c);
       pp.query("crown.h_crown_BTU_lb",  cr.h_crown_BTU_lb);
       pp.query("crown.wind_10m_kmh",    cr.wind_10m_kmh);
       pp.query("crown.ros_model",       cr.ros_model);
       pp.query("crown.use_passive_blend", cr.use_passive_blend);

       // Phase 9: Flame temperature diagnostic
       pp.query("flame_temp_method", flame_temp_method);
       pp.query("flame_temp_T_amb",  flame_temp_T_amb);
         
       // Phase 9: Flame tilt angle diagnostic
       pp.query("compute_flame_tilt", compute_flame_tilt);
       pp.query("flame_tilt_T_amb",   flame_tilt_T_amb);
       pp.query("flame_tilt_rho_air", flame_tilt_rho_air);

       // Phase 10: fuel map
       auto& fm = fuel_map;
       pp.query("fuel_map.file",               fm.fuel_map_file);
       pp.query("fuel_map.format",             fm.fuel_map_format);
       pp.query("fuel_map.blending_fraction",  fm.blending_fraction);

       // Phase 11: multi-ignition schedule and polygon ignition
       auto& ig = ignition;
       pp.query("ignition.schedule_file", ig.ignition_schedule_file);
       pp.query("ignition.polygon_file",  ig.polygon_file);
       pp.query("ignition.polygon_type",  ig.polygon_type);
       pp.query("ignition.polyline_width", ig.polyline_width);

       // Phase 12: fire acceleration
       auto& ac = accel;
       pp.query("accel.enable",          ac.enable);
       pp.query("accel.use_temporal",    ac.use_temporal);
       pp.query("accel.L_acc",           ac.L_acc);
       pp.query("accel.A_point",         ac.A_point);
       pp.query("accel.A_line",          ac.A_line);
       pp.query("accel.perim_limit",     ac.perim_limit);
       pp.query("accel.enable_wind_lag", ac.enable_wind_lag);
       pp.query("accel.tau_wind",        ac.tau_wind);

       // Phase 13B: ROS model selection
       pp.query("ros_model", ros_model);

       // Balbi parameters
       auto& bl = balbi;
       pp.query("balbi.T_a",     bl.T_a);
       pp.query("balbi.T_f",     bl.T_f);
       pp.query("balbi.T_i",     bl.T_i);
       pp.query("balbi.delta_H", bl.delta_H);
       pp.query("balbi.C_pf",    bl.C_pf);
       pp.query("balbi.r_00",    bl.r_00);
       pp.query("balbi.tau_0",   bl.tau_0);

       // Cheney-Gould parameters
       auto& cg = cheney_gould;
       pp.query("cheney_gould.moisture", cg.moisture);
       pp.query("cheney_gould.curing",   cg.curing);

       // Phase 15: BEHAVE parameters
       auto& bh = behave;
       pp.query("behave.dynamic_transfer_lo", bh.dynamic_transfer_lo);
       pp.query("behave.dynamic_transfer_hi", bh.dynamic_transfer_hi);

       // Phase 10: firebreaks (indexed, stop at first missing)
       for (int nb = 0; nb < 10; ++nb) {
           std::string type_key = "firebreak." + std::to_string(nb) + ".type";
           std::string type_val;
           if (!pp.query(type_key.c_str(), type_val)) break;
           FireBreak fb;
           fb.type = type_val;
           std::string pfx = "firebreak." + std::to_string(nb) + ".";
           pp.query((pfx + "x_lo").c_str(),   fb.x_lo);
           pp.query((pfx + "y_lo").c_str(),   fb.y_lo);
           pp.query((pfx + "x_hi").c_str(),   fb.x_hi);
           pp.query((pfx + "y_hi").c_str(),   fb.y_hi);
           pp.query((pfx + "cx").c_str(),     fb.cx);
           pp.query((pfx + "cy").c_str(),     fb.cy);
           pp.query((pfx + "radius").c_str(), fb.radius);
           firebreaks.push_back(fb);
       }
       n_firebreaks = static_cast<int>(firebreaks.size());
    }

    // Convenience accessors for coupling mode
    bool is_passive()     const { return coupling_type == "passive"; }
    bool is_lagged()      const { return coupling_type == "lagged"; }
    bool is_synchronous() const { return coupling_type == "synchronous"; }
    bool injects_flux()   const { return is_lagged() || is_synchronous(); }
};

#endif
