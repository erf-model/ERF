/**
 * @file ERF_DustParams.H
 * @brief Input parameter container for the ERF-Dust 2D surface emission layer.
 *
 * Parameters are read from the "erf.dust" ParmParse prefix at runtime.
 * Default values correspond to bare quartz soil at an arid site.
 *
 * Phase 1 creates this struct and registers it with CMake and GNUmake.
 * Phase 2 adds grid_ratio and uses this struct in ERF_DustGrid and ERF_DustLayer.
 *
 * References:
 *   Bagnold, R. A. (1941). The Physics of Blown Sand and Desert Dunes.
 *     Methuen, London.
 *   Marticorena, B., & Bergametti, G. (1995). Modeling the atmospheric dust
 *     cycle: 1. Design of a soil-derived dust emission scheme.
 *     J. Geophys. Res., 100, 16415-16430.
 *     https://doi.org/10.1029/95JD00690
 *   Shao, Y., & Lu, H. (2000). A simple expression for wind erosion threshold
 *     friction velocity.
 *     J. Geophys. Res., 105, 22437-22443.
 *     https://doi.org/10.1029/2000JD900304
 */

#ifndef ERF_DUST_PARAMS_H
#define ERF_DUST_PARAMS_H

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <string>

/**
 * @struct DustParams
 * @brief Container for ERF-Dust module parameters.
 *
 * All parameters are read from the "erf.dust" ParmParse prefix at runtime.
 */
struct DustParams {

    bool        enable              = false;    ///< Enable dust module
    bool        dust_debug          = false;    ///< Write per-timestep debug fields to plotfile
    int         n_size_bins         = 3;        ///< Number of particle size bins; each bin maps to
                                                 ///<   one RhoQ*_comp scalar component in ERF.
                                                 ///<   Bin 0: coarse  (2.5 um < d <= 10 um)
                                                 ///<   Bin 1: PM10    (1.0 um < d <= 2.5 um)
                                                 ///<   Bin 2: PM2.5   (d <= 1.0 um)
    int         grid_ratio          = 1;        ///< Refinement factor of dust grid relative to ERF
                                                 ///<   level-0 horizontal grid [dimensionless].
                                                 ///<   grid_ratio = 1 places one dust cell per
                                                 ///<   atmospheric cell. Values > 1 subdivide each
                                                 ///<   atmospheric cell into grid_ratio^2 dust cells.
                                                 ///<   All atmospheric box x,y sizes must be divisible
                                                 ///<   by grid_ratio. Set erf.dust.grid_ratio to match
                                                 ///<   the desired surface resolution.
    amrex::Vector<amrex::Real> bin_diameter_um  = {7.0, 3.5, 0.7}; ///< Representative aerodynamic diameter per bin [um]
    amrex::Real particle_density    = 2650.0;  ///< Bulk particle density [kg/m^3]; quartz default
    amrex::Real z0_dust             = 0.01;    ///< Aerodynamic roughness length over bare tailings [m]
    amrex::Real silt_fraction       = 0.10;    ///< Default surface silt mass fraction [dimensionless]
                                                 ///<   Used when no surface map file is provided.
                                                 ///<   Typical range for mine tailings: 0.05-0.30.
                                                 ///<   Reference: Marticorena & Bergametti (1995),
                                                 ///<   J. Geophys. Res., 100, 16415-16430.
                                                 ///<   https://doi.org/10.1029/95JD00690
    amrex::Real threshold_A_coeff  = 0.0123;  ///< Bagnold threshold coefficient A [dimensionless]
                                                 ///<   u*_t = A * sqrt(rho_p * g * d / rho_a)
                                                 ///<   Reference: Bagnold (1941), The Physics of Blown
                                                 ///<   Sand and Desert Dunes, Methuen, London.
    amrex::Real crust_index         = 0.0;     ///< Default surface crust strength index [dimensionless]
                                                 ///<   Range [0,1]: 0 = uncrusted, 1 = fully crusted.
                                                 ///<   Updated by PHREEQC reader in Phase 4.
    amrex::Real phreeqc_update_interval_s = 86400.0; ///< Interval at which PHREEQC output files
                                                 ///<   are re-read [s]. Default = 1 day.
                                                 ///<   Geochemical processes (crust formation,
                                                 ///<   efflorescence) evolve on day-to-week timescales,
                                                 ///<   much slower than the atmospheric timestep.
    amrex::Real alpha_crust     = 0.5;  ///< Crust strength reduction coefficient for u*_t [-].
                                         ///<   u*_t_new = u*_t_base * (1 - alpha_crust * crust_index).
                                         ///<   Reference: Marticorena & Bergametti (1995),
                                         ///<     https://doi.org/10.1029/95JD00690
    amrex::Real alpha_efflor    = 0.3;  ///< Efflorescence reduction coefficient for u*_t [-].
    amrex::Real ustar_t_base    = -1.0; ///< Base u*_t before chemistry [m/s].
                                         ///<   Negative means compute from Bagnold at init.
    amrex::Real test_ustar      = 0.0;   ///< Placeholder u* [m/s] used in Phase 6 before
                                         ///<   Phase 9 extracts u* from the MRF surface layer.
                                         ///<   Set > 0 to test emission physics.
                                         ///<   0 = no emission (all cells below threshold).
    amrex::Real rho_air         = 1.225; ///< Air density used in emission flux [kg/m^3].
                                         ///<   Standard atmosphere default.
    std::string blast_schedule_file = "";  ///< Path to blast schedule CSV.
                                           ///<   Empty = no blast events.
                                           ///<   Format: time_s cx cy radius mass_kg_m2
                                           ///<     [mineral_type] [priority]
                                           ///<   Reference: ERF_DustBlastSchedule.H
    amrex::Real blast_reactivity    = 2.0; ///< Reactivity multiplier for fresh blast
                                           ///<   surfaces relative to weathered tailings [-].
                                           ///<   Applied to injected mass flux.
    amrex::Real supp_tau_base_s     = 3600.0; ///< Base suppression half-life [s].
                                           ///<   Water: ~3600 s. MgCl2: ~43200 s.
                                           ///<   Reference: Gillies et al. (2005),
                                           ///<     https://doi.org/10.1080/10473289.2005.10464701
    amrex::Real test_surf_temp_K    = 293.15;  ///< Surface temperature placeholder [K].
                                           ///<   Phase 9 replaces with MRF T_sfc.
    amrex::Real test_wind_speed     = 5.0;     ///< 10 m wind speed placeholder [m/s].
                                           ///<   Phase 9 replaces with MRF 10m wind.
    amrex::Real zref = 10.0; ///< Wind extraction height [m]. Set equal to erf.most.zref.
                          ///<   Hong & Pan (1996):
                          ///<   https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
    std::string phreeqc_crust_var    = "crust_index";    ///< Variable/column name for crust index.
    std::string phreeqc_silt_var     = "silt_fraction";  ///< Variable/column name for silt fraction.
    std::string phreeqc_efflor_var   = "efflorescence";  ///< Variable/column name for efflorescence.
    std::string phreeqc_supp_var     = "suppression_mod";///< Variable/column name for suppression modifier.
    std::string phreeqc_metal_var    = "metal_as_bin0";  ///< Variable/column name for As mass fraction bin 0.
    std::string surface_map_file    = "";      ///< Path to NetCDF surface property map; empty = uniform
    std::string phreeqc_output_file = "";      ///< Path to PHREEQC output NetCDF; empty = disabled
    std::string soil_type_file      = "";      ///< Path to soil type raster (.asc or .nc).
                                                 ///<   Code 0 = undefined; 1-16 = STATSGO natural soils;
                                                 ///<   100 = mine tailings; 101 = Li brine pond;
                                                 ///<   102 = REE tailings; 103 = Cu tailings;
                                                 ///<   104 = unpaved haul road.
                                                 ///<   STATSGO: https://www.nrcs.usda.gov/resources/
                                                 ///<     data-and-reports/ssurgo/statsgo2-data
                                                 ///<   AP-42 road emission: https://www.epa.gov/
                                                 ///<     air-emissions-factors-and-quantification/
                                                 ///<     ap-42-compilation-air-emission-factors
    std::string silt_fraction_file  = "";      ///< Path to silt fraction raster. Empty = uniform silt_fraction.
    std::string crust_index_file    = "";      ///< Path to crust index raster. Empty = uniform crust_index.
    std::string moisture_flag_file  = "";      ///< Path to moisture inhibition flag raster. Empty = uniform 0.
    std::string suppression_file    = "";      ///< Path to suppression coverage raster. Empty = uniform 0.
                                                 ///<   Supported formats: ESRI ASCII (.asc), NetCDF (.nc).
                                                 ///<   NetCDF requires ERF_ENABLE_NETCDF=ON.

    // Phase 10: Atmosphere coupling parameters
    bool        transport_bins_separately = false; ///< When false (default), use one 3D
                                                   ///<   scalar for total dust mass density.
                                                   ///<   When true, use n_size_bins scalars
                                                   ///<   (Phase 11). Default false reduces
                                                   ///<   3D advection cost.
    amrex::Real atm_feedback = 1.0;                ///< Coupling strength [0,1]. 0 = coupling
                                                   ///<   disabled. 1 = full injection.

    // Phase 11: Stokes settling parameters
    std::vector<amrex::Real> bin_diameters = {7.0e-6, 2.5e-6, 50.0e-6}; ///< Per-bin particle
                                                   ///<   diameters [m]. Used for Stokes settling.
                                                   ///<   Default: PM10 bin (7 um), PM2.5 bin (2.5 um),
                                                   ///<   coarse bin (50 um).
                                                   ///<   Reference: Seinfeld & Pandis (2006), Ch. 9.

    // Phase 12: Dry deposition parameters
    amrex::Real deposition_E0 = 3.0e-3; ///< Surface collection efficiency [-].
                                        ///<   Zhang et al. (2001):
                                        ///<   https://doi.org/10.1016/S1352-2310(00)00326-5
                                        ///<   Bare mine surface: 3e-3
                                        ///<   Paved road: 1e-4
                                        ///<   Vegetated buffer: 1e-2

    // Phase 13: Return fields and feedbacks
    amrex::Real loading_feedback_coeff = 0.0; ///< Shao (2001) loading feedback coefficient.
                                               ///<   Effective threshold: u*_t *= (1 + coeff * C_sfc).
                                               ///<   0 = disabled (default).
                                               ///<   Shao (2001): https://doi.org/10.1029/2001JD900171
    bool use_dynamic_moisture = false;         ///< Pass Q1fx3 to fill_dust_moist_from_atm.
                                               ///<   Requires moisture_type != None in ERF.
                                               ///<   Safe to set true with no moisture: Q1fx3 will
                                               ///<   be null and function returns immediately.
                                               ///<   Fecan et al. (1999):
                                               ///<   https://doi.org/10.1007/s00585-999-0149-7

    // Phase 14: MRF nonlocal countergradient PBL diffusion
    // MRF turbulent Schmidt number for the dust scalar is controlled via
    // erf.dust.mrf_Sc_t in the PBL turbulence parameters (ERF_TurbStruct.H).
    // Default: same diffusivity as heat (Sc_t = Pr_t).
    // Seinfeld & Pandis (2006), Ch. 16: https://doi.org/10.1002/9781118227572

    // Phase 16: Dust diagnostics and plotfile output
    int         dust_plot_int      = -1;           ///< Steps between dust plotfiles. -1 = disabled.
    std::string dust_plot_prefix   = "plt_dust_";  ///< Plotfile name prefix.
    std::string dust_diag_file     = "dust_diag.dat"; ///< CSV diagnostics file path.

    // Phase 17: EPA NAAQS compliance output
    std::string dust_naaqs_file    = "dust_naaqs.csv"; ///< NAAQS compliance CSV output.
                                                        ///< 40 CFR Part 50.

    // Phase 18: MSHA worker exposure tracking
    amrex::Real msha_pel_mg_m3       = 5.0;            ///< PEL threshold [mg/m³], 30 CFR 56.5001
    amrex::Real msha_shift_duration_s = 28800.0;       ///< Shift duration [s]. Default 8 hours.
    std::string msha_exposure_file   = "msha_exposure.csv";        ///< Per-step CSV output
    std::string msha_shift_file      = "msha_shift_summary.csv";   ///< End-of-shift summary
    amrex::Vector<std::string> msha_receptor_names;    ///< Receptor point names
    amrex::Vector<amrex::Real> msha_receptor_x;        ///< Receptor x-coordinates [m]
    amrex::Vector<amrex::Real> msha_receptor_y;        ///< Receptor y-coordinates [m]

    // Phase 19: Lagrangian super-particle source-receptor attribution
    bool        enable_particles        = false; ///< Lagrangian particle tracking.
                                                   ///< Requires ERF_USE_PARTICLES.
    int         particle_release_interval = 1;  ///< Release every N steps.

    // Phase 20: multi-mine site PHREEQC tables.
    // Each site has a name, PHREEQC file, and bounding box in domain metres.
    // Sites are indexed 1..N; 0 = unassigned (global Phase 4 table used).
    // Reference: Parkhurst & Appelo (2013). https://pubs.usgs.gov/tm/06/a43/
    amrex::Vector<std::string> site_names;         ///< Site name labels.
    amrex::Vector<std::string> site_phreeqc_files; ///< Per-site PHREEQC file paths.
    amrex::Vector<amrex::Real> site_x_lo;          ///< Site bounding box x_lo [m].
    amrex::Vector<amrex::Real> site_y_lo;          ///< Site bounding box y_lo [m].
    amrex::Vector<amrex::Real> site_x_hi;          ///< Site bounding box x_hi [m].
    amrex::Vector<amrex::Real> site_y_hi;          ///< Site bounding box y_hi [m].

    // Phase 21: PHREEQC deposition feedback file writer.
    // phreeqc_feedback_interval_s = 0 disables periodic writes.
    // phreeqc_feedback_file: deposition grid file (overwritten each interval).
    // phreeqc_site_summary_file: per-site CSV (appended each interval).
    // Reference: Parkhurst & Appelo (2013). https://pubs.usgs.gov/tm/06/a43/
    amrex::Real phreeqc_feedback_interval_s  = 0.0;
    std::string phreeqc_feedback_file        = "dust_dep_feedback.dat";
    std::string phreeqc_site_summary_file    = "dust_dep_site_summary.csv";

    // Phase 22: haul road vehicle schedule.
    // CSV format: road_name, x_lo, y_lo, x_hi, y_hi, road_width_m,
    //             vehicle_weight_t, silt_pct, vmt_per_h, start_time_s, end_time_s
    // Reference: EPA AP-42 Ch. 13.2.2. https://www.epa.gov/air-emissions-factors-and-quantification/
    //            ap-42-compilation-air-emissions-factors
    std::string road_schedule_file  = "";             ///< Road schedule CSV. Empty = disabled.
    std::string road_diag_file      = "dust_road_diag.csv"; ///< Road diagnostics CSV.

    /**
     * Constructor that reads all parameters from ParmParse.
     */
    DustParams() {
        amrex::ParmParse pp("erf.dust");

        pp.query("enable", enable);
        pp.query("dust_debug", dust_debug);
        pp.query("n_size_bins", n_size_bins);
        pp.queryarr("bin_diameter_um", bin_diameter_um);
        pp.query("particle_density", particle_density);
        pp.query("z0_dust", z0_dust);
        pp.query("silt_fraction", silt_fraction);
        pp.query("threshold_A_coeff", threshold_A_coeff);
        pp.query("crust_index", crust_index);
        pp.query("grid_ratio", grid_ratio);
        pp.query("phreeqc_update_interval_s", phreeqc_update_interval_s);
        pp.query("alpha_crust", alpha_crust);
        pp.query("alpha_efflor", alpha_efflor);
        pp.query("ustar_t_base", ustar_t_base);
        pp.query("test_ustar", test_ustar);
        pp.query("rho_air", rho_air);
        pp.query("blast_schedule_file", blast_schedule_file);
        pp.query("blast_reactivity", blast_reactivity);
        pp.query("supp_tau_base_s",  supp_tau_base_s);
        pp.query("test_surf_temp_K", test_surf_temp_K);
        pp.query("test_wind_speed",  test_wind_speed);
        pp.query("zref", zref);
        pp.query("phreeqc_crust_var", phreeqc_crust_var);
        pp.query("phreeqc_silt_var", phreeqc_silt_var);
        pp.query("phreeqc_efflor_var", phreeqc_efflor_var);
        pp.query("phreeqc_supp_var", phreeqc_supp_var);
        pp.query("phreeqc_metal_var", phreeqc_metal_var);
        pp.query("surface_map_file", surface_map_file);
        pp.query("phreeqc_output_file", phreeqc_output_file);
        pp.query("soil_type_file",     soil_type_file);
        pp.query("silt_fraction_file", silt_fraction_file);
        pp.query("crust_index_file",   crust_index_file);
        pp.query("moisture_flag_file", moisture_flag_file);
        pp.query("suppression_file",   suppression_file);
        pp.query("transport_bins_separately", transport_bins_separately);
        pp.query("atm_feedback", atm_feedback);
        pp.query("deposition_E0", deposition_E0);
        pp.query("loading_feedback_coeff", loading_feedback_coeff);
        pp.query("use_dynamic_moisture",   use_dynamic_moisture);
        pp.query("dust_plot_int",    dust_plot_int);
        pp.query("dust_plot_prefix", dust_plot_prefix);
        pp.query("dust_diag_file",   dust_diag_file);
        pp.query("dust_naaqs_file",  dust_naaqs_file);
        pp.query("msha_pel_mg_m3", msha_pel_mg_m3);
        pp.query("msha_shift_duration_s", msha_shift_duration_s);
        pp.query("msha_exposure_file", msha_exposure_file);
        pp.query("msha_shift_file", msha_shift_file);
        pp.queryarr("msha_receptor_names", msha_receptor_names);
        pp.queryarr("msha_receptor_x", msha_receptor_x);
        pp.queryarr("msha_receptor_y", msha_receptor_y);
        
        pp.query("enable_particles", enable_particles);
        pp.query("particle_release_interval", particle_release_interval);
        
        // Allow per-bin diameter override from input file.
        // erf.dust.bin_diameters = d0 d1 d2 ...
        {
            amrex::Vector<amrex::Real> tmp;
            int n = pp.queryarr("bin_diameters", tmp);
            if (n > 0) {
                bin_diameters.assign(tmp.begin(), tmp.end());
            }
        }

        // Phase 20: multi-mine site configuration
        pp.queryarr("site_names",         site_names);
        pp.queryarr("site_phreeqc_files", site_phreeqc_files);
        pp.queryarr("site_x_lo",          site_x_lo);
        pp.queryarr("site_y_lo",          site_y_lo);
        pp.queryarr("site_x_hi",          site_x_hi);
        pp.queryarr("site_y_hi",          site_y_hi);

        // Phase 21: PHREEQC deposition feedback file writer
        pp.query("phreeqc_feedback_interval_s", phreeqc_feedback_interval_s);
        pp.query("phreeqc_feedback_file",        phreeqc_feedback_file);
        pp.query("phreeqc_site_summary_file",    phreeqc_site_summary_file);

        // Phase 22: road schedule CSV.
        pp.query("road_schedule_file", road_schedule_file);
        pp.query("road_diag_file",     road_diag_file);
   }  // end DustParams constructor
 };     // end struct DustParams
#endif // ERF_DUST_PARAMS_H
