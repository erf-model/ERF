#ifndef ERF_DATA_STRUCT_H_
#define ERF_DATA_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_AdvStruct.H>
#include <ERF_DiffStruct.H>
#include <ERF_SpongeStruct.H>
#include <ERF_TurbStruct.H>
#include <ERF_TurbPertStruct.H>

enum struct ABLDriverType {
    None, PressureGradient, GeostrophicWind
};

enum struct CouplingType {
    OneWay, TwoWay
};

enum struct TerrainType {
    Static, Moving
};

enum struct MoistureModelType{
    Eulerian, Lagrangian, Undefined
};

enum struct MoistureType {
    Kessler, SAM, SAM_NoIce, SAM_NoPrecip_NoIce, Kessler_NoRain, None
};

enum struct WindFarmType {
    Fitch, EWP, SimpleAD, GeneralAD, None
};

enum struct WindFarmLocType{
    lat_lon, x_y, None
};

enum struct LandSurfaceType {
    SLM, MM5, None
};

enum struct Coord {
    x, y, z
};

enum Rayleigh {
    ubar, vbar, wbar, thetabar, nvars
};

enum Sponge {
    ubar_sponge, vbar_sponge, nvars_sponge
};

enum struct PerturbationType {
    perturbSource, perturbDirect, None
};

/**
 * Container holding many of the algorithmic options and parameters
 */

struct SolverChoice {
  public:
    void init_params(int max_level)
    {
        amrex::ParmParse pp(pp_prefix);

        // Do we have terrain (or grid stretching)?
        pp.query("use_terrain", use_terrain);

        pp.query("grid_stretching_ratio", grid_stretching_ratio);
        if (grid_stretching_ratio != 0) {
            AMREX_ASSERT_WITH_MESSAGE((grid_stretching_ratio >= 1.),
                                      "The grid stretching ratio must be greater than 1");
        }
        if (grid_stretching_ratio >= 1) {
            if (!use_terrain) {
                amrex::Print() << "Turning terrain on to enable grid stretching" << std::endl;
                use_terrain = true;
            }
            pp.query("zsurface", zsurf);
            if (zsurf != 0.0) {
                amrex::Print() << "Nominal zsurface height != 0, may result in unexpected behavior"
                    << std::endl;
            }
            pp.get("initial_dz", dz0);
        }

        // Do we set map scale factors to 0.5 instead of 1 for testing?
        pp.query("test_mapfactor", test_mapfactor);

        // What type of moisture model to use
        static std::string moisture_model_string = "None";
        pp.query("moisture_model", moisture_model_string);
        if (moisture_model_string == "SAM") {
            moisture_type = MoistureType::SAM;
            RhoQv_comp = RhoQ1_comp;
            RhoQr_comp = RhoQ4_comp;
        } else if (moisture_model_string == "SAM_NoIce") {
            moisture_type = MoistureType::SAM_NoIce;
            RhoQv_comp = RhoQ1_comp;
            RhoQr_comp = RhoQ4_comp;
        } else if (moisture_model_string == "SAM_NoPrecip_NoIce") {
            moisture_type = MoistureType::SAM_NoPrecip_NoIce;
            RhoQv_comp = RhoQ1_comp;
        } else if (moisture_model_string == "Kessler") {
            moisture_type = MoistureType::Kessler;
            RhoQv_comp = RhoQ1_comp;
            RhoQr_comp = RhoQ3_comp;
        }else if (moisture_model_string == "Kessler_NoRain") {
            moisture_type = MoistureType::Kessler_NoRain;
            RhoQv_comp = RhoQ1_comp;
        } else {
            moisture_type = MoistureType::None;
        }

        // TODO: should we set default for dry??
        // Set a different default for moist vs dry
        if (moisture_type != MoistureType::None) {
            if (moisture_type == MoistureType::Kessler_NoRain ||
                moisture_type == MoistureType::SAM ||
                moisture_type == MoistureType::SAM_NoIce ||
                moisture_type == MoistureType::SAM_NoPrecip_NoIce) {
                buoyancy_type = 1; // asserted in make buoyancy
            } else {
                buoyancy_type = 2; // uses Tprime
            }
        }

        // Which expression (1,2 or 3) to use for buoyancy
        pp.query("buoyancy_type", buoyancy_type);

        // What type of land surface model to use
        static std::string lsm_model_string = "None";
        pp.query("land_surface_model", lsm_model_string);
        if (lsm_model_string == "SLM") {
            lsm_type = LandSurfaceType::SLM;
        } else if (lsm_model_string == "MM5") {
            lsm_type = LandSurfaceType::MM5;
        } else {
            lsm_type = LandSurfaceType::None;
        }

        // Is the terrain static or moving?
        static std::string terrain_type_string = "Static";
        pp.query("terrain_type",terrain_type_string);
        if (terrain_type_string == "Moving" || terrain_type_string == "moving") {
            terrain_type = TerrainType::Moving;
        } else if (terrain_type_string == "Static" || terrain_type_string == "static") {
            terrain_type = TerrainType::Static;
        } else {
            amrex::Abort("terrain_type can be either Moving/moving or Static/static");
        }

        // Use lagged_delta_rt in the fast integrator?
        pp.query("use_lagged_delta_rt", use_lagged_delta_rt);

        // These default to true but are used for unit testing
        pp.query("use_gravity", use_gravity);
        gravity = use_gravity? CONST_GRAV: 0.0;

        pp.query("c_p", c_p);
        rdOcp = R_d / c_p;

        read_int_string(max_level, "no_substepping", no_substepping, 0);

#if defined(ERF_USE_POISSON_SOLVE)
        // Should we project the initial velocity field to make it divergence-free?
        pp.query("project_initial_velocity", project_initial_velocity);

        read_int_string(max_level, "anelastic", anelastic, 0);

        bool any_anelastic = false;
        for (int i = 0; i <= max_level; ++i) {
            if (anelastic[i] == 1) any_anelastic = true;
        }

        // If anelastic at all, we do not advect rho -- it is always == rho0
        if (any_anelastic == 1) {
            constant_density = true;
            buoyancy_type = 2; // uses Tprime
        } else {
            constant_density = false;  // We default to false but allow the user to set it
            pp.query("constant_density", constant_density);
        }

        pp.query("project_every_stage", project_every_stage);
        pp.query("ncorr", ncorr);
        pp.query("poisson_abstol", poisson_abstol);
        pp.query("poisson_reltol", poisson_reltol);

        for (int lev = 0; lev <= max_level; lev++) {
            if (anelastic[lev] != 0)
            {
                no_substepping[lev] = 1;
            }
        }
#else
        anelastic.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) anelastic[i] = 0;
#endif

        pp.query("force_stage1_single_substep", force_stage1_single_substep);

        // Include Coriolis forcing?
        pp.query("use_coriolis", use_coriolis);

        // Include Rayleigh damping (separate flags for each variable)
        pp.query("rayleigh_damp_U", rayleigh_damp_U);
        pp.query("rayleigh_damp_V", rayleigh_damp_V);
        pp.query("rayleigh_damp_W", rayleigh_damp_W);
        pp.query("rayleigh_damp_T", rayleigh_damp_T);
        pp.query("rayleigh_dampcoef", rayleigh_dampcoef);
        pp.query("rayleigh_zdamp", rayleigh_zdamp);

        // Flag to do explicit MOST formulation
        pp.query("use_explicit_most",use_explicit_most);

        // Flag to do MOST rotations with terrain
        pp.query("use_rotate_most",use_rotate_most);
        if (use_rotate_most) {
            AMREX_ASSERT_WITH_MESSAGE(use_terrain,"MOST stress rotations are only valid with terrain!");
            AMREX_ASSERT_WITH_MESSAGE(use_explicit_most, "MOST Stress rotations are only valid with explicit MOST!");
        }

        // Which external forcings?
        static std::string abl_driver_type_string = "None";
        pp.query("abl_driver_type",abl_driver_type_string);

        if (!abl_driver_type_string.compare("PressureGradient")) {
            abl_driver_type = ABLDriverType::PressureGradient;
        } else if (!abl_driver_type_string.compare("GeostrophicWind")) {
            abl_driver_type = ABLDriverType::GeostrophicWind;
        } else if (!abl_driver_type_string.compare("None")) {
            abl_driver_type = ABLDriverType::None; // No ABL driver for simulating classical fluid dynamics problems
        } else {
            amrex::Error("Don't know this abl_driver_type");
        }

        // Which type of inflow turbulent generation
        static std::string turb_pert_type_string = "None";
        pp.query("perturbation_type",turb_pert_type_string);
        if (!turb_pert_type_string.compare("source")) {
            pert_type = PerturbationType::perturbSource;
            amrex::Print() << "Perturb temperature via source term\n";
        } else if (!turb_pert_type_string.compare("direct")) {
            pert_type = PerturbationType::perturbDirect;
            amrex::Print() << "Perturb temperature via direct field add\n";
        } else if (!turb_pert_type_string.compare("None")) {
            pert_type = PerturbationType::None;
            amrex::Print() << "No inflow perturbation method selected\n";
        } else {
            amrex::Abort("Dont know this perturbation_type");
        }

        amrex::Vector<amrex::Real> abl_pressure_grad_in = {0.0, 0.0, 0.0};
        pp.queryarr("abl_pressure_grad",abl_pressure_grad_in);
        for(int i = 0; i < AMREX_SPACEDIM; ++i) abl_pressure_grad[i] = abl_pressure_grad_in[i];

        amrex::Vector<amrex::Real> abl_geo_forcing_in = {0.0, 0.0, 0.0};
        if(pp.queryarr("abl_geo_forcing",abl_geo_forcing_in)) {
            amrex::Print() << "Specified abl_geo_forcing: (";
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                abl_geo_forcing[i] = abl_geo_forcing_in[i];
                amrex::Print() << abl_geo_forcing[i] << " ";
            }
            amrex::Print() << ")" << std::endl;
        }

        if (use_coriolis)
        {
            build_coriolis_forcings();
        }

        pp.query("add_custom_rhotheta_forcing", custom_rhotheta_forcing);
        pp.query("add_custom_moisture_forcing", custom_moisture_forcing);
        pp.query("add_custom_w_subsidence", custom_w_subsidence);
        pp.query("add_custom_geostrophic_profile", custom_geostrophic_profile);
        pp.query("custom_forcing_uses_primitive_vars", custom_forcing_prim_vars);

        pp.query("nudging_from_input_sounding", nudging_from_input_sounding);

        have_geo_wind_profile = (!abl_geo_wind_table.empty() || custom_geostrophic_profile);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!(!abl_geo_wind_table.empty() && custom_geostrophic_profile),
            "Should not have both abl_geo_wind_table and custom_geostrophic_profile set.");

        pp.query("Ave_Plane", ave_plane);

        pp.query("mp_clouds", do_cloud);
        pp.query("mp_precip", do_precip);
        pp.query("use_moist_background", use_moist_background);

        // Use numerical diffusion?
        pp.query("use_NumDiff",use_NumDiff);
        if(use_NumDiff) {
            pp.query("NumDiffCoeff",NumDiffCoeff);
            AMREX_ASSERT_WITH_MESSAGE(( (NumDiffCoeff >= 0.) && (NumDiffCoeff <= 1.) ),
                                      "Numerical diffusion coefficient must be between 0 & 1.");
            NumDiffCoeff *= std::pow(2.0,-6);
        }

        // Use monotonic advection?
        pp.query("use_mono_adv",use_mono_adv);

           advChoice.init_params();
          diffChoice.init_params();
        spongeChoice.init_params();

        turbChoice.resize(max_level+1);
        for (int lev = 0; lev <= max_level; lev++) {
            turbChoice[lev].init_params(lev,max_level);
        }

        // YSU PBL: use consistent coriolis frequency
        for (int lev = 0; lev <= max_level; lev++) {
            if (turbChoice[lev].pbl_ysu_use_consistent_coriolis) {
                if (use_coriolis) {
                    turbChoice[lev].pbl_ysu_coriolis_freq = coriolis_factor * sinphi;
                    if (lev == 0) {
                        amrex::Print() << "YSU PBL using ERF coriolis frequency: " << turbChoice[lev].pbl_ysu_coriolis_freq << std::endl;
                    }
                } else {
                    amrex::Abort("YSU cannot use ERF coriolis frequency if not using coriolis");
                }
            }
        }


        // Which type of refinement
        static std::string coupling_type_string = "TwoWay";
        pp.query("coupling_type",coupling_type_string);
        if (coupling_type_string == "TwoWay") {
            coupling_type = CouplingType::TwoWay;
        } else if (coupling_type_string == "OneWay") {
            coupling_type = CouplingType::OneWay;
        } else {
            amrex::Abort("Dont know this coupling_type");
        }

        // Which type of windfarm model
        static std::string windfarm_type_string = "None";
        pp.query("windfarm_type", windfarm_type_string);
        if (windfarm_type_string == "Fitch") {
            windfarm_type = WindFarmType::Fitch;
        }
        else if (windfarm_type_string == "EWP") {
            windfarm_type = WindFarmType::EWP;
        }
        else if (windfarm_type_string == "SimpleActuatorDisk") {
            windfarm_type = WindFarmType::SimpleAD;
        }
        else if (windfarm_type_string == "GeneralActuatorDisk") {
            windfarm_type = WindFarmType::GeneralAD;
        }
        else if (windfarm_type_string != "None") {
            amrex::Abort("Are you using windfarms? Dont know this windfarm_type. windfarm_type"
                         " has to be Fitch, EWP, SimpleActuatorDisk, GeneralActuatorDisk or None.");
        }

        static std::string windfarm_loc_type_string = "None";
        windfarm_loc_type = WindFarmLocType::None;
        pp.query("windfarm_loc_type", windfarm_loc_type_string);
        if (windfarm_loc_type_string == "lat_lon") {
            windfarm_loc_type = WindFarmLocType::lat_lon;
        }
        else if (windfarm_loc_type_string == "x_y") {
            windfarm_loc_type = WindFarmLocType::x_y;
        }
        else if (windfarm_loc_type_string != "None") {
            amrex::Abort("Are you using windfarms? Dont know this windfarm_loc_type."
                         " windfarm_loc_type has to be specified as lat_lon or x_y.");
        }

        pp.query("windfarm_loc_table",  windfarm_loc_table);
        pp.query("windfarm_spec_table", windfarm_spec_table);
        pp.query("windfarm_blade_table", windfarm_blade_table);
        pp.query("windfarm_airofil_tables", windfarm_airfoil_tables);

        // Sampling distance upstream of the turbine to find the
        // incoming free stream velocity as a factor of the diameter of the
        // turbine. ie. the sampling distance will be this number multiplied
        // by the diameter of the turbine
        pp.query("sampling_distance_by_D", sampling_distance_by_D);
        pp.query("turb_disk_angle_from_x", turb_disk_angle);

        pp.query("windfarm_x_shift",windfarm_x_shift);
        pp.query("windfarm_y_shift",windfarm_y_shift);
        // Test if time averaged data is to be output
        pp.query("time_avg_vel",time_avg_vel);

        check_params(max_level);
    }

    void check_params(int max_level)
    {
        // Warn for PBL models and moisture - these may not yet be compatible
        for (int lev = 0; lev <= max_level; lev++) {
            if ((moisture_type != MoistureType::None) && (turbChoice[lev].pbl_type != PBLType::None)) {
                amrex::Warning("\n*** WARNING: Moisture may not yet be compatible with PBL models, \n    proceed with caution ***");
            }
        }
        //
        // Buoyancy type check
        //
        if (buoyancy_type != 1 && buoyancy_type != 2 && buoyancy_type != 3 && buoyancy_type != 4) {
            amrex::Abort("buoyancy_type must be 1, 2, 3 or 4");
        }

        if (!use_lagged_delta_rt && !(terrain_type == TerrainType::Moving)) {
            amrex::Error("Can't turn off lagged_delta_rt when terrain not moving");
        }

        //
        // Wind farm checks
        //
        if (windfarm_type==WindFarmType::SimpleAD and sampling_distance_by_D < 0.0) {
             amrex::Abort("To use simplified actuator disks, you need to provide a variable"
                           " erf.sampling_distance_by_D in the inputs which specifies the upstream"
                           " distance as a factor of the turbine diameter at which the incoming free stream"
                           " velocity will be computed at.");
        }
        pp.query("turb_disk_angle_from_x", turb_disk_angle);
        if((windfarm_type==WindFarmType::SimpleAD or windfarm_type==WindFarmType::GeneralAD) and turb_disk_angle < 0.0) {
            amrex::Abort("To use simplified actuator disks, you need to provide a variable"
                          " erf.turb_disk_angle_from_x in the inputs which is the angle of the face of the"
                          " turbine disk from the x-axis. A turbine facing an oncoming flow in the x-direction"
                          " will have turb_disk_angle value of 90 deg.");
        }
        if (windfarm_loc_type == WindFarmLocType::lat_lon and (windfarm_x_shift < 0.0 or windfarm_y_shift < 0.0)) {
            amrex::Abort("You are using windfarms with latitude-logitude option to position the turbines."
                         " For this you should provide the inputs erf.windfarm_x_shift and"
                         " erf.windfarm_y_shift which are the values by which the bounding box of the"
                         " windfarm is shifted from the x and the y axes.");
        }
    }

    void display(int max_level)
    {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "force_stage1_single_substep : "  << force_stage1_single_substep << std::endl;
        for (int lev = 0; lev <= max_level; lev++) {
            amrex::Print() << "anelastic      at level : " << lev << " is " <<     anelastic[lev] << std::endl;
            amrex::Print() << "no_substepping at level : " << lev << " is " << no_substepping[lev] << std::endl;
        }
        amrex::Print() << "use_coriolis                : " << use_coriolis << std::endl;
        amrex::Print() << "use_gravity                 : " << use_gravity << std::endl;

        if (coupling_type == CouplingType::TwoWay) {
            amrex::Print() << "Using two-way coupling " << std::endl;
        } else if (coupling_type == CouplingType::OneWay) {
            amrex::Print() << "Using one-way coupling " << std::endl;
        }

        if (terrain_type == TerrainType::Static) {
            amrex::Print() << "Using static terrain " << std::endl;
        } else {
            amrex::Print() << "Using moving terrain " << std::endl;
        }

        if (abl_driver_type == ABLDriverType::None) {
            amrex::Print() << "ABL Driver Type: " << "None" << std::endl;
            amrex::Print() << "No ABL driver selected " << std::endl;
        } else if (abl_driver_type == ABLDriverType::PressureGradient) {
            amrex::Print() << "ABL Driver Type: " << "PressureGradient" << std::endl;
            amrex::Print() << "Driving abl_pressure_grad: (";
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
                amrex::Print() << abl_pressure_grad[i] << " ";
            amrex::Print() << ")" << std::endl;
        } else if (abl_driver_type == ABLDriverType::GeostrophicWind) {
            amrex::Print() << "ABL Driver Type: " << "GeostrophicWind" << std::endl;
            amrex::Print() << "Driving abl_geo_forcing: (";
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
                amrex::Print() << abl_geo_forcing[i] << " ";
            amrex::Print() << ")" << std::endl;
        }

        amrex::Print() << "Buoyancy_type               : " << buoyancy_type << std::endl;

           advChoice.display();
          diffChoice.display();
        spongeChoice.display();

        for (int lev = 0; lev <= max_level; lev++) {
            turbChoice[lev].display(lev);
        }
    }

    void build_coriolis_forcings()
    {
        amrex::ParmParse pp(pp_prefix);

        // Read the rotational time period (in seconds)
        amrex::Real rot_time_period = 86400.0;
        pp.query("rotational_time_period", rot_time_period);

        coriolis_factor = 2.0 * 2.0 * PI / rot_time_period;
        amrex::Print() << "Coriolis factor = " << coriolis_factor << std::endl;

        amrex::Real latitude = 90.0;
        pp.query("latitude", latitude);

        pp.query("coriolis_3d", coriolis_3d);

        // Convert to radians
        latitude *= (PI/180.);
        sinphi = std::sin(latitude);
        if (coriolis_3d) {
            cosphi = std::cos(latitude);
        }

        amrex::Print() << "Coriolis frequency, f = " << coriolis_factor * sinphi << " 1/s" << std::endl;

        if (abl_driver_type == ABLDriverType::GeostrophicWind) {
            // Read in the geostrophic wind -- we only use this to construct
            //     the forcing term so no need to keep it
            amrex::Vector<amrex::Real> abl_geo_wind(3);
            pp.queryarr("abl_geo_wind",abl_geo_wind);

            if(!pp.query("abl_geo_wind_table",abl_geo_wind_table)) {
                abl_geo_forcing = {
                    -coriolis_factor * (abl_geo_wind[1]*sinphi - abl_geo_wind[2]*cosphi),
                     coriolis_factor *  abl_geo_wind[0]*sinphi,
                    -coriolis_factor *  abl_geo_wind[0]*cosphi
                };
            } else {
                amrex::Print() << "NOTE: abl_geo_wind_table provided, ignoring input abl_geo_wind" << std::endl;
            }
        }
    }

    void read_int_string (int max_level, const char* string_to_read,
                          amrex::Vector<int>& vec_to_fill, int default_int)
    {
        amrex::ParmParse pp("erf");
        int nvals = pp.countval(string_to_read);
        AMREX_ALWAYS_ASSERT(nvals == 0 || nvals == 1 || nvals >= max_level+1);
        amrex::Vector<int> temp; temp.resize(nvals);
        pp.queryarr(string_to_read,temp);
        if (nvals == 0) {
            for (int i = 0; i <= max_level; ++i) vec_to_fill.push_back(default_int);
        } else if (nvals == 1) {
            for (int i = 0; i <= max_level; ++i) vec_to_fill.push_back(temp[0]);
        } else {
            for (int i = 0; i <= max_level; ++i) vec_to_fill.push_back(temp[i]);
        }
    }


    AdvChoice    advChoice;
    DiffChoice   diffChoice;
    SpongeChoice spongeChoice;
    amrex::Vector<TurbChoice>  turbChoice;

    // Default prefix
    std::string pp_prefix {"erf"};

    int         force_stage1_single_substep = 1;

    amrex::Vector<int> no_substepping;
    amrex::Vector<int> anelastic;

    int         constant_density    = 0;
    int         project_every_stage = 1;
    int         ncorr               = 1;
    amrex::Real poisson_abstol      = 1e-10;
    amrex::Real poisson_reltol      = 1e-10;

    bool        test_mapfactor         = false;

    bool        use_terrain            = false;

    int         buoyancy_type          = 1; // uses rhoprime directly

    // Specify what additional physics/forcing modules we use
    bool        use_gravity            = false;
    bool        use_coriolis           = false;
    bool        coriolis_3d            = true;

    bool        rayleigh_damp_U        = false;
    bool        rayleigh_damp_V        = false;
    bool        rayleigh_damp_W        = false;
    bool        rayleigh_damp_T        = false;
    amrex::Real rayleigh_dampcoef      = 0.2; // inverse time scale [1/s]
    amrex::Real rayleigh_zdamp         = 500.0; // damping layer depth [m]
    amrex::Real rayleigh_ztop;

    // This defaults to true but can be set to false for moving terrain cases only
    bool        use_lagged_delta_rt    = true;

    // Constants
    amrex::Real gravity;
    amrex::Real c_p = Cp_d; // specific heat at constant pressure for dry air [J/(kg-K)]
    amrex::Real rdOcp;

    // Staggered z levels for vertical grid stretching
    amrex::Real grid_stretching_ratio   = 0;
    amrex::Real zsurf                   = 0.0;
    amrex::Real dz0;

#if defined(ERF_USE_POISSON_SOLVE)
    int project_initial_velocity = 1;
#endif

    // Coriolis forcing
    amrex::Real coriolis_factor = 0.0;
    amrex::Real cosphi          = 0.0;
    amrex::Real sinphi          = 0.0;

    // User-specified forcings in problem definition
    bool custom_rhotheta_forcing  = false;
    bool custom_moisture_forcing  = false;
    bool custom_w_subsidence      = false;
    bool custom_geostrophic_profile = false;
    bool custom_forcing_prim_vars = false;

    // Do we use source terms to nudge the solution towards
    //    the time-varying data provided in input sounding files?
    bool nudging_from_input_sounding = false;

    // User specified MOST BC type
    bool use_explicit_most = false;

    // MOST stress rotations
    bool use_rotate_most = false;

    // User wishes to output time averaged velocity fields
    bool time_avg_vel = false;

    // Type of perturbation
    PerturbationType pert_type;

    // Numerical diffusion
    bool use_NumDiff{false};
    amrex::Real NumDiffCoeff{0.};

    // Monotonic advection limiter
    bool use_mono_adv{false};

    CouplingType coupling_type;
    TerrainType  terrain_type;
    MoistureType moisture_type;
    WindFarmType windfarm_type;
    WindFarmLocType windfarm_loc_type;
    LandSurfaceType lsm_type;

    ABLDriverType abl_driver_type;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_forcing;
    std::string abl_geo_wind_table;
    bool have_geo_wind_profile {false};

    int ave_plane {2};
    // Microphysics params
    bool do_cloud {true};
    bool do_precip {true};
    bool use_moist_background {false};
    int RhoQv_comp {-1};

    // This component will be model-dependent:
    // if a model with no rain, this will stay -1
    // if Kessler, then it will be set to RhoQ3
    // if     SAM, then it will be set to RhoQ4
    int RhoQr_comp {-1};

    std::string windfarm_loc_table, windfarm_spec_table;
    std::string windfarm_blade_table, windfarm_airfoil_tables;
    amrex::Real sampling_distance_by_D = -1.0;
    amrex::Real turb_disk_angle = -1.0;
    amrex::Real windfarm_x_shift = -1.0;
    amrex::Real windfarm_y_shift = -1.0;
};
#endif
