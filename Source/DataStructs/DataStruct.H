#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <IndexDefines.H>
#include <AdvStruct.H>
#include <DiffStruct.H>
#include <SpongeStruct.H>
#include <TurbStruct.H>

enum struct ABLDriverType {
    None, PressureGradient, GeostrophicWind
};

enum struct CouplingType {
    OneWay, TwoWay
};

enum struct TerrainType {
    Static, Moving
};

enum struct MoistureType {
    Kessler, SAM, FastEddy, Kessler_NoRain, None
};

enum struct WindFarmType {
    Fitch
};

enum struct LandSurfaceType {
    SLM, MM5, None
};

enum struct Coord {
    x, y, z
};

enum Rayleigh {
    tau, ubar, vbar, wbar, thetabar, nvars
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
        } else if (moisture_model_string == "Kessler") {
            moisture_type = MoistureType::Kessler;
        }else if (moisture_model_string == "Kessler_NoRain") {
            moisture_type = MoistureType::Kessler_NoRain;
        } else if (moisture_model_string == "FastEddy") {
            moisture_type = MoistureType::FastEddy;
        } else {
            moisture_type = MoistureType::None;
        }

        // TODO: should we set default for dry??
        // Set a different default for moist vs dry
        if (moisture_type != MoistureType::None) {
            if (moisture_type == MoistureType::FastEddy ||
                moisture_type == MoistureType::SAM) {
                buoyancy_type = 1; // asserted in make buoyancy
            } else {
                buoyancy_type = 2; // uses Tprime
            }
        }

        // Which expression (1,2 or 3) to use for buoyancy
        pp.query("buoyancy_type", buoyancy_type);
        if (buoyancy_type != 1 && buoyancy_type != 2 && buoyancy_type != 3 && buoyancy_type != 4) {
            amrex::Abort("buoyancy_type must be 1, 2, 3 or 4");
        }

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

        if (!use_lagged_delta_rt && !(terrain_type == TerrainType::Moving)) {
            amrex::Error("Can't turn off lagged_delta_rt when terrain not moving");
        }

        // These default to true but are used for unit testing
        pp.query("use_gravity", use_gravity);
        gravity = use_gravity? CONST_GRAV: 0.0;

        pp.query("c_p", c_p);
        rdOcp = R_d / c_p;


#if defined(ERF_USE_POISSON_SOLVE)
        // Should we project the initial velocity field to make it divergence-free?
        pp.query("project_initial_velocity", project_initial_velocity);
#endif

        // Turn off acoustic substepping?
        pp.query("no_substepping", no_substepping);

        pp.query("force_stage1_single_substep", force_stage1_single_substep);
        pp.query("incompressible", incompressible);

        // If this is set, it must be even
        if (incompressible != 0 && no_substepping == 0)
        {
            amrex::Abort("If you specify incompressible, you must specific no_substepping");
        }

        // Include Coriolis forcing?
        pp.query("use_coriolis", use_coriolis);

        // Include Rayleigh damping?
        pp.query("use_rayleigh_damping", use_rayleigh_damping);
        pp.query("rayleigh_damp_U", rayleigh_damp_U);
        pp.query("rayleigh_damp_V", rayleigh_damp_V);
        pp.query("rayleigh_damp_W", rayleigh_damp_W);
        pp.query("rayleigh_damp_T", rayleigh_damp_T);

        // Which external forcings?
        static std::string abl_driver_type_string = "None";
        pp.query("abl_driver_type",abl_driver_type_string);

        if (!abl_driver_type_string.compare("PressureGradient")) {
            abl_driver_type = ABLDriverType::PressureGradient;
        } else if (!abl_driver_type_string.compare("GeostrophicWind")) {
            abl_driver_type = ABLDriverType::GeostrophicWind;
        } else if (!abl_driver_type_string.compare("None")){
            abl_driver_type = ABLDriverType::None; // No ABL driver for simulating classical fluid dynamics problems
        } else {
            amrex::Error("Don't know this abl_driver_type");
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
        pp.query("custom_forcing_uses_primitive_vars", custom_forcing_prim_vars);

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

           advChoice.init_params();
          diffChoice.init_params();
        spongeChoice.init_params();

        turbChoice.resize(max_level+1);
        for (int lev = 0; lev <= max_level; lev++) {
            turbChoice[lev].init_params(lev,max_level);
        }

        // Which type of refinement
        static std::string coupling_type_string = "OneWay";
        pp.query("coupling_type",coupling_type_string);
        if (coupling_type_string == "TwoWay") {
            coupling_type = CouplingType::TwoWay;
        } else if (coupling_type_string == "OneWay") {
            coupling_type = CouplingType::OneWay;
        } else {
            amrex::Abort("Dont know this coupling_type");
        }

        // Which type of windfarm model
        static std::string windfarm_type_string = "Fitch";
        pp.query("windfarm_type", windfarm_type_string);
        if (windfarm_type_string == "Fitch") {
            windfarm_type = WindFarmType::Fitch;
        }
        else {
            amrex::Abort("Dont know this windfarm_type. Currently the only option is Fitch.");
        }
    }

    void display()
    {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "no_substepping              : " << no_substepping << std::endl;
        amrex::Print() << "force_stage1_single_substep : "  << force_stage1_single_substep << std::endl;
        amrex::Print() << "incompressible              : "  << incompressible << std::endl;
        amrex::Print() << "use_coriolis                : " << use_coriolis << std::endl;
        amrex::Print() << "use_rayleigh_damping        : " << use_rayleigh_damping << std::endl;
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

           advChoice.display();
          diffChoice.display();
        spongeChoice.display();

        int max_level = turbChoice.size()-1;
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

        // Convert to radians
        latitude *= (PI/180.);
        sinphi = std::sin(latitude);
        cosphi = std::cos(latitude);

        if (abl_driver_type == ABLDriverType::GeostrophicWind)
        {
            // Read in the geostrophic wind -- we only use this to construct
            //     the forcing term so no need to keep it
            amrex::Vector<amrex::Real> abl_geo_wind(3);
            pp.queryarr("abl_geo_wind",abl_geo_wind);

            abl_geo_forcing = {
                -coriolis_factor * (abl_geo_wind[1]*sinphi - abl_geo_wind[2]*cosphi),
                 coriolis_factor *  abl_geo_wind[0]*sinphi,
                -coriolis_factor *  abl_geo_wind[0]*cosphi
            };
        }
    }

    AdvChoice    advChoice;
    DiffChoice   diffChoice;
    SpongeChoice spongeChoice;
    amrex::Vector<TurbChoice>  turbChoice;

    // Default prefix
    std::string pp_prefix {"erf"};

    int         no_substepping              = 0;
    int         force_stage1_single_substep = 1;
    int         incompressible              = 0;

    bool        test_mapfactor         = false;

    bool        use_terrain            = false;
    int         buoyancy_type          = 1; // uses rhoprime directly

    // Specify what additional physics/forcing modules we use
    bool        use_gravity            = false;
    bool        use_coriolis           = false;
    bool        use_rayleigh_damping   = false;

    bool        rayleigh_damp_U        = false;
    bool        rayleigh_damp_V        = false;
    bool        rayleigh_damp_W        = true;
    bool        rayleigh_damp_T        = false;

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
    amrex::Real cosphi          = 0.0 ;
    amrex::Real sinphi          = 0.0;

    // User-specified forcings in problem definition
    bool custom_rhotheta_forcing = false;
    bool custom_moisture_forcing = false;
    bool custom_forcing_prim_vars = false;

    // Numerical diffusion
    bool use_NumDiff{false};
    amrex::Real NumDiffCoeff{0.};

    CouplingType coupling_type;
    TerrainType  terrain_type;
    MoistureType moisture_type;
    WindFarmType windfarm_type;
    LandSurfaceType lsm_type;

    ABLDriverType abl_driver_type;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_forcing;

    int ave_plane {2};
    // Microphysics params
    bool do_cloud {true};
    bool do_precip {true};
    bool use_moist_background {false};
};
#endif
