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
#include <ERF_DampingStruct.H>
#include <ERF_DiffStruct.H>
#include <ERF_EBStruct.H>
#include <ERF_SpongeStruct.H>
#include <ERF_TurbStruct.H>
#include <ERF_TurbPertStruct.H>

enum MapFacType {
//  This version assumes isotropic
    m_x, u_x, v_x, num,
    m_y = 0, u_y = 1, v_y = 2
//  This version allows for non-isotropic
//  m_x, u_x, v_x,
//  m_y, u_y, v_y, num
};

enum TauType {
    tau11, tau22, tau33, tau12, tau13, tau23, tau21, tau31, tau32
};

AMREX_ENUM(InitType,
    None, Input_Sounding, NCFile, WRFInput, Metgrid, Uniform, HindCast
);

AMREX_ENUM(SoundingType,
    ConstantDensity, Ideal, Isentropic, DryIsentropic
);

AMREX_ENUM(ABLDriverType,
    None, PressureGradient, GeostrophicWind
);

AMREX_ENUM(CouplingType,
    OneWay, TwoWay
);

AMREX_ENUM(SubsteppingType,
    None, Implicit
);

AMREX_ENUM(MeshType,
    ConstantDz, StretchedDz, VariableDz
);

AMREX_ENUM(TerrainType,
    None, StaticFittedMesh, MovingFittedMesh, EB, ImmersedForcing
);

AMREX_ENUM(BuildingsType,
    None, ImmersedForcing
);

AMREX_ENUM(MoistureModelType,
    Eulerian, Lagrangian, Undefined
);

AMREX_ENUM(MoistureType,
    SAM, SAM_NoIce, SAM_NoPrecip_NoIce, Kessler, Kessler_NoRain, SatAdj, Morrison, Morrison_NoIce, None
);

AMREX_ENUM(WindFarmType,
    Fitch, EWP, SimpleAD, GeneralAD, None
);

AMREX_ENUM(WindFarmLocType,
    lat_lon, x_y, None
);

AMREX_ENUM(LandSurfaceType,
    SLM, MM5, None, NOAHMP
);

AMREX_ENUM(RadiationType,
    None, RRTMGP
);

enum struct Coord {
    x, y, z
};

// These are used as integers so must be enum not enum struct
enum Rayleigh {
    ubar, vbar, wbar, thetabar, nvars
};

// These are used as integers so must be enum not enum struct
enum Sponge {
    ubar_sponge, vbar_sponge, nvars_sponge
};

struct MoistureComponentIndices {
    int qv = -1;  // Water vapor
    int qc = -1;  // Cloud liquid water
    int qi = -1;  // Cloud ice
    int qr = -1;  // Rain
    int qs = -1;  // Snow
    int qg = -1;  // Graupel

    // Constructor for easy initialization
    MoistureComponentIndices (int qv_comp, int qc_comp,
                              int qi_comp=-1,
                              int qr_comp=-1,
                              int qs_comp=-1,
                              int qg_comp=-1)
        : qv(qv_comp), qc(qc_comp), qi(qi_comp), qr(qr_comp), qs(qs_comp), qg(qg_comp) {}

    // Default constructor
    MoistureComponentIndices() = default;
};

/**
 * Container holding many of the algorithmic options and parameters
 */

struct SolverChoice {
  public:
    void init_params (int max_level, std::string pp_prefix)
    {
        amrex::ParmParse pp(pp_prefix);

        bool bogus;
        if (pp.query("use_terrain",bogus) > 0)  {
            amrex::Abort("The input use_terrain is deprecated.  Set terrain_type instead.");
        }

        // Do we set map scale factors to 0.5 instead of 1 for testing?
        pp.query("test_mapfactor", test_mapfactor);

        // Which horizontal pressure gradient formulation to use with terrain fitted coords?
        // 0: dp/dx with dp/dz correction (default)
        // 1: gradient of vertically interpolated p, see Klemp 2011
        pp.query("gradp_type", gradp_type);
        AMREX_ALWAYS_ASSERT(gradp_type == 0 || gradp_type == 1);

        // What type of moisture model to use?
        moisture_type = MoistureType::None; // Default
        if (pp.query("moisture_type",moisture_type) > 0) {
            amrex::Abort("The input moisture_type is deprecated.  Set moisture_model instead.");
        }
        pp.query_enum_case_insensitive("moisture_model",moisture_type);
        if ( (moisture_type == MoistureType::Morrison) ||
             (moisture_type == MoistureType::SAM) ) {
            moisture_indices = MoistureComponentIndices(
                RhoQ1_comp,  // water vapor
                RhoQ2_comp,  // cloud water
                RhoQ3_comp,  // cloud ice
                RhoQ4_comp,  // rain
                RhoQ5_comp,  // snow
                RhoQ6_comp   // graupel
            );
        } else if ( (moisture_type == MoistureType::Morrison_NoIce) ||
                    (moisture_type == MoistureType::SAM_NoIce) ) {
            moisture_indices = MoistureComponentIndices(
                RhoQ1_comp,  // water vapor
                RhoQ2_comp,  // cloud water
                -1,          // cloud ice
                RhoQ4_comp   // rain
            );
        } else if ( (moisture_type == MoistureType::SAM_NoPrecip_NoIce) ||
                    (moisture_type == MoistureType::Kessler_NoRain) ||
                    (moisture_type == MoistureType::SatAdj) ) {
            moisture_indices = MoistureComponentIndices(
                RhoQ1_comp,  // water vapor
                RhoQ2_comp   // cloud water
            );
        } else if (moisture_type == MoistureType::Kessler) {
            moisture_indices = MoistureComponentIndices(
                RhoQ1_comp,  // water vapor
                RhoQ2_comp,  // cloud water
                -1,          // cloud ice
                RhoQ3_comp   // rain
            );
        }

        // Set a default for dry
        buoyancy_type.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) {
            buoyancy_type[i] = 1; // uses Rhoprime
        }

        // Set a default for moist
        if (moisture_type != MoistureType::None) {
            if ( !(moisture_type == MoistureType::Kessler_NoRain     ||
                   moisture_type == MoistureType::SAM                ||
                   moisture_type == MoistureType::SAM_NoIce          ||
                   moisture_type == MoistureType::SAM_NoPrecip_NoIce ||
                   moisture_type == MoistureType::Morrison ||
                   moisture_type == MoistureType::Morrison_NoIce ||
                   moisture_type == MoistureType::SatAdj) )
            {
                for (int i = 0; i <= max_level; ++i) {
                    buoyancy_type[i] = 2; // uses Tprime
                }
            }

            pp.query("moisture_tight_coupling",moisture_tight_coupling);
        }

        // Which expression (1,2/3 or 4) to use for buoyancy
        int default_buoyancy_type = buoyancy_type[0];
        read_int_string(max_level, "buoyancy_type", buoyancy_type, default_buoyancy_type);

        // What type of land surface model to use
        lsm_type = LandSurfaceType::None; // Default
        pp.query_enum_case_insensitive("land_surface_model",lsm_type);

        // What type of radiation model to use
        rad_type = RadiationType::None; // Default
        pp.query_enum_case_insensitive("radiation_model", rad_type);

        // Verify that radiation model cannot be RRTMGP if ERF was not compiled with RRTMGP
#ifndef ERF_USE_RRTMGP
        if (rad_type == RadiationType::RRTMGP)
        {
            amrex::Abort("ERF was not compiled with RRTMGP enabled!");
        }
#endif

        // Is the terrain none, static or moving?
        std::string terrain_type_temp = "";
        pp.query("terrain_type", terrain_type_temp);
        if (terrain_type_temp == "Moving") {
            amrex::Warning("erf.terrain_type = Moving is deprecated; please replace Moving by MovingFittedMesh");
            terrain_type = TerrainType::MovingFittedMesh;
        } else if (terrain_type_temp == "Static") {
            amrex::Warning("erf.terrain_type = Static is deprecated; please replace Static by StaticFittedMesh");
            terrain_type = TerrainType::StaticFittedMesh;
        } else {
            pp.query_enum_case_insensitive("terrain_type",terrain_type);
        }

        // Get buildings type
        std::string buildings_type_temp = "";
        pp.query("buildings_type", buildings_type_temp);
        if (buildings_type_temp == "ImmersedForcing") {
            buildings_type = BuildingsType::ImmersedForcing;
        }

        //
        // Read the init_type here to make sure we correctly set the mesh and terrain types
        //
        std::string init_type_temp_string;
        pp.query("init_type",init_type_temp_string);
        if ( (init_type_temp_string == "Real") || (init_type_temp_string == "real") ) {
            amrex::Abort("erf.init_type = Real is deprecated; please replace Real by WRFInput");
        } else if ( (init_type_temp_string == "Ideal") || (init_type_temp_string == "ideal") ) {
            amrex::Abort("erf.init_type = Ideal is deprecated; please replace Ideal by WRFInput");
        } else if (init_type_temp_string == "hindcast") {
            init_type = InitType::HindCast;
        }
        else {
            pp.query_enum_case_insensitive("init_type",init_type);
            use_real_bcs = ( (init_type == InitType::WRFInput) || (init_type == InitType::Metgrid) );
        }

        if ( (init_type == InitType::WRFInput) || (init_type == InitType::Metgrid) ) {
            if (terrain_type != TerrainType::StaticFittedMesh) {
                amrex::Abort("Only terrain_type = StaticFittedMesh are allowed with init_type = WRFInput or Metgrid");
            }
        }

        if (init_type == InitType::WRFInput) {
            if (moisture_type == MoistureType::None) {
                amrex::Abort("Can't have moisture_type = None with init_type = WRFInput");
            }

            // NetCDF wrfbdy lateral boundary file
            std::string nc_bdy_file_temp_string;
            bool has_bdy = pp.query("nc_bdy_file", nc_bdy_file_temp_string);
            if (!has_bdy) use_real_bcs = false;

            bool use_real_bcs_temp = use_real_bcs;
            pp.query("use_real_bcs", use_real_bcs_temp);
            if (use_real_bcs && !use_real_bcs_temp) {
                use_real_bcs = false;
            }
        }

        if (use_real_bcs) {
            pp.query("upwind_real_bcs",upwind_real_bcs);
        }

        // Check for rebalancing with wrfinput
        if (init_type == InitType::WRFInput) {
            pp.query("rebalance_wrfinput",rebalance_wrfinput);
        }

        // How to interpret input_sounding
        if (init_type == InitType::Input_Sounding) {
            pp.query_enum_case_insensitive("sounding_type",sounding_type);
        }

        if (terrain_type == TerrainType::StaticFittedMesh ||
            terrain_type == TerrainType::MovingFittedMesh) {
            mesh_type = MeshType::VariableDz;
        }

        pp.query("grid_stretching_ratio", grid_stretching_ratio);
        if (grid_stretching_ratio != 0) {
            AMREX_ASSERT_WITH_MESSAGE((grid_stretching_ratio >= 1.),
                                      "The grid stretching ratio must be greater than 1");
        }
        if (grid_stretching_ratio >= 1) {
            if (terrain_type == TerrainType::None) {
                terrain_type = TerrainType::StaticFittedMesh;
            }
            if (mesh_type == MeshType::ConstantDz) {
                mesh_type = MeshType::StretchedDz;
            }
            pp.query("zsurface", zsurf);
            if (zsurf != 0.0) {
                amrex::Print() << "Nominal zsurface height != 0, may result in unexpected behavior"
                    << std::endl;
            }
            pp.get("initial_dz", dz0);
        }

        int n_zlevels = pp.countval("terrain_z_levels");
        if (n_zlevels > 0)
        {
            if (terrain_type == TerrainType::None) {
                terrain_type = TerrainType::StaticFittedMesh;
            }
            if (mesh_type == MeshType::ConstantDz) {
                mesh_type = MeshType::StretchedDz;
            }
        }

        // Use lagged_delta_rt in the fast integrator?
        pp.query("use_lagged_delta_rt", use_lagged_delta_rt);

        // These default to true but are used for unit testing
        pp.query("use_gravity", use_gravity);
        gravity = use_gravity? CONST_GRAV: 0.0;

        pp.query("c_p", c_p);
        rdOcp = R_d / c_p;

        // *******************************************************************************
        // Read anelastic etc and over-ride if necessary
        // *******************************************************************************

        read_int_string(max_level, "anelastic", anelastic, 0);
        read_int_string(max_level, "fixed_density", fixed_density, 0);
        read_int_string(max_level, "project_initial_velocity", project_initial_velocity, 0);

        for (int i = 0; i <= max_level; ++i) {
            if (anelastic[i] == 1) {
                project_initial_velocity[i] = 1;
                fixed_density[i]            = 1; // We default to true but are allowed to override below
                buoyancy_type[i]            = 3; // (This isn't actually used when anelastic is set)
            }
        }

        // *******************************************************************************

        if (!upwind_real_bcs && anelastic[0]) {
            amrex::Print() << "Setting upwind_real_bcs to true because we are doing anelastic at level 0" << std::endl;
            upwind_real_bcs = true;
        }

        // *******************************************************************************
        // Read substepping_type and allow for different values at each level
        // *******************************************************************************
        substepping_type.resize(max_level+1);

        for (int i = 0; i <= max_level; i++) {
            substepping_type[i] = SubsteppingType::Implicit;
        }

        int nvals = pp.countval("substepping_type");
        AMREX_ALWAYS_ASSERT(nvals == 0 || nvals == 1 || nvals >= max_level+1);

        if (nvals == 1) {
            pp.query_enum_case_insensitive("substepping_type",substepping_type[0]);
            for (int i = 1; i <= max_level; i++) {
                substepping_type[i] = substepping_type[0];
            }
        } else if (nvals > 1) { // in this case we have asserted nvals >= max_level+1
            for (int i = 0; i <= max_level; i++) {
                pp.query_enum_case_insensitive("substepping_type",substepping_type[i],i);
            }
        }

        pp.query("substepping_diag", substepping_diag);

        pp.query("beta_s", beta_s);



        // *******************************************************************************
        // Error check on deprecated input
        // *******************************************************************************
        int nvals_old = pp.countval("no_substepping");
        if (nvals_old > 0) {
            amrex::Abort("The no_substepping flag is deprecated -- set substepping_type instead");
        }

        // *******************************************************************************

        pp.query("ncorr", ncorr);
        pp.query("poisson_abstol", poisson_abstol);
        pp.query("poisson_reltol", poisson_reltol);

        for (int lev = 0; lev <= max_level; lev++) {
            if (anelastic[lev] != 0)
            {
                substepping_type[lev] = SubsteppingType::None;
            }
        }

        pp.query("force_stage1_single_substep", force_stage1_single_substep);

        // Include Coriolis forcing?
        pp.query("use_coriolis", use_coriolis);
        pp.query("has_lat_lon", has_lat_lon);
        pp.query("variable_coriolis", variable_coriolis);

        // Include four stream radiation approximation
        pp.query("four_stream_radiation", four_stream_radiation);

        // flags for whether to apply other source terms in substep only
        pp.query("immersed_forcing_substep", immersed_forcing_substep); // apply immersed forcing source terms in substep only
        pp.query("forest_substep", forest_substep);                     // apply canopy-related source terms in substep only

        // immersed forcing parameters
        pp.query("if_Cd_scalar", if_Cd_scalar);
        pp.query("if_Cd_momentum", if_Cd_momentum);
        pp.query("if_z0", if_z0);
        pp.query("if_surf_temp_flux", if_surf_temp_flux);
        pp.query("if_init_surf_temp", if_init_surf_temp);
        pp.query("if_surf_heating_rate", if_surf_heating_rate);
        pp.query("if_Olen", if_Olen_in);
        pp.query("if_use_most",if_use_most);

        if ((if_init_surf_temp > 0.0 && if_surf_temp_flux != 1e-8) ||
            (if_init_surf_temp > 0.0 && if_Olen_in != 1e-8) ||
            (if_Olen_in != 1e-8 && if_surf_temp_flux != 1e-8))
        {
            amrex::Abort("Can only specify one of init_surf_temp, surf_temp_flux, or Olen");
        }

        if (if_use_most && buildings_type == BuildingsType::ImmersedForcing)
        {
            amrex::Abort("MOST wall-model with immersed forcing for buildings is not currently supported");
        }

        if (if_surf_temp_flux != 1e-8 && buildings_type == BuildingsType::ImmersedForcing)
        {
            amrex::Abort("Specifying surf_temp_flux with immersed forcing for buildings is not currently supported");
        }

        if (if_Olen_in != 1e-8  && buildings_type == BuildingsType::ImmersedForcing)
        {
            amrex::Abort("Specifying Olen with immersed forcing for buildings is not currently supported");
        }

        // Flag to do MOST rotations with terrain
        pp.query("use_rotate_surface_flux",use_rotate_surface_flux);
        if (use_rotate_surface_flux) {
            AMREX_ASSERT_WITH_MESSAGE(terrain_type != TerrainType::None,"MOST stress rotations are only valid with terrain!");
        }

        // Which external forcings?
        abl_driver_type = ABLDriverType::None; // Default: no ABL driver for simulating classical fluid dynamics problems
        pp.query_enum_case_insensitive("abl_driver_type",abl_driver_type);
        pp.query("const_massflux_u", const_massflux_u);
        pp.query("const_massflux_v", const_massflux_v);
        pp.query("const_massflux_tau", const_massflux_tau);
        pp.query("const_massflux_layer_lo", const_massflux_layer_lo);
        pp.query("const_massflux_layer_hi", const_massflux_layer_hi);

        // Which type of inflow turbulent generation
        pert_type = PerturbationType::None; // Default
        pp.query_enum_case_insensitive("perturbation_type",pert_type);

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
            build_coriolis_forcings_const_lat(pp_prefix);
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

        pp.query("use_moist_background", use_moist_background);

        // Use numerical diffusion?
        pp.query("num_diff_coeff",num_diff_coeff);
        AMREX_ASSERT_WITH_MESSAGE(( (num_diff_coeff >= 0.) && (num_diff_coeff <= 1.) ),
                                  "Numerical diffusion coefficient must be between 0 & 1.");
        use_num_diff = (num_diff_coeff > 0);
        if (use_num_diff) {
            amrex::Print() << "6th-order numerical diffusion turned on with coefficient = "
                << num_diff_coeff << std::endl;
            num_diff_coeff *= std::pow(2.0,-6);
        }

            advChoice.init_params(pp_prefix);
           diffChoice.init_params(pp_prefix);
        dampingChoice.init_params(pp_prefix);
         spongeChoice.init_params(pp_prefix);
             ebChoice.init_params(pp_prefix);

        turbChoice.resize(max_level+1);
        for (int lev = 0; lev <= max_level; lev++) {
            turbChoice[lev].init_params(lev,max_level,pp_prefix);
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
        // MRF
        for (int lev = 0; lev <= max_level; lev++) {
            if (turbChoice[lev].pbl_ysu_use_consistent_coriolis) {
                if (use_coriolis) {
                    turbChoice[lev].pbl_ysu_coriolis_freq = coriolis_factor * sinphi;
                    if (lev == 0) {
                        amrex::Print() << "MRF PBL using ERF coriolis frequency: " << turbChoice[lev].pbl_ysu_coriolis_freq << std::endl;
                    }
                } else {
                    amrex::Abort("MRF cannot use ERF coriolis frequency if not using coriolis");
                }
            }
        }

        // Are we using SHOC? (test on compilation done in turb struct)
        if (turbChoice[0].pbl_type == PBLType::SHOC) { use_shoc = true; }

        // Implicit vertical diffusion (not available with Shoc)
        if (!use_shoc) {
            // This controls the time-centering of the vertical differences in the diffusive term
             bool do_vert_implicit = false;
             if (pp.query("vert_implicit", do_vert_implicit) && do_vert_implicit) {
                 // set to default here
                 vert_implicit_fac[0] = 1.0;
                 vert_implicit_fac[1] = 1.0;
                 vert_implicit_fac[2] = 0.0;
             }

             // This may be one value for all RK stages or a different value in each stage
             int n_impfac = pp.countval("vert_implicit_fac");
             AMREX_ALWAYS_ASSERT(n_impfac == 0 || n_impfac == 1 || n_impfac==3);
             if (n_impfac > 0 && do_vert_implicit) {
                 amrex::Print() << "Overriding defaults with specified implicit factor(s)" << std::endl;
             }

             if (n_impfac == 1) {
                 amrex::Real fac_in;
                 pp.get("vert_implicit_fac", fac_in);
                 for (int i=0; i<3; ++i) { vert_implicit_fac[i] = fac_in; }
             } else if (n_impfac == 3) {
                 pp.getarr("vert_implicit_fac", vert_implicit_fac);
             }

             // If true (default), include implicit contributions to vertical
             // thermal diffusion
             pp.query("implicit_thermal_diffusion", implicit_thermal_diffusion);

             // If true (default), include implicit contributions in tau13, tau23,
             // (and if ERF_IMPLICIT_W is set, tau33) to correct u, v, (and w).
             pp.query("implicit_momentum_diffusion", implicit_momentum_diffusion);

             if (!implicit_thermal_diffusion && !implicit_momentum_diffusion) {
                 amrex::Print() << "Thermal and momentum diffusion are both turned off -- turning off vertical implicit solve" << std::endl;
                 vert_implicit_fac[0] = 0.0;
                 vert_implicit_fac[1] = 0.0;
                 vert_implicit_fac[2] = 0.0;
             }

             // This controls when the vertical implicit solve for the diffusive terms will happen relative to
             //      the acoustic substepping (if it happens, i.e. if any of the implicit_fac > 0.)
             // The default is true (i.e. that it happens before the acoustic substepping).
             pp.query("implicit_before_substep", implicit_before_substep);
         }

        // Which type of multilevel coupling
        coupling_type = CouplingType::TwoWay; // Default
        pp.query_enum_case_insensitive("coupling_type",coupling_type);

        // Test for hybrid (compressible + anelastic) -- in this case we must use one-way coupling
        bool any_anelastic = false;
        bool any_compress  = false;
        for (int lev = 0; lev <= max_level; lev++) {
            if (anelastic[lev] == 0) {
                any_compress = true;
            } else {
                any_anelastic = true;
            }
        }
        if (any_anelastic && any_compress) {
            coupling_type = CouplingType::OneWay;
        }

        // Which type of windfarm model
        windfarm_type = WindFarmType::None; // Default
        pp.query_enum_case_insensitive("windfarm_type",windfarm_type);

        static std::string windfarm_loc_type_string = "None";
        windfarm_loc_type = WindFarmLocType::None;
        pp.query_enum_case_insensitive("windfarm_loc_type",windfarm_loc_type);

        pp.query("windfarm_loc_table",  windfarm_loc_table);
        pp.query("windfarm_spec_table", windfarm_spec_table);
        pp.query("windfarm_blade_table", windfarm_blade_table);
        pp.query("windfarm_airfoil_tables", windfarm_airfoil_tables);
        pp.query("windfarm_spec_table_extra", windfarm_spec_table_extra);

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

        pp.query("hindcast_lateral_forcing", hindcast_lateral_forcing);

        if (hindcast_lateral_forcing) {
            pp.query("hindcast_boundary_data_dir", hindcast_boundary_data_dir);

            if(hindcast_boundary_data_dir.empty()) {
                amrex::Abort("ERROR: Missing input parameter 'erf.hindcast_boundary_data_dir' for boundary data for lateral forcing");
            }
            pp.query("hindcast_data_interval_in_hrs", hindcast_data_interval_in_hrs);
            if(hindcast_data_interval_in_hrs < 0.0) {
                amrex::Abort("ERROR: Input parameter 'erf.hindcast_data_interval_in_hrs' which is the time interval between the "
                             "data files is either missing or set to less than zero");
            }
            pp.query("hindcast_lateral_sponge_strength", hindcast_lateral_sponge_strength);
            pp.query("hindcast_lateral_sponge_length", hindcast_lateral_sponge_length);

            pp.query("hindcast_zhi_sponge_length", hindcast_zhi_sponge_length);
            pp.query("hindcast_zhi_sponge_strength", hindcast_zhi_sponge_strength);

            pp.query("hindcast_zhi_sponge_damping", hindcast_zhi_sponge_damping);

            if (hindcast_lateral_forcing and hindcast_lateral_sponge_strength < 0.0) {
                amrex::Abort("ERROR: Missing input parameter 'erf.hindcast_lateral_sponge_strength' or it is specified to be less than zero");
            }

            if (hindcast_lateral_forcing and hindcast_lateral_sponge_length < 0.0) {
                amrex::Abort("ERROR: Missing input parameter 'erf.hindcast_lateral_sponge_length' or it is specified to be less than zero");
            }

            if (hindcast_zhi_sponge_damping and hindcast_zhi_sponge_strength < 0.0) {
                amrex::Abort("ERROR: Missing input parameter 'erf.hindcast_zhi_sponge_strength' or it is specified to be less than zero");
            }

            if (hindcast_zhi_sponge_damping and hindcast_zhi_sponge_strength < 0.0) {
                amrex::Abort("ERROR: Missing input parameter 'erf.hindcast_zhi_sponge_strength' or it is specified to be less than zero");
            }
        }

        pp.query("io_hurricane_eye_tracker", io_hurricane_eye_tracker);
        if (io_hurricane_eye_tracker) {
            pp.query("hurricane_eye_latitude", hurricane_eye_latitude);
            pp.query("hurricane_eye_longitude", hurricane_eye_longitude);
            if(hurricane_eye_latitude == -1e10 or hurricane_eye_longitude == -1e10) {
                amrex::Abort("ERROR: You are using 'erf.io_hurricane_eye_tracker' to write out the files that track the eye of the hurricane"
                              " but have not provided the initial location of the eye of the hurricane to be tracked. There has to be two"
                              " options in the inputs - erf.hurricane_eye_latitude and erf.hurricane_eye_longitude that gives an approximate"
                              " location of the eye in the initial condition");
            }
        }
    }

    void check_params (int max_level, const amrex::Vector<amrex::Geometry>& geom_vect, amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2> phys_bc_type)
    {
#if 0
        // Warn for PBL models and moisture - these may not yet be compatible
        for (int lev = 0; lev <= max_level; lev++) {
            if ((moisture_type != MoistureType::None) && (turbChoice[lev].pbl_type != PBLType::None)) {
                amrex::Warning("\n*** WARNING: Moisture may not yet be compatible with PBL models, \n    proceed with caution ***");
            }
        }
#endif
        //
        // Buoyancy type check
        //
        for (int lev = 0; lev <= max_level; lev++) {
            if (buoyancy_type[lev] != 1 && buoyancy_type[lev] != 2 && buoyancy_type[lev] != 3 && buoyancy_type[lev] != 4) {
                amrex::Abort("buoyancy_type must be 1, 2, 3 or 4");
            }
        }

        if (!use_lagged_delta_rt && !(terrain_type == TerrainType::MovingFittedMesh)) {
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
        if ( (windfarm_type==WindFarmType::SimpleAD ||
              windfarm_type==WindFarmType::GeneralAD ) && turb_disk_angle < 0.0) {
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


        if ( (const_massflux_u != 0) && !(geom_vect[0].isPeriodic(0)) ) {
            amrex::Abort("Constant mass flux (in x) should be used with periodic boundaries");
        }

        if ( (const_massflux_v != 0) && !(geom_vect[0].isPeriodic(1)) ) {
            amrex::Abort("Constant mass flux (in y) should be used with periodic boundaries");
        }

        // Mesoscale diffusion -- test if LES appropriate with this dx,dy
        for (int lev = 0; lev <= max_level; lev++) {
            if ((geom_vect[lev].CellSize(0) > 2000.) || (geom_vect[lev].CellSize(1) > 2000.))
            {
                if ( (turbChoice[lev].les_type == LESType::Smagorinsky) && !turbChoice[lev].smag2d ) {
                    amrex::Warning("Should use 2-D Smagorinsky for mesoscale resolution");
                } else if (turbChoice[lev].les_type == LESType::Deardorff) {
                    amrex::Warning("Should not use Deardorff LES for mesoscale resolution");
                }
            }
        }

        // Turn off implicit solve if we have no diffusion
        bool l_use_kturb = turbChoice[0].use_kturb;
        for (int lev = 1; lev <= max_level; lev++) {
            l_use_kturb = (l_use_kturb || turbChoice[lev].use_kturb);
        }
        bool l_use_diff    = ( (diffChoice.molec_diff_type != MolecDiffType::None) || l_use_kturb );
        bool l_implicit_diff = (vert_implicit_fac[0] > 0 ||
                                vert_implicit_fac[1] > 0 ||
                                vert_implicit_fac[2] > 0);
        if (l_implicit_diff && !l_use_diff) {
            amrex:: Print() << "No molecular or turbulent diffusion, turning off implicit solve" << std::endl;
            vert_implicit_fac[0] = 0;
            vert_implicit_fac[1] = 0;
            vert_implicit_fac[2] = 0;
        }

        for (int lev = 0; lev <= max_level; lev++) {
            turbChoice[lev].check_params(phys_bc_type);
        }
    }

    void display (int max_level, std::string pp_prefix)
    {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        for (int lev = 0; lev <= max_level; lev++) {
            amrex::Print() << "At level " << lev << " : " << std::endl;
            if (anelastic[lev]) {
                amrex::Print() << "   anelastic with no substepping" << std::endl;
            } else {
               if (substepping_type[lev] == SubsteppingType::None) {
                   amrex::Print() << "   compressible with no substepping" << std::endl;
               } else if (substepping_type[lev] == SubsteppingType::Implicit) {
                   amrex::Print() << "   compressible with implicit substepping" << std::endl;
               }
           }
           if (fixed_density[lev]) {
               amrex::Print() << "   and fixed density" << std::endl;
           }
        }

        amrex::Print() << "vert_implicit_fac           : " << vert_implicit_fac[0] << " "
                                                           << vert_implicit_fac[1] << " "
                                                           << vert_implicit_fac[2];
        if (vert_implicit_fac[0] > 0 ||
            vert_implicit_fac[0] > 1 ||
            vert_implicit_fac[0] > 2)
        {
            amrex::Print() << " (theta=" << implicit_thermal_diffusion
                           << ", momenta=" << implicit_momentum_diffusion;
#ifdef ERF_IMPLICIT_W
            amrex::Print() << ", including w";
#endif
            amrex::Print() << ")";
        }
        amrex::Print() << std::endl;
        amrex::Print() << "use_coriolis                : " << use_coriolis << std::endl;
        amrex::Print() << "use_gravity                 : " << use_gravity << std::endl;

        if (moisture_type == MoistureType::SAM) {
            amrex::Print() << "Moisture Model: SAM" << std::endl;
        } else if (moisture_type == MoistureType::SAM_NoIce) {
            amrex::Print() << "Moisture Model: SAM No Ice" << std::endl;
        } else if (moisture_type == MoistureType::SAM_NoPrecip_NoIce) {
            amrex::Print() << "Moisture Model: SAM No Precip No Ice" << std::endl;
        } else if (moisture_type == MoistureType::Morrison) {
            amrex::Print() << "Moisture Model: Morrison" << std::endl;
        } else if (moisture_type == MoistureType::Morrison_NoIce) {
            amrex::Print() << "Moisture Model: Morrison_NoIce" << std::endl;
        } else if (moisture_type == MoistureType::Kessler) {
            amrex::Print() << "Moisture Model: Kessler" << std::endl;
        } else if (moisture_type == MoistureType::Kessler_NoRain) {
            amrex::Print() << "Moisture Model: Kessler No Rain" << std::endl;
        } else if (moisture_type == MoistureType::SatAdj) {
            amrex::Print() << "Moisture Model: Saturation Adjustment" << std::endl;
        } else {
            amrex::Print() << "Moisture Model: None" << std::endl;
        }

        if (terrain_type == TerrainType::StaticFittedMesh) {
            amrex::Print() << "Terrain Type: StaticFittedMesh" << std::endl;
        } else if (terrain_type == TerrainType::MovingFittedMesh) {
            amrex::Print() << "Terrain Type: MovingFittedMesh" << std::endl;
        } else if (terrain_type == TerrainType::EB) {
            amrex::Print() << "Terrain Type: EB" << std::endl;
            ebChoice.display();
        } else if (terrain_type == TerrainType::ImmersedForcing) {
            amrex::Print() << "Terrain Type: ImmersedForcing" << std::endl;
        } else {
            amrex::Print() << "Terrain Type: None" << std::endl;
        }

        if (buildings_type == BuildingsType::ImmersedForcing) {
            amrex::Print() << "Buildings Type: ImmersedForcing" << std::endl;
        } else {
            amrex::Print() << "Buildings Type: None" << std::endl;
        }

        if (mesh_type == MeshType::ConstantDz) {
            amrex::Print() << "   Mesh Type: ConstantDz" << std::endl;
        } else if (mesh_type == MeshType::StretchedDz) {
            amrex::Print() << "   Mesh Type: StretchedDz" << std::endl;
        } else if (mesh_type == MeshType::VariableDz) {
            amrex::Print() << "   Mesh Type: VariableDz" << std::endl;
        } else {
            amrex::Abort("No mesh_type set!");
        }

        amrex::Print() << "ABL Driver Type: " << std::endl;
        if (abl_driver_type == ABLDriverType::None) {
            amrex::Print() << "    None" << std::endl;
        } else if (abl_driver_type == ABLDriverType::PressureGradient) {
            amrex::Print() << "    Pressure Gradient "
                           << amrex::RealVect(abl_pressure_grad[0],abl_pressure_grad[1],abl_pressure_grad[2])
                           << std::endl;
        } else if (abl_driver_type == ABLDriverType::GeostrophicWind) {
            amrex::Print() << "    Geostrophic Wind "
                           << amrex::RealVect(abl_geo_forcing[0],abl_geo_forcing[1],abl_geo_forcing[2])
                           << std::endl;
        }

        if (max_level > 0) {
            amrex::Print() << "Coupling Type: " << std::endl;
            if (coupling_type == CouplingType::TwoWay) {
                amrex::Print() << "    Two-way" << std::endl;
            } else if (coupling_type == CouplingType::OneWay) {
                amrex::Print() << "    One-way" << std::endl;
            }
        }

        if (rad_type == RadiationType::RRTMGP) {
            amrex::Print() << "Radiation Model: RRTMGP" << std::endl;
        } else {
            amrex::Print() << "Radiation Model: None" << std::endl;
        }

        amrex::Print() << "Gradp_type                  : " << gradp_type << std::endl;

        for (int lev = 0; lev <= max_level; lev++) {
            amrex::Print() << "Buoyancy_type at level " << lev << " : " << buoyancy_type[lev] << std::endl;
        }

            advChoice.display(pp_prefix);
           diffChoice.display();
         spongeChoice.display();
        dampingChoice.display();

        for (int lev = 0; lev <= max_level; lev++) {
            turbChoice[lev].display(lev);
        }
    }

    void build_coriolis_forcings_const_lat (std::string pp_prefix)
    {
        amrex::ParmParse pp(pp_prefix);

        // Read the rotational time period (in seconds)
        amrex::Real rot_time_period = 86400.0;
        pp.query("rotational_time_period", rot_time_period);

        coriolis_factor = 2.0 * 2.0 * PI / rot_time_period;

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

        if (vec_to_fill.size() < max_level+1) {
            vec_to_fill.resize(max_level+1);
        }

        if (nvals == 0) {
            for (int i = 0; i <= max_level; ++i) vec_to_fill[i] = default_int;
        } else if (nvals == 1) {
            for (int i = 0; i <= max_level; ++i) vec_to_fill[i] = temp[0];
        } else {
            for (int i = 0; i <= max_level; ++i) vec_to_fill[i] = temp[i];
        }
    }

    inline static
    InitType init_type = InitType::None;

    inline static
    SoundingType sounding_type = SoundingType::Ideal;

    inline static
    TerrainType terrain_type = TerrainType::None;

    inline static
    BuildingsType buildings_type = BuildingsType::None;

    inline static
    bool use_real_bcs = false;

    inline static
    bool upwind_real_bcs = false;

    inline static
    MeshType mesh_type = MeshType::ConstantDz;

    static
    void set_mesh_type (MeshType new_mesh_type)
    {
        mesh_type = new_mesh_type;
    }

    AdvChoice     advChoice;
    DiffChoice    diffChoice;
    DampingChoice dampingChoice;
    SpongeChoice  spongeChoice;
    amrex::Vector<TurbChoice>  turbChoice;
    EBChoice      ebChoice;

    int         force_stage1_single_substep = 1;

    amrex::Vector<SubsteppingType> substepping_type;
    amrex::Vector<int> anelastic;
    amrex::Vector<int> fixed_density;
    amrex::Vector<int> project_initial_velocity;
    amrex::Vector<int> buoyancy_type;

    // do some extra CFL diagnostics for compressible with substepping
    bool substepping_diag = false;

    // time off-centering coefficient, > 0 for forward weighting (i.e., bias
    // towards the future time step)
    amrex::Real beta_s = 0.1;

    // This controls the time-centering of the *vertical* differences in the diffusive term for
    // theta, u, v (and w if ERF_IMPLICIT_W is set)
    // 0: fully explicit
    // 1: fully implicit
    amrex::Vector<amrex::Real> vert_implicit_fac = {0.0, 0.0, 0.0}; // one value per RK stage
    // if any vert_implicit_fac > 0, then the following apply:
    bool implicit_thermal_diffusion  = true;
    bool implicit_momentum_diffusion = true;
    bool implicit_before_substep     = true;

    int         ncorr               = 1;
    amrex::Real poisson_abstol      = 1e-8;
    amrex::Real poisson_reltol      = 1e-8;

    bool        test_mapfactor         = false;

    int         gradp_type             = 0;

    // Specify what additional physics/forcing modules we use
    bool        use_gravity            = false;
    bool        use_coriolis           = false;
    bool        coriolis_3d            = true;

    bool        four_stream_radiation  = false;

    // Specify whether to apply other various source terms on substep only
    bool        immersed_forcing_substep  = false;
    bool        forest_substep            = false;

    // immersed forcing parameters
    amrex::Real if_Cd_scalar          = 10.0;
    amrex::Real if_Cd_momentum        = 10.0;
    // immersed forcing MOST parameters.
    amrex::Real if_z0                 = 0.1;   // [m]
    amrex::Real if_surf_temp_flux     = 1e-8;  // [K m/s]
    amrex::Real if_init_surf_temp     = 0.0;   // [K]
    amrex::Real if_surf_heating_rate  = 0.0;   // [K/hr]
    amrex::Real if_Olen_in            = 1e-8;  // [m]
    bool        if_use_most           = false;

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

    // MOST stress rotations
    bool use_rotate_surface_flux = false;

    // Should we use SHOC?
    bool use_shoc = false;

    // User wishes to output time averaged velocity fields
    bool time_avg_vel = false;

    // Type of perturbation
    PerturbationType pert_type;

    // Numerical diffusion
    bool use_num_diff{false};
    amrex::Real num_diff_coeff{0.};

    // Rebalance wrfinput
    bool rebalance_wrfinput{false};

    CouplingType coupling_type;
    MoistureType moisture_type;
    WindFarmType windfarm_type;
    WindFarmLocType windfarm_loc_type;
    LandSurfaceType lsm_type;
    RadiationType rad_type;

    ABLDriverType abl_driver_type;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_forcing;
    std::string abl_geo_wind_table;
    bool have_geo_wind_profile {false};

    bool has_lat_lon{false};
    bool variable_coriolis{false};

    int ave_plane {2};

    // Microphysics params
    bool use_moist_background {false};
    MoistureComponentIndices moisture_indices;

    bool moisture_tight_coupling {false};

    std::string windfarm_loc_table, windfarm_spec_table, windfarm_spec_table_extra;
    std::string windfarm_blade_table, windfarm_airfoil_tables;
    amrex::Real sampling_distance_by_D = -1.0;
    amrex::Real turb_disk_angle = -1.0;
    amrex::Real windfarm_x_shift = -1.0;
    amrex::Real windfarm_y_shift = -1.0;

    // Use forest canopy model?
    bool do_forest_drag {false};

    // Enforce constant mass flux?
    amrex::Real const_massflux_u = 0.0;
    amrex::Real const_massflux_v = 0.0;
    amrex::Real const_massflux_tau = 1.0;
    amrex::Real const_massflux_layer_lo = -1e34;
    amrex::Real const_massflux_layer_hi =  1e34;
    int massflux_klo {0}; // these are updated in ERF.cpp
    int massflux_khi {0};

    std::string hindcast_boundary_data_dir;
    amrex::Real hindcast_data_interval_in_hrs = -1.0;
    bool hindcast_lateral_forcing = false;
    amrex::Real hindcast_lateral_sponge_strength = -1.0, hindcast_lateral_sponge_length = -1.0;
    amrex::Real hindcast_zhi_sponge_strength = -1.0, hindcast_zhi_sponge_length = -1.0;
    bool hindcast_zhi_sponge_damping = false;

    bool io_hurricane_eye_tracker = false;
    amrex::Real hurricane_eye_latitude = -1e10, hurricane_eye_longitude = -1e10;
};
#endif
