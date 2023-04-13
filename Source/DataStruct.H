#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <ERF_Math.H>

enum class ABLDriverType {
    None, PressureGradient, GeostrophicWind
};

enum class Coord {
    x, y, z
};

enum class LESType {
    None, Smagorinsky, Deardorff
};

enum class PBLType {
    None, MYNN25
};

enum class MolecDiffType {
    None, Constant, ConstantAlpha
};

struct SolverChoice {
  public:
    void init_params()
    {
        amrex::ParmParse pp(pp_prefix);

        // Do we have terrain (or grid stretching)?
        pp.query("use_terrain", use_terrain);

        // Do we use map scale factors?
        pp.query("test_mapfactor", test_mapfactor);

        // Which expression (1,2 or 3) to use for buoyancy
        pp.query("buoyancy_type", buoyancy_type);
        if (buoyancy_type != 1 && buoyancy_type != 2 && buoyancy_type != 3) {
            amrex::Abort("buoyancy_type must be 1, 2 or 3");
        }

        // Is the terrain static or moving?
        pp.query("terrain_type", terrain_type);

        // Use lagged_delta_rt in the fast integrator?
        pp.query("use_lagged_delta_rt", use_lagged_delta_rt);
        if (!use_lagged_delta_rt && !terrain_type == 1) {
            amrex::Error("Can't turn off lagged_delta_rt when terrain not moving");
        }

        // These default to true but are used for unit testing
        pp.query("use_gravity", use_gravity);
        gravity = use_gravity? CONST_GRAV: 0.0;

        pp.query("c_p", c_p);
        rdOcp = R_d / c_p;

#if defined(ERF_USE_WARM_NO_PRECIP)
        pp.query("tau_cond", tau_cond);
#endif

#if defined(ERF_USE_POISSON_SOLVE)
        // Should we project the initial velocity field to make it divergence-free?
        pp.query("project_initial_velocity", project_initial_velocity);
#endif

        // Which LES closure?
        static std::string les_type_string = "None";
        pp.query("les_type",les_type_string);

        if (!les_type_string.compare("Smagorinsky")) {
            les_type = LESType::Smagorinsky;
        } else if (!les_type_string.compare("Deardorff")) {
            les_type = LESType::Deardorff;
        } else if (!les_type_string.compare("None")) {
            les_type = LESType::None; // Means DNS
        } else {
            amrex::Error("Don't know this les_type");
        }

        // Include molecular transport? Default to no.
        static std::string molec_diff_type_string = "None";

        // If not running LES then default to "Constant"
        if (les_type == LESType::None) { // DNS
            molec_diff_type_string = "Constant"; // For DNS, use constant molecular viscosity
        }

        // But allow the user to reset to None
        pp.query("molec_diff_type",molec_diff_type_string);

        if (!molec_diff_type_string.compare("Constant")) {
            molec_diff_type = MolecDiffType::Constant;
        } else if (!molec_diff_type_string.compare("ConstantAlpha")) {
            molec_diff_type = MolecDiffType::ConstantAlpha;
        } else if (!molec_diff_type_string.compare("None")) {
            molec_diff_type = MolecDiffType::None;
        } else {
            amrex::Error("Don't know this molec_diff_type");
        }

        if ( (molec_diff_type == MolecDiffType::ConstantAlpha) &&
             (les_type        != LESType::None) ) {
            amrex::Error("We don't allow LES with MolecDiffType::ConstantAlpha");
        }

        // Which PBL Closure
        static std::string pbl_type_string = "None";
        pp.query("pbl_type",pbl_type_string);
        if (pbl_type_string == "MYNN2.5") {
            pbl_type = PBLType::MYNN25;
        } else if (pbl_type_string == "None") {
            pbl_type = PBLType::None;
        } else {
            amrex::Error("Don't know this pbl_type");
        }

        // Do some more stuff for PBL Modeling
        if (pbl_type != PBLType::None) {
            // Check for compatibility between PBL, LES, Molec Transport
            if (les_type != LESType::None) {
                amrex::Print() << "Selected a PBL model and an LES model: " <<
                "Using PBL for vertical transport, LES for horizontal" << std::endl;
            } else if (les_type == LESType::Deardorff) {
                amrex::Error("It is not recommended to use Deardorff LES and a PBL model");
            }
            pp.query("pbl_A1", pbl_A1);
            pp.query("pbl_A2", pbl_A2);
            pp.query("pbl_B1", pbl_B1);
            pp.query("pbl_B2", pbl_B2);
            pp.query("pbl_C1", pbl_C1);
            pp.query("pbl_C2", pbl_C2);
            pp.query("pbl_C3", pbl_C3);
            pp.query("pbl_C4", pbl_C4);
            pp.query("pbl_C5", pbl_C5);
        }

        // Right now, solving the QKE equation is only supported when MYNN PBL is turned on
        if (pbl_type == PBLType::MYNN25) use_QKE = true;
        if (use_QKE) {
            pp.query("diffuse_QKE_3D", diffuse_QKE_3D);
            pp.query("advect_QKE", advect_QKE);
        }

        // Diffusive/viscous/LES constants...
        pp.query("alpha_T", alpha_T);
        pp.query("alpha_C", alpha_C);
        pp.query("dynamicViscosity", dynamicViscosity);
        pp.query("Cs" , Cs);
        pp.query("CI" , CI);
        pp.query("Pr_t", Pr_t);
        pp.query("Sc_t", Sc_t);
        pp.query("rho0_trans", rho0_trans);

        if (molec_diff_type == MolecDiffType::ConstantAlpha) {
            amrex::Print() << "Using constant kinematic diffusion coefficients" << std::endl;
            amrex::Print() << "  momentum : " << dynamicViscosity/rho0_trans << " m^2/s" << std::endl;
            amrex::Print() << "  temperature : " << alpha_T << " m^2/s" << std::endl;
            amrex::Print() << "  scalar : " << alpha_C << " m^2/s" << std::endl;
        }

        pp.query("Ce" , Ce);
        pp.query("Ce_wall" , Ce_wall);
        pp.query("sigma_k" , sigma_k);

        if (les_type == LESType::Deardorff) {
            pp.query("Ck" , Ck);
        }

        pp.query("theta_ref", theta_ref);

        // Compute relevant forms of diffusion parameters
        Pr_t_inv = 1.0 / Pr_t;
        Sc_t_inv = 1.0 / Sc_t;
        rhoAlpha_T = rho0_trans * alpha_T;
        rhoAlpha_C = rho0_trans * alpha_C;

        // Order of spatial discretization
        pp.query("horiz_spatial_order", horiz_spatial_order);
        pp.query("vert_spatial_order",   vert_spatial_order);

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
        pp.queryarr("abl_geo_forcing",abl_geo_forcing_in);
        for(int i = 0; i < AMREX_SPACEDIM; ++i) abl_geo_forcing[i] = abl_geo_forcing_in[i];

        if (use_coriolis)
        {
            build_coriolis_forcings();
        }

        static std::string ic_bc_type_string = "Ideal";
        pp.query("ic_bc_type", ic_bc_type_string);

        pp.query("Ave_Plane", ave_plane);

#ifdef ERF_USE_MOISTURE
        pp.query("mp_clouds", do_cloud);
        pp.query("mp_precip", do_precip);
        moist_use_WENO = true;
#endif

        // Use WENO for advection?
        pp.query("all_use_WENO", all_use_WENO);
        if (all_use_WENO || moist_use_WENO) {
            pp.query("spatial_order_WENO", spatial_order_WENO);
            AMREX_ASSERT_WITH_MESSAGE(( (spatial_order_WENO==3) || (spatial_order_WENO==5) ),
                                      "WENO advection only supports orders 3 & 5.");
        }

        // Use numerical diffusion?
        pp.query("use_NumDiff",use_NumDiff);
        if(use_NumDiff) {
            pp.query("NumDiffCoeff",NumDiffCoeff);
            AMREX_ASSERT_WITH_MESSAGE(( (NumDiffCoeff >= 0.) && (NumDiffCoeff <= 1.) ),
                                      "Numerical diffusion coefficient must be between 0 & 1.");
            NumDiffCoeff *= std::pow(2.0,-6);
        }

    }

    void display()
    {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "use_coriolis          : " << use_coriolis << std::endl;
        amrex::Print() << "use_rayleigh_damping  : " << use_rayleigh_damping << std::endl;
        amrex::Print() << "use_gravity           : " << use_gravity << std::endl;
        amrex::Print() << "rho0_trans            : " << rho0_trans << std::endl;
        amrex::Print() << "alpha_T               : " << alpha_T << std::endl;
        amrex::Print() << "alpha_C               : " << alpha_C << std::endl;
        amrex::Print() << "dynamicViscosity      : " << dynamicViscosity << std::endl;
        amrex::Print() << "Cs                    : " << Cs << std::endl;
        amrex::Print() << "CI                    : " << CI << std::endl;
        amrex::Print() << "Ce                    : " << Ce << std::endl;
        amrex::Print() << "Ce at wall            : " << Ce_wall << std::endl;
        amrex::Print() << "Ck                    : " << Ck << std::endl;
        amrex::Print() << "reference theta       : " << theta_ref << std::endl;
        amrex::Print() << "sigma_k               : " << sigma_k << std::endl;
        amrex::Print() << "Pr_t                  : " << Pr_t << std::endl;
        amrex::Print() << "Sc_t                  : " << Sc_t << std::endl;
        amrex::Print() << "horiz spatial_order   : " << horiz_spatial_order << std::endl;
        amrex::Print() << "vert  spatial_order   : " << vert_spatial_order << std::endl;

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

        if (les_type == LESType::None) {
            amrex::Print() << "Using DNS model " << std::endl;
        } else if  (les_type == LESType::Smagorinsky) {
            amrex::Print() << "Using Smagorinsky LES model " << std::endl;
        } else if  (les_type == LESType::Deardorff) {
            amrex::Print() << "Using Deardorff LES model " << std::endl;
        }

        if (molec_diff_type == MolecDiffType::Constant) {
            amrex::Print() << "Using constant molecular diffusivity (relevant for DNS)" << std::endl;
        } else if (molec_diff_type == MolecDiffType::None) {
            amrex::Print() << "Not using any molecular diffusivity, i.e. using the modeled turbulent diffusivity"
            << std::endl;
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

    // Default prefix
    std::string pp_prefix {"erf"};

    bool        use_terrain            = false;
    bool        test_mapfactor         = false;
    int         terrain_type           = 0;
#ifdef ERF_USE_MOISTURE
    int         buoyancy_type          = 2; // uses Tprime
#else
    int         buoyancy_type          = 1; // uses rhoprime directly
#endif

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
    amrex::Real c_p = 1004.5; // specific heat at constant pressure for dry air [J/(kg-K)]
    amrex::Real rdOcp;

#if defined(ERF_USE_WARM_NO_PRECIP)
    amrex::Real tau_cond = 1.0; // Default time of 1 sec -- this is somewhat arbitray
#endif

#if defined(ERF_USE_POISSON_SOLVE)
    int project_initial_velocity = 1;
#endif

    // Molecular transport model
    MolecDiffType molec_diff_type;
    // Diffusive/viscous coefficients [m2/s]
    amrex::Real alpha_T = 0.0;
    amrex::Real alpha_C = 0.0;
    // Density for computation of rhoAlpha (which is assumed constant) [kg/m3]
    amrex::Real rho0_trans = 1.0;
    // Dynamic diffusion coefficients [kg/(m-s)]
    amrex::Real rhoAlpha_T = 0.0;
    amrex::Real rhoAlpha_C = 0.0;
    amrex::Real dynamicViscosity = 0.0;

    // LES model
    LESType les_type;
    // Smagorinsky Cs coefficient
    amrex::Real Cs = 0.0;
    // Smagorinsky CI coefficient
    amrex::Real CI = 0.0;
    // Smagorinsky Turbulent Prandtl Number
    amrex::Real Pr_t = 1. / 3.;
    amrex::Real Pr_t_inv = 3.0;
    // Smagorinsky Turbulent Schmidt Number
    amrex::Real Sc_t = 1.0;
    amrex::Real Sc_t_inv = 1.0;

    // Deardorff Ce coefficient
    amrex::Real Ce = 0.93;
    amrex::Real Ce_wall = 0.0; // if > 0, then set Ce to this at k=0
    // Deardorff Ck coefficient
    amrex::Real Ck = 0.1;
    // Deardorff sigma_k coefficient
    amrex::Real sigma_k = 0.5;

    amrex::Real theta_ref = 300.0;

    // PBL model
    PBLType pbl_type;
    // Model coefficients
    amrex::Real pbl_A1 = 1.18;
    amrex::Real pbl_A2 = 0.665;
    amrex::Real pbl_B1 = 24.0;
    amrex::Real pbl_B2 = 15.0;
    amrex::Real pbl_C1 = 0.137;
    amrex::Real pbl_C2 = 0.75;
    amrex::Real pbl_C3 = 0.352;
    amrex::Real pbl_C4 = 0.0;
    amrex::Real pbl_C5 = 0.2;
    // QKE stuff - default is not to use it, if MYNN2.5 PBL is used default is turb transport in Z-direction only
    bool use_QKE = false;
    bool diffuse_QKE_3D = false;
    bool advect_QKE = true;

    // Coriolis forcing
    amrex::Real coriolis_factor = 0.0;
    amrex::Real cosphi          = 0.0 ;
    amrex::Real sinphi          = 0.0;

    // Spatial discretization
    int   horiz_spatial_order = 2;
    int    vert_spatial_order = 2;

    // Positive definite advection scheme (3/5 order)
    bool all_use_WENO{false};
    bool moist_use_WENO{false};
    int  spatial_order_WENO{3};

    // Numerical diffusion
    bool use_NumDiff{false};
    amrex::Real NumDiffCoeff{0.};

    ABLDriverType abl_driver_type;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_forcing;

    int ave_plane {2};
#ifdef ERF_USE_MOISTURE
    // Microphysics params
    bool do_cloud {true};
    bool do_precip {true};
#endif
};
#endif
