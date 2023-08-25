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

/**
 * Container holding many of the algorithmic options and parameters
 */

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
        if (!use_lagged_delta_rt && !(terrain_type == 1)) {
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

        // Turn off acoustic substepping?
        pp.query("no_substepping", no_substepping);

        pp.query("force_stage1_single_substep", force_stage1_single_substep);
        pp.query("incompressible", incompressible);

        // If this is set, it must be even
        if (incompressible != 0 && no_substepping == 0)
        {
            amrex::Abort("If you specify incompressible, you must specific no_substepping");
        }

        // Order and type of spatial discretizations used in advection
        pp.query("use_efficient_advection", use_efficient_advection);
        std::string dycore_horiz_adv_string    = "" ; std::string dycore_vert_adv_string   = "";
        std::string dryscal_horiz_adv_string   = "" ; std::string dryscal_vert_adv_string  = "";
        pp.query("dycore_horiz_adv_type"   , dycore_horiz_adv_string);
        pp.query("dycore_vert_adv_type"    , dycore_vert_adv_string);
        pp.query("dryscal_horiz_adv_type"  , dryscal_horiz_adv_string);
        pp.query("dryscal_vert_adv_type"   , dryscal_vert_adv_string);

#if defined(ERF_USE_MOISTURE) or defined(ERF_USE_WARM_NO_PRECIP)
        std::string moistscal_horiz_adv_string = ""; std::string moistscal_vert_adv_string = "";
        pp.query("moistscal_horiz_adv_type", moistscal_horiz_adv_string);
        pp.query("moistscal_vert_adv_type" , moistscal_vert_adv_string);
#endif

        if (use_efficient_advection){
           amrex::Print() << "Using efficient advection scheme" << std::endl;;
        }

        if ( (dycore_horiz_adv_string == "Centered_2nd") ||
             (dycore_horiz_adv_string == "Upwind_3rd"  ) ||
             (dycore_horiz_adv_string == "Centered_4th") ||
             (dycore_horiz_adv_string == "Upwind_5th"  ) ||
             (dycore_horiz_adv_string == "Centered_6th") )
        {
            dycore_horiz_adv_type = adv_type_convert_string_to_advtype(dycore_horiz_adv_string);
            amrex::Print() << "Using dycore_horiz_adv_type: " << dycore_horiz_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dycore_horiz_adv_type" << std::endl;;
        }

        if ( (dycore_vert_adv_string == "Centered_2nd") ||
             (dycore_vert_adv_string == "Upwind_3rd"  ) ||
             (dycore_vert_adv_string == "Centered_4th") ||
             (dycore_vert_adv_string == "Upwind_5th"  ) ||
             (dycore_vert_adv_string == "Centered_6th") )
        {
            dycore_vert_adv_type = adv_type_convert_string_to_advtype(dycore_vert_adv_string);
            amrex::Print() << "Using dycore_vert_adv_type: " << dycore_vert_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dycore_vert_adv_type" << std::endl;;
        }

        if ( (dryscal_horiz_adv_string == "Centered_2nd") ||
             (dryscal_horiz_adv_string == "Upwind_3rd"  ) ||
             (dryscal_horiz_adv_string == "Centered_4th") ||
             (dryscal_horiz_adv_string == "Upwind_5th"  ) ||
             (dryscal_horiz_adv_string == "Centered_6th") ||
             (dryscal_horiz_adv_string == "WENO3"       ) ||
             (dryscal_horiz_adv_string == "WENOZ3"      ) ||
             (dryscal_horiz_adv_string == "WENOMZQ3"    ) ||
             (dryscal_horiz_adv_string == "WENO5"       ) ||
             (dryscal_horiz_adv_string == "WENOZ5"      ) )
        {
            dryscal_horiz_adv_type = adv_type_convert_string_to_advtype(dryscal_horiz_adv_string);
            amrex::Print() << "Using dryscal_horiz_adv_type: " << dryscal_horiz_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dryscal_horiz_adv_type" << std::endl;;
        }

        if ( (dryscal_vert_adv_string == "Centered_2nd") ||
             (dryscal_vert_adv_string == "Upwind_3rd"  ) ||
             (dryscal_vert_adv_string == "Centered_4th") ||
             (dryscal_vert_adv_string == "Upwind_5th"  ) ||
             (dryscal_vert_adv_string == "Centered_6th") ||
             (dryscal_vert_adv_string == "WENO3"       ) ||
             (dryscal_vert_adv_string == "WENOZ3"      ) ||
             (dryscal_vert_adv_string == "WENOMZQ3"     ) ||
             (dryscal_vert_adv_string == "WENO5"       ) ||
             (dryscal_vert_adv_string == "WENOZ5"      ) )
        {
            dryscal_vert_adv_type = adv_type_convert_string_to_advtype(dryscal_vert_adv_string);
            amrex::Print() << "Using dryscal_vert_adv_type: " << dryscal_vert_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dryscal_vert_adv_type" << std::endl;;
        }

#if defined(ERF_USE_MOISTURE) or defined(ERF_USE_WARM_NO_PRECIP)
        if ( (moistscal_horiz_adv_string == "Centered_2nd") ||
             (moistscal_horiz_adv_string == "Upwind_3rd"  ) ||
             (moistscal_horiz_adv_string == "Centered_4th") ||
             (moistscal_horiz_adv_string == "Upwind_5th"  ) ||
             (moistscal_horiz_adv_string == "Centered_6th") ||
             (moistscal_horiz_adv_string == "WENO3"       ) ||
             (moistscal_horiz_adv_string == "WENOZ3"      ) ||
             (moistscal_horiz_adv_string == "WENOMZQ3"    ) ||
             (moistscal_horiz_adv_string == "WENO5"       ) ||
             (moistscal_horiz_adv_string == "WENOZ5"      ) )
        {
            moistscal_horiz_adv_type = adv_type_convert_string_to_advtype(moistscal_horiz_adv_string);
            amrex::Print() << "Using moistscal_horiz_adv_type: " << moistscal_horiz_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default moistscal_horiz_adv_type" << std::endl;;
        }

        if ( (moistscal_vert_adv_string == "Centered_2nd") ||
             (moistscal_vert_adv_string == "Upwind_3rd"  ) ||
             (moistscal_vert_adv_string == "Centered_4th") ||
             (moistscal_vert_adv_string == "Upwind_5th"  ) ||
             (moistscal_vert_adv_string == "Centered_6th") ||
             (moistscal_vert_adv_string == "WENO3"       ) ||
             (moistscal_vert_adv_string == "WENOZ3"      ) ||
             (moistscal_vert_adv_string == "WENOMZQ3"    ) ||
             (moistscal_vert_adv_string == "WENO5"       ) ||
             (moistscal_vert_adv_string == "WENOZ5"      ) )
        {
            moistscal_vert_adv_type = adv_type_convert_string_to_advtype(moistscal_vert_adv_string);
            amrex::Print() << "Using moistscal_vert_adv_type: " << moistscal_vert_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default moistscal_vert_adv_type" << std::endl;;
        }
#endif

        // Include Coriolis forcing?
        pp.query("use_coriolis", use_coriolis);

        // Include Rayleigh damping?
        pp.query("use_rayleigh_damping", use_rayleigh_damping);
        pp.query("rayleigh_damp_U", rayleigh_damp_U);
        pp.query("rayleigh_damp_V", rayleigh_damp_V);
        pp.query("rayleigh_damp_W", rayleigh_damp_W);
        pp.query("rayleigh_damp_T", rayleigh_damp_T);

        // Include sponge damping
        pp.query("use_xlo_sponge_damping", use_xlo_sponge_damping);
        pp.query("use_xhi_sponge_damping", use_xhi_sponge_damping);
        pp.query("use_ylo_sponge_damping", use_ylo_sponge_damping);
        pp.query("use_yhi_sponge_damping", use_yhi_sponge_damping);
        pp.query("use_zlo_sponge_damping", use_zlo_sponge_damping);
        pp.query("use_zhi_sponge_damping", use_zhi_sponge_damping);

        pp.query("sponge_strength" , sponge_strength);
        pp.query("xlo_sponge_end"  , xlo_sponge_end);
        pp.query("xhi_sponge_start", xhi_sponge_start);
        pp.query("ylo_sponge_end"  , ylo_sponge_end);
        pp.query("yhi_sponge_start", yhi_sponge_start);
        pp.query("zlo_sponge_end"  , zlo_sponge_end);
        pp.query("zhi_sponge_start", zhi_sponge_start);

        pp.query("sponge_density" , sponge_density);
        pp.query("sponge_x_velocity" , sponge_x_velocity);
        pp.query("sponge_y_velocity" , sponge_y_velocity);
        pp.query("sponge_z_velocity" , sponge_z_velocity);

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

        pp.query("Ave_Plane", ave_plane);

#ifdef ERF_USE_MOISTURE
        pp.query("mp_clouds", do_cloud);
        pp.query("mp_precip", do_precip);
#endif

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
        amrex::Print() << "no_substepping              : " << no_substepping << std::endl;
        amrex::Print() << "force_stage1_single_substep : "  << force_stage1_single_substep << std::endl;
        amrex::Print() << "incompressible              : "  << incompressible << std::endl;
        amrex::Print() << "use_coriolis                : " << use_coriolis << std::endl;
        amrex::Print() << "use_rayleigh_damping        : " << use_rayleigh_damping << std::endl;
        amrex::Print() << "use_gravity                 : " << use_gravity << std::endl;
        amrex::Print() << "rho0_trans                  : " << rho0_trans << std::endl;
        amrex::Print() << "alpha_T                     : " << alpha_T << std::endl;
        amrex::Print() << "alpha_C                     : " << alpha_C << std::endl;
        amrex::Print() << "dynamicViscosity            : " << dynamicViscosity << std::endl;
        amrex::Print() << "Cs                          : " << Cs << std::endl;
        amrex::Print() << "CI                          : " << CI << std::endl;
        amrex::Print() << "Ce                          : " << Ce << std::endl;
        amrex::Print() << "Ce at wall                  : " << Ce_wall << std::endl;
        amrex::Print() << "Ck                          : " << Ck << std::endl;
        amrex::Print() << "reference theta             : " << theta_ref << std::endl;
        amrex::Print() << "sigma_k                     : " << sigma_k << std::endl;
        amrex::Print() << "Pr_t                        : " << Pr_t << std::endl;
        amrex::Print() << "Sc_t                        : " << Sc_t << std::endl;
        amrex::Print() << "dycore_horiz_adv_type       : " << adv_type_convert_int_to_string(dycore_horiz_adv_type) << std::endl;
        amrex::Print() << "dycore_vert_adv_type        : " << adv_type_convert_int_to_string(dycore_vert_adv_type) << std::endl;
        amrex::Print() << "dryscal_horiz_adv_type      : " << adv_type_convert_int_to_string(dryscal_horiz_adv_type) << std::endl;
        amrex::Print() << "dryscal_vert_adv_type       : " << adv_type_convert_int_to_string(dryscal_vert_adv_type) << std::endl;
#if defined(ERF_USE_MOISTURE) or defined(ERF_USE_WARM_NO_PRECIP)
        amrex::Print() << "moistscal_horiz_adv_type    : " << adv_type_convert_int_to_string(moistscal_horiz_adv_type) << std::endl;
        amrex::Print() << "moistscal_vert_adv_type     : " << adv_type_convert_int_to_string(moistscal_vert_adv_type) << std::endl;
#endif

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

    std::string
    adv_type_convert_int_to_string (AdvType adv_int)
    {
        if (adv_int == AdvType::Centered_2nd) {
            return "Centered_2nd";
        } else if (adv_int == AdvType::Upwind_3rd) {
            return "Upwind_3rd";
        } else if (adv_int == AdvType::Centered_4th) {
            return "Centered_4th";
        } else if (adv_int == AdvType::Upwind_5th) {
            return "Upwind_5th";
        } else if (adv_int == AdvType::Centered_6th) {
            return "Centered_6th";
        } else if (adv_int == AdvType::Weno_3) {
            return "WENO3";
        } else if (adv_int == AdvType::Weno_3Z) {
            return "WENOZ3";
        } else if (adv_int == AdvType::Weno_5) {
            return "WENO5";
        } else if (adv_int == AdvType::Weno_5Z) {
            return "WENOZ5";
        } else if (adv_int == AdvType::Weno_3MZQ) {
            return "WENOMZQ3";
        } else {
            return "Unknown";
        }
    }


    AdvType adv_type_convert_string_to_advtype (std::string adv_string)
    {
        if (adv_string == "Centered_2nd") {
            return AdvType::Centered_2nd;
        } else if (adv_string == "Upwind_3rd") {
            return AdvType::Upwind_3rd;
        } else if (adv_string == "Centered_4th") {
            return AdvType::Centered_4th;
        } else if (adv_string == "Upwind_5th") {
            return AdvType::Upwind_5th;
        } else if (adv_string == "Centered_6th") {
            return AdvType::Centered_6th;
        } else if (adv_string == "WENO3") {
            return AdvType::Weno_3;
        } else if (adv_string == "WENOZ3") {
            return AdvType::Weno_3Z;
        } else if (adv_string == "WENO5") {
            return AdvType::Weno_5;
        } else if (adv_string == "WENOZ5") {
            return AdvType::Weno_5Z;
        } else if (adv_string == "WENOMZQ3") {
            return AdvType::Weno_3MZQ;
        } else {
            return AdvType::Unknown;
        }
    }

    // Default prefix
    std::string pp_prefix {"erf"};

    int         no_substepping              = 0;
    int         force_stage1_single_substep = 1;
    int         incompressible              = 0;

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

    bool    use_xlo_sponge_damping     = false;
    bool    use_xhi_sponge_damping     = false;
    bool    use_ylo_sponge_damping     = false;
    bool    use_yhi_sponge_damping     = false;
    bool    use_zlo_sponge_damping     = false;
    bool    use_zhi_sponge_damping     = false;
    amrex::Real sponge_strength;
    amrex::Real xlo_sponge_end, xhi_sponge_start;
    amrex::Real ylo_sponge_end, yhi_sponge_start;
    amrex::Real zlo_sponge_end, zhi_sponge_start;
    amrex::Real sponge_density, sponge_x_velocity, sponge_y_velocity, sponge_z_velocity;

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
    // Order and type of spatial discretizations used in advection
    bool use_efficient_advection = false;
    AdvType dycore_horiz_adv_type    = AdvType::Centered_2nd;
    AdvType dycore_vert_adv_type     = AdvType::Centered_2nd;
    AdvType dryscal_horiz_adv_type   = AdvType::Centered_2nd;
    AdvType dryscal_vert_adv_type    = AdvType::Centered_2nd;
#if defined(ERF_USE_MOISTURE) or defined(ERF_USE_WARM_NO_PRECIP)
    AdvType moistscal_horiz_adv_type = AdvType::Weno_3;
    AdvType moistscal_vert_adv_type  = AdvType::Weno_3;
#endif

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
