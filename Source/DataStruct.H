#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

#include <ERF_Constants.H>

enum class ABLDriverType {
    None, PressureGradient, GeostrophicWind
};

enum class Coord {
    x, y, z
};

enum class Stagger {
    None, x, y, z
};

enum class AdvectedQuantity {
    unity, u, v, w, theta, scalar
};

enum class AdvectingQuantity {
    rho_u, rho_v, rho_w
};

enum class AdvectionDir {
    x, y, z
};

enum class DiffusionDir {
    x, y, z
};

enum class MomentumEqn {
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

enum class IC_BC_Type {
    Ideal, Real
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
        pp.query("Ck" , Ck);
        pp.query("sigma_k" , sigma_k);

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

        if (!ic_bc_type_string.compare("Ideal")) {
            icBcType = IC_BC_Type::Ideal;
        } else if (!ic_bc_type_string.compare("Real")) {
            icBcType = IC_BC_Type::Real;
        } else {
            amrex::Error("Don't know this ic_bc_type");
        }

        pp.query("Ave_Plane", ave_plane);

#ifdef ERF_USE_MOISTURE
        pp.query("mp_clouds", do_cloud);
        pp.query("mp_precip", do_precip);
#endif
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
        amrex::Print() << "Ck                    : " << Ck << std::endl;
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

        if (icBcType == IC_BC_Type::Ideal) {
            amrex::Print() << "Using idealized initial and boundary conditions (No mesoscale forcing!)" << std::endl;
        }
        else if (icBcType == IC_BC_Type::Real) {
            amrex::Print() << "Using REAL initial and boundary conditions (Mesoscale forcing!)" << std::endl;
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
    amrex::Real Pr_t = 1.0;
    amrex::Real Pr_t_inv = 1.0;
    // Smagorinsky Turbulent Schmidt Number
    amrex::Real Sc_t = 1.0;
    amrex::Real Sc_t_inv = 1.0;

    // Deardorff Ce coefficient
    amrex::Real Ce = 0.0;
    // Deardorff Ck coefficient
    amrex::Real Ck = 0.0;
    // Deardorff sigma_k coefficient
    amrex::Real sigma_k = 1.0;

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
    // QKE stuff - default is not to use it, if MYNN2.5 PBL is used default is turb transport in Z-directioon only
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

    ABLDriverType abl_driver_type;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_forcing;

    // IC and BC Type: "Ideal" or "Real"
    IC_BC_Type icBcType;

    int ave_plane {2};
#ifdef ERF_USE_MOISTURE
    // Microphysics params
    bool do_cloud {true};
    bool do_precip {true};
#endif
};

struct InputSoundingData {
public:
    InputSoundingData() {}

    void read_from_file(const std::string input_sounding_file,
                        const amrex::Real ztop)
    {
        z_inp_sound.resize(0);
        theta_inp_sound.resize(0);
        qv_inp_sound.resize(0);
        U_inp_sound.resize(0);
        V_inp_sound.resize(0);

        pm_integ.resize(0);
        rhom_integ.resize(0);

        pd_integ.resize(0);
        rhod_integ.resize(0);

        // Read the input_sounding file
        amrex::Print() << "input_sounding file location : " << input_sounding_file << std::endl;
        std::ifstream input_sounding_reader(input_sounding_file);
        if(!input_sounding_reader.is_open()) {
            amrex::Error("Error opening the input_sounding file\n");
        }
        else {
            // Read the contents of the input_sounding file
            amrex::Print() << "Successfully opened the input_sounding file. Now reading... " << std::endl;
            std::string  line;

            // Read the first line
            std::getline(input_sounding_reader, line);
            std::istringstream iss(line);
            iss >> press_ref_inp_sound >> theta_ref_inp_sound >> qv_ref_inp_sound;
            qv_ref_inp_sound *= 0.001;

            // Add surface
            z_inp_sound.push_back(0);
            theta_inp_sound.push_back(theta_ref_inp_sound);
            qv_inp_sound.push_back(qv_ref_inp_sound);
            U_inp_sound.push_back(0);
            V_inp_sound.push_back(0);

            // Read the vertical profile at each given height
            amrex::Real z, theta, qv, U, V;
            while(std::getline(input_sounding_reader, line)) {
                std::istringstream iss_z(line);
                iss_z >> z >> theta >> qv >> U >> V;
                z_inp_sound.push_back(z);
                theta_inp_sound.push_back(theta);
                qv_inp_sound.push_back(qv*0.001);
                U_inp_sound.push_back(U);
                V_inp_sound.push_back(V);
                if (z >= ztop) break;
            }

            // At this point, we have an input_sounding with
            // z_inp_sound[N-1] <= ztop
            if (z != ztop)
            {
                int N = size();
                amrex::Real z1  =     z_inp_sound[N-2];
                amrex::Real th1 = theta_inp_sound[N-2];
                amrex::Real qv1 =    qv_inp_sound[N-2];
                amrex::Real U1  =     U_inp_sound[N-2];
                amrex::Real V1  =     V_inp_sound[N-2];
                amrex::Real z2  =     z_inp_sound[N-1];
                amrex::Real th2 = theta_inp_sound[N-1];
                amrex::Real qv2 =    qv_inp_sound[N-1];
                amrex::Real U2  =     U_inp_sound[N-1];
                amrex::Real V2  =     V_inp_sound[N-1];

                amrex::Real zfrac = (ztop - z1) / (z2 - z1);
                amrex::Real th_top = th1 + zfrac*(th2 - th1);
                amrex::Real qv_top = qv1 + zfrac*(qv2 - qv1);
                amrex::Real  U_top =  U1 + zfrac*( U2 -  U1);
                amrex::Real  V_top =  V1 + zfrac*( V2 -  V1);

                if (z < ztop)
                {
                    // Add domain top with extrapolated values
                    z_inp_sound.push_back(ztop);
                    theta_inp_sound.push_back(th_top);
                    qv_inp_sound.push_back(qv_top);
                    U_inp_sound.push_back(U_top);
                    V_inp_sound.push_back(V_top);
                }
                else // z > ztop
                {
                    // Overwrite top value with interpolated values
                    z_inp_sound[N-1] = ztop;
                    theta_inp_sound[N-1] = th_top;
                    qv_inp_sound[N-1] = qv_top;
                    U_inp_sound[N-1] = U_top;
                    V_inp_sound[N-1] = V_top;
                }
            }
        }
        amrex::Print() << "Successfully read the input_sounding file..." << std::endl;
        input_sounding_reader.close();

        host_to_device();
    }

    void calc_rho_p(const amrex::Real ztop)
    {
        /* Calculate density and pressure, roughly following the procedure in
         * WRF dyn_em/module_initialize_ideal.F. We integrate hydrostatically
         * from the surface up through the air column to get the total (moist)
         * density and pressure, then back down to get the dry density and
         * pressure.
         *
         * This deviates from that code slightly because we calculate the
         * correct moist theta, equal to virtual potential temperature, and
         * also remove the factor (1+qv) for consistency with the the surface
         * density.
         */
        const int maxiter = 10;
        const int Ninp = size();
        pm_integ.resize(Ninp);
        pd_integ.resize(Ninp);
        rhom_integ.resize(Ninp);
        rhod_integ.resize(Ninp);

        // evaluate surface quantities (k=0)
        amrex::Real thm_surf = theta_ref_inp_sound *
            (1 + (R_v/R_d-1) * qv_ref_inp_sound); // _moist_ theta == virtual potential temperature
        pm_integ[0] = 100. * press_ref_inp_sound; // _total_ pressure (incl moisture)
        rhom_integ[0] = 1 / (
            R_d/p_0 * thm_surf * std::pow(pm_integ[0]/p_0, -iGamma)); // density of _moist_ air

        amrex::Print() << "ideal sounding init: surface density of moist air = "
            << rhom_integ[0] << std::endl;
        amrex::Print() << "integrating moist air column up to " << ztop << std::endl;

        // Note:
        //   p_dry = rho_d R_d T
        //   p_tot = rho_m R_d T_v
        //         = rho_d(1 + q_v) R_d T_v

        // integrate from surface to domain top
        amrex::Real qvf, dz;
        for (int k=1; k < size(); ++k)
        {
            qvf = 1 + (R_v/R_d-1) * qv_inp_sound[k];
            dz = z_inp_sound[k] - z_inp_sound[k-1];
            rhom_integ[k] = rhom_integ[k-1]; // guess
            for (int it=0; it < maxiter; ++it)
            {
                pm_integ[k] = pm_integ[k-1]
                    - 0.5*dz*(rhom_integ[k] + rhom_integ[k-1])*CONST_GRAV;
                AMREX_ALWAYS_ASSERT(pm_integ[k] > 0);
                rhom_integ[k] = 1 / (
                    R_d/p_0 * theta_inp_sound[k]*qvf * std::pow(pm_integ[k]/p_0, -iGamma));
            }
        }

        // now, integrate from the top of the sounding (where it's dry) back
        // down to get the dry air column properties
        pd_integ[Ninp-1] = pm_integ[Ninp-1];
        rhod_integ[Ninp-1] = 1 / (
                    R_d/p_0 * theta_inp_sound[Ninp-1] * std::pow(pd_integ[Ninp-1]/p_0, -iGamma));
        amrex::Print() << "z  p_d  rho_d  theta" << std::endl;
        amrex::Print() << z_inp_sound[Ninp-1] << " " << pd_integ[Ninp-1] << " " << rhod_integ[Ninp-1] << " " << theta_inp_sound[Ninp-1] << std::endl; // DEBUG
        for (int k=Ninp-2; k >= 0; --k)
        {
            dz = z_inp_sound[k+1] - z_inp_sound[k];
            rhod_integ[k] = rhod_integ[k+1]; // guess
            for (int it=0; it < maxiter; ++it)
            {
                pd_integ[k] = pd_integ[k+1]
                    + 0.5*dz*(rhod_integ[k] + rhod_integ[k+1])*CONST_GRAV;
                AMREX_ALWAYS_ASSERT(pd_integ[k] > 0);
                rhod_integ[k] = 1 / (
                    R_d/p_0 * theta_inp_sound[k] * std::pow(pd_integ[k]/p_0, -iGamma));
            }
            amrex::Print() << z_inp_sound[k] << " " << pd_integ[k] << " " << rhod_integ[k] << " " << theta_inp_sound[k] << std::endl; // DEBUG
        }
        // Note: at this point, the surface pressure, density of the dry air
        // column is stored in pd_integ[0], rhod_integ[0]

        // update
        host_to_device();
    }

    void host_to_device()
    {
        z_inp_sound_d.resize(size());
        theta_inp_sound_d.resize(size());
        qv_inp_sound_d.resize(size());
        U_inp_sound_d.resize(size());
        V_inp_sound_d.resize(size());

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         z_inp_sound.begin(), z_inp_sound.end(),
                         z_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         theta_inp_sound.begin(), theta_inp_sound.end(),
                         theta_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         qv_inp_sound.begin(), qv_inp_sound.end(),
                         qv_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         U_inp_sound.begin(), U_inp_sound.end(),
                         U_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         V_inp_sound.begin(), V_inp_sound.end(),
                         V_inp_sound_d.begin());

        if (rhod_integ.size() > 0)
        {
            //amrex::Print() << "Copying rho_d, p_d to device" << std::endl;
            rho_inp_sound_d.resize(size()+2);
            p_inp_sound_d.resize(size()+2);
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             rhod_integ.begin(), rhod_integ.end(),
                             rho_inp_sound_d.begin());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             pd_integ.begin(), pd_integ.end(),
                             p_inp_sound_d.begin());
        }
    }

    int size() const
    {
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == theta_inp_sound.size());
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == qv_inp_sound.size());
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == U_inp_sound.size());
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == V_inp_sound.size());
        return z_inp_sound.size();
    }

    // Members
    // - read from file
    amrex::Real press_ref_inp_sound, theta_ref_inp_sound, qv_ref_inp_sound; // input
    amrex::Vector<amrex::Real> z_inp_sound, theta_inp_sound, qv_inp_sound, U_inp_sound, V_inp_sound;
    amrex::Gpu::DeviceVector<amrex::Real> z_inp_sound_d, theta_inp_sound_d, qv_inp_sound_d, U_inp_sound_d, V_inp_sound_d;
    // - moist profiles
    amrex::Vector<amrex::Real> pm_integ, rhom_integ; // from integrating up air column
    // - dry profiles
    amrex::Vector<amrex::Real> pd_integ, rhod_integ; // from integrating down air column
    // - to set solution fields
    amrex::Gpu::DeviceVector<amrex::Real> p_inp_sound_d, rho_inp_sound_d;
};
#endif
