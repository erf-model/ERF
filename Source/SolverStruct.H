#ifndef _SOLVER_STRUCT_H_
#define _SOLVER_STRUCT_H_

#include <ERF_Constants.H>
#include <AMReX_ParmParse.H>
#include <DataStruct.H>

struct SolverChoice {
  public:
    void init_params()
    {
        amrex::ParmParse pp("erf");

        // These default to true but are used for unit testing
        pp.query("use_state_advection", use_state_advection);
        pp.query("use_momentum_advection", use_momentum_advection);
        pp.query("use_thermal_diffusion", use_thermal_diffusion);
        pp.query("use_scalar_diffusion", use_scalar_diffusion);
        pp.query("use_momentum_diffusion", use_momentum_diffusion);
        pp.query("use_pressure", use_pressure);
        pp.query("use_gravity", use_gravity);

        gravity = use_gravity? CONST_GRAV: 0.0;

        // Which LES closure?
        static std::string les_type_string = "None";
        pp.query("les_type",les_type_string);
        if (!les_type_string.compare("Smagorinsky")) {
            les_type = LESType::Smagorinsky;
        } else if (!les_type_string.compare("Deardorff")) {
            les_type = LESType::Deardorff;
        } else if (!les_type_string.compare("None")) {
            les_type = LESType::None;
        } else {
            amrex::Error("Don't know this les_type");
        }

        if (les_type == LESType::Deardorff)
        {
            amrex::Error("Deardorff LES model not yet implemented");
        }

        // Diffusive/viscous/LES constants...
        pp.query("alpha_T", alpha_T);
        pp.query("alpha_C", alpha_C);
        pp.query("dynamicViscosity", dynamicViscosity);
        pp.query("Cs" , Cs);
        pp.query("CI" , CI);
        pp.query("PrT", PrT);

        // Order of spatial discretization
        pp.query("spatial_order", spatial_order);

        // Include Coriolis forcing?
        pp.query("use_coriolis", use_coriolis);

        // Which external forcings?
        static std::string abl_driver_type_string = "None";
        pp.query("abl_driver_type",abl_driver_type_string);
        if (!abl_driver_type_string.compare("PressureGradient")) {
            abl_driver_type = ABLDriverType::PressureGradient;
        } else if (!abl_driver_type_string.compare("GeostrophicWind")) {
            abl_driver_type = ABLDriverType::GeostrophicWind;
        } else {
            abl_driver_type = ABLDriverType::None;
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
    }

    void display()
    {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "use_state_advection   : " << use_state_advection << std::endl;
        amrex::Print() << "use_momentum_advection: " << use_momentum_advection << std::endl;
        amrex::Print() << "use_thermal_diffusion : " << use_thermal_diffusion << std::endl;
        amrex::Print() << "use_scalar_diffusion  : " << use_scalar_diffusion << std::endl;
        amrex::Print() << "use_momentum_diffusion: " << use_momentum_diffusion << std::endl;
        amrex::Print() << "use_pressure          : " << use_pressure << std::endl;
        amrex::Print() << "use_coriolis          : " << use_coriolis << std::endl;
        amrex::Print() << "use_gravity           : " << use_gravity << std::endl;
        amrex::Print() << "alpha_T               : " << alpha_T << std::endl;
        amrex::Print() << "alpha_C               : " << alpha_C << std::endl;
        amrex::Print() << "dynamicViscosity      : " << dynamicViscosity << std::endl;
        amrex::Print() << "Cs                    : " << Cs << std::endl;
        amrex::Print() << "CI                    : " << CI << std::endl;
        amrex::Print() << "PrT                   : " << PrT << std::endl;
        amrex::Print() << "spatial_order         : " << spatial_order << std::endl;

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
            amrex::Print() << "Using Smagorinsky model " << std::endl;
        } else if  (les_type == LESType::Deardorff) {
            amrex::Print() << "Using Deardorff model " << std::endl;
        }
    }

    void build_coriolis_forcings()
    {
        amrex::ParmParse pp("erf");

        // Read the rotational time period (in seconds)
        amrex::Real rot_time_period = 86400.0;
        pp.query("rotational_time_period", rot_time_period);

        coriolis_factor = 2.0 * 2.0 * PI / rot_time_period;
        amrex::Print() << "Coriolis factor = " << coriolis_factor << std::endl;

        amrex::Real latitude = 90.0;
        pp.query("latitude", latitude);
        AMREX_ALWAYS_ASSERT(amrex::Math::abs(latitude - 90.0) < 1.e-12);

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

    // These default to true but are used for unit testing
    bool        use_state_advection    = true;
    bool        use_momentum_advection = true;

    bool        use_thermal_diffusion  = true;
    bool        use_scalar_diffusion   = true;
    bool        use_momentum_diffusion = true;
    bool        use_pressure           = true;

    bool        use_gravity            = false;
    bool        use_coriolis           = false;

    // Diffusive/viscous coefficients
    amrex::Real alpha_T = 0.0;
    amrex::Real alpha_C= 0.0;
    amrex::Real dynamicViscosity = 0.0;

    // This will be set to 0 or CONST_GRAV, depending on the use_gravity flag
    amrex::Real gravity;

    // LES model
    LESType les_type;
    // Smagorinsky Cs coefficient
    amrex::Real Cs = 0.0;
    // Smagorinsky CI coefficient
    amrex::Real CI = 0.0;
    // Smagorinsky PrT coefficient
    amrex::Real PrT = 0.0;

    // Coriolis forcing
    amrex::Real coriolis_factor = 0.0;
    amrex::Real cosphi          = 0.0 ;
    amrex::Real sinphi          = 0.0;

    // Spatial discretization
    int         spatial_order = 2;

    ABLDriverType abl_driver_type;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_forcing;
};

#endif
