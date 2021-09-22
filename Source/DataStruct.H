#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

struct SolverChoice {
public:
    enum class abl_driver {
      None, PressureGradient, GeostrophicWind
    };

    //SolverChoice() {}; //Default
    SolverChoice(const bool use_state_advection_in, const bool use_momentum_advection_in,
                 bool  use_thermal_diffusion_in, const amrex::Real alpha_T_in,
                 const bool use_scalar_diffusion_in, const amrex::Real alpha_S_in,
                 const bool use_momentum_diffusion_in, const amrex::Real dynamicViscosity_in,
                 const bool use_smagorinsky_in, const amrex::Real Cs_in,
                 const bool use_pressure_in, const bool use_gravity_in, const int spatial_order_in,
                 const std::string &abl_driver_type_in,
                 const std::vector<amrex::Real> &abl_pressure_grad_in,
                 const std::vector<amrex::Real> &abl_geo_wind_in):
            use_state_advection(use_state_advection_in), use_momentum_advection(use_momentum_advection_in),
            use_thermal_diffusion(use_thermal_diffusion_in), alpha_T(alpha_T_in),
            use_scalar_diffusion(use_scalar_diffusion_in), alpha_S(alpha_S_in),
            use_momentum_diffusion(use_momentum_diffusion_in), dynamicViscosity(dynamicViscosity_in),
            use_smagorinsky(use_smagorinsky_in), Cs(Cs_in),
            use_pressure(use_pressure_in), use_gravity(use_gravity_in), spatial_order(spatial_order_in),
            abl_driver_type((!abl_driver_type_in.compare("PressureGradient")) ? abl_driver::PressureGradient 
                     : (!abl_driver_type_in.compare("GeostrophicWind"))  ? abl_driver::GeostrophicWind 
                     : abl_driver::None),

            abl_pressure_grad(get_gpu_array(abl_pressure_grad_in)),
            abl_geo_wind(get_gpu_array(abl_geo_wind_in)){
        // All members are constant and are initialized at construction
    }
    void display() {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "use_state_advection: " << use_state_advection << std::endl;
        amrex::Print() << "use_momentum_advection: " << use_momentum_advection << std::endl;
        amrex::Print() << "use_thermal_diffusion: " << use_thermal_diffusion << std::endl;
        amrex::Print() << "alpha_T: " << alpha_T << std::endl;
        amrex::Print() << "use_scalar_diffusion: " << use_scalar_diffusion << std::endl;
        amrex::Print() << "alpha_S: " << alpha_S << std::endl;
        amrex::Print() << "use_momentum_diffusion: " << use_momentum_diffusion << std::endl;
        amrex::Print() << "dynamicViscosity: " << dynamicViscosity << std::endl;
        amrex::Print() << "use_smagorinsky: " << use_smagorinsky << std::endl;
        amrex::Print() << "Cs: " << Cs << std::endl;
        amrex::Print() << "use_pressure: " << use_pressure << std::endl;
        amrex::Print() << "use_gravity: " << use_gravity << std::endl;
        amrex::Print() << "spatial_order: " << spatial_order << std::endl;
        if (abl_driver_type == abl_driver::None)
            amrex::Print() << "abl_driver_type: " << "None" << std::endl;
            amrex::Print() << "No ABL driver selected " << std::endl;
        if (abl_driver_type == abl_driver::PressureGradient) {
            amrex::Print() << "abl_driver_type: " << "PressureGradient" << std::endl;
            amrex::Print() << "Driving abl_pressure_grad: (";
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
                amrex::Print() << abl_pressure_grad[i] << " ";
            amrex::Print() << ")" << std::endl;
        }
        if (abl_driver_type == abl_driver::GeostrophicWind) {
            amrex::Print() << "abl_driver_type: " << "GeostrophicWind" << std::endl;
            amrex::Print() << "Driving abl_geo_wind: (";
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
                amrex::Print() << abl_geo_wind[i] << " ";
            amrex::Print() << ")" << std::endl;
        }
    }

    // Members
    const bool        use_state_advection;
    const bool        use_momentum_advection;
    const bool        use_thermal_diffusion;
    const amrex::Real alpha_T;
    const bool        use_scalar_diffusion;
    const amrex::Real alpha_S;
    const bool        use_momentum_diffusion;
    const amrex::Real dynamicViscosity;
    const bool        use_smagorinsky;
    const amrex::Real Cs;
    const bool        use_pressure;
    const bool        use_gravity;
    const int         spatial_order;
    const abl_driver  abl_driver_type;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_wind;

private:
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> get_gpu_array(const std::vector<amrex::Real> & vec) {
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> arr {0.0, 0.0, 0.0};
      for(int i = 0; i < AMREX_SPACEDIM; ++i) arr[i] = vec[i];
      return arr;
    }

};

enum class Coord {
    x, y, z
};

enum class NextOrPrev {
    prev = -1,
    local = 0,
    next = 1
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

enum class TurbulenceModel {
    DNS, Smagorinsky
};

enum class ABLDriverType {
    None, GeostrophicWind, PressureGradient
};

#endif
