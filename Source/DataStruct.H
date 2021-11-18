#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

enum class ABLDriverType {
    None, PressureGradient, GeostrophicWind
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

struct SolverChoice {
  public:
    SolverChoice(const bool use_state_advection_in, const bool use_momentum_advection_in,
                 const bool  use_thermal_diffusion_in, const amrex::Real alpha_T_in,
                 const bool use_scalar_diffusion_in, const amrex::Real alpha_C_in,
                 const bool use_momentum_diffusion_in, const amrex::Real dynamicViscosity_in,
                 const bool use_smagorinsky_in, const amrex::Real Cs_in,
         const amrex::Real Pr_t_in, const amrex::Real Sc_t_in,
         const amrex::Real rho_0_in,
                 const bool use_pressure_in, const bool use_gravity_in,  const bool use_coriolis_in,
                 const amrex::Real coriolis_factor_in, const amrex::Real sinphi_in,
                 const amrex::Real cosphi_in, const int spatial_order_in,
                 const std::string &abl_driver_type_in,
                 const amrex::Vector<amrex::Real> &abl_pressure_grad_in,
                 const amrex::Vector<amrex::Real> &abl_geo_forcing_in):
            use_state_advection(use_state_advection_in), use_momentum_advection(use_momentum_advection_in),
            use_thermal_diffusion(use_thermal_diffusion_in), alpha_T(alpha_T_in), rhoAlpha_T(alpha_T_in*rho_0_in),
            use_scalar_diffusion(use_scalar_diffusion_in), alpha_C(alpha_C_in), rhoAlpha_C(alpha_C_in*rho_0_in),
            use_momentum_diffusion(use_momentum_diffusion_in), dynamicViscosity(dynamicViscosity_in),
            use_smagorinsky(use_smagorinsky_in), Cs(Cs_in), Pr_t_inv(1.0/Pr_t_in), Sc_t_inv(1.0/Sc_t_in),
            use_pressure(use_pressure_in), use_gravity(use_gravity_in),  use_coriolis(use_coriolis_in),
            coriolis_factor(coriolis_factor_in), cosphi(cosphi_in),
            sinphi(sinphi_in), spatial_order(spatial_order_in),
            abl_driver_type((!abl_driver_type_in.compare("PressureGradient")) ? ABLDriverType::PressureGradient
                     : (!abl_driver_type_in.compare("GeostrophicWind"))  ? ABLDriverType::GeostrophicWind
                     : ABLDriverType::None),

            abl_pressure_grad(get_gpu_array(abl_pressure_grad_in)),
            abl_geo_forcing(get_gpu_array(abl_geo_forcing_in)){
        // All members are constant and are initialized at construction
    }

    void display() {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "use_state_advection: " << use_state_advection << std::endl;
        amrex::Print() << "use_momentum_advection: " << use_momentum_advection << std::endl;
        amrex::Print() << "use_thermal_diffusion: " << use_thermal_diffusion << std::endl;
        amrex::Print() << "alpha_T: " << alpha_T << std::endl;
        amrex::Print() << "rhoAlpha_T: " << rhoAlpha_T << std::endl;
        amrex::Print() << "use_scalar_diffusion: " << use_scalar_diffusion << std::endl;
        amrex::Print() << "alpha_C: " << alpha_C << std::endl;
        amrex::Print() << "rhoAlpha_C: " << rhoAlpha_C << std::endl;
        amrex::Print() << "use_momentum_diffusion: " << use_momentum_diffusion << std::endl;
        amrex::Print() << "dynamicViscosity: " << dynamicViscosity << std::endl;
        amrex::Print() << "use_smagorinsky: " << use_smagorinsky << std::endl;
        amrex::Print() << "Cs: " << Cs << std::endl;
        amrex::Print() << "Pr_t: " << 1.0/Pr_t_inv << std::endl;
        amrex::Print() << "Sc_t: " << 1.0/Sc_t_inv << std::endl;
        amrex::Print() << "use_pressure: " << use_pressure << std::endl;
        amrex::Print() << "use_coriolis " << use_coriolis << std::endl;
        amrex::Print() << "use_gravity: " << use_gravity << std::endl;
        amrex::Print() << "spatial_order: " << spatial_order << std::endl;
        if (abl_driver_type == ABLDriverType::None) {
            amrex::Print() << "ABL Driver Type: " << "None" << std::endl;
            amrex::Print() << "No ABL driver selected " << std::endl;
        }
        if (abl_driver_type == ABLDriverType::PressureGradient) {
            amrex::Print() << "ABL Driver Type: " << "PressureGradient" << std::endl;
            amrex::Print() << "Driving abl_pressure_grad: (";
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
                amrex::Print() << abl_pressure_grad[i] << " ";
            amrex::Print() << ")" << std::endl;
        }
        if (abl_driver_type == ABLDriverType::GeostrophicWind) {
            amrex::Print() << "ABL Driver Type: " << "GeostrophicWind" << std::endl;
            amrex::Print() << "Driving abl_geo_forcing: (";
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
                amrex::Print() << abl_geo_forcing[i] << " ";
            amrex::Print() << ")" << std::endl;
        }
    }

    // Members
    const bool        use_state_advection;
    const bool        use_momentum_advection;
    const bool        use_thermal_diffusion;
    const amrex::Real alpha_T;
    const amrex::Real rhoAlpha_T = 0.0;
    const bool        use_scalar_diffusion;
    const amrex::Real alpha_C;
    const amrex::Real rhoAlpha_C = 0.0;
    const bool        use_momentum_diffusion;
    const amrex::Real dynamicViscosity;
    const bool        use_smagorinsky;
    const amrex::Real Cs;
    const amrex::Real Pr_t_inv;
    const amrex::Real Sc_t_inv;
    const bool        use_pressure;
    const bool        use_gravity;
    const bool        use_coriolis;
    const amrex::Real coriolis_factor;
    const amrex::Real cosphi;
    const amrex::Real sinphi;
    const int         spatial_order;
    const ABLDriverType  abl_driver_type;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_pressure_grad;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> abl_geo_forcing;

  private:
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> get_gpu_array(const std::vector<amrex::Real> & vec) {
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> arr {0.0, 0.0, 0.0};
      for(int i = 0; i < AMREX_SPACEDIM; ++i) arr[i] = vec[i];
      return arr;
    }

};

#endif
