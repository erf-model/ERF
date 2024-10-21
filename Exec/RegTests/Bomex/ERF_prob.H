#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_prob_common.H"

struct ProbParm : ProbParmDefaults {
    amrex::Real rho_0 = 0.0;
    amrex::Real T_0   = 0.0;
    amrex::Real A_0   = 1.0;
    amrex::Real KE_0 = 0.1;

    amrex::Real U_0 = 0.0;
    amrex::Real V_0 = 0.0;
    amrex::Real W_0 = 0.0;

    // random initial perturbations (legacy code)
    amrex::Real U_0_Pert_Mag = 0.0;
    amrex::Real V_0_Pert_Mag = 0.0;
    amrex::Real W_0_Pert_Mag = 0.0;
    amrex::Real T_0_Pert_Mag = 0.0; // perturbation to rho*Theta
    amrex::Real qv_0_Pert_Mag = 0.0;

    // divergence-free initial perturbations
    amrex::Real pert_deltaU = 0.0;
    amrex::Real pert_deltaV = 0.0;
    amrex::Real pert_periods_U = 5.0;
    amrex::Real pert_periods_V = 5.0;
    amrex::Real pert_ref_height = 100.0;

    //==============================================
    // USER-DEFINED INPUTS
    // source terms
    amrex::Real advection_heating_rate = 0.0;
    amrex::Real restart_time = 9e9;
    amrex::Real cutoff = 500.0;
    amrex::Real cutoff_transition = 50.0;

    amrex::Real advection_moisture_rate = 0.0;
    amrex::Real moisture_cutoff = 500.0;
    amrex::Real moisture_cutoff_transition = 50.0;

    amrex::Real wbar_sub_max    = -0.65;
    amrex::Real wbar_cutoff_max = 1500.0;
    amrex::Real wbar_cutoff_min = 2100.0;

    bool custom_TKE = false;

    //==============================================

    // helper vars
    amrex::Real aval;
    amrex::Real bval;
    amrex::Real ufac;
    amrex::Real vfac;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem (const amrex::Real* problo, const amrex::Real* probhi);

#include "Prob/ERF_init_constant_density_hse.H"
#include "Prob/ERF_init_rayleigh_damping.H"

    void init_custom_pert (
        const amrex::Box&  bx,
        const amrex::Box& xbx,
        const amrex::Box& ybx,
        const amrex::Box& zbx,
        amrex::Array4<amrex::Real const> const& state,
        amrex::Array4<amrex::Real      > const& state_pert,
        amrex::Array4<amrex::Real      > const& x_vel_pert,
        amrex::Array4<amrex::Real      > const& y_vel_pert,
        amrex::Array4<amrex::Real      > const& z_vel_pert,
        amrex::Array4<amrex::Real      > const& r_hse,
        amrex::Array4<amrex::Real      > const& p_hse,
        amrex::Array4<amrex::Real const> const& z_nd,
        amrex::Array4<amrex::Real const> const& z_cc,
        amrex::GeometryData const& geomdata,
        amrex::Array4<amrex::Real const> const& mf_m,
        amrex::Array4<amrex::Real const> const& mf_u,
        amrex::Array4<amrex::Real const> const& mf_v,
        const SolverChoice& sc) override;

    void update_rhotheta_sources (
        const amrex::Real& /*time*/,
        amrex::Vector<amrex::Real>& src,
        amrex::Gpu::DeviceVector<amrex::Real>& d_src,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

    void update_rhoqt_sources (
        const amrex::Real& /*time*/,
        amrex::Vector<amrex::Real>& qsrc,
        amrex::Gpu::DeviceVector<amrex::Real>& d_qsrc,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

    void update_w_subsidence (
        const amrex::Real& /*time*/,
        amrex::Vector<amrex::Real>& wbar,
        amrex::Gpu::DeviceVector<amrex::Real>& d_wbar,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

    void update_geostrophic_profile (
        const amrex::Real& /*time*/,
        amrex::Vector<amrex::Real>& u_geos,
        amrex::Gpu::DeviceVector<amrex::Real>& d_u_geos,
        amrex::Vector<amrex::Real>& v_geos,
        amrex::Gpu::DeviceVector<amrex::Real>& d_v_geos,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

protected:
    std::string name () override { return "BOMEX"; }

private:
    ProbParm parms;
    amrex::Vector<amrex::Real> zlevels;
};

#endif
