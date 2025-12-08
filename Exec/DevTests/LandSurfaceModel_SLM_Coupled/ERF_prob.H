#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

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
    amrex::Real pert_periods_U = 0.0;
    amrex::Real pert_periods_V = 0.0;
    amrex::Real pert_ref_height = 200.0;

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


    std::string lsf_file = "";
    std::string snd_file = "";
    amrex::Real lsf_scale = 1.0;

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

#include "Prob/ERF_InitConstantDensityHSE.H"
#include "Prob/ERF_InitRayleighDamping.H"

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
        const SolverChoice& sc,
        const int lev) override;

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

    void read_forcing_file(
        const std::string lsf_file,
        int& nt,
        int& nz,
        amrex::Vector<amrex::Real> &times,
        amrex::Vector<amrex::Vector<amrex::Real>>& z_lsf,
        amrex::Vector<amrex::Vector<amrex::Real>>& p_lsf,
        amrex::Vector<amrex::Vector<amrex::Real>>& t_lsf,
        amrex::Vector<amrex::Vector<amrex::Real>>& q_lsf,
        amrex::Vector<amrex::Vector<amrex::Real>>& u_lsf,
        amrex::Vector<amrex::Vector<amrex::Real>>& v_lsf,
        amrex::Vector<amrex::Vector<amrex::Real>>& w_lsf);

    void interp_forcing(const amrex::GeometryData &geom,
                        const amrex::Vector<amrex::Real>& zlevels_stag,
                        const amrex::Vector<amrex::Real>& times,
                        amrex::Vector<amrex::Vector<amrex::Real>>& z_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& t_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& q_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& u_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& v_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& w_lsf);

    void get_forcing_time_coeffs(const amrex::Real& time,
                                 int& curr,
                                 int& next,
                                 amrex::Real& coeff_curr,
                                 amrex::Real& coeff_next);

  protected:
    std::string name () override { return "BOMEX"; }

private:
    ProbParm parms;
    amrex::Vector<amrex::Real> zlevels;

    // inputs from large scale forcing file
    int nz_lsf = -1; // number of z levels in the LSF file
    int nt_lsf = -1; // number of times found in the LSF file
    amrex::Vector<amrex::Real> times_lsf;
    // original data from the file
    amrex::Vector<amrex::Vector<amrex::Real>> z_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> p_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> t_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> q_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> u_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> v_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> w_lsf;

    // interpolated data from file
    amrex::Vector<amrex::Vector<amrex::Real>> z_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> p_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> t_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> q_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> u_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> v_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> w_int_lsf;
};

#endif
