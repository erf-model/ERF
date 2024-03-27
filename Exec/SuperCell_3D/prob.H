#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

#include "prob_common.H"
#include "EOS.H"
#include "IndexDefines.H"
#include "TileNoZ.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real z_tr = 12000.0;
  amrex::Real height = 1200.0;
  amrex::Real theta_0 = 300.0;
  amrex::Real theta_tr = 343.0;
  amrex::Real T_tr = 213.0;

  amrex::Real x_c = 0.0;
  amrex::Real z_c = 1.5e3;
  amrex::Real x_r = 4.0e3;
  amrex::Real z_r = 1.5e3;
  amrex::Real theta_c = 3.0;

  amrex::Real C_p = 1004.0;
  amrex::Real Theta_0 = 300.0;
}; // namespace ProbParm


class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/init_density_hse_dry.H"

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


    void erf_init_dens_hse_moist (amrex::MultiFab& rho_hse,
                                  std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                                  amrex::Geometry const& geom) override;

    void erf_init_rayleigh (
        amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
        amrex::Geometry      const& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

    amrex::Real compute_theta (amrex::Real z);

    amrex::Real compute_p_k (amrex::Real& p_k,
                             const amrex::Real p_k_minus_1,
                             amrex::Real& theta_k,
                             amrex::Real& rho_k,
                             amrex::Real& q_v_k,
                             amrex::Real& T_dp,
                             amrex::Real& T_b,
                             const amrex::Real dz,
                             const amrex::Real z,
                             const amrex::Real rho_k_minus_1);

    amrex::Real compute_F (const amrex::Real& p_k,
                           const amrex::Real& p_k_minus_1,
                           amrex::Real &theta_k,
                           amrex::Real& rho_k,
                           amrex::Real& q_v_k,
                           amrex::Real& T_dp,
                           amrex::Real& T_b,
                           const amrex::Real& dz,
                           const amrex::Real& z,
                           const amrex::Real& rho_k_minus_1);

    void compute_rho (const amrex::Real& z,
                      const amrex::Real& pressure,
                      amrex::Real &theta,
                      amrex::Real& rho,
                      amrex::Real& q_v,
                      amrex::Real& T_dp,
                      amrex::Real& T_b);

    void init_isentropic_hse_no_terrain(amrex::Real *theta,
                                        amrex::Real* r,
                                        amrex::Real* p,
                                        amrex::Real *q_v,
                                        const amrex::Real& dz,
                                        const amrex::Real&  prob_lo_z,
                                        const int& khi);

protected:
    std::string name () override { return "Supercell"; }

private:
    ProbParm parms;
};

#endif
