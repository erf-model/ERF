#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_prob_common.H"

struct ProbParm : ProbParmDefaults {
  // perturbations to initial conditions
  amrex::Real rho_0 = 0.0;
  amrex::Real T_0 = 0.0;
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;
  amrex::Real W_0 = 0.0; // needed for rayleigh damping

  // Rayleigh parameters
  amrex::Real dampcoef = 0.003;
  amrex::Real zdamp = 1.0;

  // hill parameters
  amrex::Real hmax = 0.0;   // full hill height, can be negative to simulate valleys
  amrex::Real L = 0.0;      // hill length at half-height
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_init_density_hse_dry.H"
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

    void init_custom_terrain (
        const amrex::Geometry& geom,
        amrex::MultiFab& z_phys_nd,
        const amrex::Real& time) override;

protected:
    std::string name() override { return "Witch of Agnesi Hill"; }

private:
    ProbParm parms;
};

#endif
