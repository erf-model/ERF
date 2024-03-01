#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "prob_common.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = 1.2;
  amrex::Real T_0 = 300.0; // surface temperature == mean potential temperature
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;

  // random initial perturbations (legacy code)
  amrex::Real U_0_Pert_Mag = 0.0;
  amrex::Real V_0_Pert_Mag = 0.0;
  amrex::Real W_0_Pert_Mag = 0.0;

  // divergence-free initial perturbations
  amrex::Real pert_deltaU = 0.0;
  amrex::Real pert_deltaV = 0.0;
  amrex::Real pert_periods_U  = 5.0;
  amrex::Real pert_periods_V  = 5.0;
  amrex::Real pert_ref_height = 1.0;

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

#include "Prob/init_density_hse_dry_terrain.H"

    void init_custom_pert (
        const amrex::Box&  bx,
        const amrex::Box& xbx,
        const amrex::Box& ybx,
        const amrex::Box& zbx,
        amrex::Array4<amrex::Real      > const& state,
        amrex::Array4<amrex::Real      > const& x_vel,
        amrex::Array4<amrex::Real      > const& y_vel,
        amrex::Array4<amrex::Real      > const& z_vel,
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
    std::string name() override { return "3D Terrain - Hemisphere"; }

private:
    ProbParm parms;
};

#endif
