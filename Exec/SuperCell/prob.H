#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "prob_common.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real T_0 = 300.0; // surface temperature == mean potential temperature
  amrex::Real U_0 = 10.0;
  amrex::Real V_0 = 0.0;
  amrex::Real x_c = 0.0; // center of thermal perturbation
  amrex::Real z_c = 3200.0;
  amrex::Real x_r = 1000.0;
  amrex::Real z_r = 1000.0;
  amrex::Real T_pert = 5.0; // perturbation temperature
  // overridden physical constants
  amrex::Real C_p = 1004.0;
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

    void erf_init_rayleigh (
        amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
        amrex::Geometry      const& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

protected:
    std::string name() override { return "Supercell"; }

private:
    ProbParm parms;
};

#endif
