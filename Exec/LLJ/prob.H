#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "prob_common.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = 1.0;
  amrex::Real T_0 = 300.0;
  amrex::Real V_0 = 1.0;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/init_constant_density_hse.H"

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
        amrex::Array4<amrex::Real      > const& qv,
        amrex::Array4<amrex::Real      > const& qc,
        amrex::Array4<amrex::Real      > const& qi,
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
    std::string name() override { return "Low-Level Jet"; }

private:
    ProbParm parms;
};

#endif
