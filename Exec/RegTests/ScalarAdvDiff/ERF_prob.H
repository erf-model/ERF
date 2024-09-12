#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_prob_common.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = 0.0;
  amrex::Real T_0   = 1.0;
  amrex::Real A_0   = 1.0;
  amrex::Real B_0   = 0.0;
  amrex::Real u_0   = 0.0;
  amrex::Real v_0   = 0.0;
  amrex::Real w_0   = 0.0;
  amrex::Real rad_0 = 0.0;
  amrex::Real z0    = 0.1;  // Surface Roughness
  amrex::Real zRef  = 80.0; // Reference Height
  amrex::Real uRef  = 0.0;  // Reference Wind Speed

  amrex::Real xc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)
  amrex::Real yc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)
  amrex::Real zc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)

  int prob_type = -1;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_init_constant_density_hse.H"

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
    std::string name() override { return "Scalar Advection/Diffusion"; }

private:
    ProbParm parms;
};

#endif
