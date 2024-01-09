#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "prob_common.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real p_inf = p_0; //freestream pressure [Pa]
  amrex::Real T_inf = 300.0; //freestream temperature [K]
  amrex::Real M_inf = 0.2; //freestream Mach number [-]
  amrex::Real alpha = 0.0; //inflow angle, 0 --> x-aligned [rad]
  amrex::Real gamma = Gamma; //specific heat ratio [-]
  amrex::Real beta = 0.01; //non-dimensional max perturbation strength [-]
  amrex::Real sigma = 1.0; //Gaussian standard deviation, i.e., spreading parameter [-]
  amrex::Real R = 2.0; //characteristic length scale for grid [m]
  amrex::Real xc = 0.5; //normalized x-location of vortex center [-]
  amrex::Real yc = 0.5; //normalized y-location of vortex center [-]
  // calculated quantiites
  amrex::Real rho_0; //characteristic density [m/s]
  amrex::Real a_inf; //speed of sound [m/s]
  amrex::Real inv_gm1; //1/(gamma - 1) [-]
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem(const amrex_real* problo, const amrex_real* probhi);

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
        amrex::GeometryData const& geomdata,
        amrex::Array4<amrex::Real const> const& mf_m,
        amrex::Array4<amrex::Real const> const& mf_u,
        amrex::Array4<amrex::Real const> const& mf_v,
        const SolverChoice& sc) override;

protected:
    std::string name() override { return "Isentropic Vortex"; }

private:
    ProbParm parms;
};

#endif
