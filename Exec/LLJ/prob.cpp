#include "prob.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const amrex::Box&  /*bx*/,
    const amrex::Box& /*xbx*/,
    const amrex::Box& /*ybx*/,
    const amrex::Box& /*zbx*/,
    amrex::Array4<amrex::Real const> const& /*state*/,
    amrex::Array4<amrex::Real      > const& /*state_pert*/,
    amrex::Array4<amrex::Real      > const& /*x_vel_pert*/,
    amrex::Array4<amrex::Real      > const& /*y_vel_pert*/,
    amrex::Array4<amrex::Real      > const& /*z_vel_pert*/,
    amrex::Array4<amrex::Real      > const& /*r_hse*/,
    amrex::Array4<amrex::Real      > const& /*p_hse*/,
    amrex::Array4<amrex::Real const> const& /*z_nd*/,
    amrex::Array4<amrex::Real const> const& /*z_cc*/,
    amrex::GeometryData const& /*geomdata*/,
    amrex::Array4<amrex::Real const> const& /*mf_m*/,
    amrex::Array4<amrex::Real const> const& /*mf_u*/,
    amrex::Array4<amrex::Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/)
{
  amrex::Print() << "Dummy function..Needed for linking" << std::endl;
}
