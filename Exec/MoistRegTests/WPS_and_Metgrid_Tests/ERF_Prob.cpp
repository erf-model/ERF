#include "ERF_Prob.H"

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
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp_erf("erf");
    std::string my_prob_name; pp_erf.get("prob_name",my_prob_name);

    if (my_prob_name == "WPS") {
#include "Prob/ERF_InitCustomPert_KE.H"
    }

    amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_pert_vels (
    const Box& /*xbx*/,
    const Box& /*ybx*/,
    const Box& /*zbx*/,
    Array4<Real      > const& /*x_vel_pert*/,
    Array4<Real      > const& /*y_vel_pert*/,
    Array4<Real      > const& /*z_vel_pert*/,
    Array4<Real const> const& /*z_nd*/,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp("erf");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);

    amrex::Gpu::streamSynchronize();
}
