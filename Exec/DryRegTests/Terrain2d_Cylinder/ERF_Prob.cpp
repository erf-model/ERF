#include "ERF_Prob.H"
#include "ERF_EOS.H"
#include "ERF_TerrainMetrics.H"

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
  Real rho_0 = 1.2; pp.query("rho_0", rho_0);
  Real T_0 = 300.0; pp.query("T_0", T_0);

  init_base_parms(rho_0, T_0);
}

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& /*state_pert*/,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
}

void
Problem::init_custom_pert_vels (
    const Box& xbx,
    const Box& /*ybx*/,
    const Box& /*zbx*/,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& /*y_vel_pert*/,
    Array4<Real      > const& /*z_vel_pert*/,
    Array4<Real const> const& /*z_nd*/,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp("prob");
    Real U_0 = 10.0; pp.query("U_0", U_0);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = U_0;
    });

    amrex::Gpu::streamSynchronize();
}
