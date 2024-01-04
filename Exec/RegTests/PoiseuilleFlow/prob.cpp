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

  pp.query("prob_type", parms.prob_type);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    amrex::Array4<Real> const& state,
    amrex::Array4<Real> const& x_vel,
    amrex::Array4<Real> const& y_vel,
    amrex::Array4<Real> const& z_vel,
    amrex::Array4<Real> const&,
    amrex::Array4<Real> const&,
    amrex::Array4<Real const> const&,
    amrex::Array4<Real const> const&,
    amrex::GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Arbitrarily choose the radius of the bubble to be 0.05 times the length of the domain

    // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
    state(i, j, k, RhoScalar_comp) = 0.0;

    if (use_moisture) {
        state(i, j, k, RhoQ1_comp) = 0.0;
        state(i, j, k, RhoQ2_comp) = 0.0;
    }

  });

  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx      = geomdata.CellSize();
      const Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

      // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
      if (parms.prob_type == 10)
          x_vel(i, j, k) = 1.0 - z_h * z_h;
      else
          x_vel(i, j, k) = 0.0;
  });

  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx      = geomdata.CellSize();
      const Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

      // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
      if (parms.prob_type == 11)
         y_vel(i, j, k) = 1.0 - z_h * z_h;
      else
         y_vel(i, j, k) = 0.0;
  });

  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    z_vel(i, j, k) = 0.0;
  });
}
