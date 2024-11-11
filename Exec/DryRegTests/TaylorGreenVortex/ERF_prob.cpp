#include "ERF_prob.H"
#include "ERF_EOS.H"

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
  pp.query("M_0", parms.M_0);
  pp.query("V_0", parms.V_0);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

  ParallelFor(bx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Initial potential temperature (actually rho*theta) perturbation
    const Real p = parms_d.rho_0 * parms_d.V_0*parms_d.V_0*
                          (
                             1.0 / (Gamma * parms_d.M_0 * parms_d.M_0)
                          + (1.0 / 16.0) * (cos(2 * x) + cos(2 * y)) * (cos(2 * z) + 2)
                          );
    state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p) - getRhoThetagivenP(p_hse(i,j,k));

    // Set scalar = 0 everywhere
    state_pert(i, j, k, RhoScalar_comp) = 1.0 * parms_d.rho_0;

    if (use_moisture) {
        state_pert(i, j, k, RhoQ1_comp) = 0.0;
        state_pert(i, j, k, RhoQ2_comp) = 0.0;
    }
  });

  // Set the x-velocity
  ParallelFor(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.0) * dx[0];
      const Real y = prob_lo[1] + (j + 0.5) * dx[1];
      const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the x-velocity
      x_vel_pert(i, j, k) = parms_d.V_0 * sin(x) * cos(y) * cos(z);
  });

  // Set the y-velocity
  ParallelFor(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const Real y = prob_lo[1] + (j + 0.0) * dx[1];
      const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the y-velocity
      y_vel_pert(i, j, k) = - parms_d.V_0 * cos(x) * sin(y) * cos(z);
  });

  // Set the z-velocity
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel_pert(i, j, k) = 0.0;
  });
}
