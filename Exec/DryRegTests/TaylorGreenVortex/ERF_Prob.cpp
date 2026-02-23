#include "ERF_Prob.H"
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
  ParmParse pp("prob");
  Real rho_0 =   1.0; pp.query("rho_0", rho_0);
  Real   T_0 = 300.0; pp.query("T_0", T_0);

  init_base_parms(rho_0, T_0);
}

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp("prob");
    Real rho_0 =   1.0; pp.query("rho_0", rho_0);
    Real   T_0 = 300.0; pp.query("T_0", T_0);
    Real M_0   =   1.0; pp.query("M_0", M_0);
    Real V_0   =   1.0; pp.query("V_0", V_0);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Initial potential temperature (actually rho*theta) perturbation
        const Real p = rho_0 * V_0*V_0*
                              ( 1.0 / (Gamma * M_0 * M_0)
                              + (1.0 / 16.0) * (cos(2 * x) + cos(2 * y)) * (cos(2 * z) + 2));
        state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p) - getRhoThetagivenP(p_hse(i,j,k));

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 1.0 * rho_0;
    });
}

void
Problem::init_custom_pert_vels (
    const Box& xbx,
    const Box& ybx,
    const Box& /*zbx*/,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& /*z_vel_pert*/,
    Array4<Real const> const& /*z_nd*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp("prob");
    Real V_0   =   1.0; pp.query("V_0", V_0);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.0) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = V_0 * sin(x) * cos(y) * cos(z);
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.0) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = - V_0 * cos(x) * sin(y) * cos(z);
    });
}
