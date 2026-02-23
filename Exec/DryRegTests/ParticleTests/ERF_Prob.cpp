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
  Real T_0 = 300.0; pp.query("T_0", T_0);
  Real rho_0 = 1.0; pp.query("rho_0", T_0);

  init_base_parms(rho_0, T_0);
}

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
  ParmParse pp("prob");
  Real x_c =    0.0; pp.query("x_c", x_c);
  Real z_c = 3000.0; pp.query("z_c", z_c);
  Real x_r = 4000.0; pp.query("x_r", x_r);
  Real z_r = 2000.0; pp.query("z_r", z_r);
  Real T_pert = -15.0; pp.query("T_pert", T_pert);

  // Overridden physical constants
  Real C_p = 1004.0;

  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = z_cc(i,j,k);

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

    Real L = std::sqrt(
        std::pow((x - x_c)/x_r, 2) +
        std::pow((z - z_c)/z_r, 2)
    );
    if (L <= 1.0) {
        Real dT = T_pert * (std::cos(PI*L) + 1.0)/2.0;

        // Note: dT is a temperature perturbation, theta_perturbed is base state + perturbation in theta
        Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), R_d/C_p);

        // This version perturbs rho but not p
        state_pert(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
    }
  });
  amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_pert_vels (
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real const> const& z_nd,
    GeometryData const& geomdata,
    Array4<Real const> const& mf_u,
    Array4<Real const> const& mf_v,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
  ParmParse pp("prob");
  Real U_0 = 0.0; pp.query("U_0",U_0);

  const int klo = geomdata.Domain().smallEnd()[2];
  const int khi = geomdata.Domain().bigEnd()[2];

  // Set the x-velocity
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      Real ztop = z_nd(i,j,khi+1);
      Real zht  = z_nd(i,j,klo);
      x_vel_pert(i, j, k) = U_0 * ztop / (ztop - zht);
  });

  const auto dx = geomdata.CellSize();
  amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
  dxInv[0] = 1. / dx[0];
  dxInv[1] = 1. / dx[1];
  dxInv[2] = 1. / dx[2];

  // Set the z-velocity from impenetrable condition
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel_pert(i, j, k) = WFromOmega(i, j, k, 0.0,
                                       x_vel_pert, y_vel_pert,
                                       mf_u, mf_v, z_nd, dxInv);
  });

  amrex::Gpu::streamSynchronize();
}

