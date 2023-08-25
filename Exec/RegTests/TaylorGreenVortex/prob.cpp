#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "ERF_Constants.H"
#include "IndexDefines.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

using namespace amrex;

ProbParm parms;

void
erf_init_rayleigh(Vector<Real>& /*tau*/,
                  Vector<Real>& /*ubar*/,
                  Vector<Real>& /*vbar*/,
                  Vector<Real>& /*wbar*/,
                  Vector<Real>& /*thetabar*/,
                  Geometry const& /*geom*/)
{
   amrex::Error("Should never get here for TaylorGreenVortex problem");
}

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>&,
                  std::unique_ptr<MultiFab>&,
                  amrex::Geometry const&)
{
    Real rho_0 = parms.rho_0;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(1);
        const Array4<Real> rho_hse_arr = rho_hse[mfi].array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            rho_hse_arr(i,j,k) = rho_0;
        });
    }
}

void
init_custom_prob(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& state,
    Array4<Real      > const& x_vel,
    Array4<Real      > const& y_vel,
    Array4<Real      > const& z_vel,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real const> const&,
    Array4<Real const> const&,
#if defined(ERF_USE_MOISTURE)
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
#elif defined(ERF_USE_WARM_NO_PRECIP)
    Array4<Real      > const&,
    Array4<Real      > const&,
#endif
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice&)
{
  ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial potential temperature (Actually rho*theta)
    const Real p = parms.rho_0 * parms.V_0*parms.V_0*
                          (
                             1.0 / (Gamma * parms.M_0 * parms.M_0)
                          + (1.0 / 16.0) * (cos(2 * x) + cos(2 * y)) * (cos(2 * z) + 2)
                          );
    state(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p);

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 1.0 * state(i,j,k,Rho_comp);

#if defined(ERF_USE_MOISTURE)
    state(i, j, k, RhoQt_comp) = 0.0;
    state(i, j, k, RhoQp_comp) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Set the x-velocity
  ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.0) * dx[0];
      const Real y = prob_lo[1] + (j + 0.5) * dx[1];
      const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the x-velocity
      x_vel(i, j, k) = parms.V_0 * sin(x) * cos(y) * cos(z);
  });

  // Set the y-velocity
  ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const Real y = prob_lo[1] + (j + 0.0) * dx[1];
      const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the y-velocity
      y_vel(i, j, k) = - parms.V_0 * cos(x) * sin(y) * cos(z);
  });

  // Set the z-velocity
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });
}

void
init_custom_terrain (const Geometry& /*geom*/,
                           MultiFab& z_phys_nd,
                     const Real& /*time*/)
{
    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<amrex::Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,0) = 0.0;
        });
    }
}

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.Theta_0);
  pp.query("M_0", parms.M_0);
  pp.query("V_0", parms.V_0);
}
