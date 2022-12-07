#include "prob.H"
#include "prob_common.H"

#include "IndexDefines.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

using namespace amrex;

ProbParm parms;

void
erf_init_rayleigh(Vector<Real>& /*tau*/,
                  Vector<Real>& /*ubar*/,
                  Vector<Real>& /*vbar*/,
                  Vector<Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for PoiseuilleFlow problem");
}

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>&,
                  std::unique_ptr<MultiFab>&,
                  amrex::Geometry const& /*geom*/)
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
  amrex::Array4<Real> const& state,
  amrex::Array4<Real> const& x_vel,
  amrex::Array4<Real> const& y_vel,
  amrex::Array4<Real> const& z_vel,
  amrex::Array4<Real> const&,
  amrex::Array4<Real> const&,
  amrex::Array4<Real const> const&,
  amrex::Array4<Real const> const&,
  amrex::GeometryData const& geomdata,
  Array4<Real const> const& mf_m,
  Array4<Real const> const& mf_u,
  Array4<Real const> const& mf_v)
{
  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Arbitrarily choose the radius of the bubble to be 0.05 times the length of the domain

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial potential temperature
    state(i, j, k, RhoTheta_comp) = parms.rho_0 * parms.Theta_0;

    // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
    state(i, j, k, RhoScalar_comp) = 0.0;

#ifdef ERF_USE_MOISTURE
    state(i, j, k, RhoQt_comp) = 0.0;
    state(i, j, k, RhoQp_comp) = 0.0;
#endif
  });

  const Box& xbx = amrex::surroundingNodes(bx,0);
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

  const Box& ybx = amrex::surroundingNodes(bx,1);
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

  const Box& zbx = amrex::surroundingNodes(bx,2);
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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

        Array4<Real> z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,0) = 0.;
        });
    }
}

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.Theta_0);

  pp.query("prob_type", parms.prob_type);
}
