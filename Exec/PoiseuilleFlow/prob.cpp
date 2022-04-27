#include "prob.H"
#include "prob_common.H"

#include "IndexDefines.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

ProbParm parms;

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& /*tau*/,
                  amrex::Vector<amrex::Real>& /*ubar*/,
                  amrex::Vector<amrex::Real>& /*vbar*/,
                  amrex::Vector<amrex::Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for PoiseuilleFlow problem");
}

#ifdef ERF_USE_TERRAIN
void
erf_init_dens_hse(amrex::MultiFab& rho_hse,
                  amrex::MultiFab const& z_phys_nd,
                  amrex::MultiFab const& z_phys_cc,
                  amrex::Geometry const& geom)
{
    amrex::Error("This problem not set up to use with terrain");
}

#else

void
erf_init_dens_hse(amrex::Real* dens_hse_ptr,
                  amrex::Geometry const& geom,
                  const int ng_dens_hse)
{
  const int khi = geom.Domain().bigEnd()[2];
  for (int k = -ng_dens_hse; k <= khi+ng_dens_hse; k++)
  {
      dens_hse_ptr[k] = parms.rho_0;
  }
}

#endif

void
erf_init_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
#ifdef ERF_USE_TERRAIN
  amrex::Array4<amrex::Real> const& r_hse,
  amrex::Array4<amrex::Real> const& p_hse,
  amrex::Array4<amrex::Real const> const& z_nd,
  amrex::Array4<amrex::Real const> const& z_cc,
#endif
  amrex::GeometryData const& geomdata)
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
  });

  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* dx      = geomdata.CellSize();
      const amrex::Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

      // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
      if (parms.prob_type == 10)
          x_vel(i, j, k) = 1.0 - z_h * z_h;
      else
          x_vel(i, j, k) = 0.0;
  });

  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* dx      = geomdata.CellSize();
      const amrex::Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

      // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
      if (parms.prob_type == 11)
         y_vel(i, j, k) = 1.0 - z_h * z_h;
      else
         y_vel(i, j, k) = 0.0;
  });

  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    z_vel(i, j, k) = 0.0;
  });
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
