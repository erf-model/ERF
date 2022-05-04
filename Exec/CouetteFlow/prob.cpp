#include "prob.H"
#include "prob_common.H"

#include "IndexDefines.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"

ProbParm parms;

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
erf_init_rayleigh(amrex::Vector<amrex::Real>& /*tau*/,
                  amrex::Vector<amrex::Real>& /*ubar*/,
                  amrex::Vector<amrex::Real>& /*vbar*/,
                  amrex::Vector<amrex::Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for CouetteFlow problem");
}

void
init_custom_prob(
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
    amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial potential temperature
    state(i, j, k, RhoTheta_comp) = parms.rho_0 * parms.Theta_0;

    // Set scalar to 0
    state(i, j, k, RhoScalar_comp) = 0.0;
  });

  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const auto prob_hi  = geomdata.ProbHi();
    const auto dx       = geomdata.CellSize();
    const amrex::Real z = (k + 0.5) * dx[2];
    x_vel(i, j, k) = parms.u_0 * z / prob_hi[2];
  });

  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const auto prob_hi  = geomdata.ProbHi();
    const auto dx       = geomdata.CellSize();
    const amrex::Real z = (k + 0.5) * dx[2];
    y_vel(i, j, k) = parms.v_0 * z / prob_hi[2];
  });

  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = parms.w_0;
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
  pp.query("u_0", parms.u_0);
  pp.query("v_0", parms.v_0);
  pp.query("w_0", parms.w_0);
}
