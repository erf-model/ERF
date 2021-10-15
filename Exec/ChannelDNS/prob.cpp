#include "prob.H"

ProbParm parms;

void
erf_init_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata)
{
  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Geometry
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Define a point (xc,yc,zc) at the center of the domain
    const amrex::Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
    const amrex::Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
    const amrex::Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

    const amrex::Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

    // Arbitrarily choose the radius of the bubble to be 0.05 times the length of the domain

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial potential temperature
    state(i, j, k, RhoTheta_comp) = parms.rho_0 * parms.T_0;

    // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
    state(i, j, k, RhoScalar_comp) = parms.A_0 * exp(-10.*r*r);
  });

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    x_vel(i, j, k) = 10.0;
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });
}

AMREX_GPU_DEVICE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}

void
erf_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);
  pp.query("A_0", parms.A_0);
}
}
