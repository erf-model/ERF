#include "prob.H"
#include "prob_common.H"

#include "AMReX_Random.H"
#include "AMReX_ParmParse.H"
#include "IndexDefines.H"

ProbParm parms;

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

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& tau,
                  amrex::Vector<amrex::Real>& ubar,
                  amrex::Vector<amrex::Real>& vbar,
                  amrex::Vector<amrex::Real>& thetabar,
                  amrex::Geometry      const& geom)
{
  const amrex::Real* prob_lo = geom.ProbLo();
  const auto dx              = geom.CellSize();
  const int khi              = geom.Domain().bigEnd()[2];

  for (int k = 0; k <= khi; k++)
  {
      const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

      tau[k]      = 1.0;

      ubar[k]     = parms.U_0;
      vbar[k]     = parms.V_0;
      thetabar[k] = parms.Theta_0;
  }
}

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

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial Rho0*Theta0
    state(i, j, k, RhoTheta_comp) = parms.rho_0 * parms.Theta_0;

    // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
    state(i, j, k, RhoScalar_comp) = parms.A_0 * exp(-10.*r*r);
  });

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelForRNG(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    // Note that this is called on a box of x-faces
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the x-velocity
    amrex::Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
    amrex::Real x_vel_prime = (rand_double*2.0 - 1.0)*parms.U_0_Pert_Mag;

    x_vel(i, j, k) = parms.U_0;
    if(z <= 100.0) {
        x_vel(i, j, k) += x_vel_prime;
    }
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelForRNG(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    // Note that this is called on a box of y-faces
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the y-velocity
    amrex::Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
    amrex::Real y_vel_prime = (rand_double*2.0 - 1.0)*parms.V_0_Pert_Mag;

    y_vel(i, j, k) = parms.V_0;
    if(z <= 100.) {
        y_vel(i, j, k) += y_vel_prime;
    }
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity
  amrex::ParallelForRNG(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
  // Note that this is called on a box of z-faces
//    const amrex::Real* dx = geomdata.CellSize();
//    const amrex::Real* prob_lo = geomdata.ProbLo();
    const int dom_lo_z = geomdata.Domain().smallEnd()[2];
    const int dom_hi_z = geomdata.Domain().bigEnd()[2];

    // Set the z-velocity
    amrex::Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
    amrex::Real z_vel_prime = (rand_double*2.0 - 1.0)*parms.W_0_Pert_Mag;

    if (k == dom_lo_z || k == dom_hi_z+1) {
        z_vel(i, j, k) = 0.0;
    } else {
        z_vel(i, j, k) = parms.W_0 + z_vel_prime;
    }
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
  pp.query("A_0", parms.A_0);

  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
  pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
  pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
}
