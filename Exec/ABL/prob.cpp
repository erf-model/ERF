#include "prob.H"
#include "prob_common.H"

#include "AMReX_Random.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "IndexDefines.H"

using namespace amrex;

ProbParm parms;

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>&,
                  std::unique_ptr<MultiFab>&,
                  Geometry const&)
{
    Real R0 = parms.rho_0;
    for ( MFIter mfi(rho_hse, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
       Array4<Real> rho_arr = rho_hse.array(mfi);
       const Box& gbx = mfi.growntilebox(1);
       ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
           rho_arr(i,j,k) = R0;
       });
    }
}

void
erf_init_rayleigh(Vector<Real>& tau,
                  Vector<Real>& ubar,
                  Vector<Real>& vbar,
                  Vector<Real>& wbar,
                  Vector<Real>& thetabar,
                  Geometry      const& geom)
{
  const auto ztop = geom.ProbHi()[2];
  amrex::Print() << "Rayleigh damping (coef="<<parms.dampcoef<<") between "
    << ztop-parms.zdamp << " and " << ztop << std::endl;

  const Real* prob_lo = geom.ProbLo();
  const auto dx       = geom.CellSize();
  const int khi       = geom.Domain().bigEnd()[2];

  for (int k = 0; k <= khi; k++)
  {
      // WRF's vertical velocity damping layer structure, which is based
      // on Durran and Klemp 1983
      const Real z = prob_lo[2] + (k + 0.5) * dx[2];
      const Real zfrac = 1 - (ztop - z) / parms.zdamp;
      if (zfrac >= 0)
      {
          const Real sinefac = std::sin(PIoTwo*zfrac);
          tau[k]      = parms.dampcoef * sinefac * sinefac;
          ubar[k]     = parms.U_0;
          vbar[k]     = parms.V_0;
          wbar[k]     = parms.W_0;
          thetabar[k] = parms.Theta_0;
      }
      else
      {
          tau[k]      = 0.0;
          ubar[k]     = 0.0;
          vbar[k]     = 0.0;
          wbar[k]     = 0.0;
          thetabar[k] = 0.0;
      }
  }
}

void
init_custom_prob(
  const Box& bx,
  Array4<Real> const& state,
  Array4<Real> const& x_vel,
  Array4<Real> const& y_vel,
  Array4<Real> const& z_vel,
  Array4<Real> const&,
  Array4<Real> const&,
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
  ParallelForRNG(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Define a point (xc,yc,zc) at the center of the domain
    const Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
    const Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
    const Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

    const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial Rho0*Theta0
    Real RhoTheta = parms.rho_0 * parms.Theta_0;
    if ((z <= parms.pert_ref_height) && (parms.T_0_Pert_Mag != 0.0)) {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        RhoTheta += (rand_double*2.0 - 1.0)*parms.T_0_Pert_Mag;
    }
    state(i, j, k, RhoTheta_comp) = RhoTheta;

    // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
    state(i, j, k, RhoScalar_comp) = parms.A_0 * exp(-10.*r*r);

    // Set an initial value for QKE
    state(i, j, k, RhoQKE_comp) = parms.QKE_0;

#if defined(ERF_USE_MOISTURE)
    state(i, j, k, RhoQt_comp) = 0.0;
    state(i, j, k, RhoQp_comp) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Construct a box that is on x-faces
  const Box& xbx = surroundingNodes(bx,0);
  // Set the x-velocity
  ParallelForRNG(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    // Note that this is called on a box of x-faces
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the x-velocity
    x_vel(i, j, k) = parms.U_0;
    if ((z <= parms.pert_ref_height) && (parms.U_0_Pert_Mag != 0.0))
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real x_vel_prime = (rand_double*2.0 - 1.0)*parms.U_0_Pert_Mag;
        x_vel(i, j, k) += x_vel_prime;
    }
    if (parms.pert_deltaU != 0.0)
    {
        const amrex::Real yl = y - prob_lo[1];
        const amrex::Real zl = z / parms.pert_ref_height;
        const amrex::Real damp = std::exp(-0.5 * zl * zl);
        x_vel(i, j, k) += parms.ufac * damp * z * std::cos(parms.aval * yl);
    }
  });

  // Construct a box that is on y-faces
  const Box& ybx = surroundingNodes(bx,1);
  // Set the y-velocity
  ParallelForRNG(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    // Note that this is called on a box of y-faces
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the y-velocity
    y_vel(i, j, k) = parms.V_0;
    if ((z <= parms.pert_ref_height) && (parms.V_0_Pert_Mag != 0.0))
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real y_vel_prime = (rand_double*2.0 - 1.0)*parms.V_0_Pert_Mag;
        y_vel(i, j, k) += y_vel_prime;
    }
    if (parms.pert_deltaV != 0.0)
    {
        const amrex::Real xl = x - prob_lo[0];
        const amrex::Real zl = z / parms.pert_ref_height;
        const amrex::Real damp = std::exp(-0.5 * zl * zl);
        y_vel(i, j, k) += parms.vfac * damp * z * std::cos(parms.bval * xl);
    }
  });

  // Construct a box that is on z-faces
  const Box& zbx = surroundingNodes(bx,2);
  // Set the z-velocity
  ParallelForRNG(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
  // Note that this is called on a box of z-faces
//    const Real* dx = geomdata.CellSize();
//    const Real* prob_lo = geomdata.ProbLo();
    const int dom_lo_z = geomdata.Domain().smallEnd()[2];
    const int dom_hi_z = geomdata.Domain().bigEnd()[2];

    // Set the z-velocity
    if (k == dom_lo_z || k == dom_hi_z+1)
    {
        z_vel(i, j, k) = 0.0;
    }
    else if (parms.W_0_Pert_Mag != 0.0)
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real z_vel_prime = (rand_double*2.0 - 1.0)*parms.W_0_Pert_Mag;
        z_vel(i, j, k) = parms.W_0 + z_vel_prime;
    }
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
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.Theta_0);
  pp.query("A_0", parms.A_0);
  pp.query("QKE_0", parms.QKE_0);

  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
  pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
  pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
  pp.query("T_0_Pert_Mag", parms.T_0_Pert_Mag);

  pp.query("pert_deltaU", parms.pert_deltaU);
  pp.query("pert_deltaV", parms.pert_deltaV);
  pp.query("pert_periods_U", parms.pert_periods_U);
  pp.query("pert_periods_V", parms.pert_periods_V);
  pp.query("pert_ref_height", parms.pert_ref_height);
  parms.aval = parms.pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
  parms.bval = parms.pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
  parms.ufac = parms.pert_deltaU * std::exp(0.5) / parms.pert_ref_height;
  parms.vfac = parms.pert_deltaV * std::exp(0.5) / parms.pert_ref_height;

  pp.query("dampcoef", parms.dampcoef);
  pp.query("zdamp", parms.zdamp);
}
