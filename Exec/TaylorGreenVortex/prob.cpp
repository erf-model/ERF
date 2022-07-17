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
init_isentropic_hse(const Real& r_sfc, const Real& theta,
                          Real* r,           Real* p,
                    const Real& dz,    const Real&  prob_lo_z,
                    const int& khi)
{
  // r_sfc / p_0 are the density / pressure at the surface
  int k0 = 0;

  // Initial guess
  r[k0] = r_sfc;
  p[k0] = p_0 - (0.5*dz) * r[k0] * CONST_GRAV;

  int MAX_ITER = 10;
  Real TOL = 1.e-8;
  {
      // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
      bool converged_hse = false;
      Real p_hse;
      Real p_eos;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p_0 -  (0.5*dz) * r[k0] * CONST_GRAV;
          p_eos = getPgivenRTh(r[k0]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k0],theta);

          Real drho = A / (dpdr + 0.5 * dz * CONST_GRAV);

          r[k0] = r[k0] + drho;
          p[k0] = getPgivenRTh(r[k0]*theta);

          if (std::abs(drho) < TOL)
          {
              converged_hse = true;
              break;
          }
      }

      if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k0 << std::endl;
      if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
  }

  // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
  for (int k = 1; k <= khi; k++)
  {
      // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
      // to discretely satisfy HSE -- here we assume spatial_order = 2 -- we can generalize this later if needed
      bool converged_hse = false;

      r[k] = r[k-1];

      Real p_eos = getPgivenRTh(r[k]*theta);
      Real p_hse;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p[k-1] -  dz * 0.5 * (r[k-1]+r[k]) * CONST_GRAV;
          p_eos = getPgivenRTh(r[k]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k],theta);
          // Gamma * p_0 * std::pow( (R_d * theta / p_0), Gamma) * std::pow(r[k], Gamma-1.0) ;

          Real drho = A / (dpdr + 0.5 * dz * CONST_GRAV);

          r[k] = std::max(0.9*r[k-1], std::min(r[k] + drho, 1.1*r[k-1]));
          p[k] = getPgivenRTh(r[k]*theta);

          if (std::abs(drho) < TOL * r[k-1])
          {
              converged_hse = true;
              //amrex::Print() << " converged " << std::endl;
              break;
          }
      }

      if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
      if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
  }
}

void
erf_init_rayleigh(Vector<Real>& /*tau*/,
                  Vector<Real>& /*ubar*/,
                  Vector<Real>& /*vbar*/,
                  Vector<Real>& /*thetabar*/,
                  Geometry const& /*geom*/)
{
   amrex::Error("Should never get here for TaylorGreenVortex problem");
}

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>&,
                  std::unique_ptr<MultiFab>&,
                  Geometry const& geom)
{
  const Real prob_lo_z = geom.ProbLo()[2];
  const Real dz        = geom.CellSize()[2];
  const int khi        = geom.Domain().bigEnd()[2];

  const Real& T_sfc    = 300.;
  const Real& rho_sfc  = p_0 / (R_d*T_sfc);
  const Real& Thetabar = T_sfc;

  // These are at cell centers (unstaggered)
  Vector<Real> h_r(khi+1);
  Vector<Real> h_p(khi+1);

  amrex::Gpu::DeviceVector<amrex::Real> d_r(khi+1);
  amrex::Gpu::DeviceVector<amrex::Real> d_p(khi+1);

  init_isentropic_hse(rho_sfc,Thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

  amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
  amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

  Real* r = d_r.data();
  Real* p = d_p.data();


#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(1);
        const Array4<Real> rho_hse_arr = rho_hse[mfi].array();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int kk = std::max(k,0);
            rho_hse_arr(i,j,k) = r[kk];
        });
    }
}

void
init_custom_prob(
        const Box& bx,
        Array4<Real      > const& state,
        Array4<Real      > const& x_vel,
        Array4<Real      > const& y_vel,
        Array4<Real      > const& z_vel,
        Array4<Real      > const& r_hse,
        Array4<Real      > const& p_hse,
        Array4<Real const> const&,
        Array4<Real const> const&,
        GeometryData const& geomdata)
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

#ifdef ERF_USE_MOISTURE
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Construct a box that is on x-faces
  const Box& xbx = surroundingNodes(bx,0);
  // Set the x-velocity
  ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.0) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the x-velocity
    x_vel(i, j, k) = parms.V_0 * sin(x) * cos(y) * cos(z);
  });

  // Construct a box that is on y-faces
  const Box& ybx = surroundingNodes(bx,1);
  // Set the y-velocity
  ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.0) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the y-velocity
    y_vel(i, j, k) = - parms.V_0 * cos(x) * sin(y) * cos(z);
  });

  // Construct a box that is on z-faces
  const Box& zbx = surroundingNodes(bx,2);
  // Set the z-velocity
  ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });
}

void
init_custom_terrain(const Geometry& geom, MultiFab& z_phys_nd)
{
    auto dx = geom.CellSizeArray();

    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<amrex::Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
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
