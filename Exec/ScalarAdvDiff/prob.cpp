#include "prob.H"
#include "prob_common.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "ERF_Constants.H"
#include "IndexDefines.H"

using namespace amrex;

ProbParm parms;

void
erf_init_rayleigh(amrex::Vector<Real>& tau,
                  amrex::Vector<Real>& ubar,
                  amrex::Vector<Real>& vbar,
                  amrex::Vector<Real>& wbar,
                  amrex::Vector<Real>& thetabar,
                  amrex::Geometry      const& geom)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      tau[k]  = 1.0;
      ubar[k] = 2.0;
      vbar[k] = 1.0;
      wbar[k] = 0.0;
      thetabar[k] = parms.Theta_0;
  }
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
    const amrex::Box& bx,
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
  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Define a point (xc,yc,zc) at the center of the domain
    const Real xc = parms.xc_frac * (prob_lo[0] + prob_hi[0]);
    const Real yc = parms.yc_frac * (prob_lo[1] + prob_hi[1]);
    const Real zc = parms.zc_frac * (prob_lo[2] + prob_hi[2]);

    // Define ellipse parameters
    const Real r0   = parms.rad_0 * (prob_hi[0] - prob_lo[0]);

    const Real r3d    = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
    const Real r2d_xy = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
    const Real r2d_xz = std::sqrt((x-xc)*(x-xc) + (z-zc)*(z-zc));

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial potential temperature
    state(i, j, k, RhoTheta_comp) = parms.rho_0 * parms.Theta_0;

    if (parms.prob_type == 10)
    {
        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain,
        //            + B_0*sin(x)
        state(i, j, k, RhoScalar_comp) = parms.A_0 * exp(-10.*r3d*r3d) + parms.B_0*sin(x);

    } else if (parms.prob_type == 11) {
        state(i, j, k, RhoScalar_comp) = parms.A_0 * 0.25 * (1.0 + std::cos(PI * std::min(r2d_xy, r0) / r0));
    } else if (parms.prob_type == 12) {
        state(i, j, k, RhoScalar_comp) = parms.A_0 * 0.25 * (1.0 + std::cos(PI * std::min(r2d_xz, r0) / r0));
    } else if (parms.prob_type == 13) {
        const Real r0_z = parms.rad_0 * (prob_hi[2] - prob_lo[2]);
        const Real r2d_xz = std::sqrt((x-xc)*(x-xc)/(r0*r0) + (z-zc)*(z-zc)/(r0_z*r0_z)); //ellipse for mapfac shear validation
        state(i, j, k, RhoScalar_comp) = parms.A_0 * 0.25 * (1.0 + std::cos(PI * std::min(r2d_xz, r0_z) / r0_z));
    } else {
        // Set scalar = A_0 in a ball of radius r0 and 0 elsewhere
        if (r3d < r0) {
           state(i, j, k, RhoScalar_comp) = parms.A_0;
        } else {
           state(i, j, k, RhoScalar_comp) = 0.0;
        }
    }

    state(i, j, k, RhoScalar_comp) *= parms.rho_0;

#if defined(ERF_USE_MOISTURE)
    state(i, j, k, RhoQt_comp) = 0.0;
    state(i, j, k, RhoQp_comp) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif

  });

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel(i, j, k) = parms.u_0;

      const Real* prob_lo = geomdata.ProbLo();
      const Real*      dx = geomdata.CellSize();
      const Real        z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the x-velocity
      x_vel(i, j, k) = parms.u_0 + parms.uRef *
                       std::log((z + parms.z0)/parms.z0)/
                       std::log((parms.zRef +parms.z0)/parms.z0);
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel(i, j, k) = parms.v_0;
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });
}

void
init_custom_terrain(const Geometry& geom, MultiFab& z_phys_nd,
                    const Real& /*time*/)
{
    auto dx = geom.CellSizeArray();

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& gbx = mfi.growntilebox(1);
        Array4<Real> z_arr = z_phys_nd.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real z = k * dx[2];

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k) = z;
        });
    }
    z_phys_nd.FillBoundary(geom.periodicity());
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
    pp.query("B_0", parms.B_0);
    pp.query("u_0", parms.u_0);
    pp.query("v_0", parms.v_0);
    pp.query("rad_0", parms.rad_0);
    pp.query("z0", parms.z0);
    pp.query("zRef", parms.zRef);
    pp.query("uRef", parms.uRef);

    pp.query("xc_frac", parms.xc_frac);
    pp.query("yc_frac", parms.yc_frac);
    pp.query("zc_frac", parms.zc_frac);

   pp.query("prob_type", parms.prob_type);
}

AMREX_GPU_DEVICE
Real
dhdt(int /*i*/, int /*j*/,
     const GpuArray<Real,AMREX_SPACEDIM> /*dx*/,
     const Real /*time_mt*/, const Real /*delta_t*/)
{
    return 0.;
}
