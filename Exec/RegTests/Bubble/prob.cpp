#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "TileNoZ.H"

using namespace amrex;

ProbParm parms;

void
init_isentropic_hse_no_terrain(const Real& r_sfc, const Real& theta,
                               Real* r,           Real* p,
                               const Real& dz,    const Real&  /*prob_lo_z*/,
                               const int& khi)
{
  // r_sfc / p_0 are the density / pressure at the surface
  int k0 = 0;

  // Initial guess
  Real half_dz = 0.5*dz;
  r[k0] = r_sfc;
  p[k0] = p_0 - half_dz * r[k0] * CONST_GRAV;

  int MAX_ITER = 10;
  Real TOL = 1.e-8;
  {
      // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
      bool converged_hse = false;
      Real p_hse;
      Real p_eos;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p_0 -  half_dz * r[k0] * CONST_GRAV;
          p_eos = getPgivenRTh(r[k0]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k0],theta);

          Real drho = A / (dpdr + half_dz * CONST_GRAV);

          r[k0] = r[k0] + drho;
          p[k0] = getPgivenRTh(r[k0]*theta);

          if (std::abs(drho) < TOL)
          {
              converged_hse = true;
              break;
          }
      }

      // if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k0 << std::endl;
      // if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
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
          Real r_avg = 0.5 * (r[k-1]+r[k]);
          p_hse = p[k-1] -  dz * r_avg * CONST_GRAV;
          p_eos = getPgivenRTh(r[k]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k],theta);
          // Gamma * p_0 * std::pow( (R_d * theta / p_0), Gamma) * std::pow(r[k], Gamma-1.0) ;

          Real drho = A / (dpdr + dz * CONST_GRAV);

          r[k] = r[k] + drho;
          p[k] = getPgivenRTh(r[k]*theta);

          if (std::abs(drho) < TOL * r[k-1])
          {
              converged_hse = true;
              //amrex::Print() << " converged " << std::endl;
              break;
          }
      }

      // if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
      // if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
  }
  r[khi+1] = r[khi];
}

AMREX_GPU_DEVICE
static
void
init_isentropic_hse_terrain(int i, int j,
                            const Real& r_sfc, const Real& theta,
                                  Real* r,           Real* p,
                            const Array4<Real const> z_cc,
                            const int& khi)
{
  // r_sfc / p_0 are the density / pressure at the surface
  int k0  = 0;
  Real half_dz = z_cc(i,j,k0);

  // Initial guess
  r[k0] = r_sfc;
  p[k0] = p_0 - half_dz * r[k0] * CONST_GRAV;

  int MAX_ITER = 10;
  Real TOL = 1.e-8;

  {
      // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
      bool converged_hse = false;
      Real p_hse;
      Real p_eos;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p_0 - half_dz * r[k0] * CONST_GRAV;
          p_eos = getPgivenRTh(r[k0]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k0],theta);

          Real drho = A / (dpdr + half_dz * CONST_GRAV);

          r[k0] = r[k0] + drho;
          p[k0] = getPgivenRTh(r[k0]*theta);

          if (std::abs(drho) < TOL)
          {
              converged_hse = true;
              break;
          }
      }

      //if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k0 << std::endl;
      //if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
  }

  // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
  for (int k = 1; k <= khi; k++)
  {
      // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
      // to discretely satisfy HSE -- here we assume spatial_order = 2 -- we can generalize this later if needed
      bool converged_hse = false;

      Real dz_loc = (z_cc(i,j,k) - z_cc(i,j,k-1));

      r[k] = r[k-1];

      Real p_eos = getPgivenRTh(r[k]*theta);
      Real p_hse;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p[k-1] -  dz_loc * 0.5 * (r[k-1]+r[k]) * CONST_GRAV;
          p_eos = getPgivenRTh(r[k]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k],theta);

          Real drho = A / (dpdr + dz_loc * CONST_GRAV);

          r[k] = r[k] + drho;
          p[k] = getPgivenRTh(r[k]*theta);

          if (std::abs(drho) < TOL * r[k-1])
          {
              converged_hse = true;
              //amrex::Print() << " converged " << std::endl;
              break;
          }
      }

      //if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
      //if (!converged_hse) amrex::Error("Didn't converge the iterations in init");
  }
}

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>& z_phys_nd,
                  std::unique_ptr<MultiFab>& z_phys_cc,
                  Geometry const& geom)
{
    // if erf.init_type is specified, then the base density should
    // already have been calculated. In that case, these assumed (dry)
    // surface conditions are probably incorrect.
    AMREX_ALWAYS_ASSERT(parms.T_0 > 0);

    const Real prob_lo_z = geom.ProbLo()[2];
    const Real dz        = geom.CellSize()[2];
    const int khi        = geom.Domain().bigEnd()[2];

    const Real T_sfc    = 300.;
    const Real rho_sfc  = p_0 / (R_d*T_sfc);
    const Real Thetabar = T_sfc;

    // use_terrain = 1
    if (z_phys_nd) {

        if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

        for ( MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > rho_arr  = rho_hse.array(mfi);
            Array4<Real const> z_cc_arr = z_phys_cc->const_array(mfi);

            // Create a flat box with same horizontal extent but only one cell in vertical
            const Box& tbz = mfi.nodaltilebox(2);
            Box b2d = tbz; // Copy constructor
            b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
            b2d.setRange(2,0);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
              Array1D<Real,0,255> r;;
              Array1D<Real,0,255> p;;

              init_isentropic_hse_terrain(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

              for (int k = 0; k <= khi; k++) {
                 rho_arr(i,j,k) = r(k);
              }
              rho_arr(i,j,   -1) = rho_arr(i,j,0);
              rho_arr(i,j,khi+1) = rho_arr(i,j,khi);
            });
        } // mfi

    } else { // use_terrain = 0

        // These are at cell centers (unstaggered)
        Vector<Real> h_r(khi+2);
        Vector<Real> h_p(khi+2);

        amrex::Gpu::DeviceVector<Real> d_r(khi+2);
        amrex::Gpu::DeviceVector<Real> d_p(khi+2);

        init_isentropic_hse_no_terrain(rho_sfc,Thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

        Real* r = d_r.data();

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
          } // mfi
    } // no terrain
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
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
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
    const SolverChoice& sc)
{
    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    const Real rho_sfc   = p_0 / (R_d*parms.T_0);
    const Real thetabar  = parms.T_0;
    const Real dz        = geomdata.CellSize()[2];
    const Real prob_lo_z = geomdata.ProbLo()[2];
    const Real rdOcp     = sc.rdOcp;

    // These are at cell centers (unstaggered)
    Vector<Real> h_r(khi+1);
    Vector<Real> h_p(khi+1);

    amrex::Gpu::DeviceVector<Real> d_r(khi+1);
    amrex::Gpu::DeviceVector<Real> d_p(khi+1);

    amrex::Print() << "Bubble delta T = " << parms.T_pert << " K" << std::endl;
    amrex::Print() << "  centered at ("
        << parms.x_c << " " << parms.y_c << " " << parms.z_c << ")" << std::endl;
    amrex::Print() << "  with extent ("
        << parms.x_r << " " << parms.y_r << " " << parms.z_r << ")" << std::endl;

    if (parms.T_0 <= 0)
    {
        amrex::Print() << "Ignoring parms.T_0 = " << parms.T_0
            << ", background fields should have been initialized with erf.init_type"
            << std::endl;
    }

    if (z_cc) { // nonflat terrain

#if 0
        if (parms.T_0 > 0)
        {
            // Create a flat box with same horizontal extent but only one cell in vertical
            Box b2d = surroundingNodes(bx); // Copy constructor
            b2d.setRange(2,0);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                Array1D<Real,0,255> r;;
                Array1D<Real,0,255> p;;

                init_isentropic_hse_terrain(i,j,rho_sfc,thetabar,&(r(0)),&(p(0)),z_cc,khi);

                for (int k = 0; k <= khi; k++) {
                   r_hse(i,j,k) = r(k);
                   p_hse(i,j,k) = p(k);
                }
                r_hse(i,j,   -1) = r_hse(i,j,0);
                r_hse(i,j,khi+1) = r_hse(i,j,khi);
            });
        }
#endif

        amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const Real z = z_cc(i,j,k);

            perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k),
                              parms, rdOcp,
                              state(i, j, k, Rho_comp),
                              state(i, j, k, RhoTheta_comp));

            state(i, j, k, RhoScalar_comp) = 0.0;

#if defined(ERF_USE_MOISTURE)
            state(i, j, k, RhoQt_comp) = 0.0;
            state(i, j, k, RhoQp_comp) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
            state(i, j, k, RhoQv_comp) = 0.0;
            state(i, j, k, RhoQc_comp) = 0.0;
#endif

        });
    } else {

#if 0
        if (parms.T_0 > 0)
        {
            init_isentropic_hse_no_terrain(rho_sfc,thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

            Real* r = d_r.data();
            Real* p = d_p.data();

            amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                for (int k = 0; k <= khi; k++) {
                   r_hse(i,j,k) = r[k];
                   p_hse(i,j,k) = p[k];
                }
                r_hse(i,j,   -1) = r_hse(i,j,0);
                r_hse(i,j,khi+1) = r_hse(i,j,khi);
            });
        }
#endif

        amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const Real z = prob_lo[2] + (k + 0.5) * dx[2];

            perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k),
                              parms, rdOcp,
                              state(i, j, k, Rho_comp),
                              state(i, j, k, RhoTheta_comp));

            state(i, j, k, RhoScalar_comp) = 0.0;

#ifdef ERF_USE_MOISTURE
            state(i, j, k, RhoQt_comp) = 0.0;
            state(i, j, k, RhoQp_comp) = 0.0;
#endif
        });
    }

    const Real u0 = parms.U_0;
    const Real v0 = parms.V_0;

    // Set the x-velocity
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel(i, j, k) = u0;
    });

    // Set the y-velocity
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel(i, j, k) = v0;
    });

    // Set the z-velocity
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel(i, j, k) = 0.0;
    });

    amrex::Gpu::streamSynchronize();
}

void
init_custom_terrain (const Geometry& /*geom*/,
                           MultiFab& z_phys_nd,
                     const Real& /*time*/)
{
    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}

void
erf_init_rayleigh(Vector<Real>& tau,
                  Vector<Real>& ubar,
                  Vector<Real>& vbar,
                  Vector<Real>& wbar,
                  Vector<Real>& thetabar,
                  amrex::Geometry const& geom)
{
  const auto ztop = geom.ProbHi()[2];
  amrex::Print() << "Rayleigh damping (coef="<<parms.dampcoef<<") between "
    << ztop-parms.zdamp << " and " << ztop << std::endl;

  const int khi = geom.Domain().bigEnd()[2];
  const auto prob_lo = geom.ProbLo();
  const auto dx = geom.CellSize();

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
          wbar[k]     = 0.0;
          thetabar[k] = parms.T_0;
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
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("x_c", parms.x_c);
  pp.query("y_c", parms.y_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("y_r", parms.y_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
  pp.query("T_pert_is_airtemp", parms.T_pert_is_airtemp);
  pp.query("perturb_rho", parms.perturb_rho);
  pp.query("dampcoef", parms.dampcoef);
  pp.query("zdamp", parms.zdamp);
}
