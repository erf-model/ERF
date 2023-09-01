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
static void
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
  }
}

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>& z_phys_nd,
                  std::unique_ptr<MultiFab>& z_phys_cc,
                  Geometry const& geom)
{
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_temperature(amrex::Real z, amrex::Real z_0, amrex::Real z_trop, amrex::Real z_top,
                           amrex::Real T_0, amrex::Real T_trop, amrex::Real T_top)
{
    if (z <= z_trop) {
      amrex::Real lapse = - (T_trop - T_0) / (z_trop - z_0);
      return T_0 - lapse * (z - z_0);
    } else {
      amrex::Real lapse = - (T_top - T_trop) / (z_top - z_trop);
      return T_trop - lapse * (z - z_trop);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_pressure(amrex::Real z, amrex::Real z_0,
                        amrex::Real z_trop, amrex::Real z_top,
                        amrex::Real T_0, amrex::Real T_trop, amrex::Real T_top)
{
    if (z <= z_trop) {
      amrex::Real lapse = - (T_trop - T_0) / (z_trop - z_0);
      amrex::Real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
      return p_0 * std::pow( T / T_0 , CONST_GRAV/(R_d*lapse) );
    } else {
      // Get pressure at the tropopause
      amrex::Real lapse = - (T_trop - T_0) / (z_trop - z_0);
      amrex::Real p_trop = p_0 * std::pow( T_trop / T_0 , CONST_GRAV/(R_d*lapse) );
      // Get pressure at requested height
      lapse = - (T_top - T_trop) / (z_top - z_trop);
      if (lapse != 0) {
        amrex::Real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
        return p_trop * std::pow( T / T_trop , CONST_GRAV/(R_d*lapse) );
      } else {
        return p_trop * std::exp(-CONST_GRAV*(z-z_trop)/(R_d*T_trop));
      }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_sat_mix(amrex::Real press, amrex::Real T ) {
    return 380./(press) * std::exp(17.27*(T-273.)/(T-36.));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_relhum(amrex::Real z, amrex::Real z_trop)
{
  if (z <= z_trop) {
      return 1. - 0.75 * pow(z / z_trop , 1.25);
  } else {
      return 0.25;
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
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
#if defined(ERF_USE_MOISTURE)
    Array4<Real      > const& qv,
    Array4<Real      > const& qc,
    Array4<Real      > const& qi,
#elif defined(ERF_USE_WARM_NO_PRECIP)
    Array4<Real      > const&   ,
    Array4<Real      > const&   ,
#endif
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];
  const amrex::Real& prob_hi_z = geomdata.ProbHi()[2];

  const amrex::Real rdOcp   = sc.rdOcp;

  // const amrex::Real thetatr = 343.0;
  // const amrex::Real theta0  = 300.0;
  const amrex::Real ztr     = 12000.0;
  const amrex::Real Ttr     = 213.0;
  const amrex::Real Ttop    = 213.0;
  const amrex::Real deltaz  = 1000.*0.0;
  const amrex::Real zs      = 5000.;
  const amrex::Real us      = 30.;
  const amrex::Real uc      = 15.;

  amrex::ParallelForRNG(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo         = geomdata.ProbLo();
    const auto dx              = geomdata.CellSize();
    const amrex::Real z        = prob_lo[2] + (k + 0.5) * dx[2];

    amrex::Real relhum = init_supercell_relhum(z, ztr);
    amrex::Real temp   = init_supercell_temperature(z, prob_lo_z, ztr, prob_hi_z, parms.T_0, Ttr, Ttop);
    amrex::Real press  = init_supercell_pressure(z, prob_lo_z, ztr, prob_hi_z, parms.T_0, Ttr, Ttop);
    amrex::Real qvs    = init_supercell_sat_mix(press, temp);

#if 0
    amrex::Real thetaeq;
    amrex::Real faceq;
    if (z <= ztr) {
       thetaeq = theta0 + (thetatr - theta0)*std::pow(z/ztr, 5./4.);
    }
    else {
       thetaeq = thetatr*std::exp(CONST_GRAV*(z-ztr)/c_p*Ttr);
    }

    faceq = 1.;
    for (int km = 0; km < k; ++km) {
       amrex::Real zloc = prob_lo[2]+(km+0.5)*dx[2];
       amrex::Real theq;
      if (zloc <= ztr) {
       theq = theta0 + (thetatr - theta0)*std::pow(zloc/ztr, 5./4.);
      } else {
       theq = thetatr*std::exp(CONST_GRAV*(zloc-ztr)/c_p*Ttr);
      }
       faceq -= CONST_GRAV/(c_p*theq)*dz;
    }

    temp = faceq*thetaeq;
#endif

    if (relhum*qvs > 0.014) relhum = 0.014/qvs;
    amrex::Real qvapor  = std::min(0.014, qvs*relhum);
    amrex::Real rho     = press/(R_d+qvapor*R_v)/temp;

    // perturb theta
    amrex::Real rand_double = amrex::Random(engine) - 1.0;        // Random number in [-1,1]
    amrex::Real scaling = (khi-static_cast<amrex::Real>(k))/khi;  // Less effect at higher levels
    amrex::Real deltaT = parms.T_pert*scaling*rand_double;

    amrex::Real theta = getThgivenRandT(rho, temp+deltaT, rdOcp);

    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = rho*theta;
    state(i, j, k, Rho_comp)      = rho;

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;

    // mean states
#if defined(ERF_USE_MOISTURE)
    state(i, j, k, RhoQt_comp) = rho*qvapor;
    state(i, j, k, RhoQp_comp) = 0.0;
    qv(i, j, k) = qvapor;
    qc(i, j, k) = 0.0;
    qi(i, j, k) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    state(i, j, k, RhoQv_comp) = rho*qvapor;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real z = prob_lo_z + (k+0.5) * dz;
    if (z < zs-deltaz) {
      x_vel(i, j, k) = us*(z/zs) - uc;
    } else if (std::abs(z-zs) < deltaz) {
      x_vel(i, j, k) = (-0.8+3.*(z/zs)-1.25*(z/zs)*(z/zs))*us-uc;
    } else {
      x_vel(i, j, k) = us-uc;
    }
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
erf_init_rayleigh(amrex::Vector<amrex::Real>& tau,
                  amrex::Vector<amrex::Real>& ubar,
                  amrex::Vector<amrex::Real>& vbar,
                  amrex::Vector<amrex::Real>& wbar,
                  amrex::Vector<amrex::Real>& thetabar,
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
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
}
