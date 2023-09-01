#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "TerrainMetrics.H"
#include "TileNoZ.H"

using namespace amrex;

ProbParm parms;

AMREX_GPU_DEVICE
static
void
init_isentropic_hse(int i, int j,
                    const Real& r_sfc, const Real& theta,
                          Real* r,           Real* p,
                    const Array4<Real const> z_cc,
                    const int& khi)
{
  // r_sfc / p_0 are the density / pressure at the surface
  int k0  = 0;
  Real hz = z_cc(i,j,k0);

  // Initial guess
  r[k0] = r_sfc;
  p[k0] = p_0 - hz * r[k0] * CONST_GRAV;

  int MAX_ITER = 10;
  Real TOL = 1.e-8;

  {
      // We do a Newton iteration to satisfy the EOS & HSE (with constant theta)
      bool converged_hse = false;
      Real p_hse;
      Real p_eos;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p_0 - hz * r[k0] * CONST_GRAV;
          p_eos = getPgivenRTh(r[k0]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k0],theta);

          Real drho = A / (dpdr + hz * CONST_GRAV);

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

      r[k] = r[k-1];

      Real p_eos = getPgivenRTh(r[k]*theta);
      Real p_hse;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          Real dz_loc = (z_cc(i,j,k) - z_cc(i,j,k-1));
          p_hse = p[k-1] -  dz_loc * 0.5 * (r[k-1]+r[k]) * CONST_GRAV;
          p_eos = getPgivenRTh(r[k]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k],theta);

          Real drho = A / (dpdr + dz_loc * CONST_GRAV);

          r[k] = r[k] + drho;
          p[k] = getPgivenRTh(r[k]*theta);

          if (std::abs(drho) < TOL)
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
                  std::unique_ptr<MultiFab>& /*z_phys_nd*/,
                  std::unique_ptr<MultiFab>& z_phys_cc,
                  Geometry const& geom)
{
  //const Real prob_lo_z = geom.ProbLo()[2];
  const int khi        = geom.Domain().bigEnd()[2];

  const Real T_sfc    = parms.T_0;
  const Real rho_sfc  = p_0 / (R_d*T_sfc);
  const Real Thetabar = T_sfc;

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

         init_isentropic_hse(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

         for (int k = 0; k <= khi; k++) {
            rho_arr(i,j,k) = r(k);
         }
         rho_arr(i,j,   -1) = rho_arr(i,j,0);
         rho_arr(i,j,khi+1) = rho_arr(i,j,khi);
       });
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
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& z_nd,
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
    const SolverChoice&)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
//const Real z0 = (0.5) * dx[2] + prob_lo[2];
//const Real Tbar = parms.T_0 - z0 * CONST_GRAV / parms.C_p;
//const Real pbar = p_0 * std::pow(Tbar/parms.T_0, R_d/parms.C_p); // from Straka1993
//const Real pbar = p_0 * std::pow(Tbar/parms.T_0, parms.C_p/R_d); // isentropic relation, consistent with exner pressure def
//const Real rhobar = pbar / (R_d*Tbar); // UNUSED

  const Real rho_sfc   = p_0 / (R_d*parms.T_0);
  const Real thetabar  = parms.T_0;

  // Create a flat box with same horizontal extent but only one cell in vertical
  Box b2d = surroundingNodes(bx); // Copy constructor
  b2d.setRange(2,0);

  ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
  {
     Array1D<Real,0,255> r;;
     Array1D<Real,0,255> p;;

     init_isentropic_hse(i,j,rho_sfc,thetabar,&(r(0)),&(p(0)),z_cc,khi);

     for (int k = 0; k <= khi; k++) {
        r_hse(i,j,k) = r(k);
        p_hse(i,j,k) = p(k);
     }
     r_hse(i,j,   -1) = r_hse(i,j,0);
     r_hse(i,j,khi+1) = r_hse(i,j,khi);
  });

  // Geometry (note we must include these here to get the data on device)
  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = z_cc(i,j,k);

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

    Real L = std::sqrt(
        std::pow((x - parms.x_c)/parms.x_r, 2) +
        std::pow((z - parms.z_c)/parms.z_r, 2)
    );
    Real dT;
    if (L > 1.0) {
        dT = 0.0;
    }
    else {
        dT = parms.T_pert * (std::cos(PI*L) + 1.0)/2.0;
    }

    // Note: dT is T perturbation, theta_perturbed is theta PLUS perturbation in theta
    Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), R_d/parms.C_p);

    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = std::pow(p_hse(i,j,k)/p_0,1.0/Gamma) * p_0 / R_d;
    state(i, j, k, Rho_comp) = state(i, j, k, RhoTheta_comp) / theta_perturbed;

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;

#if defined(ERF_USE_MOISTURE)
    state(i, j, k, RhoQt_comp) = 0.0;
    state(i, j, k, RhoQp_comp) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel(i, j, k) = parms.U_0;
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel(i, j, k) = 0.0;
  });

  const auto dx = geomdata.CellSize();
  amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
  dxInv[0] = 1. / dx[0];
  dxInv[1] = 1. / dx[1];
  dxInv[2] = 1. / dx[2];

  // Set the z-velocity from impenetrable condition
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = WFromOmega(i, j, k, 0.0, x_vel, y_vel, z_nd, dxInv);
  });

  amrex::Gpu::streamSynchronize();

}

void
erf_init_rayleigh(amrex::Vector<Real>& /*tau*/,
                  amrex::Vector<Real>& /*ubar*/,
                  amrex::Vector<Real>& /*vbar*/,
                  amrex::Vector<Real>& /*wbar*/,
                  amrex::Vector<Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for WitchOfAgnesi problem");
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

void
init_custom_terrain (const Geometry& geom,
                           MultiFab& z_phys_nd,
                     const Real& /*time*/)
{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    auto ProbHiArr = geom.ProbHiArray();

    // Domain valid box (z_nd is nodal)
    const amrex::Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    // int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2);

    // User function parameters
    Real a    = 0.5;
    Real num  = 8 * a * a * a;
    Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
    // Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Populate bottom plane
    int k0 = domlo_z;

    for ( amrex::MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        amrex::Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Clip indices for ghost-cells
            int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
            // int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

            // Location of nodes
            Real x = (ii  * dx[0] - xcen);
            // Real y = (jj  * dx[1] - ycen);

            // WoA Hill in x-direction
            Real height = num / (x*x + 4 * a * a);

            // Populate terrain height
            z_arr(i,j,k0) = height;
        });
    }
}
