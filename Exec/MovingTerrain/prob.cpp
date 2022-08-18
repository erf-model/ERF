#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "IndexDefines.H"
#include "TerrainMetrics.H"

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
  const Real thetabar = T_sfc;

  if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

  for ( MFIter mfi(rho_hse, TilingIfNotGPU()); mfi.isValid(); ++mfi )
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

         init_isentropic_hse(i,j,rho_sfc,thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

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
  const amrex::Box& bx,
  Array4<Real      > const& state,
  Array4<Real      > const& x_vel,
  Array4<Real      > const& y_vel,
  Array4<Real      > const& z_vel,
  Array4<Real      > const& r_hse,
  Array4<Real      > const& p_hse,
  Array4<Real const> const& z_nd,
  Array4<Real const> const& z_cc,
  GeometryData const& geomdata)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  const Real T_sfc    = parms.T_0;
  const Real rho_sfc  = p_0 / (R_d*T_sfc);
  const Real thetabar = T_sfc;

  // These are at cell centers (unstaggered)
  amrex::Vector<Real> h_r(khi+1);
  amrex::Vector<Real> h_p(khi+1);

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

  Real H           = geomdata.ProbHi()[2];
  Real Ampl        = parms.Ampl;
  Real wavelength  = 100.;
  Real kp          = 2. * 3.1415926535 / wavelength;
  Real g           = 9.8;
  Real omega       = std::sqrt(g * kp);

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = z_cc(i,j,k);

    state(i, j, k, Rho_comp) = r_hse(i,j,k);

    Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
    Real p_prime = -(omega * omega / kp) * fac * Ampl * std::sin(kp * x);
    Real p_total = p_hse(i,j,k) + p_prime;

    // Define (rho theta) given pprime
    state(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p_total);

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;

#ifdef ERF_USE_MOISTURE
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const auto prob_lo  = geomdata.ProbLo();
      const auto dx       = geomdata.CellSize();

      const Real x = prob_lo[0] + i * dx[0];
      const Real z = 0.5 * (z_nd(i,j,k) + z_nd(i,j,k+1));

      Real fac    = std::cosh( kp * (z - H) ) / std::sinh(kp * H);

      x_vel(i, j, k) = -omega * fac * Ampl * std::sin(k * x);
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on z-faces
  const auto dx         = geomdata.CellSize();
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
  dxInv[0] = 1. / dx[0];
  dxInv[1] = 1. / dx[1];
  dxInv[2] = 1. / dx[2];

  // Set the z-velocity from impenetrable condition
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      Real x   = (i + 0.5) * dx[0];
      Real z   = z_cc(i,j,k);
      Real fac = std::sinh( kp * (z - H) ) / std::sinh(kp * H);
      z_vel(i, j, k) = omega * fac * Ampl * std::cos(kp * x);
  });

  amrex::Gpu::streamSynchronize();

}

void
erf_init_rayleigh(amrex::Vector<Real>& /*tau*/,
                  amrex::Vector<Real>& /*ubar*/,
                  amrex::Vector<Real>& /*vbar*/,
                  amrex::Vector<Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for DensityCurrent problem");
}

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("Ampl", parms.Ampl);
  pp.query("T_0" , parms.T_0);
}

void
init_custom_terrain (const Geometry& geom,
                           MultiFab& z_phys_nd,
                     const Real& time)
{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();

    // Domain valid box (z_nd is nodal)
    const amrex::Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    int domlo_z = domain.smallEnd(2);

    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Populate bottom plane
    int k0 = domlo_z;

    Real Ampl        = parms.Ampl;
    Real wavelength  = 100.;
    Real kp          = 2. * 3.1415926535 / wavelength;
    Real g           = 9.8;
    Real omega       = std::sqrt(g * kp);

    for ( amrex::MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        amrex::Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Clip indices for ghost-cells
            int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

            // Location of nodes
            Real x = ii  * dx[0];

            // Wave heigght
            Real height = Ampl * std::sin(kp * x - omega * time);

            // Populate terrain height
            z_arr(i,j,k0) = height;
        });
    }
}

AMREX_GPU_DEVICE
Real
dhdt(int i, int j,
     const GpuArray<Real,AMREX_SPACEDIM> dx,
     const Real time_mt, const Real delta_t)
{
    Real Ampl        = parms.Ampl;
    Real wavelength  = 100.;
    Real kp          = 2. * 3.1415926535 / wavelength;
    Real g           = 9.8;
    Real omega       = std::sqrt(g * kp);

    // Location of w-faces
    const Real x  = (i + 0.5) * dx[0];

    Real h_old = Ampl * std::sin(kp * x - omega *  time_mt);
    Real h_new = Ampl * std::sin(kp * x - omega * (time_mt + delta_t));

    return ( (h_new - h_old) / delta_t );
}
