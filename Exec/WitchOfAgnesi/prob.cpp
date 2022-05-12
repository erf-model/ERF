#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "IndexDefines.H"

using namespace amrex;

#ifdef ERF_USE_TERRAIN

ProbParm parms;

AMREX_GPU_DEVICE
void
init_isentropic_hse(int i, int j,
                    const Real& r_sfc, const Real& theta,
                          Real* r,           Real* p,
                    const Array4<Real const> z_nd,
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
          // Gamma * p_0 * std::pow( (R_d * theta / p_0), Gamma) * std::pow(r[k], Gamma-1.0) ;

          Real drho = A / (dpdr + 0.5 * dz_loc * CONST_GRAV);

          r[k] = std::max(0.9*r[k-1], std::min(r[k] + drho, 1.1*r[k-1]));
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
                  MultiFab const& z_phys_nd,
                  MultiFab const& z_phys_cc,
                  Geometry const& geom)
{
  const Real prob_lo_z = geom.ProbLo()[2];
  const int khi        = geom.Domain().bigEnd()[2];

  const Real& T_sfc    = parms.T_0;
  const Real& rho_sfc  = p_0 / (R_d*T_sfc);
  const Real& Thetabar = T_sfc;

  if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

  for ( MFIter mfi(rho_hse, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
       Array4<Real      > rho_arr  = rho_hse.array(mfi);
       Array4<Real const> z_nd_arr = z_phys_nd.const_array(mfi);
       Array4<Real const> z_cc_arr = z_phys_cc.const_array(mfi);

       // Create a flat box with same horizontal extent but only one cell in vertical
       const Box& tbz = mfi.nodaltilebox(2);
       Box b2d = tbz; // Copy constructor
       b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
       b2d.setRange(2,0);

       ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
         Array1D<Real,0,255> r;;
         Array1D<Real,0,255> p;;

         init_isentropic_hse(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_nd_arr,z_cc_arr,khi);

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
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::Array4<amrex::Real> const& r_hse,
  amrex::Array4<amrex::Real> const& p_hse,
  amrex::Array4<amrex::Real const> const& z_nd,
  amrex::Array4<amrex::Real const> const& z_cc,
  amrex::GeometryData const& geomdata)
{
  const int khi              = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
//const amrex::Real z0 = (0.5) * dx[2] + prob_lo[2];
//const amrex::Real Tbar = parms.T_0 - z0 * CONST_GRAV / parms.C_p;
//const amrex::Real pbar = p_0 * std::pow(Tbar/parms.T_0, R_d/parms.C_p); // from Straka1993
//const amrex::Real pbar = p_0 * std::pow(Tbar/parms.T_0, parms.C_p/R_d); // isentropic relation, consistent with exner pressure def
//const amrex::Real rhobar = pbar / (R_d*Tbar); // UNUSED

  const amrex::Real& rho_sfc   = p_0 / (R_d*parms.T_0);
  const amrex::Real& thetabar  = parms.T_0;
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];

  // These are at cell centers (unstaggered)
  amrex::Vector<amrex::Real> h_r(khi+1);
  amrex::Vector<amrex::Real> h_p(khi+1);

  // Create a flat box with same horizontal extent but only one cell in vertical
  Box b2d = surroundingNodes(bx); // Copy constructor
  b2d.setRange(2,0);

  ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
  {
     Array1D<Real,0,255> r;;
     Array1D<Real,0,255> p;;

     init_isentropic_hse(i,j,rho_sfc,thetabar,&(r(0)),&(p(0)),z_nd,z_cc,khi);

     for (int k = 0; k <= khi; k++) {
        r_hse(i,j,k) = r(k);
        p_hse(i,j,k) = p(k);
     }
     r_hse(i,j,   -1) = r_hse(i,j,0);
     r_hse(i,j,khi+1) = r_hse(i,j,khi);
  });

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo         = geomdata.ProbLo();
    const auto dx              = geomdata.CellSize();

    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real z = z_cc(i,j,k);

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const amrex::Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

    amrex::Real L = std::sqrt(
        std::pow((x - parms.x_c)/parms.x_r, 2) +
        std::pow((z - parms.z_c)/parms.z_r, 2)
    );
    amrex::Real dT;
    if (L > 1.0) {
        dT = 0.0;
    }
    else {
        dT = parms.T_pert * (std::cos(PI*L) + 1.0)/2.0;
    }

    // Note: dT is a perturbation in temperature, theta_perturbed is theta PLUS perturbation in theta
    amrex::Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), R_d/parms.C_p);

    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = std::pow(p_hse(i,j,k)/p_0,1.0/Gamma) * p_0 / R_d;
    state(i, j, k, Rho_comp) = state(i, j, k, RhoTheta_comp) / theta_perturbed;

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
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    x_vel(i, j, k) = parms.U_0;
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
}

#else

ProbParm parms;

void
init_isentropic_hse(const Real& r_sfc, const Real& theta,
                          Real* r,           Real* p,
                    const Real& dz,    const Real&  prob_lo_z,
                    const int& khi)
{
  // r_sfc / p_0 are the density / pressure at the surface
  r[0] = r_sfc;
  p[0] = p_0 - (0.5*dz) * r[0] * CONST_GRAV;

  int MAX_ITER = 10;
  Real TOL = 1.e-8;

  int k0 = 0;
  {
      // We do a Newton iteration to satisfy the EOS (with constant theta) and to discretely satisfy HSE
      bool converged_hse = false;

      // Initial guess
      r[k0] = r_sfc;
      p[k0] = getPgivenRTh(r[k0]*theta);

      Real p_eos = getPgivenRTh(r[k0]*theta);
      Real p_hse;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p_0 -  (0.5*dz) * r_sfc * CONST_GRAV;
          p_eos = getPgivenRTh(r[k0]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k0],theta);

          Real drho = A / (dpdr + 0.5 * dz * CONST_GRAV);

          r[k0] = std::max(0.9*r_sfc, std::min(r[k0] + drho, 1.1*r_sfc));
          p[k0] = getPgivenRTh(r[k0]*theta);

          if (std::abs(drho) < TOL * r_sfc)
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
erf_init_dens_hse(Real* dens_hse_ptr,
                  Geometry const& geom,
                  const int /*ng_dens_hse*/)
{
  const Real prob_lo_z = geom.ProbLo()[2];
  const Real dz        = geom.CellSize()[2];
  const int khi        = geom.Domain().bigEnd()[2];

  const Real& T_sfc    = 300.;
  const Real& rho_sfc  = p_0 / (R_d*T_sfc);
  const Real& Thetabar = T_sfc;

  Vector<Real> r(khi+1);
  Vector<Real> p(khi+1);

  init_isentropic_hse(rho_sfc,Thetabar,r.data(),p.data(),dz,prob_lo_z,khi);

  for (int k = 0; k <= khi; k++)
  {
      dens_hse_ptr[k] = r[k];
  }

  dens_hse_ptr[   -1] = dens_hse_ptr[  0];
  dens_hse_ptr[khi+1] = dens_hse_ptr[khi];
}

void
init_custom_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata)
{
  const int khi              = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
//const amrex::Real z0 = (0.5) * dx[2] + prob_lo[2];
//const amrex::Real Tbar = parms.T_0 - z0 * CONST_GRAV / parms.C_p;
//const amrex::Real pbar = p_0 * std::pow(Tbar/parms.T_0, R_d/parms.C_p); // from Straka1993
//const amrex::Real pbar = p_0 * std::pow(Tbar/parms.T_0, parms.C_p/R_d); // isentropic relation, consistent with exner pressure def
//const amrex::Real rhobar = pbar / (R_d*Tbar); // UNUSED

  const amrex::Real& rho_sfc   = p_0 / (R_d*parms.T_0);
  const amrex::Real& thetabar  = parms.T_0;
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];

  // These are at cell centers (unstaggered)
  amrex::Vector<amrex::Real> h_r(khi+1);
  amrex::Vector<amrex::Real> h_p(khi+1);

  amrex::Gpu::DeviceVector<amrex::Real> d_r(khi+1);
  amrex::Gpu::DeviceVector<amrex::Real> d_p(khi+1);

  init_isentropic_hse(rho_sfc,thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

  amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
  amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

  amrex::Real* r = d_r.data();
  amrex::Real* p = d_p.data();

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo         = geomdata.ProbLo();
    const auto dx              = geomdata.CellSize();

    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const amrex::Real Tbar_hse = p[k] / (R_d * r[k]);

    amrex::Real L = std::sqrt(
        std::pow((x - parms.x_c)/parms.x_r, 2) +
        std::pow((z - parms.z_c)/parms.z_r, 2)
    );
    amrex::Real dT;
    if (L > 1.0) {
        dT = 0.0;
    }
    else {
        dT = parms.T_pert * (std::cos(PI*L) + 1.0)/2.0;
    }

    // Note: dT is a perturbation in temperature, theta_perturbed is theta PLUS perturbation in theta
    amrex::Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p[k], R_d/parms.C_p);

    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = std::pow(p[k]/p_0,1.0/Gamma) * p_0 / R_d;
    state(i, j, k, Rho_comp) = state(i, j, k, RhoTheta_comp) / theta_perturbed;

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;
  });

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    x_vel(i, j, k) = parms.U_0;
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
}

#endif

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& /*tau*/,
                  amrex::Vector<amrex::Real>& /*ubar*/,
                  amrex::Vector<amrex::Real>& /*vbar*/,
                  amrex::Vector<amrex::Real>& /*thetabar*/,
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
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
}
