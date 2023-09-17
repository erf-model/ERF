#include "prob.H"
#include "prob_common.H"

#include "IndexDefines.H"
#include "ERF_Constants.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

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
          wbar[k]     = parms.W_0; //0.0;
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
init_custom_prob(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
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
    amrex::GeometryData const& geomdata,
    Array4<Real const> const&,
    Array4<Real const> const&,
    Array4<Real const> const&,
    const SolverChoice& sc)
{

// QKE for PBL
/*
  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      // Set initial value for QKE
      state(i, j, k, RhoQKE_comp) = 0.1; //parms.QKE_0;
  });
*/

// Initialize vortex here

// u-velocity component
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();
       
      const Real x = prob_lo[0] + i * dx[0]; // face center
      const Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center
      const Real z = prob_lo[2] + (k + 0.5) * dx[2]; // cell center
      // const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);

                     //(parms.M_inf * std::cos(parms.alpha)
                     //- (y - parms.yc)/parms.R * Omg) * parms.a_inf;

      // Zero-out the velocity
      x_vel(i, j, k) = 0;
        
      // Get vortex location
      const Real Xc = parms.Xc_0;
      const Real Yc = parms.Yc_0;
      // Calculate u-velocity for vortex
      const Real rr = std::pow((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc),0.5); // Radius from the center
      const Real v_max = parms.VMAX; // 15; // Maximum horizontal velocity in vortex
      const Real R_max = parms.RMAX; // 100; // Radius of maximum winds
      const Real R_0 = parms.RZERO; //800; // Radius of zero wind speed
      const Real z_0 = parms.ZZERO; //2000; // Height of zero wind speed

      if (z > z_0) {
          x_vel(i, j, k) = 0.0;
      } else {
            if (rr > R_0) {
                  x_vel(i, j, k) = 0.0;
            } else {
                  const Real II = (z_0-z)/z_0;
                  const Real term1 = (v_max*v_max)*(rr/R_max)*(rr/R_max);
                  const Real term2 = std::pow(2*R_max/(rr+R_max),3) - std::pow(2*R_max/(R_0+R_max),3);
                  const Real term3 = std::pow(2*7.2921*std::pow(10,-5)*std::sin(20*3.1415/180),2)*(rr*rr)/4;
                  const Real term4 = 2*7.2921*std::pow(10,-5)*std::sin(20*3.1415/180)*rr/2 ;
                  const Real v_tang = II*(std::pow(term1*term2 + term3,0.5) - term4);
                  const Real thet_angl = std::atan2(y-Yc,x-Xc);
                  x_vel(i, j, k) = -1*std::abs(v_tang)*std::sin(thet_angl);
            }
      }

      state(i, j, k, RhoQKE_comp) = parms.QKE_0;
      amrex::Print() <<"QKE="<<state(i, j, k, RhoQKE_comp)<< std::endl;
  });

// v-velocity component
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();

      const Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
      const Real y = prob_lo[1] + j * dx[1]; // face center
      const Real z = prob_lo[2] + (k + 0.5) * dx[2]; // cell center
      // Zero-out the velocity
      y_vel(i, j, k) = 0;

      // Get vortex location
      const Real Xc = parms.Xc_0;
      const Real Yc = parms.Yc_0;
      // Calculate v-velocity for vortex
      const Real rr = std::pow((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc),0.5); // Radius from the center
      const Real v_max = parms.VMAX; // 15; // Maximum horizontal velocity in vortex
      const Real R_max = parms.RMAX; // 100; // Radius of maximum winds
      const Real R_0 = parms.RZERO; //800; // Radius of zero wind speed
      const Real z_0 = parms.ZZERO; //2000; // Height of zero wind speed

      if (z > z_0) {
          y_vel(i, j, k) = 0.0;
      } else {
            if (rr > R_0) {
                  y_vel(i, j, k) = 0.0;
            } else {
                  const Real II = (z_0-z)/z_0;
                  const Real term1 = (v_max*v_max)*(rr/R_max)*(rr/R_max);
                  const Real term2 = std::pow(2*R_max/(rr+R_max),3) - std::pow(2*R_max/(R_0+R_max),3);
                  const Real term3 = std::pow(2*7.2921*std::pow(10,-5)*std::sin(20*3.1415/180),2)*(rr*rr)/4;
                  const Real term4 = 2*7.2921*std::pow(10,-5)*std::sin(20*3.1415/180)*rr/2 ;
                  const Real v_tang = II*(std::pow(term1*term2 + term3,0.5) - term4);
                  const Real thet_angl = std::atan2(y-Yc,x-Xc);
                  y_vel(i, j, k) = std::abs(v_tang)*std::cos(thet_angl);
            }
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

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

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

  pp.query("Xc_0", parms.Xc_0);
  pp.query("Yc_0", parms.Yc_0);
  pp.query("VMAX", parms.VMAX);
  pp.query("RMAX", parms.RMAX);
  pp.query("RZERO", parms.RZERO);
  pp.query("ZZERO", parms.ZZERO);


  pp.query("QKE_0", parms.QKE_0);
//  pp.query("U_0", parms.U_0);
//  pp.query("V_0", parms.V_0);
//  pp.query("W_0", parms.W_0);
  pp.query("dampcoef", parms.dampcoef);
  pp.query("zdamp", parms.zdamp);
}
