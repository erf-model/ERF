#include "prob.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);

  pp.query("Xc_0", parms.Xc_0);
  pp.query("Yc_0", parms.Yc_0);
  pp.query("VMAX", parms.VMAX);
  pp.query("RMAX", parms.RMAX);
  pp.query("RZERO", parms.RZERO);
  pp.query("ZZERO", parms.ZZERO);

  pp.query("QKE_0", parms.QKE_0);
  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("dampcoef", parms.dampcoef);
  pp.query("zdamp", parms.zdamp);
  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
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

  Real QKE_0 = parms.QKE_0;
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      state(i, j, k, RhoQKE_comp) = QKE_0;
      amrex::Print() <<"QKE="<<state(i, j, k, RhoQKE_comp)<< std::endl;
  });

// Initialize vortex here

  // Get vortex location
  const Real Xc = parms.Xc_0;
  const Real Yc = parms.Yc_0;
  const Real v_max = parms.VMAX; // 15; // Maximum horizontal velocity in vortex
  const Real R_max = parms.RMAX; // 100; // Radius of maximum winds
  const Real R_0 = parms.RZERO; //800; // Radius of zero wind speed
  const Real z_0 = parms.ZZERO; //2000; // Height of zero wind speed

// u-velocity component
  amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();

      const Real x = prob_lo[0] + i * dx[0]; // face center
      const Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center
      const Real z = prob_lo[2] + (k + 0.5) * dx[2]; // cell center

      // Zero-out the velocity
      x_vel(i, j, k) = 0;

      if (z > z_0) {
          x_vel(i, j, k) = 0.0;
      } else {
            const Real rr = std::pow((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc),0.5); // Radius from the center
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
  });

// v-velocity component
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();

      const Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
      const Real y = prob_lo[1] + j * dx[1]; // face center
      const Real z = prob_lo[2] + (k + 0.5) * dx[2]; // cell center

      // Zero-out the velocity
      y_vel(i, j, k) = 0;

      if (z > z_0) {
          y_vel(i, j, k) = 0.0;
      } else {
            const Real rr = std::pow((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc),0.5); // Radius from the center
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
Problem::init_custom_terrain(
    const Geometry& /*geom*/,
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
