#include "prob.H"
#include "prob_common.H"

#include "IndexDefines.H"
#include "ERF_Constants.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

ProbParm parms;

#ifdef ERF_USE_TERRAIN
void
erf_init_dens_hse(amrex::MultiFab& rho_hse,
                  const amrex::MultiFab& z_phys_nd,
                  const amrex::MultiFab& z_phys_cc,
                  amrex::Geometry const& geom)
{
    amrex::Real R0 = parms.rho_0;
    for ( amrex::MFIter mfi(rho_hse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
       amrex::Array4<amrex::Real> rho_arr = rho_hse.array(mfi);
       const amrex::Box& gbx = mfi.growntilebox(1);
       amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
           rho_arr(i,j,k) = R0;
       });
    }
}
#else
void
erf_init_dens_hse(amrex::Real* dens_hse_ptr,
                  amrex::Geometry const& geom,
                  const int ng_dens_hse)
{
  const int khi = geom.Domain().bigEnd()[2];
  for (int k = -ng_dens_hse; k <= khi+ng_dens_hse; k++)
  {
      dens_hse_ptr[k] = parms.rho_0;
  }
}
#endif

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& /*tau*/,
                  amrex::Vector<amrex::Real>& /*ubar*/,
                  amrex::Vector<amrex::Real>& /*vbar*/,
                  amrex::Vector<amrex::Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for Ekman Spiral problem");
}

void
erf_init_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
#ifdef ERF_USE_TERRAIN
  amrex::Array4<amrex::Real> const& r_hse,
  amrex::Array4<amrex::Real> const& p_hse,
  amrex::Array4<amrex::Real const> const& z_nd,
  amrex::Array4<amrex::Real const> const& z_cc,
#endif
  amrex::GeometryData const& geomdata)
{
  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Geometry
    //const amrex::Real* prob_lo = geomdata.ProbLo();
    //const auto dx = geomdata.CellSize();
    //const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    //const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    //const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial potential temperature (Actually rho*theta)
    state(i, j, k, RhoTheta_comp) = parms.rho_0 * parms.Theta_0;



    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;
  });

  amrex::ParmParse pp("erf");
  amrex::Real rot_time_period;
  pp.get("rotational_time_period", rot_time_period);
  amrex::Real coriolis_factor = 4.0 * PI / rot_time_period;

  amrex::Real Az;
  pp.get("dynamicViscosity", Az); // dynamic viscosity [kg-m/s]
  Az = Az / parms.rho_0; // kinematic viscosity [m^2/s]

  amrex::Vector<amrex::Real> abl_geo_wind(3);
  pp.queryarr("abl_geo_wind",abl_geo_wind);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    (amrex::Math::abs(abl_geo_wind[1]) < 1.0e-15) &&
    (amrex::Math::abs(abl_geo_wind[2]) < 1.0e-15),
    "Ekman Spiral uses geostrophic forcing of the form (V_0, 0, 0)");
  const amrex::Real u_0 = abl_geo_wind[0];

  const amrex::Real DE = std::sqrt(2.0 * Az / coriolis_factor);
  //amrex::Print() << "Ekman depth = " << DE << " m" << std::endl;
  const amrex::Real a = 1.0 / DE;

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    const auto dx = geomdata.CellSize();
    const amrex::Real z = (k + 0.5) * dx[2];

    // Set the x-velocity
    x_vel(i, j, k) = u_0 * (1.0 - std::exp(-a * z) * std::cos(-a * z));
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    const auto dx = geomdata.CellSize();
    const amrex::Real z = (k + 0.5) * dx[2];

    // Set the y-velocity
    y_vel(i, j, k) = -u_0 * std::exp(-a * z) * std::sin(-a * z);
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });
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
}
