#include "prob.H"
#include "prob_common.H"

#include "ERF_Constants.H"

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
  pp.query("T_0", parms.Theta_0);
}

void
Problem::init_custom_prob(
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
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice&)
{
  ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Geometry
    //const Real* prob_lo = geomdata.ProbLo();
    //const auto dx = geomdata.CellSize();
    //const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    //const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    //const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Set the density
    state(i, j, k, Rho_comp) = parms.rho_0;

    // Initial potential temperature (Actually rho*theta)
    state(i, j, k, RhoTheta_comp) = parms.rho_0 * parms.Theta_0;

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

  ParmParse pp("erf");
  Real rot_time_period;
  pp.get("rotational_time_period", rot_time_period);
  Real coriolis_factor = 4.0 * PI / rot_time_period;

  Real Az;
  pp.get("dynamicViscosity", Az); // dynamic viscosity [kg-m/s]
  Az = Az / parms.rho_0; // kinematic viscosity [m^2/s]

  Vector<Real> abl_geo_wind(3);
  pp.queryarr("abl_geo_wind",abl_geo_wind);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    (amrex::Math::abs(abl_geo_wind[1]) < 1.0e-15) &&
    (amrex::Math::abs(abl_geo_wind[2]) < 1.0e-15),
    "Ekman Spiral uses geostrophic forcing of the form (V_0, 0, 0)");
  const Real u_0 = abl_geo_wind[0];

  const Real DE = std::sqrt(2.0 * Az / coriolis_factor);
  //amrex::Print() << "Ekman depth = " << DE << " m" << std::endl;
  const Real a = 1.0 / DE;

  // Set the x-velocity
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    const auto dx = geomdata.CellSize();
    const Real z = (k + 0.5) * dx[2];

    // Set the x-velocity
    x_vel(i, j, k) = u_0 * (1.0 - std::exp(-a * z) * std::cos(-a * z));
  });

  // Set the y-velocity
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    const auto dx = geomdata.CellSize();
    const Real z = (k + 0.5) * dx[2];

    // Set the y-velocity
    y_vel(i, j, k) = -u_0 * std::exp(-a * z) * std::sin(-a * z);
  });

  // Set the z-velocity
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
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
