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

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    //
    // NOTE: this routine is only called when doing custom initialization!
    //
    // In the case of initializng from an input sounding or from a wrfinput file,
    // this routine will not be called.
    //

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
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
    const Real a = 1.0 / DE;

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const auto dx = geomdata.CellSize();
        const Real z = (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = u_0 * (1.0 - std::exp(-a * z) * std::cos(-a * z));
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const auto dx = geomdata.CellSize();
        const Real z = (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = -u_0 * std::exp(-a * z) * std::sin(-a * z);
    });

    // Set the z-velocity
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        z_vel_pert(i, j, k) = 0.0;
    });
}
