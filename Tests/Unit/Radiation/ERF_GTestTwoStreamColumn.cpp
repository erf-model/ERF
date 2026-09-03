#include <cmath>
#include <vector>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_RadStruct.H>
#include <ERF_TwoStreamColumn.H>

// Two-stream column unit-test contract
// ------------------------------------
// These tests run vertical_two_stream_sweep() on a single column and pin
// the properties that a correctly oriented, correctly signed sweep must have:
//   1. Layer temperature is the absolute temperature from the equation of
//      state (Exner function), not the potential temperature.
//   2. k = kmin is the surface layer and k = kmax the top layer: SW heating
//      is positive everywhere and decreases monotonically from the top layer
//      to the surface layer in a uniform-density column.
//   3. The absorbed SW surface flux equals the Beer-Lambert direct beam
//      through the whole column times (1 - albedo).
//   4. In an isothermal column over a black surface at the same temperature
//      the LW heating is non-positive everywhere, strongest at the top layer
//      (cooling to space), and the surface net LW equals the analytic
//      sigma*T^4 * exp(-tau_column).
//   5. A gray surface reflects (1 - eps) of the downwelling LW, so the
//      surface net LW equals eps * sigma*T^4 * exp(-tau_column) rather than
//      the large negative value obtained without the reflected term.
//   6. isothermal_test zeroes the LW heating; night zeroes the SW heating.
//
// Portability: the sweep is a device function, so it is launched through
// amrex::ParallelFor on a one-cell horizontal box and results are copied
// back to the host before any GTest assertion.

namespace {

constexpr int kNz = 8;
constexpr amrex::Real kDz = 100.0;
constexpr amrex::Real kSigma = 5.670374419e-8;

struct ColumnResult {
    std::vector<amrex::Real> q_sw;
    std::vector<amrex::Real> q_lw;
    amrex::Real max_heating = 0.0;
    amrex::Real sw_surface = 0.0;
    amrex::Real lw_net_surface = 0.0;
};

// Run the sweep for a column with uniform density `rho` and uniform
// absolute temperature `T_air` (converted to rho*theta through the EOS).
ColumnResult run_uniform_column (const RadChoice& rad_choice, amrex::Real rho, amrex::Real T_air)
{
    const amrex::Box bx(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, kNz - 1));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {kDz, kDz, kNz * kDz});
    const int is_periodic[3] = {1, 1, 0};
    const amrex::Geometry geom(bx, &real_box, 0, is_periodic);

    amrex::FArrayBox state(bx, 2);
    amrex::FArrayBox qheating(bx, 2);
    const amrex::Real theta = getThgivenRandT(rho, T_air, RdoCp);
    state.setVal<amrex::RunOn::Device>(rho, bx, Rho_comp, 1);
    state.setVal<amrex::RunOn::Device>(rho * theta, bx, RhoTheta_comp, 1);
    qheating.setVal<amrex::RunOn::Device>(0.0);

    amrex::Gpu::DeviceVector<amrex::Real> scalars(3, 0.0);
    amrex::Real* scalar_ptr = scalars.data();
    const auto state_arr = state.const_array();
    const auto qheating_arr = qheating.array();
    const amrex::Array4<const amrex::Real> no_z_phys{};

    const amrex::Box xy_box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
    {
        amrex::Real max_heating = 0.0;
        amrex::Real sw_surface = 0.0;
        amrex::Real lw_net = 0.0;
        vertical_two_stream_sweep(i, j, bx, geom, state_arr, rad_choice, /*cloudy=*/false,
                                  qheating_arr, max_heating, sw_surface, lw_net, no_z_phys);
        scalar_ptr[0] = max_heating;
        scalar_ptr[1] = sw_surface;
        scalar_ptr[2] = lw_net;
    });
    amrex::Gpu::streamSynchronize();

    ColumnResult result;
    std::vector<amrex::Real> host_scalars(3);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, scalars.begin(), scalars.end(), host_scalars.begin());
    result.max_heating = host_scalars[0];
    result.sw_surface = host_scalars[1];
    result.lw_net_surface = host_scalars[2];

    amrex::FArrayBox host_q(bx, 2, amrex::The_Pinned_Arena());
    host_q.copy<amrex::RunOn::Device>(qheating);
    amrex::Gpu::streamSynchronize();
    const auto hq = host_q.const_array();
    for (int k = 0; k < kNz; ++k) {
        result.q_sw.push_back(hq(0, 0, k, 0));
        result.q_lw.push_back(hq(0, 0, k, 1));
    }
    return result;
}

RadChoice base_choice ()
{
    RadChoice rc;
    rc.rad_type = RadType::TwoStream;
    rc.sw_enabled = true;
    rc.lw_enabled = true;
    rc.tau_per_layer = 0.05;
    rc.tau_lw_per_layer = 1.0;
    rc.solar_zenith_deg = 60.0;
    rc.S0 = 1361.0;
    rc.surface_albedo_sw = 0.3;
    rc.surface_emissivity_lw = 1.0;
    rc.surface_temp_k = 290.0;
    rc.isothermal_test = false;
    return rc;
}

} // namespace

TEST(TwoStreamColumn, TemperatureComesFromExnerFunction)
{
    const amrex::Real rho = 1.0;
    const amrex::Real theta = 300.0;
    const amrex::Real T = get_temperature_from_rhotheta(rho * theta, rho);

    // p = p_0 (R_d rho theta / p_0)^gamma is below p_0 for this state, so
    // T = theta (p/p_0)^(R_d/c_p) must be well below theta.
    EXPECT_LT(T, theta - 5.0);
    EXPECT_NEAR(T, getTgivenRandRTh(rho, rho * theta), 1.0e-9);

    // Round trip through the EOS.
    const amrex::Real theta_back = getThgivenRandT(rho, T, RdoCp);
    EXPECT_NEAR(theta_back, theta, 1.0e-8 * theta);

    // Defensive fallbacks for unphysical input.
    EXPECT_EQ(get_temperature_from_rhotheta(-1.0, rho), 288.15);
    EXPECT_EQ(get_temperature_from_rhotheta(rho * theta, 0.0), 288.15);
}

TEST(TwoStreamColumn, ShortwaveHeatingDecreasesFromTopToSurface)
{
    const RadChoice rc = base_choice();
    const ColumnResult r = run_uniform_column(rc, 1.0, 290.0);

    ASSERT_EQ(static_cast<int>(r.q_sw.size()), kNz);
    for (int k = 0; k < kNz; ++k) {
        EXPECT_GT(r.q_sw[k], 0.0) << "k = " << k;
    }
    // Orientation: the beam enters at k = kmax and is weaker at every lower
    // layer, so with uniform density the heating must decrease toward k = 0.
    for (int k = kNz - 1; k > 0; --k) {
        EXPECT_GT(r.q_sw[k], r.q_sw[k - 1]) << "k = " << k;
    }

    // Absorbed surface flux: Beer-Lambert through kNz layers, times (1 - albedo).
    const amrex::Real mu0 = std::cos(rc.solar_zenith_deg * M_PI / 180.0);
    const amrex::Real expected = rc.S0 * mu0 * std::exp(-kNz * rc.tau_per_layer / mu0)
                                 * (1.0 - rc.surface_albedo_sw);
    EXPECT_NEAR(r.sw_surface, expected, 1.0e-9 * expected);
    EXPECT_GT(r.max_heating, 0.0);
}

TEST(TwoStreamColumn, LongwaveCoolsToSpaceFromTheTopLayer)
{
    const RadChoice rc = base_choice();
    const amrex::Real T = rc.surface_temp_k;   // air and surface at the same T
    const ColumnResult r = run_uniform_column(rc, 1.0, T);

    // An isothermal column over a black surface at the same temperature can
    // only lose energy to space, so no layer warms and the top layer cools most.
    for (int k = 0; k < kNz; ++k) {
        EXPECT_LE(r.q_lw[k], 1.0e-12) << "k = " << k;
    }
    for (int k = 0; k < kNz - 1; ++k) {
        EXPECT_LT(r.q_lw[kNz - 1], r.q_lw[k]) << "k = " << k;
    }
    EXPECT_LT(r.q_lw[kNz - 1], 0.0);

    // Surface: F_up(0) = sigma T^4, F_down(0) = sigma T^4 (1 - exp(-tau_col)).
    const amrex::Real B = kSigma * T * T * T * T;
    const amrex::Real tau_col = kNz * rc.tau_lw_per_layer;
    EXPECT_NEAR(r.lw_net_surface, B * std::exp(-tau_col), 1.0e-9 * B);
}

TEST(TwoStreamColumn, GraySurfaceReflectsDownwellingLongwave)
{
    RadChoice rc = base_choice();
    rc.surface_emissivity_lw = 0.5;
    const amrex::Real T = rc.surface_temp_k;
    const ColumnResult r = run_uniform_column(rc, 1.0, T);

    // F_up(0) = eps B + (1 - eps) F_down(0), F_down(0) = B (1 - exp(-tau_col))
    // => F_up(0) - F_down(0) = eps B exp(-tau_col).
    const amrex::Real B = kSigma * T * T * T * T;
    const amrex::Real tau_col = kNz * rc.tau_lw_per_layer;
    const amrex::Real expected = rc.surface_emissivity_lw * B * std::exp(-tau_col);
    EXPECT_NEAR(r.lw_net_surface, expected, 1.0e-9 * B);

    // Without the reflected term the surface would appear to absorb roughly
    // (1 - eps) B, i.e. a net flux of order -0.5 B. Guard against that.
    EXPECT_GT(r.lw_net_surface, -1.0e-6 * B);
}

TEST(TwoStreamColumn, IsothermalTestModeZeroesLongwaveHeating)
{
    RadChoice rc = base_choice();
    rc.isothermal_test = true;
    const ColumnResult r = run_uniform_column(rc, 1.0, 290.0);
    for (int k = 0; k < kNz; ++k) {
        EXPECT_EQ(r.q_lw[k], 0.0) << "k = " << k;
    }
    EXPECT_EQ(r.lw_net_surface, 0.0);
}

TEST(TwoStreamColumn, NightHasNoShortwave)
{
    RadChoice rc = base_choice();
    rc.solar_zenith_deg = 120.0;   // sun below the horizon
    const ColumnResult r = run_uniform_column(rc, 1.0, 290.0);
    for (int k = 0; k < kNz; ++k) {
        EXPECT_EQ(r.q_sw[k], 0.0) << "k = " << k;
    }
    EXPECT_EQ(r.sw_surface, 0.0);
}

TEST(TwoStreamColumn, DisabledBandsWriteZeroHeating)
{
    RadChoice rc = base_choice();
    rc.sw_enabled = false;
    rc.lw_enabled = false;
    const ColumnResult r = run_uniform_column(rc, 1.0, 290.0);
    for (int k = 0; k < kNz; ++k) {
        EXPECT_EQ(r.q_sw[k], 0.0);
        EXPECT_EQ(r.q_lw[k], 0.0);
    }
    EXPECT_EQ(r.sw_surface, 0.0);
    EXPECT_EQ(r.lw_net_surface, 0.0);
    EXPECT_EQ(r.max_heating, 0.0);
}
