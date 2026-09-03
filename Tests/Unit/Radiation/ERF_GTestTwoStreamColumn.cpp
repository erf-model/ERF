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
//   2b. The stored heating rates are potential-temperature tendencies,
//      (dT/dt) / pi, matching the RhoTheta source term and RRTMGP.
//   3. The absorbed SW surface flux equals the Beer-Lambert direct beam
//      through the whole column times (1 - albedo), the reflected beam
//      leaves the top attenuated by the diffuse transmittance, and the
//      column SW energy budget closes (absorbed by air + surface = incident
//      - reflected at the top).
//   4. In an isothermal column over a black surface at the same temperature
//      the LW heating is non-positive everywhere, strongest at the top layer
//      (cooling to space), and the surface net LW equals the analytic
//      sigma*T^4 * exp(-tau_column).
//   5. A gray surface reflects (1 - eps) of the downwelling LW, so the
//      surface net LW equals eps * sigma*T^4 * exp(-tau_column) rather than
//      the large negative value obtained without the reflected term.
//   6. Night zeroes the SW heating; disabled bands write zero.
//   8. tau_model = mass: the SW column optical depth is set by the mass
//      path (resolution independent), layer optics are extinction-weighted
//      mixtures of the constituents, Rayleigh-only columns absorb nothing,
//      cloud water brightens the column, and the LW band uses the mass path.
//   7. The mass-based gray LW option makes the column optical depth
//      independent of the vertical resolution and reproduces the Stephens
//      (1978) cloud emissivity; separate direct/diffuse albedo and the
//      Earth-Sun distance factor behave as documented.
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
    amrex::Real sw_up_toa = 0.0;
    amrex::Real lw_net_surface = 0.0;
    amrex::Real lw_up_toa = 0.0;
};

// Run the sweep for a column with uniform density `rho` and uniform
// absolute temperature `T_air` (converted to rho*theta through the EOS).
ColumnResult run_uniform_column (const RadChoice& rad_choice_in, amrex::Real rho, amrex::Real T_air,
                                 int nz = kNz, amrex::Real dz = kDz, amrex::Real qv = 0.0,
                                 amrex::Real qc = 0.0)
{
    const TwoStreamParams rad_choice = make_two_stream_params(rad_choice_in);
    const amrex::Box bx(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, nz - 1));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {dz, dz, nz * dz});
    const int is_periodic[3] = {1, 1, 0};
    const amrex::Geometry geom(bx, &real_box, 0, is_periodic);

    // Carry the moisture components so qv can be set (zero by default).
    amrex::FArrayBox state(bx, RhoQ2_comp + 1);
    amrex::FArrayBox qheating(bx, 2);
    const amrex::Real theta = getThgivenRandT(rho, T_air, RdoCp, qv);
    state.setVal<amrex::RunOn::Device>(0.0);
    state.setVal<amrex::RunOn::Device>(rho, bx, Rho_comp, 1);
    state.setVal<amrex::RunOn::Device>(rho * theta, bx, RhoTheta_comp, 1);
    state.setVal<amrex::RunOn::Device>(rho * qv, bx, RhoQ1_comp, 1);
    state.setVal<amrex::RunOn::Device>(rho * qc, bx, RhoQ2_comp, 1);
    qheating.setVal<amrex::RunOn::Device>(0.0);

    amrex::Gpu::DeviceVector<amrex::Real> scalars(5, 0.0);
    amrex::Real* scalar_ptr = scalars.data();
    const auto state_arr = state.const_array();
    const auto qheating_arr = qheating.array();
    const amrex::Array4<const amrex::Real> no_z_phys{};

    const amrex::Box xy_box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
    {
        amrex::Real max_heating = 0.0;
        amrex::Real sw_surface = 0.0;
        amrex::Real sw_up = 0.0;
        amrex::Real lw_net = 0.0;
        amrex::Real lw_up = 0.0;
        vertical_two_stream_sweep(i, j, bx, geom, state_arr, rad_choice, /*cloudy=*/false,
                                  qheating_arr, max_heating, sw_surface, sw_up, lw_net, lw_up,
                                  no_z_phys);
        scalar_ptr[0] = max_heating;
        scalar_ptr[1] = sw_surface;
        scalar_ptr[2] = sw_up;
        scalar_ptr[3] = lw_net;
        scalar_ptr[4] = lw_up;
    });
    amrex::Gpu::streamSynchronize();

    ColumnResult result;
    std::vector<amrex::Real> host_scalars(5);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, scalars.begin(), scalars.end(), host_scalars.begin());
    result.max_heating = host_scalars[0];
    result.sw_surface = host_scalars[1];
    result.sw_up_toa = host_scalars[2];
    result.lw_net_surface = host_scalars[3];
    result.lw_up_toa = host_scalars[4];

    amrex::FArrayBox host_q(bx, 2, amrex::The_Pinned_Arena());
    host_q.copy<amrex::RunOn::Device>(qheating);
    amrex::Gpu::streamSynchronize();
    const auto hq = host_q.const_array();
    for (int k = 0; k < nz; ++k) {
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
    // Outgoing LW at the top: sigma T^4 from the (isothermal) column.
    EXPECT_NEAR(r.lw_up_toa, B, 1.0e-9 * B);
}

TEST(TwoStreamColumn, HeatingRatesArePotentialTemperatureTendencies)
{
    // Non-scattering column: the top layer's temperature tendency follows
    // from the direct beam it absorbs plus the reflected beam absorbed on
    // the way up, and the stored value must be that divided by the Exner
    // function.
    const RadChoice rc = base_choice();
    const amrex::Real rho = 1.0;
    const amrex::Real T_air = 290.0;
    const ColumnResult r = run_uniform_column(rc, rho, T_air);

    const amrex::Real mu0 = std::cos(rc.solar_zenith_deg * M_PI / 180.0);
    const amrex::Real tau = rc.tau_per_layer;
    const amrex::Real tau_col = kNz * tau;
    const amrex::Real F0 = rc.S0 * mu0;
    const amrex::Real F_dir_sfc = F0 * std::exp(-tau_col / mu0);
    // Net downward flux at the top interface and below the top layer.
    const amrex::Real u_top = rc.surface_albedo_sw * F_dir_sfc * std::exp(-2.0 * tau_col);
    const amrex::Real u_below = rc.surface_albedo_sw * F_dir_sfc * std::exp(-2.0 * (tau_col - tau));
    const amrex::Real F_net_top = F0 - u_top;
    const amrex::Real F_net_below = F0 * std::exp(-tau / mu0) - u_below;
    const amrex::Real dTdt = (F_net_top - F_net_below) / (rho * Cp_d * kDz);

    const amrex::Real theta = getThgivenRandT(rho, T_air, RdoCp);
    const amrex::Real exner = getExnergivenRTh(rho * theta, RdoCp);
    EXPECT_LT(exner, 1.0);   // the test state sits below p_0
    EXPECT_NEAR(r.q_sw[kNz - 1], dTdt / exner, 1.0e-9 * dTdt / exner);
    // The stored value must NOT be the raw temperature tendency.
    EXPECT_GT(r.q_sw[kNz - 1], dTdt * (1.0 + 1.0e-6));
}

TEST(TwoStreamColumn, ShortwaveEnergyBudgetClosesWithSurfaceReflection)
{
    const RadChoice rc = base_choice();
    const amrex::Real rho = 1.0;
    const ColumnResult r = run_uniform_column(rc, rho, 290.0);

    // Non-scattering layers: the reflected beam alpha * F_dir(0) travels up
    // as diffuse light with transmittance exp(-2 tau) per layer.
    const amrex::Real mu0 = std::cos(rc.solar_zenith_deg * M_PI / 180.0);
    const amrex::Real tau_col = kNz * rc.tau_per_layer;
    const amrex::Real F_dir_sfc = rc.S0 * mu0 * std::exp(-tau_col / mu0);
    const amrex::Real expected_up = rc.surface_albedo_sw * F_dir_sfc * std::exp(-2.0 * tau_col);
    EXPECT_NEAR(r.sw_up_toa, expected_up, 1.0e-6 * expected_up);

    // Energy budget: absorbed by the air (sum rho cp dz Q) plus absorbed by
    // the surface equals incident minus reflected at the top.
    // q_sw is a potential-temperature tendency; the energy absorbed by a
    // layer is rho cp dz dT/dt = rho cp dz (pi q_sw).
    const amrex::Real theta = getThgivenRandT(rho, 290.0, RdoCp);
    const amrex::Real exner = getExnergivenRTh(rho * theta, RdoCp);
    amrex::Real absorbed_air = 0.0;
    for (int k = 0; k < kNz; ++k) {
        absorbed_air += rho * Cp_d * kDz * exner * r.q_sw[k];
    }
    const amrex::Real incident = rc.S0 * mu0;
    EXPECT_NEAR(absorbed_air + r.sw_surface, incident - r.sw_up_toa, 1.0e-9 * incident);
}

TEST(TwoStreamColumn, ConservativeScatteringDepositsNoEnergyInTheAir)
{
    RadChoice rc = base_choice();
    rc.single_scattering_albedo = 1.0;   // every layer scatters, nothing absorbs
    rc.asymmetry_factor = 0.6;
    const ColumnResult r = run_uniform_column(rc, 1.0, 290.0);

    const amrex::Real mu0 = std::cos(rc.solar_zenith_deg * M_PI / 180.0);
    const amrex::Real incident = rc.S0 * mu0;
    for (int k = 0; k < kNz; ++k) {
        EXPECT_NEAR(r.q_sw[k], 0.0, 1.0e-7 * incident / (Cp_d * kDz)) << "k = " << k;
    }
    // Whatever is not reflected at the top is absorbed by the surface.
    EXPECT_NEAR(r.sw_surface, incident - r.sw_up_toa, 1.0e-6 * incident);
    // A scattering column reflects more than the bare surface would.
    EXPECT_GT(r.sw_up_toa, rc.surface_albedo_sw * incident * std::exp(-2.0 * kNz * rc.tau_per_layer));
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
    EXPECT_EQ(r.sw_up_toa, 0.0);
    EXPECT_EQ(r.lw_net_surface, 0.0);
    EXPECT_EQ(r.lw_up_toa, 0.0);
    EXPECT_EQ(r.max_heating, 0.0);
}

TEST(TwoStreamColumn, CloudBandIsLocatedByLayerCenterHeight)
{
    RadChoice rc = base_choice();
    rc.tau_profile_type = TauProfileType::CloudLayer;
    rc.cloud_base_height_m = 300.0;
    rc.cloud_top_height_m = 700.0;
    rc.cloud_tau_per_layer = 0.5;
    rc.tau_per_layer = 0.05;
    auto P = [&]() { return make_two_stream_params(rc); };

    // Inside the band: base + cloud enhancement for the cloudy column only.
    EXPECT_TRUE(is_cloud_level(500.0, P()));
    EXPECT_NEAR(tau_layer_value(500.0, rc.tau_per_layer, P(), true), 0.55, 1.0e-12);
    EXPECT_NEAR(tau_layer_value(500.0, rc.tau_per_layer, P(), false), 0.05, 1.0e-12);
    // Band edges are inclusive; outside the band the base value is returned.
    EXPECT_TRUE(is_cloud_level(300.0, P()));
    EXPECT_TRUE(is_cloud_level(700.0, P()));
    EXPECT_FALSE(is_cloud_level(299.9, P()));
    EXPECT_FALSE(is_cloud_level(700.1, P()));
    EXPECT_NEAR(tau_layer_value(900.0, rc.tau_per_layer, P(), true), 0.05, 1.0e-12);

    // Constant profile ignores the band entirely.
    rc.tau_profile_type = TauProfileType::Constant;
    EXPECT_NEAR(tau_layer_value(500.0, rc.tau_per_layer, P(), true), 0.05, 1.0e-12);

    // Cloud scattering properties follow the same band test.
    rc.tau_profile_type = TauProfileType::CloudLayer;
    rc.cloud_single_scattering_albedo = 0.9;
    rc.cloud_asymmetry_factor = 0.85;
    amrex::Real omega = -1.0, g = -1.0;
    select_scattering_props(500.0, P(), true, omega, g);
    EXPECT_EQ(omega, 0.9);
    EXPECT_EQ(g, 0.85);
    select_scattering_props(900.0, P(), true, omega, g);
    EXPECT_EQ(omega, rc.single_scattering_albedo);
    EXPECT_EQ(g, rc.asymmetry_factor);
}

TEST(TwoStreamColumn, DynamicOpticalDepthIsLinearInMoisture)
{
    // Zero coefficients reproduce the static value exactly.
    EXPECT_EQ(diagnose_tau_dynamic(0.05, 0.01, 0.001, 0.0, 0.0), 0.05);
    // Linear in qv and qc.
    EXPECT_NEAR(diagnose_tau_dynamic(0.05, 0.01, 0.001, 10.0, 200.0), 0.05 + 0.1 + 0.2, 1.0e-12);
    // Negative or non-finite mixing ratios contribute nothing.
    EXPECT_EQ(diagnose_tau_dynamic(0.05, -0.01, -1.0, 10.0, 200.0), 0.05);
    EXPECT_EQ(diagnose_tau_dynamic(0.05, std::nan(""), 0.0, 10.0, 200.0), 0.05);
    // Clamped to the physical range.
    EXPECT_EQ(diagnose_tau_dynamic(0.05, 1.0, 1.0, 1.0e3, 1.0e3), 100.0);
}

TEST(TwoStreamColumn, MoistureHelpersGuardMissingComponents)
{
    // A dry state with only (Rho, RhoTheta) must report zero moisture rather
    // than reading past the end of the component range.
    const amrex::Box bx(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::FArrayBox dry(bx, 2, amrex::The_Pinned_Arena());
    dry.setVal<amrex::RunOn::Host>(1.0, bx, Rho_comp, 1);
    dry.setVal<amrex::RunOn::Host>(300.0, bx, RhoTheta_comp, 1);
    EXPECT_EQ(get_qv_from_state(0, 0, 0, dry.const_array()), 0.0);
    EXPECT_EQ(get_qc_from_state(0, 0, 0, dry.const_array()), 0.0);

    // With moisture components present the mixing ratios are RhoQ / Rho.
    amrex::FArrayBox moist(bx, RhoQ2_comp + 1, amrex::The_Pinned_Arena());
    moist.setVal<amrex::RunOn::Host>(0.0);
    moist.setVal<amrex::RunOn::Host>(2.0, bx, Rho_comp, 1);
    moist.setVal<amrex::RunOn::Host>(600.0, bx, RhoTheta_comp, 1);
    moist.setVal<amrex::RunOn::Host>(0.02, bx, RhoQ1_comp, 1);
    moist.setVal<amrex::RunOn::Host>(0.004, bx, RhoQ2_comp, 1);
    EXPECT_NEAR(get_qv_from_state(0, 0, 0, moist.const_array()), 0.01, 1.0e-15);
    EXPECT_NEAR(get_qc_from_state(0, 0, 0, moist.const_array()), 0.002, 1.0e-15);

    // Negative stored values are treated as zero.
    moist.setVal<amrex::RunOn::Host>(-0.02, bx, RhoQ1_comp, 1);
    EXPECT_EQ(get_qv_from_state(0, 0, 0, moist.const_array()), 0.0);
}

TEST(TwoStreamColumn, MassBasedLongwaveIsIndependentOfVerticalResolution)
{
    RadChoice rc = base_choice();
    rc.sw_enabled = false;
    rc.lw_mass_absorption_enable = true;
    rc.lw_kabs_dry = 1.0e-4;
    rc.lw_kabs_vapor = 0.1;
    rc.surface_temp_k = 300.0;
    const amrex::Real rho = 1.0;
    const amrex::Real qv = 0.01;

    // Same 800 m column, air at 280 K over a 300 K surface, on 8 and 32 layers.
    const ColumnResult coarse = run_uniform_column(rc, rho, 280.0, 8, 100.0, qv);
    const ColumnResult fine   = run_uniform_column(rc, rho, 280.0, 32, 25.0, qv);

    // Column optical depth rho H (k_dry + k_v qv) = 0.88: partly transparent,
    // so the outgoing LW carries a surface contribution and the surface net
    // LW is far from zero; both must not depend on the layering.
    const amrex::Real B_s = kSigma * 300.0 * 300.0 * 300.0 * 300.0;
    EXPECT_GT(coarse.lw_net_surface, 0.05 * B_s);
    EXPECT_NEAR(coarse.lw_up_toa, fine.lw_up_toa, 1.0e-9 * B_s);
    EXPECT_NEAR(coarse.lw_net_surface, fine.lw_net_surface, 1.0e-9 * B_s);

    // With the fixed per-layer value instead, refining the grid quadruples
    // the column optical depth and changes the fluxes.
    rc.lw_mass_absorption_enable = false;
    rc.tau_lw_per_layer = 0.11;
    const ColumnResult coarse_fixed = run_uniform_column(rc, rho, 280.0, 8, 100.0, qv);
    const ColumnResult fine_fixed   = run_uniform_column(rc, rho, 280.0, 32, 25.0, qv);
    EXPECT_GT(std::abs(coarse_fixed.lw_up_toa - fine_fixed.lw_up_toa), 1.0e-2 * B_s);
}

TEST(TwoStreamColumn, MassBasedLongwaveOpticalDepthFollowsTheMassPath)
{
    RadChoice rc = base_choice();
    rc.lw_mass_absorption_enable = true;
    rc.lw_kabs_dry = 2.0e-4;
    rc.lw_kabs_vapor = 0.05;
    rc.lw_kabs_cloud = 158.0;
    const TwoStreamParams p = make_two_stream_params(rc);

    const amrex::Box bx(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::FArrayBox state(bx, RhoQ2_comp + 1, amrex::The_Pinned_Arena());
    const amrex::Real rho = 1.1, qv = 0.008, qc = 5.0e-4, dz = 50.0;
    state.setVal<amrex::RunOn::Host>(0.0);
    state.setVal<amrex::RunOn::Host>(rho, bx, Rho_comp, 1);
    state.setVal<amrex::RunOn::Host>(rho * 300.0, bx, RhoTheta_comp, 1);
    state.setVal<amrex::RunOn::Host>(rho * qv, bx, RhoQ1_comp, 1);
    state.setVal<amrex::RunOn::Host>(rho * qc, bx, RhoQ2_comp, 1);

    const amrex::Real tau = diagnose_layer_tau(0, 0, 0, dz, 25.0, state.const_array(),
                                               /*tau_base=*/1.0, /*is_sw=*/false, /*cloudy=*/false, p);
    const amrex::Real expected = rho * dz * (2.0e-4 + 0.05 * qv + 158.0 * qc);
    EXPECT_NEAR(tau, expected, 1.0e-12);

    // Cloud term: Stephens (1978) emissivity 1 - exp(-0.158 LWP[g/m^2]).
    const amrex::Real lwp_g = rho * qc * dz * 1.0e3;
    const amrex::Real tau_cloud = tau - rho * dz * (2.0e-4 + 0.05 * qv);
    EXPECT_NEAR(1.0 - std::exp(-tau_cloud), 1.0 - std::exp(-0.158 * lwp_g), 1.0e-12);

    // The SW band is unaffected by the LW option.
    EXPECT_EQ(diagnose_layer_tau(0, 0, 0, dz, 25.0, state.const_array(), 0.05, true, false, p), 0.05);
    // Disabled: the fixed per-layer value is returned.
    rc.lw_mass_absorption_enable = false;
    EXPECT_EQ(diagnose_layer_tau(0, 0, 0, dz, 25.0, state.const_array(), 1.0, false, false,
                                 make_two_stream_params(rc)), 1.0);
}

TEST(TwoStreamColumn, DiffuseAlbedoAppliesToTheDiffuseFluxOnly)
{
    RadChoice rc = base_choice();
    rc.single_scattering_albedo = 1.0;   // conservative scattering: diffuse flux reaches the surface
    rc.asymmetry_factor = 0.6;
    rc.surface_albedo_sw = 0.2;
    const amrex::Real mu0 = std::cos(rc.solar_zenith_deg * M_PI / 180.0);
    const amrex::Real incident = rc.S0 * mu0;

    rc.surface_albedo_sw_diffuse = -1.0;   // same as direct
    const ColumnResult same = run_uniform_column(rc, 1.0, 290.0);
    rc.surface_albedo_sw_diffuse = 0.8;    // reflect most diffuse light
    const ColumnResult bright = run_uniform_column(rc, 1.0, 290.0);

    // More diffuse reflection: less absorbed at the surface, more leaving the top,
    // and the column budget still closes (no absorption in a conservative column).
    EXPECT_LT(bright.sw_surface, same.sw_surface);
    EXPECT_GT(bright.sw_up_toa, same.sw_up_toa);
    EXPECT_NEAR(bright.sw_surface, incident - bright.sw_up_toa, 1.0e-6 * incident);

    // Without scattering there is no diffuse flux at the surface, so the
    // diffuse albedo cannot matter.
    rc.single_scattering_albedo = 0.0;
    rc.surface_albedo_sw_diffuse = 0.8;
    const ColumnResult dark = run_uniform_column(rc, 1.0, 290.0);
    rc.surface_albedo_sw_diffuse = -1.0;
    const ColumnResult ref = run_uniform_column(rc, 1.0, 290.0);
    EXPECT_NEAR(dark.sw_surface, ref.sw_surface, 1.0e-12 * incident);
}

TEST(TwoStreamColumn, EarthSunDistanceFactorScalesTheSolarConstant)
{
    // Spencer (1971): perihelion in early January, aphelion in early July.
    EXPECT_NEAR(compute_earth_sun_distance_factor(3.0), 1.034, 2.0e-3);
    EXPECT_NEAR(compute_earth_sun_distance_factor(185.0), 0.967, 2.0e-3);
    amrex::Real mean = 0.0;
    for (int d = 1; d <= 365; ++d) mean += compute_earth_sun_distance_factor(d);
    EXPECT_NEAR(mean / 365.0, 1.0, 1.0e-3);
    EXPECT_EQ(compute_earth_sun_distance_factor(std::nan("")), 1.0);

    RadChoice rc = base_choice();
    rc.day_of_year = 3.0;
    const ColumnResult off = run_uniform_column(rc, 1.0, 290.0);
    rc.earth_sun_distance_enable = true;
    const ColumnResult on = run_uniform_column(rc, 1.0, 290.0);
    const amrex::Real f = compute_earth_sun_distance_factor(3.0);
    EXPECT_NEAR(on.sw_surface, f * off.sw_surface, 1.0e-9 * off.sw_surface);
    EXPECT_NEAR(on.sw_up_toa, f * off.sw_up_toa, 1.0e-9 * off.sw_up_toa);
}

TEST(TwoStreamColumn, MassModelShortwaveIsIndependentOfVerticalResolution)
{
    RadChoice rc = base_choice();
    rc.tau_model = TauModel::Mass;
    rc.sw_kabs_dry = 4.0e-6;
    rc.sw_kscat_dry = 3.0e-6;
    rc.sw_kabs_vapor = 4.0e-3;
    const amrex::Real rho = 1.0, qv = 0.01;
    const amrex::Real mu0 = std::cos(rc.solar_zenith_deg * M_PI / 180.0);
    const amrex::Real incident = rc.S0 * mu0;

    const ColumnResult coarse = run_uniform_column(rc, rho, 290.0, 8, 100.0, qv);
    const ColumnResult fine   = run_uniform_column(rc, rho, 290.0, 32, 25.0, qv);
    EXPECT_NEAR(coarse.sw_surface, fine.sw_surface, 1.0e-9 * incident);
    EXPECT_NEAR(coarse.sw_up_toa, fine.sw_up_toa, 1.0e-9 * incident);
    EXPECT_NEAR(coarse.lw_up_toa, fine.lw_up_toa, 1.0e-9 * incident);
    // The column absorbs something (vapor) and scatters something (Rayleigh).
    EXPECT_LT(coarse.sw_surface, (1.0 - rc.surface_albedo_sw) * incident);
    EXPECT_GT(coarse.sw_up_toa, 0.0);
    EXPECT_LT(coarse.sw_up_toa, incident);

    // The per-layer model, by contrast, changes with the layering.
    rc.tau_model = TauModel::PerLayer;
    rc.tau_per_layer = 0.02;
    const ColumnResult coarse_fixed = run_uniform_column(rc, rho, 290.0, 8, 100.0, qv);
    const ColumnResult fine_fixed   = run_uniform_column(rc, rho, 290.0, 32, 25.0, qv);
    EXPECT_GT(std::abs(coarse_fixed.sw_surface - fine_fixed.sw_surface), 1.0e-2 * incident);
}

TEST(TwoStreamColumn, MassModelLayerOpticsAreExtinctionWeighted)
{
    RadChoice rc = base_choice();
    rc.tau_model = TauModel::Mass;
    rc.sw_kabs_dry = 4.0e-6;
    rc.sw_kscat_dry = 3.0e-6;
    rc.sw_kabs_vapor = 4.0e-3;
    rc.sw_kext_cloud = 150.0;
    rc.sw_cloud_omega = 0.9999;
    rc.sw_cloud_g = 0.85;
    const TwoStreamParams p = make_two_stream_params(rc);

    const amrex::Box bx(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::FArrayBox state(bx, RhoQ2_comp + 1, amrex::The_Pinned_Arena());
    const amrex::Real rho = 1.1, qv = 0.008, qc = 2.0e-4, dz = 50.0;
    state.setVal<amrex::RunOn::Host>(0.0);
    state.setVal<amrex::RunOn::Host>(rho, bx, Rho_comp, 1);
    state.setVal<amrex::RunOn::Host>(rho * 300.0, bx, RhoTheta_comp, 1);
    state.setVal<amrex::RunOn::Host>(rho * qv, bx, RhoQ1_comp, 1);
    state.setVal<amrex::RunOn::Host>(rho * qc, bx, RhoQ2_comp, 1);

    amrex::Real tau = -1.0, omega = -1.0, g = -1.0;
    diagnose_layer_optics(0, 0, 0, dz, 25.0, state.const_array(), 0.05, true, false, p, tau, omega, g);

    const amrex::Real path = rho * dz;
    const amrex::Real t_abs = path * 4.0e-6, t_ray = path * 3.0e-6, t_vap = path * 4.0e-3 * qv,
                      t_cld = path * 150.0 * qc;
    const amrex::Real ext = t_abs + t_ray + t_vap + t_cld;
    const amrex::Real sca = t_ray + 0.9999 * t_cld;
    EXPECT_NEAR(tau, ext, 1.0e-12);
    EXPECT_NEAR(omega, sca / ext, 1.0e-12);
    EXPECT_NEAR(g, 0.85 * 0.9999 * t_cld / sca, 1.0e-12);

    // The longwave band takes the LW mass path in this model.
    diagnose_layer_optics(0, 0, 0, dz, 25.0, state.const_array(), 1.0, false, false, p, tau, omega, g);
    EXPECT_NEAR(tau, path * (rc.lw_kabs_dry + rc.lw_kabs_vapor * qv + rc.lw_kabs_cloud * qc), 1.0e-12);
    EXPECT_EQ(omega, 0.0);

    // Per-layer model: unchanged behaviour (fixed tau, input scattering props).
    rc.tau_model = TauModel::PerLayer;
    rc.single_scattering_albedo = 0.3;
    rc.asymmetry_factor = 0.5;
    diagnose_layer_optics(0, 0, 0, dz, 25.0, state.const_array(), 0.05, true, false,
                          make_two_stream_params(rc), tau, omega, g);
    EXPECT_EQ(tau, 0.05);
    EXPECT_EQ(omega, 0.3);
    EXPECT_EQ(g, 0.5);
}

TEST(TwoStreamColumn, MassModelRayleighOnlyColumnAbsorbsNothing)
{
    RadChoice rc = base_choice();
    rc.tau_model = TauModel::Mass;
    rc.sw_kabs_dry = 0.0;
    rc.sw_kscat_dry = 2.0e-5;   // exaggerated Rayleigh so the effect is visible
    rc.sw_kabs_vapor = 0.0;
    const amrex::Real mu0 = std::cos(rc.solar_zenith_deg * M_PI / 180.0);
    const amrex::Real incident = rc.S0 * mu0;
    const ColumnResult r = run_uniform_column(rc, 1.0, 290.0);
    for (int k = 0; k < kNz; ++k) {
        EXPECT_NEAR(r.q_sw[k], 0.0, 1.0e-7 * incident / (Cp_d * kDz)) << "k = " << k;
    }
    EXPECT_NEAR(r.sw_surface, incident - r.sw_up_toa, 1.0e-6 * incident);
    // Rayleigh scattering sends more to space than the surface alone would.
    EXPECT_GT(r.sw_up_toa, rc.surface_albedo_sw * incident);
}

TEST(TwoStreamColumn, MassModelCloudWaterBrightensTheColumn)
{
    RadChoice rc = base_choice();
    rc.tau_model = TauModel::Mass;
    const ColumnResult clear = run_uniform_column(rc, 1.0, 290.0, kNz, kDz, 0.005, 0.0);
    const ColumnResult cloudy = run_uniform_column(rc, 1.0, 290.0, kNz, kDz, 0.005, 3.0e-4);
    // LWP = 1 * 3e-4 * 800 = 0.24 kg/m^2 -> tau_cloud = 36: a thick cloud.
    EXPECT_GT(cloudy.sw_up_toa, 2.0 * clear.sw_up_toa);
    EXPECT_LT(cloudy.sw_surface, 0.5 * clear.sw_surface);
    // The cloud also makes the column opaque in the longwave.
    EXPECT_LT(std::abs(cloudy.lw_net_surface), std::abs(clear.lw_net_surface));
}
