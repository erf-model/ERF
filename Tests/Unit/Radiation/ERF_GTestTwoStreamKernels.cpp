#include <cmath>

#include <AMReX_REAL.H>

#include <gtest/gtest.h>

#include <ERF_TwoStreamSW.H>
#include <ERF_TwoStreamLW.H>

// Two-stream kernel unit-test contract
// ------------------------------------
// These tests pin the algebraic behaviour of the per-layer SW and LW kernels
// in ERF_TwoStreamSW.H / ERF_TwoStreamLW.H:
//   1. Beer-Lambert direct beam: S0*mu0 at tau = 0, exp(-tau/mu0) attenuation,
//      zero when the sun is below the horizon.
//   2. SW heating: energy converging into a layer warms it; unphysical
//      inputs (dz, rho, cp <= 0) give exactly zero.
//   3. Diffuse SW source: identically zero for a non-scattering layer,
//      positive and bounded by the attenuated direct beam otherwise.
//   4. Gray-gas LW: an isothermal layer leaves sigma*T^4 unchanged,
//      transmittance is exp(-tau), and a transparent layer is a no-op.
//   5. LW heating sign: net upward flux increasing with height means the
//      layer loses energy, i.e. cooling (negative heating rate). This is the
//      convention that the column sweep relies on.

namespace {

constexpr amrex::Real kSigma = amrex::Real(5.670374419e-8);
constexpr amrex::Real kRelTol = amrex::Real(1.0e-12);

} // namespace

TEST(TwoStreamSWKernels, DirectBeamMatchesBeerLambert)
{
    const amrex::Real S0 = 1361.0;
    const amrex::Real mu0 = 0.5;

    EXPECT_NEAR(compute_sw_direct_flux(0.0, S0, mu0), S0 * mu0, kRelTol * S0);

    const amrex::Real tau = 0.2;
    const amrex::Real expected = S0 * mu0 * std::exp(-tau / mu0);
    EXPECT_NEAR(compute_sw_direct_flux(tau, S0, mu0), expected, kRelTol * S0);

    // Night: no direct beam regardless of optical depth.
    EXPECT_EQ(compute_sw_direct_flux(0.0, S0, 0.0), 0.0);
    EXPECT_EQ(compute_sw_direct_flux(tau, S0, -0.3), 0.0);
}

TEST(TwoStreamSWKernels, DirectBeamIsMonotoneInOpticalDepth)
{
    const amrex::Real S0 = 1361.0;
    const amrex::Real mu0 = 0.7;
    amrex::Real previous = compute_sw_direct_flux(0.0, S0, mu0);
    for (int n = 1; n <= 20; ++n) {
        const amrex::Real current = compute_sw_direct_flux(0.05 * n, S0, mu0);
        EXPECT_LT(current, previous);
        EXPECT_GT(current, 0.0);
        previous = current;
    }
}

TEST(TwoStreamSWKernels, HeatingIsFluxConvergenceOverRhoCp)
{
    const amrex::Real flux_top = 500.0;
    const amrex::Real flux_bot = 480.0;
    const amrex::Real dz = 16.0;
    const amrex::Real rho = 1.1;
    const amrex::Real cp = 1005.0;

    const amrex::Real expected = (flux_top - flux_bot) / dz / (rho * cp);
    EXPECT_NEAR(compute_sw_heating_rate(flux_top, flux_bot, dz, rho, cp), expected, kRelTol);
    EXPECT_GT(compute_sw_heating_rate(flux_top, flux_bot, dz, rho, cp), 0.0);

    // Unphysical inputs give exactly zero rather than Inf/NaN.
    EXPECT_EQ(compute_sw_heating_rate(flux_top, flux_bot, 0.0, rho, cp), 0.0);
    EXPECT_EQ(compute_sw_heating_rate(flux_top, flux_bot, dz, 0.0, cp), 0.0);
    EXPECT_EQ(compute_sw_heating_rate(flux_top, flux_bot, dz, rho, -1.0), 0.0);
}

TEST(TwoStreamSWKernels, DiffuseSourceVanishesWithoutScattering)
{
    const amrex::Real tau = 0.5;
    const amrex::Real F_dir_top = 600.0;
    const amrex::Real mu0 = 0.5;

    // omega == 0 must reduce EXACTLY to the direct-beam-only model.
    EXPECT_EQ(compute_sw_diffuse_flux(tau, F_dir_top, mu0, 0.0, 0.85), 0.0);
    // No incident beam or no layer: nothing to scatter.
    EXPECT_EQ(compute_sw_diffuse_flux(tau, 0.0, mu0, 0.9, 0.85), 0.0);
    EXPECT_EQ(compute_sw_diffuse_flux(0.0, F_dir_top, mu0, 0.9, 0.85), 0.0);
    // Night.
    EXPECT_EQ(compute_sw_diffuse_flux(tau, F_dir_top, 0.0, 0.9, 0.85), 0.0);
}

TEST(TwoStreamSWKernels, DiffuseSourceIsPositiveAndBounded)
{
    const amrex::Real tau = 0.5;
    const amrex::Real F_dir_top = 600.0;
    const amrex::Real mu0 = 0.5;
    const amrex::Real omega = 0.9999;
    const amrex::Real g = 0.85;

    const amrex::Real diffuse = compute_sw_diffuse_flux(tau, F_dir_top, mu0, omega, g);
    EXPECT_GT(diffuse, 0.0);
    // The diffuse source cannot exceed the direct beam entering the layer.
    EXPECT_LT(diffuse, F_dir_top);
    EXPECT_TRUE(std::isfinite(diffuse));

    // More scattering (larger omega) produces more diffuse flux.
    EXPECT_LT(compute_sw_diffuse_flux(tau, F_dir_top, mu0, 0.5, g), diffuse);
}

TEST(TwoStreamLWKernels, ThermalIntensityIsStefanBoltzmann)
{
    const amrex::Real T = 288.15;
    EXPECT_NEAR(compute_thermal_intensity(T, kSigma), kSigma * T * T * T * T, 1.0e-9);
    EXPECT_EQ(compute_thermal_intensity(0.0, kSigma), 0.0);
    EXPECT_EQ(compute_thermal_intensity(-5.0, kSigma), 0.0);
}

TEST(TwoStreamLWKernels, TransmittanceIsExpMinusTau)
{
    EXPECT_EQ(compute_lw_transmit(0.0), 1.0);
    EXPECT_NEAR(compute_lw_transmit(1.0), std::exp(-1.0), kRelTol);
    EXPECT_LT(compute_lw_transmit(50.0), 1.0e-20);
    // Negative (unphysical) optical depth is treated as transparent.
    EXPECT_EQ(compute_lw_transmit(-0.1), 1.0);
}

TEST(TwoStreamLWKernels, IsothermalLayerIsAFixedPoint)
{
    // A layer at temperature T that receives sigma*T^4 must re-emit exactly
    // sigma*T^4 in both directions, for any optical depth.
    const amrex::Real T = 280.0;
    const amrex::Real B = kSigma * T * T * T * T;
    for (amrex::Real tau : {0.0, 0.1, 1.0, 10.0}) {
        EXPECT_NEAR(compute_lw_flux_up(B, T, kSigma, tau), B, 1.0e-9);
        EXPECT_NEAR(compute_lw_flux_down(B, T, kSigma, tau), B, 1.0e-9);
    }
}

TEST(TwoStreamLWKernels, TransparentLayerIsANoOp)
{
    const amrex::Real T = 250.0;
    EXPECT_EQ(compute_lw_flux_up(123.0, T, kSigma, 0.0), 123.0);
    EXPECT_EQ(compute_lw_flux_down(45.0, T, kSigma, 0.0), 45.0);
}

TEST(TwoStreamLWKernels, OpaqueLayerEmitsAtItsOwnTemperature)
{
    const amrex::Real T = 250.0;
    const amrex::Real B = kSigma * T * T * T * T;
    // With tau -> infinity the incoming flux is fully absorbed and replaced
    // by the layer's own emission.
    EXPECT_NEAR(compute_lw_flux_up(1000.0, T, kSigma, 60.0), B, 1.0e-9);
    EXPECT_NEAR(compute_lw_flux_down(0.0, T, kSigma, 60.0), B, 1.0e-9);
}

TEST(TwoStreamLWKernels, NetUpwardFluxIncreasingWithHeightCools)
{
    // Layer emits more through its top than it receives through its bottom.
    const amrex::Real F_net_bot = 50.0;   // up - down at the bottom interface
    const amrex::Real F_net_top = 80.0;   // up - down at the top interface
    const amrex::Real dz = 16.0;
    const amrex::Real rho = 1.0;
    const amrex::Real cp = 1005.0;

    const amrex::Real Q = compute_lw_heating_rate(F_net_top, F_net_bot, dz, rho, cp);
    EXPECT_LT(Q, 0.0);
    EXPECT_NEAR(Q, -(F_net_top - F_net_bot) / dz / (rho * cp), kRelTol);

    // Reverse the divergence: energy converges into the layer -> warming.
    EXPECT_GT(compute_lw_heating_rate(F_net_bot, F_net_top, dz, rho, cp), 0.0);
    // No divergence -> no heating.
    EXPECT_EQ(compute_lw_heating_rate(F_net_top, F_net_top, dz, rho, cp), 0.0);
    // Unphysical inputs give exactly zero.
    EXPECT_EQ(compute_lw_heating_rate(F_net_top, F_net_bot, 0.0, rho, cp), 0.0);
}

TEST(TwoStreamLWKernels, SWAndLWHeatingShareASignConvention)
{
    // Same physical situation expressed in both conventions: 20 W/m^2 are
    // absorbed by the layer. For SW that is downwelling flux decreasing
    // downward; for LW that is net upward flux decreasing with height.
    const amrex::Real dz = 10.0, rho = 1.0, cp = 1005.0;
    const amrex::Real Q_sw = compute_sw_heating_rate(120.0, 100.0, dz, rho, cp);
    const amrex::Real Q_lw = compute_lw_heating_rate(30.0, 50.0, dz, rho, cp);
    EXPECT_NEAR(Q_sw, Q_lw, kRelTol);
    EXPECT_GT(Q_sw, 0.0);
}
