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
//   3. Two-stream layer solution: a non-scattering layer only absorbs, a
//      conservative layer (omega = 1) neither absorbs nor creates energy,
//      and reflectances/transmittances stay within physical bounds.
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

TEST(TwoStreamSWKernels, LayerWithoutScatteringOnlyAbsorbs)
{
    const amrex::Real tau = 0.5;
    const amrex::Real mu0 = 0.5;
    const TwoStreamLayerSW L = compute_sw_layer_two_stream(tau, 0.0, 0.85, mu0);
    EXPECT_EQ(L.R_dif, 0.0);
    EXPECT_EQ(L.R_dir, 0.0);
    EXPECT_EQ(L.T_dir, 0.0);
    // Diffuse light crosses a pure absorber with the diffusivity factor 2.
    EXPECT_NEAR(L.T_dif, std::exp(-2.0 * tau), 1.0e-9);
    EXPECT_NEAR(L.T_noscat, std::exp(-tau / mu0), kRelTol);

    // An empty layer is transparent to everything.
    const TwoStreamLayerSW E = compute_sw_layer_two_stream(0.0, 0.9, 0.85, mu0);
    EXPECT_EQ(E.R_dif, 0.0);
    EXPECT_EQ(E.T_dif, 1.0);
    EXPECT_EQ(E.R_dir, 0.0);
    EXPECT_EQ(E.T_dir, 0.0);
    EXPECT_EQ(E.T_noscat, 1.0);
}

TEST(TwoStreamSWKernels, ConservativeLayerNeitherAbsorbsNorCreatesEnergy)
{
    for (amrex::Real g : {0.0, 0.5, 0.85}) {
        for (amrex::Real tau : {0.01, 0.5, 5.0}) {
            for (amrex::Real mu0 : {0.2, 0.5, 1.0}) {
                SCOPED_TRACE("g=" + std::to_string(g) + " tau=" + std::to_string(tau) +
                             " mu0=" + std::to_string(mu0));
                const TwoStreamLayerSW L = compute_sw_layer_two_stream(tau, 1.0, g, mu0);
                // omega = 1: everything incident is reflected or transmitted.
                EXPECT_NEAR(L.R_dif + L.T_dif, 1.0, 1.0e-4);
                EXPECT_NEAR(L.R_dir + L.T_dir + L.T_noscat, 1.0, 1.0e-4);
                EXPECT_GE(L.R_dif, 0.0);
                EXPECT_GE(L.R_dir, 0.0);
                EXPECT_GE(L.T_dir, 0.0);
            }
        }
    }
}

TEST(TwoStreamSWKernels, PartlyAbsorbingLayerStaysWithinPhysicalBounds)
{
    for (amrex::Real omega : {0.3, 0.9, 0.9999}) {
        for (amrex::Real tau : {0.05, 1.0, 20.0}) {
            SCOPED_TRACE("omega=" + std::to_string(omega) + " tau=" + std::to_string(tau));
            const TwoStreamLayerSW L = compute_sw_layer_two_stream(tau, omega, 0.85, 0.5);
            EXPECT_GE(L.R_dif, 0.0);
            EXPECT_GE(L.T_dif, 0.0);
            EXPECT_LE(L.R_dif + L.T_dif, 1.0 + 1.0e-12);
            EXPECT_GE(L.R_dir, 0.0);
            EXPECT_GE(L.T_dir, 0.0);
            EXPECT_LE(L.R_dir + L.T_dir + L.T_noscat, 1.0 + 1.0e-12);
            // Some absorption must remain when omega < 1.
            EXPECT_LT(L.R_dif + L.T_dif, 1.0);
            EXPECT_LT(L.R_dir + L.T_dir + L.T_noscat, 1.0);
        }
    }
    // A thick, strongly scattering layer reflects most of the direct beam.
    const TwoStreamLayerSW thick = compute_sw_layer_two_stream(20.0, 0.9999, 0.85, 0.5);
    EXPECT_GT(thick.R_dir, 0.5);
    EXPECT_LT(thick.T_noscat, 1.0e-15);
}

TEST(TwoStreamSWKernels, MoreScatteringReflectsMore)
{
    const amrex::Real tau = 0.5, g = 0.85, mu0 = 0.5;
    amrex::Real previous = -1.0;
    for (amrex::Real omega : {0.2, 0.5, 0.8, 0.99}) {
        const TwoStreamLayerSW L = compute_sw_layer_two_stream(tau, omega, g, mu0);
        EXPECT_GT(L.R_dir, previous);
        previous = L.R_dir;
    }
    // Night: no direct beam and no direct-beam scattering.
    const TwoStreamLayerSW night = compute_sw_layer_two_stream(tau, 0.9, g, 0.0);
    EXPECT_EQ(night.T_noscat, 0.0);
    EXPECT_EQ(night.R_dir, 0.0);
    EXPECT_EQ(night.T_dir, 0.0);
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
