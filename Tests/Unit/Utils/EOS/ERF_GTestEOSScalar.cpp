#include <gtest/gtest.h>

#include "ERF_GTestEOSCommon.H"

using namespace eos_test;

namespace {

void expect_near_relative (const amrex::Real actual,
                           const amrex::Real expected,
                           const amrex::Real factor = kEOSRelTol)
{
    EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, factor));
}

} // namespace

// Motivation: The ERF constants encode the thermodynamic exponent contract that
// closes every inverse-property test below.
TEST(ERFEOSConstants, KappaGammaContract)
{
    expect_near_relative(R_d / Cp_d, (Gamma - amrex::Real(1.0)) / Gamma);
    expect_near_relative(iGamma, amrex::Real(1.0) / Gamma);
    expect_near_relative(iGamma, amrex::Real(1.0) - (R_d / Cp_d));
}

// Motivation: Temperature and dry potential temperature are the most basic EOS
// inverses. This catches exponent and pressure-reference mistakes immediately.
TEST(ERFEOSScalar, TemperatureAndPotentialTemperatureAreInverses)
{
    for (int index = 0; index < kNumEOSStates; ++index) {
        const EOSState state = get_test_state(index);
        SCOPED_TRACE(state_label(index, state));

        const amrex::Real theta = getThgivenTandP(state.temperature, state.pressure, kRdOcp);
        const amrex::Real recovered_temperature = getTgivenPandTh(state.pressure, theta, kRdOcp);
        const amrex::Real recovered_theta =
            getThgivenTandP(getTgivenPandTh(state.pressure, theta, kRdOcp), state.pressure, kRdOcp);

        expect_near_relative(recovered_temperature, state.temperature);
        expect_near_relative(recovered_theta, theta);
        expect_near_relative(getExnergivenP(state.pressure, kRdOcp), state.temperature / theta);
    }
}

// Motivation: The pressure/rho-theta inverse is the core compressible EOS
// closure. This catches exponent mistakes and missing qv moist factors.
TEST(ERFEOSScalar, PressureAndRhoThetaAreInverses)
{
    for (int index = 0; index < kNumEOSStates; ++index) {
        const EOSState state = get_test_state(index);
        SCOPED_TRACE(state_label(index, state));

        const amrex::Real rho = getRhogivenTandPress(state.temperature, state.pressure, state.qv);
        const amrex::Real theta = getThgivenTandP(state.temperature, state.pressure, kRdOcp);
        const amrex::Real rhotheta = rho * theta;

        expect_near_relative(getPgivenRTh(rhotheta, state.qv), state.pressure);
        expect_near_relative(getRhoThetagivenP(state.pressure, state.qv), rhotheta);
        expect_near_relative(getPgivenRTh(getRhoThetagivenP(state.pressure, state.qv), state.qv),
                             state.pressure);
    }
}

// Motivation: Density can be diagnosed from either T or theta at fixed p. The
// two closures must agree for the same thermodynamic state.
TEST(ERFEOSScalar, DensityClosuresAreConsistent)
{
    for (int index = 0; index < kNumEOSStates; ++index) {
        const EOSState state = get_test_state(index);
        SCOPED_TRACE(state_label(index, state));

        const amrex::Real rho_tp = getRhogivenTandPress(state.temperature, state.pressure, state.qv);
        const amrex::Real theta = getThgivenTandP(state.temperature, state.pressure, kRdOcp);
        const amrex::Real rho_thp = getRhogivenThetaPress(theta, state.pressure, kRdOcp, state.qv);

        expect_near_relative(rho_tp, rho_thp);
        expect_near_relative(state.pressure,
                             rho_tp * R_d * state.temperature * moisture_factor(state.qv));
    }
}

// Motivation: The rho-T and rho-rhoTheta closures should diagnose the same dry
// theta and recover the original temperature and pressure.
TEST(ERFEOSScalar, TemperatureDensityRhoThetaClosureIsConsistent)
{
    for (int index = 0; index < kNumEOSStates; ++index) {
        const EOSState state = get_test_state(index);
        SCOPED_TRACE(state_label(index, state));

        const amrex::Real rho = getRhogivenTandPress(state.temperature, state.pressure, state.qv);
        const amrex::Real theta_from_tp = getThgivenTandP(state.temperature, state.pressure, kRdOcp);
        const amrex::Real theta_from_rhot =
            getThgivenRandT(rho, state.temperature, kRdOcp, state.qv);
        const amrex::Real rhotheta = rho * theta_from_rhot;

        expect_near_relative(theta_from_rhot, theta_from_tp);
        expect_near_relative(getTgivenRandRTh(rho, rhotheta, state.qv), state.temperature);
        expect_near_relative(getPgivenRTh(rhotheta, state.qv), state.pressure);
    }
}

// Motivation: Exner is used in multiple equivalent EOS closures. This catches
// drift between the pressure and rho-theta implementations.
TEST(ERFEOSScalar, ExnerClosuresAgree)
{
    for (int index = 0; index < kNumEOSStates; ++index) {
        const EOSState state = get_test_state(index);
        SCOPED_TRACE(state_label(index, state));

        const amrex::Real rhotheta = getRhoThetagivenP(state.pressure, state.qv);

        expect_near_relative(getExnergivenRTh(rhotheta, kRdOcp, state.qv),
                             getExnergivenP(state.pressure, kRdOcp));
    }
}

// Motivation: The fixed-theta pressure derivative controls acoustic response.
// This checks the analytic EOS derivative against its equivalent closures.
TEST(ERFEOSScalar, PressureDerivativeAtConstantThetaIsCorrect)
{
    for (int index = 0; index < kNumEOSStates; ++index) {
        const EOSState state = get_test_state(index);
        SCOPED_TRACE(state_label(index, state));

        const amrex::Real rho = getRhogivenTandPress(state.temperature, state.pressure, state.qv);
        const amrex::Real theta = getThgivenTandP(state.temperature, state.pressure, kRdOcp);
        const amrex::Real rhotheta = rho * theta;
        const amrex::Real pressure = getPgivenRTh(rhotheta, state.qv);
        const amrex::Real dpdrho = getdPdRgivenConstantTheta(rho, theta, state.qv);

        expect_near_relative(dpdrho, Gamma * pressure / rho);
        expect_near_relative(dpdrho,
                             Gamma * R_d * state.temperature * moisture_factor(state.qv));
        EXPECT_GT(dpdrho, amrex::Real(0.0));
    }
}

// Motivation: RH enters the vapor-pressure helper as a fraction. This prevents
// silent percent-vs-fraction regressions in moisture setups.
TEST(ERFEOSScalar, VaporPressureScalesWithFractionalRH)
{
    const amrex::Real ps = amrex::Real(2350.0);

    expect_near_relative(compute_vapor_pressure(ps, amrex::Real(0.0)), amrex::Real(0.0));
    expect_near_relative(compute_vapor_pressure(ps, amrex::Real(1.0)), ps);
    expect_near_relative(compute_vapor_pressure(ps, amrex::Real(0.75)), amrex::Real(0.75) * ps);
}