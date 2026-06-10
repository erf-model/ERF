#include <array>

#include <gtest/gtest.h>

#include "ERF_GTestMicrophysicsCommon.H"

using namespace microphysics_test;

namespace {

void expect_near_relative (const amrex::Real actual,
                           const amrex::Real expected,
                           const amrex::Real factor = kValueRelTol)
{
    EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, factor));
}

} // namespace

// Motivation: Morrison and SAM only call this gamma wrapper with positive
// arguments. The test locks in that positive-only contract and key identities.
TEST(MicrophysicsGamma, PositiveIdentities)
{
    expect_near_relative(erf_gammafff(amrex::Real(1.0)), amrex::Real(1.0));
    expect_near_relative(erf_gammafff(amrex::Real(0.5)), kSqrtPi);
    expect_near_relative(erf_gammafff(amrex::Real(2.5)), amrex::Real(1.5) * erf_gammafff(amrex::Real(1.5)));
}

// Motivation: The `_cc` helper is ERF's Magnus-style exponential fallback for
// water saturation pressure. Its derivative should match both the analytic
// derivative of that closed-form expression and a finite-difference check.
TEST(MicrophysicsWaterCC, DerivativeMatchesAnalyticAndFiniteDifference)
{
    const std::array<amrex::Real, 4> temperatures = {
        amrex::Real(190.0), amrex::Real(240.0), amrex::Real(300.0), amrex::Real(350.0)};

    for (const amrex::Real temperature : temperatures) {
        const amrex::Real esat = erf_esatw_cc(temperature);
        const amrex::Real derivative = erf_dtesatw_cc(temperature);
        const amrex::Real expected = cc_derivative_formula(temperature);
        const amrex::Real finite_difference =
            central_difference([] (const amrex::Real value) { return erf_esatw_cc(value); }, temperature);

        EXPECT_GT(esat, amrex::Real(0.0));
        EXPECT_GT(derivative, amrex::Real(0.0));
        expect_near_relative(derivative, expected, kDerivativeRelTol);
        expect_near_relative(derivative, finite_difference, kDerivativeRelTol);
    }
}

// Motivation: The empirical branch is selectable at runtime. It must report the
// same mbar units as the default branch or the HSE initialization uses the
// wrong vapor pressure by a factor of 100.
TEST(MicrophysicsWaterSaturation, EmpiricalBranchUsesMbarUnits)
{
    expect_near_relative(erf_esatw(amrex::Real(273.15), true), kEmpiricalTriplePointMbar);

    const std::array<amrex::Real, 3> temperatures = {
        amrex::Real(273.15), amrex::Real(300.0), amrex::Real(343.15)};

    for (const amrex::Real temperature : temperatures) {
        const amrex::Real empirical = erf_esatw(temperature, true);
        const amrex::Real default_fit = erf_esatw(temperature);

        EXPECT_GT(empirical, amrex::Real(0.0));
        expect_near_relative(empirical, default_fit, kEmpiricalAgreementRelTol);
    }
}

// Motivation: The old low-temperature polynomial produced negative vapor
// pressures near 188 K. These regression points must stay positive.
TEST(MicrophysicsWaterSaturation, RegressionCasesStayPositive)
{
    const std::array<amrex::Real, 4> regression_temperatures = {
        amrex::Real(188.16), amrex::Real(188.17), amrex::Real(189.0), amrex::Real(189.17)};

    for (const amrex::Real temperature : regression_temperatures) {
        SCOPED_TRACE(std::to_string(temperature));
        EXPECT_GT(erf_esatw(temperature), amrex::Real(0.0));
    }
}

// Motivation: The water saturation helper should never go negative anywhere in
// the supported cold-to-warm range used by the current branch logic.
TEST(MicrophysicsWaterSaturation, NoNegativeValuesOnRepresentativeGrid)
{
    for (int index = 0; index <= 64; ++index) {
        const amrex::Real temperature = amrex::Real(188.16) + amrex::Real(2.5) * index;
        EXPECT_GT(erf_esatw(temperature), amrex::Real(0.0));
    }
}

// Motivation: The low and high water-branch switches do not need to be exactly
// continuous, but they must avoid the large jump created by the old cold fit.
TEST(MicrophysicsWaterSaturation, SwitchJumpsStayBounded)
{
    const amrex::Real low_below = erf_esatw(kWaterLowerSwitchK - amrex::Real(1.0e-3));
    const amrex::Real low_above = erf_esatw(kWaterLowerSwitchK + amrex::Real(1.0e-3));
    const amrex::Real high_below = erf_esatw(kWaterUpperSwitchK - amrex::Real(1.0e-3));
    const amrex::Real high_above = erf_esatw(kWaterUpperSwitchK + amrex::Real(1.0e-3));

    EXPECT_NEAR(low_above, low_below, scaled_tol(low_above, low_below, kSwitchValueRelTol));
    EXPECT_NEAR(high_above, high_below, scaled_tol(high_above, high_below, kSwitchValueRelTol));
}

// Motivation: The value and derivative use the same branch selector, but the
// two fitted forms are independent. This test characterizes the derivative jump
// directly rather than using a finite difference across the switch.
TEST(MicrophysicsWaterSaturation, SwitchDerivativeJumpsStayBounded)
{
    const amrex::Real low_below = erf_dtesatw(kWaterLowerSwitchK - amrex::Real(1.0e-3));
    const amrex::Real low_above = erf_dtesatw(kWaterLowerSwitchK + amrex::Real(1.0e-3));
    const amrex::Real high_below = erf_dtesatw(kWaterUpperSwitchK - amrex::Real(1.0e-3));
    const amrex::Real high_above = erf_dtesatw(kWaterUpperSwitchK + amrex::Real(1.0e-3));

    EXPECT_NEAR(low_above, low_below, scaled_tol(low_above, low_below, kSwitchDerivativeRelTol));
    EXPECT_NEAR(high_above, high_below, scaled_tol(high_above, high_below, kSwitchDerivativeRelTol));
}

// Motivation: The water branch selector is tied to the production switch in
// erf_use_positive_esatw_poly. The positivity clause is currently defensive:
// within the active [-70, 70] C switch interval the Flatau polynomial stays
// positive on this branch.
TEST(MicrophysicsWaterSaturation, BranchSelectorBoundariesAndPositivityGuard)
{
    EXPECT_FALSE(erf_use_positive_esatw_poly(kWaterLowerSwitchK - amrex::Real(1.0e-3)));
    EXPECT_TRUE(erf_use_positive_esatw_poly(amrex::Real(273.16)));
    EXPECT_FALSE(erf_use_positive_esatw_poly(kWaterUpperSwitchK + amrex::Real(1.0e-3)));

    for (int dtt_c = -69; dtt_c <= 69; ++dtt_c) {
        const amrex::Real dtt = amrex::Real(dtt_c);
        SCOPED_TRACE(std::to_string(static_cast<double>(dtt)));
        EXPECT_GT(erf_esatw_flatau_poly(dtt), amrex::Real(0.0));
    }
}

// Motivation: The water saturation helper and its derivative must represent the
// same piecewise function away from the explicit branch transitions.
TEST(MicrophysicsWaterSaturation, DerivativeMatchesFiniteDifferenceAwayFromSwitches)
{
    const std::array<amrex::Real, 4> temperatures = {
        amrex::Real(190.0), amrex::Real(210.0), amrex::Real(300.0), amrex::Real(345.0)};
    amrex::Real previous = -amrex::Real(1.0);

    for (const amrex::Real temperature : temperatures) {
        const amrex::Real esat = erf_esatw(temperature);
        const amrex::Real derivative = erf_dtesatw(temperature);
        const amrex::Real finite_difference =
            central_difference([] (const amrex::Real value) { return erf_esatw(value); }, temperature);

        EXPECT_GT(esat, amrex::Real(0.0));
        EXPECT_GT(derivative, amrex::Real(0.0));
        expect_near_relative(derivative, finite_difference, kDerivativeRelTol);

        if (previous > amrex::Real(0.0)) {
            EXPECT_GT(esat, previous);
        }
        previous = esat;
    }
}

// Motivation: Ice saturation intentionally clamps to zero above freezing.
// These checks document that preserved discontinuity explicitly.
TEST(MicrophysicsIceSaturation, PreservesAboveFreezingClamp)
{
    expect_near_relative(erf_esati(amrex::Real(273.16)), amrex::Real(6.11147274));
    EXPECT_EQ(erf_esati(amrex::Real(273.160001)), amrex::Real(0.0));
    expect_near_relative(erf_dtesati(amrex::Real(273.16)), amrex::Real(0.503223089));
    EXPECT_EQ(erf_dtesati(amrex::Real(273.160001)), amrex::Real(0.0));
}

// Motivation: The ice value and derivative polynomials should stay positive and
// agree with finite differences over the documented derivative domain.
TEST(MicrophysicsIceSaturation, DerivativeMatchesFiniteDifferenceWithinContract)
{
    const std::array<amrex::Real, 4> temperatures = {
        amrex::Real(188.1601), amrex::Real(210.0), amrex::Real(240.0), amrex::Real(260.0)};
    amrex::Real previous = -amrex::Real(1.0);

    for (const amrex::Real temperature : temperatures) {
        const amrex::Real esat = erf_esati(temperature);
        const amrex::Real derivative = erf_dtesati(temperature);

        EXPECT_GT(esat, amrex::Real(0.0));
        EXPECT_GT(derivative, amrex::Real(0.0));

        if (temperature > amrex::Real(188.16)) {
            const amrex::Real finite_difference =
                central_difference([] (const amrex::Real value) { return erf_esati(value); }, temperature);
            expect_near_relative(derivative, finite_difference, kDerivativeRelTol);
        }

        if (previous > amrex::Real(0.0)) {
            EXPECT_GT(esat, previous);
        }
        previous = esat;
    }
}

// Motivation: The ice saturation value fit is valid down to about 183.16 K.
// Invalid-domain behavior is enforced in production with AMREX_ALWAYS_ASSERT,
// but CI does not death-test those AMReX assertion paths because they can hang
// under the AMReX-initialized GoogleTest/CTest setup.
TEST(MicrophysicsIceSaturation, ValueColdContract)
{
    EXPECT_GT(erf_esati(amrex::Real(183.16) + amrex::Real(1.0e-3)), amrex::Real(0.0));
}

// Motivation: The derivative polynomial has a narrower cold contract, valid
// down to about 188.16 K.
TEST(MicrophysicsIceSaturation, DerivativeColdContract)
{
    EXPECT_GT(erf_dtesati(amrex::Real(188.16) + amrex::Real(1.0e-3)), amrex::Real(0.0));
}

// Motivation: Below freezing, saturation vapor pressure over ice should stay
// below the corresponding supercooled-water saturation vapor pressure.
TEST(MicrophysicsIceSaturation, IcePressureStaysBelowSupercooledWater)
{
    for (int temperature = 230; temperature <= 270; temperature += 5) {
        const amrex::Real t = amrex::Real(temperature);
        EXPECT_LT(erf_esati(t), erf_esatw(t));
    }
}

// Motivation: The qsat helpers should differentiate the actual capped
// implementation, not the uncapped algebraic form.
TEST(MicrophysicsQSat, WaterNormalBranchMatchesFiniteDifference)
{
    const std::array<std::pair<amrex::Real, amrex::Real>, 2> states = {{
        {amrex::Real(260.0), amrex::Real(800.0)},
        {amrex::Real(300.0), amrex::Real(800.0)}
    }};

    for (const auto& state : states) {
        const amrex::Real temperature = state.first;
        const amrex::Real pressure = state.second;
        const amrex::Real esat = erf_esatw(temperature);
        SCOPED_TRACE(tp_label("water-normal", temperature, pressure));
        ASSERT_TRUE(uses_normal_qsat_branch(esat, pressure));

        amrex::Real qsat;
        amrex::Real dqsat;
        erf_qsatw(temperature, pressure, qsat);
        erf_dtqsatw(temperature, pressure, dqsat);

        const amrex::Real expected = Rd_on_Rv * erf_dtesatw(temperature) * pressure /
                                     ((pressure - esat) * (pressure - esat));
        const amrex::Real finite_difference = central_difference(
            [pressure] (const amrex::Real value) {
                amrex::Real qsat_local;
                erf_qsatw(value, pressure, qsat_local);
                return qsat_local;
            },
            temperature);

        EXPECT_GE(qsat, amrex::Real(0.0));
        EXPECT_LE(qsat, Rd_on_Rv);
        EXPECT_GT(dqsat, amrex::Real(0.0));
        expect_near_relative(dqsat, expected, kDerivativeRelTol);
        expect_near_relative(dqsat, finite_difference, kDerivativeRelTol);
    }
}

// Motivation: The uncapped qsat helper should preserve the same algebraic
// information as the underlying saturation vapor pressure when inverted.
TEST(MicrophysicsQSat, WaterNormalBranchRoundTripRecoversEsat)
{
    const std::array<std::pair<amrex::Real, amrex::Real>, 2> states = {{
        {amrex::Real(260.0), amrex::Real(800.0)},
        {amrex::Real(300.0), amrex::Real(800.0)}
    }};

    for (const auto& state : states) {
        const amrex::Real temperature = state.first;
        const amrex::Real pressure = state.second;
        const amrex::Real esat = erf_esatw(temperature);
        SCOPED_TRACE(tp_label("water-roundtrip", temperature, pressure));
        ASSERT_TRUE(uses_normal_qsat_branch(esat, pressure));

        amrex::Real qsat;
        erf_qsatw(temperature, pressure, qsat);
        const amrex::Real recovered_esat = pressure * qsat / (Rd_on_Rv + qsat);

        expect_near_relative(recovered_esat, esat);
    }
}

// Motivation: Once the qsat cap is active, qsat should be fixed at Rd_on_Rv and
// its temperature derivative should collapse to zero.
TEST(MicrophysicsQSat, WaterCappedBranchIsConstant)
{
    const amrex::Real temperature = amrex::Real(300.0);
    const amrex::Real pressure = capped_branch_pressure(erf_esatw(temperature));

    amrex::Real qsat;
    amrex::Real dqsat;
    erf_qsatw(temperature, pressure, qsat);
    erf_dtqsatw(temperature, pressure, dqsat);

    expect_near_relative(qsat, Rd_on_Rv);
    expect_near_relative(dqsat, amrex::Real(0.0));
}

// Motivation: Ice qsat uses the same capped formula and must preserve the same
// normal-branch and capped-branch contracts.
TEST(MicrophysicsQSat, IceNormalAndCappedBranchesAreConsistent)
{
    const amrex::Real temperature = amrex::Real(260.0);
    const amrex::Real normal_pressure = amrex::Real(800.0);
    const amrex::Real capped_pressure = microphysics_test::capped_branch_pressure(erf_esati(temperature));

    amrex::Real qsat_normal;
    amrex::Real dqsat_normal;
    amrex::Real qsat_capped;
    amrex::Real dqsat_capped;
    erf_qsati(temperature, normal_pressure, qsat_normal);
    erf_dtqsati(temperature, normal_pressure, dqsat_normal);
    erf_qsati(temperature, capped_pressure, qsat_capped);
    erf_dtqsati(temperature, capped_pressure, dqsat_capped);

    const amrex::Real esat = erf_esati(temperature);
    ASSERT_TRUE(uses_normal_qsat_branch(esat, normal_pressure));
    EXPECT_GT(dqsat_normal, amrex::Real(0.0));
    EXPECT_GE(qsat_normal, amrex::Real(0.0));
    EXPECT_LE(qsat_normal, Rd_on_Rv);
    expect_near_relative(qsat_capped, Rd_on_Rv);
    expect_near_relative(dqsat_capped, amrex::Real(0.0));
}

// Motivation: The public ice qsat-derivative wrapper should preserve the
// above-freezing clamp instead of relying on callers to avoid warm inputs.
TEST(MicrophysicsQSat, IceDerivativeClampsAboveFreezing)
{
    amrex::Real dqsat;
    erf_dtqsati(amrex::Real(273.160001), amrex::Real(800.0), dqsat);
    EXPECT_EQ(dqsat, amrex::Real(0.0));
}

// Motivation: The ice qsat helper should agree with both the finite-difference
// derivative and the algebraic inversion of the uncapped branch.
TEST(MicrophysicsQSat, IceNormalBranchFiniteDifferenceAndRoundTrip)
{
    const std::array<std::pair<amrex::Real, amrex::Real>, 2> states = {{
        {amrex::Real(240.0), amrex::Real(800.0)},
        {amrex::Real(260.0), amrex::Real(800.0)}
    }};

    for (const auto& state : states) {
        const amrex::Real temperature = state.first;
        const amrex::Real pressure = state.second;
        const amrex::Real esat = erf_esati(temperature);
        SCOPED_TRACE(tp_label("ice-roundtrip", temperature, pressure));
        ASSERT_TRUE(uses_normal_qsat_branch(esat, pressure));

        amrex::Real qsat;
        amrex::Real dqsat;
        erf_qsati(temperature, pressure, qsat);
        erf_dtqsati(temperature, pressure, dqsat);

        const amrex::Real recovered_esat = pressure * qsat / (Rd_on_Rv + qsat);
        const amrex::Real finite_difference = central_difference(
            [pressure] (const amrex::Real value) {
                amrex::Real qsat_local;
                erf_qsati(value, pressure, qsat_local);
                return qsat_local;
            },
            temperature);

        expect_near_relative(recovered_esat, esat);
        expect_near_relative(dqsat, finite_difference, kDerivativeRelTol);
    }
}

// Motivation: These are physical property checks rather than regression values:
// qsat should stay bounded, increase with temperature, and decrease with
// pressure while it remains on the uncapped branch.
TEST(MicrophysicsQSat, BoundsAndMonotonicityHoldOnNormalBranch)
{
    const std::array<amrex::Real, 4> temperatures = {
        amrex::Real(240.0), amrex::Real(260.0), amrex::Real(280.0), amrex::Real(300.0)};
    const std::array<amrex::Real, 3> pressures = {
        amrex::Real(500.0), amrex::Real(800.0), amrex::Real(1000.0)};

    amrex::Real previous_water = -amrex::Real(1.0);
    amrex::Real previous_ice = -amrex::Real(1.0);

    for (const amrex::Real temperature : temperatures) {
        amrex::Real qsat_water;
        amrex::Real qsat_ice;
        erf_qsatw(temperature, amrex::Real(800.0), qsat_water);
        erf_qsati(std::min(temperature, amrex::Real(260.0)), amrex::Real(800.0), qsat_ice);

        EXPECT_GE(qsat_water, amrex::Real(0.0));
        EXPECT_LE(qsat_water, Rd_on_Rv);

        if (temperature <= amrex::Real(260.0)) {
            EXPECT_GE(qsat_ice, amrex::Real(0.0));
            EXPECT_LE(qsat_ice, Rd_on_Rv);
            if (previous_ice > amrex::Real(0.0)) {
                EXPECT_GT(qsat_ice, previous_ice);
            }
            previous_ice = qsat_ice;
        }

        if (previous_water > amrex::Real(0.0)) {
            EXPECT_GT(qsat_water, previous_water);
        }
        previous_water = qsat_water;
    }

    for (const amrex::Real temperature : {amrex::Real(260.0), amrex::Real(300.0)}) {
        amrex::Real previous = -amrex::Real(1.0);

        for (const amrex::Real pressure : pressures) {
            amrex::Real qsat;
            erf_qsatw(temperature, pressure, qsat);

            if (previous > amrex::Real(0.0)) {
                EXPECT_LT(qsat, previous);
            }
            previous = qsat;
        }
    }
}

// Motivation: The empirical branch is only useful if the downstream HSE helper
// converts both water-pressure fits from mbar to Pa the same way.
TEST(MicrophysicsHSE, SaturationPressureUsesMbarToPaForBothBranches)
{
    const std::array<amrex::Real, 2> temperatures = {
        amrex::Real(273.15), amrex::Real(300.0)};

    for (const amrex::Real temperature : temperatures) {
        expect_near_relative(HSEutils::compute_saturation_pressure(temperature, false),
                             amrex::Real(100.0) * erf_esatw(temperature, false));
        expect_near_relative(HSEutils::compute_saturation_pressure(temperature, true),
                             amrex::Real(100.0) * erf_esatw(temperature, true));
    }

    expect_near_relative(HSEutils::compute_saturation_pressure(amrex::Real(273.15), true),
                         kEmpiricalTriplePointPa);
}