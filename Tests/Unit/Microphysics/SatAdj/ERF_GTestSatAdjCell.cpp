#include <array>

#include <gtest/gtest.h>

#include "ERF_GTestSatAdjCommon.H"

// These tests exercise the scalar cell adjustment used by the production
// SatAdj kernel. They check physical properties of the adjustment rather than
// fixed output values.

using namespace satadj_test;

namespace {

void expect_theta_consistent (const CellState& state)
{
    const amrex::Real theta_expected =
        getThgivenTandP(state.tabs, amrex::Real(100.0) * state.pres_mbar, kRdOcp);
    EXPECT_NEAR(state.theta, theta_expected, scaled_tol(theta_expected, kStateTolFactor));
}

void expect_conservation (const CellState& initial, const CellState& final)
{
    const amrex::Real qt_initial = initial.qv + initial.qc;
    const amrex::Real qt_final = final.qv + final.qc;
    // SatAdj exchanges vapor and cloud water locally, so qt is conserved. At
    // fixed pressure, latent heating/cooling conserves T + (L/cp) qv.
    const amrex::Real moist_initial = initial.tabs + kFacCond * initial.qv;
    const amrex::Real moist_final = final.tabs + kFacCond * final.qv;

    EXPECT_NEAR(qt_final, qt_initial, scaled_tol(qt_initial, amrex::Real(20.0) * kStateTolFactor));
    EXPECT_NEAR(moist_final, moist_initial, scaled_tol(moist_initial, amrex::Real(20.0) * kStateTolFactor));
}

void expect_total_water_conservation (const CellState& initial, const CellState& final)
{
    const amrex::Real qt_initial = initial.qv + initial.qc;
    const amrex::Real qt_final = final.qv + final.qc;
    EXPECT_NEAR(qt_final, qt_initial, scaled_tol(qt_initial, amrex::Real(20.0) * kStateTolFactor));
}

void expect_complementarity (const CellState& state)
{
    const amrex::Real qsat_final = qsat(state.tabs, state.pres_mbar);
    EXPECT_GE(state.qc, -scaled_tol(state.qc, amrex::Real(20.0) * kStateTolFactor));
    EXPECT_LE(state.qv, qsat_final + scaled_tol(qsat_final, amrex::Real(20.0) * kStateTolFactor));

    if (state.qc > scaled_tol(amrex::Real(0.0), amrex::Real(100.0) * kStateTolFactor)) {
        EXPECT_NEAR(state.qv, qsat_final, scaled_tol(qsat_final, amrex::Real(50.0) * kStateTolFactor));
    }
}

} // namespace

// Motivation: Supersaturated cells should condense to a saturated cloudy
// state. This catches sign errors in latent heating and Newton residual
// updates.
TEST(SatAdjCell, SaturatedNewtonSolve)
{
    const amrex::Real tabs = amrex::Real(290.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    CellState state = make_cell_state(tabs, pres_mbar,
                                      qsat_initial + amrex::Real(8.0e-4),
                                      amrex::Real(4.0e-4));
    const CellState initial = state;

    adjust(state);

    EXPECT_NEAR(state.qv, qsat(state.tabs, state.pres_mbar),
                scaled_tol(state.qv, amrex::Real(50.0) * kStateTolFactor));
    EXPECT_GE(state.qc, amrex::Real(0.0));
    expect_conservation(initial, state);
    expect_theta_consistent(state);
}

// Motivation: If total water cannot support saturation, SatAdj should
// evaporate all cloud water and leave a clear subsaturated cell.
TEST(SatAdjCell, UnsaturatedAllEvaporationBranch)
{
    CellState state = make_cell_state(amrex::Real(290.0), amrex::Real(900.0),
                                      amrex::Real(4.0e-3), amrex::Real(1.0e-3));
    const CellState initial = state;

    adjust(state);

    EXPECT_EQ(state.qc, amrex::Real(0.0));
    EXPECT_NEAR(state.qv, initial.qv + initial.qc,
                scaled_tol(initial.qv + initial.qc, amrex::Real(5.0) * kStateTolFactor));
    EXPECT_NEAR(state.tabs, initial.tabs - kFacCond * initial.qc,
                scaled_tol(initial.tabs, amrex::Real(5.0) * kStateTolFactor));
    EXPECT_LE(state.qv, qsat(state.tabs, state.pres_mbar) +
                        scaled_tol(state.qv, amrex::Real(10.0) * kStateTolFactor));
    expect_theta_consistent(state);
}

// Motivation: Evaporating cloud water cools the cell and can lower qsat
// enough to require recondensation. This protects a subtle two-stage branch.
TEST(SatAdjCell, EvaporationThenRecondensationBranch)
{
    CellState state;
    ASSERT_TRUE(find_evaporation_then_recondensation_state(state));
    const CellState initial = state;

    adjust(state);

    EXPECT_GT(state.qc, amrex::Real(0.0));
    EXPECT_NEAR(state.qv, qsat(state.tabs, state.pres_mbar),
                scaled_tol(state.qv, amrex::Real(50.0) * kStateTolFactor));
    expect_conservation(initial, state);
    expect_theta_consistent(state);
}

// Motivation: For small supersaturation, the nonlinear adjustment should
// reduce to the linearized analytic condensate response.
TEST(SatAdjCell, SmallSupersaturationAsymptotic)
{
    // For small supersaturation s, linearizing the Newton residual gives
    // delta_qc ~= s / (1 + (L/cp) * d(qsat)/dT). The error should shrink as s
    // decreases.
    const amrex::Real tabs = amrex::Real(290.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    const amrex::Real dqsat_initial = dqsat_dt(tabs, pres_mbar);
    const std::array<amrex::Real, 3> supersaturation = {
        amrex::Real(1.0e-8), amrex::Real(1.0e-7), amrex::Real(1.0e-6)};

    std::array<amrex::Real, 3> errors{};

    for (int is = 0; is < static_cast<int>(supersaturation.size()); ++is) {
        CellState state = make_cell_state(tabs, pres_mbar, qsat_initial + supersaturation[is], amrex::Real(0.0));
        adjust(state);

        const amrex::Real expected = supersaturation[is] / (amrex::Real(1.0) + kFacCond * dqsat_initial);
        errors[is] = std::abs(state.qc - expected);

        EXPECT_NEAR(state.qc, expected,
                    std::max(std::abs(expected) * kAsymptoticRelTol,
                             scaled_tol(expected, amrex::Real(10.0) * kStateTolFactor)));
    }

    EXPECT_LE(errors[0], errors[1] + scaled_tol(errors[1], amrex::Real(10.0) * kStateTolFactor));
    EXPECT_LE(errors[1], errors[2] + scaled_tol(errors[2], amrex::Real(10.0) * kStateTolFactor));
}

// Motivation: SatAdj is a projection onto the final complementarity state. A
// second adjustment should leave an already adjusted cell unchanged.
TEST(SatAdjCell, Idempotence)
{
    const amrex::Real tabs = amrex::Real(290.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    CellState first = make_cell_state(tabs, pres_mbar,
                                      qsat_initial + amrex::Real(6.0e-4),
                                      amrex::Real(5.0e-4));

    adjust(first);
    CellState second = first;
    adjust(second);

    EXPECT_NEAR(second.tabs, first.tabs, scaled_tol(first.tabs, amrex::Real(20.0) * kStateTolFactor));
    EXPECT_NEAR(second.theta, first.theta, scaled_tol(first.theta, amrex::Real(20.0) * kStateTolFactor));
    EXPECT_NEAR(second.qv, first.qv, scaled_tol(first.qv, amrex::Real(20.0) * kStateTolFactor));
    EXPECT_NEAR(second.qc, first.qc, scaled_tol(first.qc, amrex::Real(20.0) * kStateTolFactor));
    expect_complementarity(second);
}

// Motivation: Negative qc can appear from transport overshoots. SatAdj should
// repair it while preserving total nonprecipitating water.
TEST(SatAdjCell, NegativeQcRepair)
{
    // Negative qc can arise from transport overshoots. SatAdj repairs it while
    // preserving total nonprecipitating water.
    const amrex::Real tabs = amrex::Real(288.0);
    const amrex::Real pres_mbar = amrex::Real(850.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    CellState state = make_cell_state(tabs, pres_mbar,
                                      qsat_initial + amrex::Real(4.0e-4),
                                      -amrex::Real(1.0e-4));
    const CellState initial = state;

    adjust(state);

    EXPECT_GE(state.qc, amrex::Real(0.0));
    expect_total_water_conservation(initial, state);
    expect_complementarity(state);
    expect_theta_consistent(state);
}