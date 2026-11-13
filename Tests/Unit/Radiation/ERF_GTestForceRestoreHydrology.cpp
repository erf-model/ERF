#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include <AMReX_MultiFab.H>

#include <ERF_ForceRestoreHydrology.H>
#include <ERF_ForceRestoreHydrologyState.H>
#include <ERF_RadStruct.H>

// One step of the force-restore hydrology (ERF_ForceRestoreHydrology.H) against closed forms:
// rain into the surface layer (NP89 Eq. 11), the column's water (NP89 Eq. 12), exact
// drainage toward field capacity (MN96 Eq. 4), exact relaxation of the deep temperature
// (NP89 Eq. 6), runoff above the porosity, the deficit, and the column budget.

namespace {
using amrex::Real;
constexpr Real tol = sizeof(Real) == 8 ? Real(1.0e-12) : Real(1.0e-5);

ForceRestoreHydrologyParams base ()
{
    ForceRestoreHydrologyParams p;
    p.d_s = Real(0.1);
    p.d_2 = Real(1.0);
    p.tau = Real(86400.0);
    p.w_fc = Real(0.3);
    p.w_sat = Real(0.45);
    p.c3 = Real(0.5);
    p.q_min = Real(0.0);
    p.q_max = Real(1.0);
    return p;
}
}

TEST(ForceRestoreHydrology, RainRateIsTheAccumulationChangeOverTheStep)
{
    EXPECT_NEAR(force_restore_rain_rate(Real(2.5), Real(1.0), Real(10.0)), Real(0.15), tol);
    // No step, no change, a reset accumulator, or a non-finite value: no rain.
    EXPECT_EQ(force_restore_rain_rate(Real(2.0), Real(1.0), Real(0.0)), Real(0.0));
    EXPECT_EQ(force_restore_rain_rate(Real(2.0), Real(1.0), Real(-1.0)), Real(0.0));
    EXPECT_EQ(force_restore_rain_rate(Real(1.0), Real(1.0), Real(10.0)), Real(0.0));
    EXPECT_EQ(force_restore_rain_rate(Real(0.0), Real(5.0), Real(10.0)), Real(0.0));
    EXPECT_EQ(force_restore_rain_rate(std::numeric_limits<Real>::quiet_NaN(), Real(0.0),
                                      Real(10.0)), Real(0.0));
    EXPECT_EQ(force_restore_rain_rate(std::numeric_limits<Real>::infinity(), Real(0.0),
                                      Real(10.0)), Real(0.0));
}

TEST(ForceRestoreHydrology, AllOffChangesNothing)
{
    const ForceRestoreHydrologyParams prm = base();  // rain, deep, drain all false
    Real qs = Real(0.25), qd = Real(0.35), td = Real(285.0);
    const auto f = force_restore_hydrology_step(qs, qd, td, Real(300.0), Real(1.0e-3),
                                                Real(1.0e-4), Real(60.0), prm);
    EXPECT_EQ(qs, Real(0.25));
    EXPECT_EQ(qd, Real(0.35));
    EXPECT_EQ(td, Real(285.0));
    EXPECT_EQ(f.runoff, Real(0.0));
    EXPECT_EQ(f.drainage, Real(0.0));
    EXPECT_EQ(f.deficit, Real(0.0));
}

TEST(ForceRestoreHydrology, RainEntersTheSurfaceLayer)
{
    ForceRestoreHydrologyParams prm = base();
    prm.rain = true;
    const Real p = Real(2.0e-3), dt = Real(30.0);
    Real qs = Real(0.2), qd = Real(0.3), td = Real(290.0);
    const auto f = force_restore_hydrology_step(qs, qd, td, Real(300.0), p, Real(0.0), dt, prm);
    // NP89 Eq. (11), C1 = 1: dq_s = P dt / (rho_w d_s) = 2e-3 * 30 / 100 = 6e-4.
    EXPECT_NEAR(qs, Real(0.2) + p * dt / (Real(1000.0) * prm.d_s), tol);
    EXPECT_EQ(qd, Real(0.3));       // no deep reservoir
    EXPECT_EQ(td, Real(290.0));
    EXPECT_EQ(f.runoff, Real(0.0));
}

TEST(ForceRestoreHydrology, ColumnTakesRainMinusEvaporation)
{
    ForceRestoreHydrologyParams prm = base();
    prm.rain = true;
    prm.deep = true;
    const Real p = Real(1.0e-3), e = Real(2.0e-4), dt = Real(100.0);
    Real qs = Real(0.2), qd = Real(0.25), td = Real(290.0);
    force_restore_hydrology_step(qs, qd, td, Real(290.0), p, e, dt, prm);
    // NP89 Eq. (12): dq_deep = (P - E) dt / (rho_w d_2) = 8e-4 * 100 / 1000 = 8e-5.
    EXPECT_NEAR(qd, Real(0.25) + (p - e) * dt / (Real(1000.0) * prm.d_2), tol);

    // Without seb_precipitation_source a rain rate handed in is ignored: the column loses
    // the evaporation alone and the surface layer gains nothing.
    prm.rain = false;
    Real qs2 = Real(0.2), qd2 = Real(0.25), td2 = Real(290.0);
    force_restore_hydrology_step(qs2, qd2, td2, Real(290.0), p, e, dt, prm);
    EXPECT_NEAR(qd2, Real(0.25) - e * dt / (Real(1000.0) * prm.d_2), tol);
    EXPECT_EQ(qs2, Real(0.2));
}

TEST(ForceRestoreHydrology, DrainageIsExactExponentialTowardFieldCapacity)
{
    ForceRestoreHydrologyParams prm = base();
    prm.deep = true;
    prm.drain = true;
    const Real q0 = Real(0.42), dt = Real(3600.0);
    Real qs = Real(0.2), qd = q0, td = Real(290.0);
    const auto f = force_restore_hydrology_step(qs, qd, td, Real(290.0), Real(0.0), Real(0.0),
                                                dt, prm);
    const Real expect = prm.w_fc + (q0 - prm.w_fc) * std::exp(-prm.c3 * dt / prm.tau);
    EXPECT_NEAR(qd, expect, tol);
    EXPECT_NEAR(f.drainage, Real(1000.0) * prm.d_2 * (q0 - expect), Real(1.0e3) * tol);

    // A step much longer than tau / C3 stops at field capacity and does not pass it.
    Real qd2 = q0;
    force_restore_hydrology_step(qs, qd2, td, Real(290.0), Real(0.0), Real(0.0),
                                 Real(1.0e9), prm);
    EXPECT_NEAR(qd2, prm.w_fc, tol);
    EXPECT_GE(qd2, prm.w_fc);

    // MN96 Table 1 gives C_3 for tau = 1 day, so the drainage does not depend on the restore
    // period of T_2 (seb_restore_timescale_s).
    ForceRestoreHydrologyParams prm_tau = prm;
    prm_tau.tau = Real(3600.0);
    Real qd_tau = q0;
    force_restore_hydrology_step(qs, qd_tau, td, Real(290.0), Real(0.0), Real(0.0), dt, prm_tau);
    EXPECT_NEAR(qd_tau, expect, tol);
    EXPECT_NEAR(expect, prm.w_fc + (q0 - prm.w_fc) * std::exp(-prm.c3 * dt / Real(86400.0)), tol);

    // Below field capacity there is no drainage.
    Real qd3 = Real(0.2);
    const auto f3 = force_restore_hydrology_step(qs, qd3, td, Real(290.0), Real(0.0), Real(0.0),
                                                 dt, prm);
    EXPECT_EQ(qd3, Real(0.2));
    EXPECT_EQ(f3.drainage, Real(0.0));
}

TEST(ForceRestoreHydrology, DeepTemperatureRelaxesExactly)
{
    ForceRestoreHydrologyParams prm = base();
    prm.deep = true;
    const Real ts = Real(300.0), t0 = Real(285.0), dt = Real(7200.0);
    Real qs = Real(0.2), qd = Real(0.25), td = t0;
    force_restore_hydrology_step(qs, qd, td, ts, Real(0.0), Real(0.0), dt, prm);
    // NP89 Eq. (6), exact: T2 = Ts + (T2_0 - Ts) exp(-dt / tau).
    EXPECT_NEAR(td, ts + (t0 - ts) * std::exp(-dt / prm.tau), Real(300.0) * tol);
    Real td2 = t0;
    force_restore_hydrology_step(qs, qd, td2, ts, Real(0.0), Real(0.0), Real(1.0e9), prm);
    EXPECT_NEAR(td2, ts, Real(300.0) * tol);
}

TEST(ForceRestoreHydrology, WaterAboveSaturationRunsOff)
{
    ForceRestoreHydrologyParams prm = base();
    prm.rain = true;
    prm.deep = true;
    const Real p = Real(0.05), dt = Real(100.0);   // 5 kg/m^2 into a 0.1 m layer: +0.05
    Real qs = Real(0.43), qd = Real(0.30), td = Real(290.0);
    const Real qd0 = qd;
    const auto f = force_restore_hydrology_step(qs, qd, td, Real(290.0), p, Real(0.0), dt, prm);
    // q_s would reach 0.48; 0.03 * 0.1 m * 1000 = 3 kg/m^2 runs off and leaves the column.
    EXPECT_NEAR(qs, prm.w_sat, tol);
    EXPECT_NEAR(f.runoff, Real(3.0), Real(1.0e3) * tol);
    EXPECT_NEAR(qd, qd0 + (p * dt - f.runoff) / (Real(1000.0) * prm.d_2), tol);

    // A column above saturation spills its excess.
    Real qs2 = Real(0.2), qd2 = Real(0.46);
    ForceRestoreHydrologyParams prm2 = base();
    prm2.deep = true;
    const auto f2 = force_restore_hydrology_step(qs2, qd2, td, Real(290.0), Real(0.0), Real(0.0),
                                                 dt, prm2);
    EXPECT_NEAR(qd2, prm2.w_sat, tol);
    EXPECT_NEAR(f2.runoff, Real(1000.0) * (Real(0.46) - prm2.w_sat), Real(1.0e3) * tol);
}

TEST(ForceRestoreHydrology, SurfaceRunoffIsAtMostTheStepsRain)
{
    // q_s already above the porosity (a restore overshoot, or a restart from a run without
    // the hydrology): only the rain of this step runs off and leaves the column; the rest
    // is clamped off q_s and the column keeps it.
    ForceRestoreHydrologyParams prm = base();
    prm.rain = true;
    prm.deep = true;
    const Real p = Real(1.0e-6), dt = Real(7200.0);   // 7.2e-3 kg/m^2 of rain
    Real qs = Real(0.70), qd = Real(0.40), td = Real(290.0);
    const auto f = force_restore_hydrology_step(qs, qd, td, Real(290.0), p, Real(0.0), dt, prm);
    EXPECT_NEAR(qs, prm.w_sat, tol);
    EXPECT_NEAR(f.runoff, p * dt, Real(1.0e-3) * tol);
    EXPECT_NEAR(qd, Real(0.40), tol);
    EXPECT_EQ(f.deficit, Real(0.0));
}

TEST(ForceRestoreHydrology, OnlyInfiltratedRainDrains)
{
    // q_s = 0.40, w_sat = 0.45, d_s = 0.1 m: the layer holds 5 kg/m^2 more. 10 kg/m^2 of rain
    // falls, so 5 kg/m^2 runs off and 5 kg/m^2 enters the column (q_deep 0.300 -> 0.305) before
    // it drains toward w_fc = 0.3 with C3 = 0.5 over tau_1d. The runoff never drains.
    ForceRestoreHydrologyParams prm = base();
    prm.rain = true;
    prm.deep = true;
    prm.drain = true;
    const Real p = Real(0.01), dt = Real(1000.0);
    Real qs = Real(0.40), qd = Real(0.30), td = Real(290.0);
    const auto f = force_restore_hydrology_step(qs, qd, td, Real(290.0), p, Real(0.0), dt, prm);
    const Real decay = std::exp(-Real(0.5) * dt / Real(86400.0));
    const Real qd_expected = Real(0.30) + Real(0.005) * decay;
    EXPECT_NEAR(f.runoff, Real(5.0), Real(1.0e3) * tol);
    EXPECT_NEAR(qs, prm.w_sat, tol);
    EXPECT_NEAR(qd, qd_expected, tol);
    EXPECT_NEAR(f.drainage, Real(5.0) * (Real(1.0) - decay), Real(1.0e3) * tol);
    EXPECT_EQ(f.deficit, Real(0.0));

    // A saturated surface layer over a column at field capacity: all rain runs off, nothing
    // enters, nothing drains, and no deficit is recorded.
    Real qs2 = Real(0.45), qd2 = Real(0.30), td2 = Real(290.0);
    const auto f2 = force_restore_hydrology_step(qs2, qd2, td2, Real(290.0), p, Real(0.0), dt, prm);
    EXPECT_NEAR(f2.runoff, p * dt, Real(1.0e3) * tol);
    EXPECT_NEAR(qd2, Real(0.30), tol);
    EXPECT_EQ(f2.drainage, Real(0.0));
    EXPECT_EQ(f2.deficit, Real(0.0));
}

TEST(ForceRestoreHydrology, EvaporationTheColumnCannotSupplyIsADeficit)
{
    ForceRestoreHydrologyParams prm = base();
    prm.deep = true;
    const Real e = Real(1.0e-3), dt = Real(1000.0);  // 1 kg/m^2 from a column holding 0.5
    Real qs = Real(0.0), qd = Real(5.0e-4), td = Real(290.0);
    const auto f = force_restore_hydrology_step(qs, qd, td, Real(290.0), Real(0.0), e, dt, prm);
    EXPECT_EQ(qd, Real(0.0));
    EXPECT_NEAR(f.deficit, Real(0.5), Real(1.0e3) * tol);
}

TEST(ForceRestoreHydrology, ColumnBudgetCloses)
{
    // Rain, evaporation, drainage and runoff all active in one step.
    ForceRestoreHydrologyParams prm = base();
    prm.rain = true;
    prm.deep = true;
    prm.drain = true;
    prm.d_2 = Real(0.5);
    const Real cases[][4] = {
        // q_s, q_deep, P, E
        {Real(0.44), Real(0.40), Real(0.02), Real(1.0e-4)},
        {Real(0.10), Real(0.31), Real(1.0e-4), Real(5.0e-5)},
        {Real(0.05), Real(1.0e-4), Real(0.0), Real(2.0e-3)},
        {Real(0.45), Real(0.449), Real(0.05), Real(0.0)},
        {Real(0.90), Real(0.20), Real(1.0e-5), Real(0.0)},   // q_s above w_sat before the rain
    };
    for (const auto& c : cases) {
        Real qs = c[0], qd = c[1], td = Real(290.0);
        const Real qd0 = qd;
        const Real dt = Real(300.0);
        const auto f = force_restore_hydrology_step(qs, qd, td, Real(295.0), c[2], c[3], dt, prm);
        const Real lhs = Real(1000.0) * prm.d_2 * (qd - qd0);
        const Real rhs = (c[2] - c[3]) * dt - f.drainage - f.runoff + f.deficit;
        EXPECT_NEAR(lhs, rhs, Real(1.0e3) * tol);
        EXPECT_GE(qd, Real(0.0));
        EXPECT_LE(qd, prm.w_sat);
        EXPECT_LE(qs, prm.w_sat);
    }
}

TEST(ForceRestoreHydrology, RainOfASkippedStepIsCarriedIntoTheNext)
{
    // A level whose balance skips a step (a new level's first step without its radiation
    // sweep) calls begin_step and end_step but not step(). The rain of that step must reach
    // the soil with the next step: begin_step keeps the earlier start of the accumulation,
    // and the rate is the mean over the time since that start, not over the last step alone.
    using namespace amrex;
    RadChoice rc;
    rc.seb_precipitation_source = SEBPrecipitationSource::Microphysics;
    const BoxArray ba3d(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba3d);
    MultiFab accum(ba3d, dm, 1, 0);
    SurfacePrecipAccumulationSources src;
    src.total.accum = &accum;
    const int klo = 0;

    ForceRestoreHydrology hydro;
    hydro.resize(1);
    hydro.define_level(0, rc, ba2d, dm, IntVect(0));
    ASSERT_TRUE(hydro.defined(0));

    MultiFab t_s(ba2d, dm, 1, 0), q_s(ba2d, dm, 1, 0), q_deep(ba2d, dm, 1, 0);
    MultiFab t_deep(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0);
    t_s.setVal(295.0); q_s.setVal(0.1); q_deep.setVal(0.1); t_deep.setVal(295.0); lh.setVal(0.0);
    ForceRestoreHydrologyParams prm = base();
    prm.rain = true;

    // Step 1 (dt = 1 s): 2 kg/m^2 of rain, skipped by the balance.
    accum.setVal(1.0);
    hydro.begin_step(0, src, klo);
    accum.setVal(3.0);
    hydro.end_step(0, src, klo, Real(1.0));
    // Step 2 (dt = 2 s): 1 kg/m^2 more; the soil takes in both steps' rain.
    hydro.begin_step(0, src, klo);
    accum.setVal(4.0);
    hydro.end_step(0, src, klo, Real(2.0));
    // 3 kg/m^2 over 1 + 2 s.
    EXPECT_NEAR(hydro.rain_rate(0)->min(0), Real(1.0), tol);
    EXPECT_NEAR(hydro.rain_rate(0)->max(0), Real(1.0), tol);
    const auto sums = hydro.step(0, t_s, q_s, q_deep, t_deep, lh, nullptr, Real(2.0), prm);
    EXPECT_EQ(sums.columns, Long(16));
    EXPECT_NEAR(hydro.totals(0)->min(ForceRestoreHydrology::Rain), Real(3.0), tol);
    EXPECT_NEAR(hydro.totals(0)->max(ForceRestoreHydrology::Rain), Real(3.0), tol);
    // All 3 kg/m^2 entered the 0.1 m surface layer: q_s 0.1 -> 0.13.
    EXPECT_NEAR(q_s.min(0), Real(0.13), tol);
    EXPECT_NEAR(q_s.max(0), Real(0.13), tol);

    // Step 3 (dt = 1 s), with step() run: only its own rain, 0.5 kg/m^2.
    hydro.begin_step(0, src, klo);
    accum.setVal(4.5);
    hydro.end_step(0, src, klo, Real(1.0));
    EXPECT_NEAR(hydro.rain_rate(0)->max(0), Real(0.5), tol);
    hydro.step(0, t_s, q_s, q_deep, t_deep, lh, nullptr, Real(1.0), prm);
    EXPECT_NEAR(hydro.totals(0)->max(ForceRestoreHydrology::Rain), Real(3.5), tol);

    // A level whose step() never runs: the rate stays the mean since the kept start, bounded
    // by the rain rate, instead of growing with each step.
    for (int n = 0; n < 4; ++n) {
        hydro.begin_step(0, src, klo);
        accum.setVal(Real(5.5) + Real(n));   // 1 kg/m^2 per 1 s step
        hydro.end_step(0, src, klo, Real(1.0));
    }
    EXPECT_NEAR(hydro.rain_rate(0)->max(0), Real(1.0), tol);
}
