#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include <AMReX_MultiFab.H>

#include <ERF_ForceRestoreHydrology.H>
#include <ERF_ForceRestoreHydrologyState.H>
#include <ERF_ForceRestoreSnow.H>
#include <ERF_ForceRestoreSnowState.H>
#include <ERF_MicrophysicsConstants.H>
#include <ERF_RadStruct.H>
#include <ERF_SurfacePrecipitation.H>

// One step of the force-restore snow layer (ERF_ForceRestoreSnow.H) against the closed forms
// of Douville, Royer and Mahfouf (1995, "DRM95") Sec. 2.1: the cover fraction (Eq. 7), the
// albedo (Eq. 8 with veg = 0), the mass balance (Eq. 6) and the melt (Eq. 12); then the
// liquid/frozen split of the surface precipitation, and the state class with the hydrology.

namespace {
using amrex::Real;
constexpr Real tol = sizeof(Real) == 8 ? Real(1.0e-12) : Real(1.0e-5);
constexpr Real t00 = Real(273.16);     // DRM95 Eq. (12)
constexpr Real l_f = Real(3.337e5);    // ERF's latent heat of fusion [J/kg]

ForceRestoreSnowParams base ()
{
    ForceRestoreSnowParams p;
    p.w_crn = Real(10.0);
    p.albedo = Real(0.70);
    p.c_s = Real(2.0e4);
    return p;
}
}

TEST(ForceRestoreSnow, ConstantsAreThePapersAndERFs)
{
    EXPECT_EQ(force_restore_snow_t_melt, t00);
    EXPECT_EQ(lat_ice, l_f);
    const ForceRestoreSnowParams p;
    EXPECT_EQ(p.w_crn, Real(10.0));    // DRM95: W_crn = 10 mm
    EXPECT_EQ(p.albedo, Real(0.70));   // DRM95: alpha_n = 0.70
    RadChoice rc;
    EXPECT_FALSE(rc.seb_snow_enable);
    EXPECT_EQ(rc.seb_snow_critical_swe_kg_m2, Real(10.0));
    EXPECT_EQ(rc.seb_snow_albedo, Real(0.70));
    EXPECT_EQ(rc.seb_snow_initial_swe_kg_m2, Real(0.0));
}

TEST(ForceRestoreSnow, CoverFractionIsEq7)
{
    // p_ng = W_n / (W_n + W_crn)
    EXPECT_NEAR(force_restore_snow_fraction(Real(10.0), Real(10.0)), Real(0.5), tol);
    EXPECT_NEAR(force_restore_snow_fraction(Real(30.0), Real(10.0)), Real(0.75), tol);
    EXPECT_NEAR(force_restore_snow_fraction(Real(2.0), Real(8.0)), Real(0.2), tol);
    EXPECT_EQ(force_restore_snow_fraction(Real(0.0), Real(10.0)), Real(0.0));
    EXPECT_EQ(force_restore_snow_fraction(Real(-1.0), Real(10.0)), Real(0.0));
    EXPECT_EQ(force_restore_snow_fraction(std::numeric_limits<Real>::quiet_NaN(), Real(10.0)),
              Real(0.0));
    EXPECT_EQ(force_restore_snow_fraction(std::numeric_limits<Real>::infinity(), Real(10.0)),
              Real(0.0));
}

TEST(ForceRestoreSnow, AlbedoIsEq8ForBareGround)
{
    const ForceRestoreSnowParams prm = base();
    // (1 - p_ng) alpha_g + p_ng alpha_n: p_ng = 0.5 -> 0.5 * 0.2 + 0.5 * 0.7.
    EXPECT_NEAR(force_restore_snow_albedo(Real(0.2), Real(10.0), prm), Real(0.45), tol);
    // p_ng = 0.75 -> 0.25 * 0.2 + 0.75 * 0.7.
    EXPECT_NEAR(force_restore_snow_albedo(Real(0.2), Real(30.0), prm), Real(0.575), tol);
    // No snow: the ground's albedo, bit for bit.
    EXPECT_EQ(force_restore_snow_albedo(Real(0.2), Real(0.0), prm), Real(0.2));
}

TEST(ForceRestoreSnow, SnowfallBuildsTheSnowBelowFreezing)
{
    const ForceRestoreSnowParams prm = base();
    Real w = Real(1.0), ts = Real(260.0);
    const auto f = force_restore_snow_step(w, ts, Real(1.0e-3), Real(0.0), Real(60.0), prm);
    EXPECT_NEAR(f.snowfall, Real(0.06), tol);
    EXPECT_EQ(f.sublimation, Real(0.0));
    EXPECT_EQ(f.melt, Real(0.0));
    EXPECT_NEAR(w, Real(1.06), tol);
    EXPECT_EQ(ts, Real(260.0));
}

TEST(ForceRestoreSnow, PartialMeltIsEq12AndHoldsTheMeltingPoint)
{
    const ForceRestoreSnowParams prm = base();
    // T_s* = T_00 + 1 K: F_n dt = (T_s* - T_00) / (C_t L_f) = C_s * 1 / L_f.
    Real w = Real(5.0), ts = t00 + Real(1.0);
    const auto f = force_restore_snow_step(w, ts, Real(0.0), Real(0.0), Real(60.0), prm);
    const Real expected = Real(2.0e4) / l_f;   // 0.0599341... kg/m^2
    EXPECT_NEAR(f.melt, expected, Real(100.0) * tol);
    EXPECT_NEAR(f.melt, Real(0.059934), Real(1.0e-6));
    EXPECT_NEAR(w, Real(5.0) - expected, Real(100.0) * tol);
    EXPECT_EQ(ts, t00);
}

TEST(ForceRestoreSnow, FullMeltLeavesTheRestOfTheHeatInTheSurface)
{
    const ForceRestoreSnowParams prm = base();
    // 0.01 kg/m^2 of snow needs L_f * 0.01 = 3337 J/m^2: 0.16685 K of the surface's 2 K.
    Real w = Real(0.01), ts = t00 + Real(2.0);
    const auto f = force_restore_snow_step(w, ts, Real(0.0), Real(0.0), Real(60.0), prm);
    EXPECT_NEAR(f.melt, Real(0.01), tol);
    EXPECT_EQ(w, Real(0.0));
    EXPECT_NEAR(ts, t00 + Real(2.0) - Real(0.16685), sizeof(Real) == 8 ? Real(1.0e-9) : Real(1.0e-3));
    EXPECT_GT(ts, t00);
}

TEST(ForceRestoreSnow, MeltHeatEqualsTheSurfaceCooling)
{
    const ForceRestoreSnowParams prm = base();
    for (const Real w0 : {Real(0.001), Real(0.05), Real(0.3), Real(40.0)}) {
        for (const Real excess : {Real(0.01), Real(0.5), Real(4.0)}) {
            Real w = w0, ts = t00 + excess;
            const auto f = force_restore_snow_step(w, ts, Real(0.0), Real(0.0), Real(10.0), prm);
            EXPECT_NEAR(prm.c_s * (t00 + excess - ts), l_f * f.melt,
                        Real(1.0e-6) * prm.c_s * excess + Real(1.0e3) * tol);
            EXPECT_GE(ts, t00);
            // While snow remains the surface is at T_00, up to the heat of the round-off
            // of W_n, which stays in T_s.
            EXPECT_TRUE(w == Real(0.0) ||
                        ts - t00 <= (sizeof(Real) == 8 ? Real(1.0e-10) : Real(1.0e-3)));
        }
    }
}

TEST(ForceRestoreSnow, NoMeltAtOrBelowTheMeltingPointOrWithoutSnow)
{
    const ForceRestoreSnowParams prm = base();
    Real w = Real(5.0), ts = t00;
    auto f = force_restore_snow_step(w, ts, Real(0.0), Real(0.0), Real(60.0), prm);
    EXPECT_EQ(f.melt, Real(0.0));
    EXPECT_EQ(w, Real(5.0));
    EXPECT_EQ(ts, t00);
    // No snow: the step changes nothing, at any temperature.
    w = Real(0.0); ts = Real(300.0);
    f = force_restore_snow_step(w, ts, Real(0.0), Real(1.0e-4), Real(60.0), prm);
    EXPECT_EQ(f.snowfall, Real(0.0));
    EXPECT_EQ(f.sublimation, Real(0.0));
    EXPECT_EQ(f.melt, Real(0.0));
    EXPECT_EQ(w, Real(0.0));
    EXPECT_EQ(ts, Real(300.0));
}

TEST(ForceRestoreSnow, SublimationIsTheCoveredFractionOfTheEvaporation)
{
    const ForceRestoreSnowParams prm = base();
    // p_ng = 0.5 at W_n = W_crn: E_n dt = 0.5 * 1e-4 * 100.
    Real w = Real(10.0), ts = Real(265.0);
    auto f = force_restore_snow_step(w, ts, Real(0.0), Real(1.0e-4), Real(100.0), prm);
    EXPECT_NEAR(f.sublimation, Real(0.005), tol);
    EXPECT_NEAR(w, Real(9.995), Real(10.0) * tol);
    // At most the snow there is: p_ng = 1/3, E dt = 100 kg/m^2.
    w = Real(5.0);
    f = force_restore_snow_step(w, ts, Real(0.0), Real(1.0), Real(100.0), prm);
    EXPECT_NEAR(f.sublimation, Real(5.0), tol);
    EXPECT_EQ(w, Real(0.0));
    // A downward flux (frost) stays with the soil; a non-finite one does nothing.
    w = Real(5.0);
    f = force_restore_snow_step(w, ts, Real(0.0), Real(-1.0e-4), Real(100.0), prm);
    EXPECT_EQ(f.sublimation, Real(0.0));
    EXPECT_EQ(w, Real(5.0));
    f = force_restore_snow_step(w, ts, std::numeric_limits<Real>::quiet_NaN(),
                                std::numeric_limits<Real>::quiet_NaN(), Real(100.0), prm);
    EXPECT_EQ(f.snowfall, Real(0.0));
    EXPECT_EQ(f.sublimation, Real(0.0));
    EXPECT_EQ(w, Real(5.0));
}

TEST(ForceRestoreSnow, ColumnBudgetClosesEq6)
{
    const ForceRestoreSnowParams prm = base();
    Real w = Real(0.0), ts = Real(268.0);
    Real snowfall = Real(0.0), sublimation = Real(0.0), melt = Real(0.0);
    for (int n = 0; n < 400; ++n) {
        // Snow for the first half with a cold surface, then a warm surface melting it.
        const Real p = (n < 200) ? Real(2.0e-3) * (Real(1.0) + std::sin(Real(0.1) * Real(n))) : Real(0.0);
        const Real e = Real(2.0e-5) * std::cos(Real(0.05) * Real(n));
        ts = (n < 200) ? Real(268.0) : t00 + Real(0.4);
        const Real w_old = w;
        const auto f = force_restore_snow_step(w, ts, p, e, Real(30.0), prm);
        EXPECT_NEAR(w - w_old, f.snowfall - f.sublimation - f.melt, Real(100.0) * tol);
        EXPECT_GE(w, Real(0.0));
        EXPECT_GE(f.sublimation, Real(0.0));
        EXPECT_GE(f.melt, Real(0.0));
        snowfall += f.snowfall; sublimation += f.sublimation; melt += f.melt;
    }
    EXPECT_GT(snowfall, Real(10.0));
    EXPECT_GT(melt, Real(1.0));
    EXPECT_NEAR(w, snowfall - sublimation - melt, Real(1.0e4) * tol);
}

TEST(ForceRestoreSnow, StepReportsWhatTheSnowActuallyChangedBy)
{
    const ForceRestoreSnowParams prm = base();
    // A negative W_n (a corrupt state the driver refuses) must not sublimate a negative
    // amount, which would take water from the soil.
    {
        Real w = Real(-1.0), t = Real(280.0);
        const ForceRestoreSnowFluxes f = force_restore_snow_step(w, t, Real(0.0), Real(1.0e-4),
                                                                 Real(1.0), prm);
        EXPECT_EQ(f.sublimation, Real(0.0));
        EXPECT_EQ(f.melt, Real(0.0));
        EXPECT_EQ(w, Real(-1.0));
        EXPECT_EQ(t, Real(280.0));
    }
    // A melt far below the resolution of a deep pack: W_n cannot take it, so no melt is
    // reported (no water for the soil) and its heat stays in T_s. The potential is
    // C_s (T - T_00) / L_f with T one representable step above T_00.
    {
        const Real t0 = std::nextafter(force_restore_snow_t_melt, Real(300.0));
        const Real w0 = Real(1.0e12);
        Real w = w0, t = t0;
        const ForceRestoreSnowFluxes f = force_restore_snow_step(w, t, Real(0.0), Real(0.0),
                                                                 Real(1.0), prm);
        EXPECT_EQ(w, w0);
        EXPECT_EQ(f.melt, Real(0.0));
        EXPECT_EQ(t, t0);
    }
    // The reported terms are the changes of W_n, exactly.
    {
        Real w = Real(50.0), t = Real(275.0);
        const ForceRestoreSnowFluxes f = force_restore_snow_step(w, t, Real(1.0e-3), Real(1.0e-4),
                                                                 Real(1.0), prm);
        EXPECT_GT(f.melt, Real(0.0));
        EXPECT_NEAR(w - Real(50.0), f.snowfall - f.sublimation - f.melt, Real(100.0) * tol);
        EXPECT_NEAR(t, force_restore_snow_t_melt, Real(1.0e3) * tol);
        EXPECT_GE(t, force_restore_snow_t_melt);
    }
}

TEST(ForceRestoreSnow, LiquidAndFrozenSplitTheAccumulation)
{
    using namespace amrex;
    const BoxArray ba3d(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba3d);
    MultiFab total(ba3d, dm, 1, 0), rain(ba3d, dm, 1, 0), snow(ba3d, dm, 1, 0);
    MultiFab graupel(ba3d, dm, 1, 0), hail(ba3d, dm, 1, 0);
    total.setVal(10.0); rain.setVal(4.0); snow.setVal(3.0); graupel.setVal(2.0); hail.setVal(1.0);
    MultiFab out(ba2d, dm, 1, 0);

    // A scheme whose total holds more than its species (Morrison's also counts the cloud
    // water and ice that fell out; here total 10, rain 4, frozen 5): Liquid is the total
    // minus the frozen part, so that Liquid + Frozen is the Total the hydrology takes
    // without the snow layer.
    SurfacePrecipAccumulationSources src;
    src.total.accum = &total; src.rain.accum = &rain; src.snow.accum = &snow;
    src.graupel.accum = &graupel; src.hail.accum = &hail;
    src.snow.native_to_kg_m2 = Real(2.0);   // a native unit of 2 kg/m^2
    snow.setVal(1.0);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Frozen);
    EXPECT_NEAR(out.min(0), Real(5.0), tol);   // 2 * 1 + 2 + 1
    EXPECT_NEAR(out.max(0), Real(5.0), tol);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Liquid);
    EXPECT_NEAR(out.min(0), Real(5.0), tol);   // 10 - 5, not the rain slot's 4
    EXPECT_NEAR(out.max(0), Real(5.0), tol);
    snow.setVal(3.0);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Total);
    EXPECT_NEAR(out.min(0), Real(10.0), tol);
    surface_precip_accum_k0(out, src, 0);
    EXPECT_NEAR(out.max(0), Real(10.0), tol);

    // No total (SAM): Liquid is the rain slot.
    {
        SurfacePrecipAccumulationSources species = src;
        species.total.accum = nullptr;
        surface_precip_accum_k0(out, species, 0, SurfacePrecipPhase::Liquid);
        EXPECT_EQ(out.min(0), Real(4.0));
        EXPECT_EQ(out.max(0), Real(4.0));
    }

    // A total and frozen species but no rain slot: rain = total - frozen, not below zero.
    src.rain.accum = nullptr; src.snow.native_to_kg_m2 = Real(1.0);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Liquid);
    EXPECT_NEAR(out.min(0), Real(4.0), tol);   // 10 - (3 + 2 + 1)
    total.setVal(5.0);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Liquid);
    EXPECT_EQ(out.min(0), Real(0.0));

    // A warm-rain scheme (Kessler): rain only; nothing frozen, liquid = total.
    SurfacePrecipAccumulationSources warm;
    warm.rain.accum = &rain;
    surface_precip_accum_k0(out, warm, 0, SurfacePrecipPhase::Frozen);
    EXPECT_EQ(out.max(0), Real(0.0));
    surface_precip_accum_k0(out, warm, 0, SurfacePrecipPhase::Liquid);
    EXPECT_EQ(out.min(0), Real(4.0));
    surface_precip_accum_k0(out, warm, 0, SurfacePrecipPhase::Total);
    EXPECT_EQ(out.min(0), Real(4.0));
}

TEST(ForceRestoreSnow, StateTakesSnowfallMeltsAndHandsTheSoilItsWater)
{
    using namespace amrex;
    RadChoice rc;
    rc.seb_snow_enable = true;
    rc.seb_precipitation_source = SEBPrecipitationSource::Microphysics;
    rc.seb_snow_initial_swe_kg_m2 = Real(2.0);
    rc.seb_surface_heat_capacity = Real(2.0e4);
    const BoxArray ba3d(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba3d);
    MultiFab rain(ba3d, dm, 1, 0), snowacc(ba3d, dm, 1, 0);
    SurfacePrecipAccumulationSources src;
    src.rain.accum = &rain; src.snow.accum = &snowacc;

    ForceRestoreSnow snow;
    snow.resize(1);
    snow.define_level(0, rc, ba2d, dm, IntVect(1, 1, 0));
    ASSERT_TRUE(snow.defined(0));
    EXPECT_EQ(snow.swe(0)->min(0), Real(2.0));
    const ForceRestoreSnowParams prm = ForceRestoreSnow::params(rc);
    EXPECT_EQ(prm.c_s, Real(2.0e4));

    ForceRestoreHydrology hydro;
    hydro.resize(1);
    hydro.define_level(0, rc, ba2d, dm, IntVect(0));

    MultiFab t_s(ba2d, dm, 1, 0), q_s(ba2d, dm, 1, 0), q_deep(ba2d, dm, 1, 0);
    MultiFab t_deep(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0), foot(ba2d, dm, 1, 0);
    t_s.setVal(t00 + 1.0); q_s.setVal(0.1); q_deep.setVal(0.1); t_deep.setVal(270.0);
    lh.setVal(0.0); foot.setVal(0.0);
    // One building footprint cell.
    for (MFIter mfi(foot); mfi.isValid(); ++mfi) { foot.array(mfi)(1, 2, 0) = Real(1.0); }

    // One 10 s step: 0.5 kg/m^2 of snow and 0.3 kg/m^2 of rain.
    rain.setVal(1.0); snowacc.setVal(4.0);
    snow.begin_step(0, src, 0); hydro.begin_step(0, src, 0);
    rain.setVal(1.3); snowacc.setVal(4.5);
    snow.end_step(0, src, 0, Real(10.0)); hydro.end_step(0, src, 0, Real(10.0));
    // The hydrology's rain is the liquid part only.
    EXPECT_NEAR(hydro.rain_rate(0)->min(0), Real(0.03), tol);
    EXPECT_NEAR(hydro.rain_rate(0)->max(0), Real(0.03), tol);

    Real ts_sum = 0.0, ts_max = 0.0;
    const auto ss = snow.step(0, t_s, lh, &foot, Real(10.0), prm, ts_sum, ts_max);
    const Real melt = Real(2.0e4) / l_f;   // Eq. 12 with a 1 K excess
    EXPECT_EQ(ss.columns, Long(15));
    EXPECT_NEAR(ss.snowfall, Real(15.0 * 0.5), Real(100.0) * tol);
    EXPECT_NEAR(ss.melt, Real(15.0) * melt, Real(100.0) * tol);
    EXPECT_LE(ss.max_residual, Real(100.0) * tol);
    EXPECT_NEAR(ts_sum, Real(15.0) * t00, Real(1.0e4) * tol);
    EXPECT_NEAR(ts_max, t00, tol);
    for (MFIter mfi(t_s); mfi.isValid(); ++mfi) {
        const auto ts = t_s.const_array(mfi);
        const auto w = snow.swe(0)->const_array(mfi);
        const auto out = snow.soil_water(0)->const_array(mfi);
        const auto tot = snow.totals(0)->const_array(mfi);
        // The footprint keeps its temperature, loses its initial snow and hands the soil
        // nothing.
        EXPECT_EQ(ts(1, 2, 0), t00 + Real(1.0));
        EXPECT_EQ(w(1, 2, 0), Real(0.0));
        EXPECT_EQ(out(1, 2, 0), Real(0.0));
        EXPECT_EQ(tot(1, 2, 0, ForceRestoreSnow::Snowfall), Real(0.0));
        // An open column: snow held at the melting point.
        EXPECT_NEAR(ts(0, 0, 0), t00, Real(1.0e3) * tol);
        EXPECT_GE(ts(0, 0, 0), t00);
        EXPECT_NEAR(w(0, 0, 0), Real(2.5) - melt, Real(100.0) * tol);
        EXPECT_NEAR(out(0, 0, 0), melt, Real(100.0) * tol);
        EXPECT_NEAR(tot(0, 0, 0, ForceRestoreSnow::Snowfall), Real(0.5), Real(100.0) * tol);
        EXPECT_NEAR(tot(0, 0, 0, ForceRestoreSnow::Melt), melt, Real(100.0) * tol);
    }

    // The hydrology takes the rain and the melt water; its Rain total counts the rain only.
    ForceRestoreHydrologyParams hp = ForceRestoreHydrology::params(rc);
    hp.rain = true; hp.d_s = Real(0.1); hp.w_sat = Real(0.45); hp.q_max = Real(1.0);
    hydro.step(0, t_s, q_s, q_deep, t_deep, lh, &foot, Real(10.0), hp, snow.soil_water(0));
    for (MFIter mfi(q_s); mfi.isValid(); ++mfi) {
        const auto q = q_s.const_array(mfi);
        const auto tot = hydro.totals(0)->const_array(mfi);
        // (0.3 + melt) kg/m^2 into a 0.1 m layer: 1 kg/m^2 is 1e-3 m of water.
        EXPECT_NEAR(q(0, 0, 0), Real(0.1) + (Real(0.3) + melt) * Real(1.0e-2), Real(100.0) * tol);
        EXPECT_NEAR(tot(0, 0, 0, ForceRestoreHydrology::Rain), Real(0.3), Real(100.0) * tol);
        EXPECT_EQ(q(1, 2, 0), Real(0.1));
    }

    // The albedo of a column with snow (Eq. 8), on valid cells; the ghost cells keep the
    // snow-free value.
    MultiFab alb(ba2d, dm, 1, IntVect(1, 1, 0));
    alb.setVal(0.2);
    snow.swe(0)->setVal(Real(10.0), 0, 1, IntVect(1, 1, 0));
    snow.apply_albedo(0, alb, prm);
    for (MFIter mfi(alb); mfi.isValid(); ++mfi) {
        const auto a = alb.const_array(mfi);
        const Box vbx = mfi.validbox();
        EXPECT_EQ(vbx.numPts(), Long(16));
        for (int j = vbx.smallEnd(1); j <= vbx.bigEnd(1); ++j) {
            for (int i = vbx.smallEnd(0); i <= vbx.bigEnd(0); ++i) {
                EXPECT_NEAR(a(i, j, 0), Real(0.45), tol);
            }
        }
        EXPECT_EQ(a(-1, 0, 0), Real(0.2));
        EXPECT_EQ(a(4, 4, 0), Real(0.2));
    }
}

TEST(ForceRestoreSnow, SnowfallOfASkippedStepIsCarriedIntoTheNext)
{
    using namespace amrex;
    RadChoice rc;
    rc.seb_snow_enable = true;
    const BoxArray ba3d(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba3d);
    MultiFab snowacc(ba3d, dm, 1, 0);
    SurfacePrecipAccumulationSources src;
    src.snow.accum = &snowacc;
    ForceRestoreSnow snow;
    snow.resize(1);
    snow.define_level(0, rc, ba2d, dm, IntVect(0));
    const ForceRestoreSnowParams prm = ForceRestoreSnow::params(rc);
    MultiFab t_s(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0);
    t_s.setVal(260.0); lh.setVal(0.0);

    // Step 1 (1 s): 2 kg/m^2, skipped by the balance. Step 2 (2 s): 1 kg/m^2 more.
    snowacc.setVal(1.0);
    snow.begin_step(0, src, 0);
    snowacc.setVal(3.0);
    snow.end_step(0, src, 0, Real(1.0));
    snow.begin_step(0, src, 0);
    snowacc.setVal(4.0);
    snow.end_step(0, src, 0, Real(2.0));
    Real ts_sum = 0.0, ts_max = 0.0;
    snow.step(0, t_s, lh, nullptr, Real(2.0), prm, ts_sum, ts_max);
    EXPECT_NEAR(snow.swe(0)->min(0), Real(3.0), tol);
    EXPECT_NEAR(snow.swe(0)->max(0), Real(3.0), tol);
    // Step 3 (1 s), used at once: its own 0.5 kg/m^2.
    snow.begin_step(0, src, 0);
    snowacc.setVal(4.5);
    snow.end_step(0, src, 0, Real(1.0));
    snow.step(0, t_s, lh, nullptr, Real(1.0), prm, ts_sum, ts_max);
    EXPECT_NEAR(snow.swe(0)->max(0), Real(3.5), tol);
    EXPECT_NEAR(snow.totals(0)->min(ForceRestoreSnow::Snowfall), Real(3.5), tol);
    // A step without begin_step (the first after a restart): no snowfall.
    snow.end_step(0, src, 0, Real(1.0));
    snow.step(0, t_s, lh, nullptr, Real(1.0), prm, ts_sum, ts_max);
    EXPECT_NEAR(snow.swe(0)->max(0), Real(3.5), tol);
}
