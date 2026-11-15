#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

#include <gtest/gtest.h>

#include <AMReX.H>
#include <AMReX_FileSystem.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

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

// amrex::Abort throws instead of ending the process, for the life of the object.
class ScopedThrowOnAbort
{
public:
    ScopedThrowOnAbort () : m_saved(amrex::system::throw_exception)
    {
        amrex::system::throw_exception = true;
    }
    ~ScopedThrowOnAbort () { amrex::system::throw_exception = m_saved; }
    ScopedThrowOnAbort (const ScopedThrowOnAbort&) = delete;
    ScopedThrowOnAbort& operator= (const ScopedThrowOnAbort&) = delete;
private:
    bool m_saved;
};

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
    const ForceRestoreSnowParams p{};
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

    // A total and frozen species but no rain slot: rain = total - frozen.
    src.rain.accum = nullptr; src.snow.native_to_kg_m2 = Real(1.0);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Liquid);
    EXPECT_NEAR(out.min(0), Real(4.0), tol);   // 10 - (3 + 2 + 1)
    // The accumulation is not bounded at zero: a bound there would change the difference
    // of two accumulations, which is what a step takes (see the test below).
    total.setVal(5.0);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Liquid);
    EXPECT_NEAR(out.min(0), Real(-1.0), tol);
    EXPECT_NEAR(out.max(0), Real(-1.0), tol);

    // Both parts from one pass: the same values, in components 0 and 1.
    MultiFab both(ba2d, dm, 2, 0);
    surface_precip_liquid_frozen_accum_k0(both, src, 0);
    EXPECT_NEAR(both.min(0), Real(-1.0), tol);
    EXPECT_NEAR(both.max(0), Real(-1.0), tol);
    EXPECT_NEAR(both.min(1), Real(6.0), tol);
    EXPECT_NEAR(both.max(1), Real(6.0), tol);

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
    ForceRestorePrecipRate precip;
    precip.resize(1);
    precip.define_level(0, rc, ba2d, dm);
    ASSERT_TRUE(precip.split());
    EXPECT_EQ(precip.rate(0)->nComp(), 2);

    MultiFab t_s(ba2d, dm, 1, 0), q_s(ba2d, dm, 1, 0), q_deep(ba2d, dm, 1, 0);
    MultiFab t_deep(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0), foot(ba2d, dm, 1, 0);
    t_s.setVal(t00 + 1.0); q_s.setVal(0.1); q_deep.setVal(0.1); t_deep.setVal(270.0);
    lh.setVal(0.0); foot.setVal(0.0);
    // One building footprint cell.
    for (MFIter mfi(foot); mfi.isValid(); ++mfi) { foot.array(mfi)(1, 2, 0) = Real(1.0); }

    // One 10 s step: 0.5 kg/m^2 of snow and 0.3 kg/m^2 of rain.
    rain.setVal(1.0); snowacc.setVal(4.0);
    precip.begin_step(0, src, 0);
    rain.setVal(1.3); snowacc.setVal(4.5);
    precip.end_step(0, src, 0, Real(10.0));
    // The hydrology's rain is the liquid part only; the snow layer's is the frozen part.
    EXPECT_NEAR(precip.rate(0)->min(ForceRestorePrecipRate::Rain), Real(0.03), tol);
    EXPECT_NEAR(precip.rate(0)->max(ForceRestorePrecipRate::Rain), Real(0.03), tol);
    EXPECT_NEAR(precip.rate(0)->min(ForceRestorePrecipRate::Frozen), Real(0.05), tol);
    EXPECT_NEAR(precip.rate(0)->max(ForceRestorePrecipRate::Frozen), Real(0.05), tol);

    // Before the first step the footprint still holds the initial snow, and the albedo
    // must not see it: no snow lies on a building.
    {
        MultiFab alb0(ba2d, dm, 1, 0);
        alb0.setVal(0.2);
        snow.apply_albedo(0, alb0, prm, &foot);
        for (MFIter mfi(alb0); mfi.isValid(); ++mfi) {
            const auto a = alb0.const_array(mfi);
            EXPECT_EQ(snow.swe(0)->const_array(mfi)(1, 2, 0), Real(2.0));
            EXPECT_EQ(a(1, 2, 0), Real(0.2));
            // An open column: p_ng = 2 / 12.
            EXPECT_NEAR(a(0, 0, 0), Real(0.2) + (Real(0.7) - Real(0.2)) / Real(6.0), tol);
        }
    }

    Real ts_sum = 0.0, ts_max = 0.0;
    const auto ss = snow.step(0, t_s, lh, &foot, Real(10.0), prm, precip.frozen(0, Real(10.0)),
                              ts_sum, ts_max);
    EXPECT_EQ(ss.bad_columns, Long(0));
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
    // The soil takes the snow layer's water only with prm.rain: without it the step stops
    // instead of dropping the water.
    {
        ForceRestoreHydrologyParams dry = hp;
        dry.rain = false;
        ScopedThrowOnAbort throw_on_abort;
        EXPECT_THROW(hydro.step(0, t_s, q_s, q_deep, t_deep, lh, &foot, Real(10.0), dry,
                                precip.rain(0, Real(10.0)), snow.soil_water(0)),
                     amrex::RuntimeError);
    }
    hydro.step(0, t_s, q_s, q_deep, t_deep, lh, &foot, Real(10.0), hp, precip.rain(0, Real(10.0)),
               snow.soil_water(0));
    precip.mark_used(0);
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

// Land where i < 2, with a building's column (2, which is land) at (0,0); water where
// i >= 2. Outside the TEST body: a device lambda cannot be launched from a private member.
namespace {
void fill_water_east_mask (amrex::iMultiFab& mask)
{
    for (amrex::MFIter mfi(mask); mfi.isValid(); ++mfi) {
        const auto m = mask.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            m(i, j, k) = (i >= 2) ? 0 : ((i == 0 && j == 0) ? 2 : 1);
        });
    }
    amrex::Gpu::streamSynchronize();
}
}

// With a land mask the snow layer is kept on land columns only. A water column (mask 0)
// keeps its temperature, loses its initial snow at the first step, hands the soil nothing,
// books no snowfall, is not counted, and the albedo there is the snow-free one before and
// after that step. 4 x 4 columns, water where i >= 2; a building's column (mask 2) is land.
// 2 kg/m^2 of snow at 1 K above the melting point with C_s = 2e4 J/m^2/K: Eq. 12 melts
// C_s / L_f, and 0.5 kg/m^2 of snow falls in the 10 s step.
TEST(ForceRestoreSnow, WaterColumnsAreLeftOut)
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
    const ForceRestoreSnowParams prm = ForceRestoreSnow::params(rc);
    ForceRestorePrecipRate precip;
    precip.resize(1);
    precip.define_level(0, rc, ba2d, dm);

    MultiFab t_s(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0);
    t_s.setVal(t00 + 1.0); lh.setVal(0.0);
    iMultiFab mask(ba2d, dm, 1, 0);
    fill_water_east_mask(mask);
    const Box land(IntVect(0, 0, 0), IntVect(1, 3, 0)), water(IntVect(2, 0, 0), IntVect(3, 3, 0));

    rain.setVal(0.0); snowacc.setVal(4.0);
    precip.begin_step(0, src, 0);
    snowacc.setVal(4.5);
    precip.end_step(0, src, 0, Real(10.0));

    // Before the first step the water columns still hold the initial snow; the albedo does
    // not see it. On land: p_ng = 2 / 12.
    const Real albedo_snow = Real(0.2) + (Real(0.7) - Real(0.2)) / Real(6.0);
    {
        MultiFab alb(ba2d, dm, 1, 0);
        alb.setVal(0.2);
        snow.apply_albedo(0, alb, prm, nullptr, &mask);
        EXPECT_EQ(snow.swe(0)->min(water, 0), Real(2.0));
        EXPECT_EQ(alb.min(water, 0), Real(0.2));
        EXPECT_EQ(alb.max(water, 0), Real(0.2));
        EXPECT_NEAR(alb.min(land, 0), albedo_snow, tol);
        EXPECT_NEAR(alb.max(land, 0), albedo_snow, tol);
    }

    Real ts_sum = 0.0, ts_max = 0.0;
    const auto ss = snow.step(0, t_s, lh, nullptr, Real(10.0), prm, precip.frozen(0, Real(10.0)),
                              ts_sum, ts_max, &mask);
    precip.mark_used(0);
    const Real melt = Real(2.0e4) / l_f;
    EXPECT_EQ(ss.columns, Long(8));
    EXPECT_EQ(ss.bad_columns, Long(0));
    EXPECT_NEAR(ss.snowfall, Real(8.0 * 0.5), Real(100.0) * tol);
    EXPECT_NEAR(ss.melt, Real(8.0) * melt, Real(100.0) * tol);
    EXPECT_NEAR(ts_sum, Real(8.0) * t00, Real(1.0e4) * tol);
    // Water: the temperature as it was, no snow, nothing for the soil, nothing booked.
    EXPECT_EQ(t_s.min(water, 0), t00 + Real(1.0));
    EXPECT_EQ(t_s.max(water, 0), t00 + Real(1.0));
    EXPECT_EQ(snow.swe(0)->max(water, 0), Real(0.0));
    EXPECT_EQ(snow.soil_water(0)->max(water, 0), Real(0.0));
    EXPECT_EQ(snow.totals(0)->max(water, ForceRestoreSnow::Snowfall), Real(0.0));
    EXPECT_EQ(snow.totals(0)->max(water, ForceRestoreSnow::Melt), Real(0.0));
    // Land, the building's column with it: snow held at the melting point.
    EXPECT_NEAR(t_s.max(land, 0), t00, Real(1.0e3) * tol);
    EXPECT_NEAR(snow.swe(0)->min(land, 0), Real(2.5) - melt, Real(100.0) * tol);
    EXPECT_NEAR(snow.swe(0)->max(land, 0), Real(2.5) - melt, Real(100.0) * tol);
    EXPECT_NEAR(snow.totals(0)->min(land, ForceRestoreSnow::Snowfall), Real(0.5), Real(100.0) * tol);

    // After the step the albedo over water is still the snow-free one.
    MultiFab alb(ba2d, dm, 1, 0);
    alb.setVal(0.2);
    snow.apply_albedo(0, alb, prm, nullptr, &mask);
    EXPECT_EQ(alb.max(water, 0), Real(0.2));
    EXPECT_GT(alb.min(land, 0), Real(0.2) + Real(0.01));
}

TEST(ForceRestoreSnow, SnowfallOfASkippedStepIsCarriedIntoTheNext)
{
    using namespace amrex;
    RadChoice rc;
    rc.seb_snow_enable = true;
    rc.seb_precipitation_source = SEBPrecipitationSource::Microphysics;
    const BoxArray ba3d(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba3d);
    MultiFab snowacc(ba3d, dm, 1, 0);
    SurfacePrecipAccumulationSources src;
    src.snow.accum = &snowacc;
    ForceRestoreSnow snow;
    snow.resize(1);
    snow.define_level(0, rc, ba2d, dm, IntVect(0));
    ForceRestorePrecipRate precip;
    precip.resize(1);
    precip.define_level(0, rc, ba2d, dm);
    const ForceRestoreSnowParams prm = ForceRestoreSnow::params(rc);
    MultiFab t_s(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0);
    t_s.setVal(260.0); lh.setVal(0.0);

    // Step 1 (1 s): 2 kg/m^2, skipped by the balance. Step 2 (2 s): 1 kg/m^2 more.
    snowacc.setVal(1.0);
    precip.begin_step(0, src, 0);
    snowacc.setVal(3.0);
    precip.end_step(0, src, 0, Real(1.0));
    precip.begin_step(0, src, 0);
    snowacc.setVal(4.0);
    precip.end_step(0, src, 0, Real(2.0));
    Real ts_sum = 0.0, ts_max = 0.0;
    snow.step(0, t_s, lh, nullptr, Real(2.0), prm, precip.frozen(0, Real(2.0)), ts_sum, ts_max);
    precip.mark_used(0);
    EXPECT_NEAR(snow.swe(0)->min(0), Real(3.0), tol);
    EXPECT_NEAR(snow.swe(0)->max(0), Real(3.0), tol);
    // Step 3 (1 s), used at once: its own 0.5 kg/m^2.
    precip.begin_step(0, src, 0);
    snowacc.setVal(4.5);
    precip.end_step(0, src, 0, Real(1.0));
    snow.step(0, t_s, lh, nullptr, Real(1.0), prm, precip.frozen(0, Real(1.0)), ts_sum, ts_max);
    precip.mark_used(0);
    EXPECT_NEAR(snow.swe(0)->max(0), Real(3.5), tol);
    EXPECT_NEAR(snow.totals(0)->min(ForceRestoreSnow::Snowfall), Real(3.5), tol);
    // A step without begin_step (the first after a restart): no snowfall.
    precip.end_step(0, src, 0, Real(1.0));
    snow.step(0, t_s, lh, nullptr, Real(1.0), prm, precip.frozen(0, Real(1.0)), ts_sum, ts_max);
    precip.mark_used(0);
    EXPECT_NEAR(snow.swe(0)->max(0), Real(3.5), tol);
}

TEST(ForceRestoreSnow, ARateIsNotAppliedAgainWithoutANewEndStep)
{
    // A level whose radiation is interpolated has neither begin_step nor end_step, but its
    // balance and snow step run. The rate of the last end_step is then still in the field,
    // and must not add its snowfall to W_n and the Snowfall total a second time.
    using namespace amrex;
    RadChoice rc;
    rc.seb_snow_enable = true;
    rc.seb_precipitation_source = SEBPrecipitationSource::Microphysics;
    const BoxArray ba3d(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba3d);
    MultiFab snowacc(ba3d, dm, 1, 0);
    SurfacePrecipAccumulationSources src;
    src.snow.accum = &snowacc;
    ForceRestoreSnow snow;
    snow.resize(1);
    snow.define_level(0, rc, ba2d, dm, IntVect(0));
    ForceRestorePrecipRate precip;
    precip.resize(1);
    precip.define_level(0, rc, ba2d, dm);
    const ForceRestoreSnowParams prm = ForceRestoreSnow::params(rc);
    MultiFab t_s(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0);
    t_s.setVal(260.0); lh.setVal(0.0);

    snowacc.setVal(1.0);
    precip.begin_step(0, src, 0);
    snowacc.setVal(1.5);
    precip.end_step(0, src, 0, Real(1.0));
    Real ts_sum = 0.0, ts_max = 0.0;
    snow.step(0, t_s, lh, nullptr, Real(1.0), prm, precip.frozen(0, Real(1.0)), ts_sum, ts_max);
    precip.mark_used(0);
    EXPECT_NEAR(snow.swe(0)->max(0), Real(0.5), tol);

    // Three more steps with no end_step between them.
    EXPECT_NEAR(precip.rate(0)->max(ForceRestorePrecipRate::Frozen), Real(0.5), tol);
    for (int n = 0; n < 3; ++n) {
        EXPECT_EQ(precip.frozen(0, Real(1.0)).scale, Real(0.0));
        const auto ss = snow.step(0, t_s, lh, nullptr, Real(1.0), prm, precip.frozen(0, Real(1.0)),
                                  ts_sum, ts_max);
        precip.mark_used(0);
        EXPECT_EQ(ss.snowfall, Real(0.0));
    }
    EXPECT_NEAR(snow.swe(0)->max(0), Real(0.5), tol);
    EXPECT_NEAR(snow.totals(0)->max(ForceRestoreSnow::Snowfall), Real(0.5), tol);
}

TEST(ForceRestoreSnow, RainOfAStepIsTheChangeOfTheTotalMinusTheFrozenPart)
{
    // Early in a run the frozen accumulators can exceed the total by round-off (a unit
    // conversion, single precision). The rain of a step is the change of (total - frozen)
    // over the step: with the accumulation bounded at zero, a start value below zero would
    // be raised to zero and the step would lose that much rain.
    using namespace amrex;
    RadChoice rc;
    rc.seb_snow_enable = true;
    rc.seb_precipitation_source = SEBPrecipitationSource::Microphysics;
    const BoxArray ba3d(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba3d);
    MultiFab total(ba3d, dm, 1, 0), snowacc(ba3d, dm, 1, 0);
    SurfacePrecipAccumulationSources src;
    src.total.accum = &total; src.snow.accum = &snowacc;
    ForceRestorePrecipRate precip;
    precip.resize(1);
    precip.define_level(0, rc, ba2d, dm);

    // Start: total 5, frozen 5.25 (liquid -0.25). End: total 6.5, frozen 5.75 (liquid 0.75).
    // Over 2 s: 1.5 kg/m^2 in all, 0.5 frozen and 1.0 rain.
    total.setVal(5.0); snowacc.setVal(5.25);
    precip.begin_step(0, src, 0);
    total.setVal(6.5); snowacc.setVal(5.75);
    precip.end_step(0, src, 0, Real(2.0));
    EXPECT_NEAR(precip.rate(0)->min(ForceRestorePrecipRate::Rain), Real(0.5), tol);
    EXPECT_NEAR(precip.rate(0)->max(ForceRestorePrecipRate::Rain), Real(0.5), tol);
    EXPECT_NEAR(precip.rate(0)->min(ForceRestorePrecipRate::Frozen), Real(0.25), tol);
    EXPECT_NEAR(precip.rate(0)->max(ForceRestorePrecipRate::Frozen), Real(0.25), tol);
    precip.mark_used(0);

    // A step in which the liquid part falls (frozen grows faster than the total): no rain,
    // never a negative rate.
    precip.begin_step(0, src, 0);
    total.setVal(6.75); snowacc.setVal(6.5);
    precip.end_step(0, src, 0, Real(1.0));
    EXPECT_EQ(precip.rate(0)->min(ForceRestorePrecipRate::Rain), Real(0.0));
    EXPECT_NEAR(precip.rate(0)->max(ForceRestorePrecipRate::Frozen), Real(0.75), tol);
}

TEST(ForceRestoreSnow, HeatAboveTheUpperBoundMeltsSnowBeforeTheBoundApplies)
{
    // The balance bounds T_s at erf.radiation.seb_prognostic_t_max_k. Under snow the bound
    // waits for the melt: the estimate comes in unbounded, the heat above the melting point
    // melts snow (Eq. 12), and the bound applies to what is left.
    ForceRestoreSnowParams p = base();
    p.t_max = Real(274.0);
    // A deep pack: all the heat of a 10 K excess goes to melt, not only the 0.84 K below
    // the bound.
    Real w = Real(100.0), t = t00 + Real(10.0);
    ForceRestoreSnowFluxes f = force_restore_snow_step(w, t, Real(0.0), Real(0.0), Real(1.0), p);
    EXPECT_NEAR(f.melt, p.c_s * Real(10.0) / l_f, Real(100.0) * tol);
    EXPECT_NEAR(t, t00, Real(1.0e3) * tol);
    // A thin pack melts in full; the rest of the heat would leave T_s above the bound,
    // which then applies.
    w = Real(0.01); t = t00 + Real(10.0);
    f = force_restore_snow_step(w, t, Real(0.0), Real(0.0), Real(1.0), p);
    EXPECT_NEAR(f.melt, Real(0.01), tol);
    EXPECT_EQ(w, Real(0.0));
    EXPECT_EQ(t, Real(274.0));
    // No bound by default: the kernel's answers without one are unchanged.
    EXPECT_EQ(ForceRestoreSnowParams().t_max, std::numeric_limits<Real>::max());
    // The state class takes the bound from the balance's input.
    RadChoice rc;
    rc.seb_prognostic_t_max_k = Real(300.0);
    EXPECT_EQ(ForceRestoreSnow::params(rc).t_max, Real(300.0));
}

TEST(ForceRestoreSnow, StepCountsTheColumnsItLeavesInARefusedState)
{
    using namespace amrex;
    RadChoice rc;
    rc.seb_snow_enable = true;
    rc.seb_snow_initial_swe_kg_m2 = Real(1.0);
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba2d);
    ForceRestoreSnow snow;
    snow.resize(1);
    snow.define_level(0, rc, ba2d, dm, IntVect(0));
    const ForceRestoreSnowParams prm = ForceRestoreSnow::params(rc);
    MultiFab t_s(ba2d, dm, 1, 0), lh(ba2d, dm, 1, 0), foot(ba2d, dm, 1, 0);
    t_s.setVal(260.0); lh.setVal(0.0); foot.setVal(0.0);
    Real ts_sum = 0.0, ts_max = 0.0;
    auto ss = snow.step(0, t_s, lh, &foot, Real(1.0), prm, ForceRestorePrecipInput{}, ts_sum, ts_max);
    EXPECT_EQ(ss.columns, Long(16));
    EXPECT_EQ(ss.bad_columns, Long(0));

    // A non-finite W_n, a negative W_n and a non-finite T_s: three columns. The same
    // values on a footprint are not counted: the step does not use them.
    for (MFIter mfi(t_s); mfi.isValid(); ++mfi) {
        snow.swe(0)->array(mfi)(0, 0, 0) = std::numeric_limits<Real>::quiet_NaN();
        snow.swe(0)->array(mfi)(1, 0, 0) = Real(-1.0);
        t_s.array(mfi)(2, 0, 0) = std::numeric_limits<Real>::infinity();
        t_s.array(mfi)(3, 3, 0) = std::numeric_limits<Real>::quiet_NaN();
        foot.array(mfi)(3, 3, 0) = Real(1.0);
    }
    ss = snow.step(0, t_s, lh, &foot, Real(1.0), prm, ForceRestorePrecipInput{}, ts_sum, ts_max);
    EXPECT_EQ(ss.columns, Long(15));
    EXPECT_EQ(ss.bad_columns, Long(3));
}

TEST(ForceRestoreSnow, RestartRefusesHalfOfTheSnowFiles)
{
    using namespace amrex;
    RadChoice rc;
    rc.seb_snow_enable = true;
    rc.seb_snow_initial_swe_kg_m2 = Real(1.0);
    const BoxArray ba2d(Box(IntVect(0, 0, 0), IntVect(3, 3, 0)));
    const DistributionMapping dm(ba2d);
    const Periodicity period = Periodicity::NonPeriodic();
    const std::string chk = "unit_force_restore_snow_chk";
    amrex::UtilCreateCleanDirectory(chk, true);
    amrex::UtilCreateCleanDirectory(chk + "/Level_0", true);
    amrex::UtilCreateCleanDirectory(chk + "/Level_1", true);

    ForceRestoreSnow written;
    written.resize(2);
    written.define_level(0, rc, ba2d, dm, IntVect(0));
    written.swe(0)->setVal(7.0);
    written.totals(0)->setVal(3.0);

    // Neither file: a checkpoint without the snow layer. The snow starts from the input,
    // and the note is printed once, not per level.
    {
        ForceRestoreSnow snow;
        snow.resize(2);
        snow.define_level(0, rc, ba2d, dm, IntVect(0));
        snow.define_level(1, rc, ba2d, dm, IntVect(0));
        testing::internal::CaptureStdout();
        snow.read_checkpoint(0, chk, period);
        snow.read_checkpoint(1, chk, period);
        const std::string report = testing::internal::GetCapturedStdout();
        const std::string note = "NOTE: the checkpoint has no force-restore snow";
        const auto first = report.find(note);
        ASSERT_NE(first, std::string::npos) << report;
        EXPECT_EQ(report.find(note, first + 1), std::string::npos) << report;
        EXPECT_EQ(snow.swe(0)->min(0), Real(1.0));
    }

    // Both files: the snow of the checkpoint.
    written.write_checkpoint(0, chk);
    {
        ForceRestoreSnow snow;
        snow.resize(1);
        snow.define_level(0, rc, ba2d, dm, IntVect(0));
        snow.read_checkpoint(0, chk, period);
        EXPECT_EQ(snow.swe(0)->min(0), Real(7.0));
        EXPECT_EQ(snow.swe(0)->max(0), Real(7.0));
        EXPECT_EQ(snow.totals(0)->min(ForceRestoreSnow::Melt), Real(3.0));
    }

    // One file without the other: the run stops and names the missing file, instead of
    // starting the snow from the input.
    for (const char* missing : {"TwoStream_SnowTotals", "TwoStream_SnowSWE"}) {
        SCOPED_TRACE(missing);
        amrex::UtilCreateCleanDirectory(chk + "/Level_0", true);
        written.write_checkpoint(0, chk);
        const std::string name = MultiFabFileFullPrefix(0, chk, "Level_", missing);
        ASSERT_EQ(std::remove((name + "_H").c_str()), 0);
        ForceRestoreSnow snow;
        snow.resize(1);
        snow.define_level(0, rc, ba2d, dm, IntVect(0));
        ScopedThrowOnAbort throw_on_abort;
        try {
            snow.read_checkpoint(0, chk, period);
            ADD_FAILURE() << "read_checkpoint did not stop";
        } catch (const amrex::RuntimeError& e) {
            const std::string what = e.what();
            EXPECT_NE(what.find(std::string("but not ") + name), std::string::npos) << what;
        }
    }
    amrex::FileSystem::RemoveAll(chk);
}
