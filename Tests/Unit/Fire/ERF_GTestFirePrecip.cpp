#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "ERF_FirePrecip.H"
#include "ERF_FireGrid.H"

/**
 * @file ERF_GTestFirePrecip.cpp
 * @brief The rain per column handed to an engine that asks for it
 *        (FireEngine::wants_surface_precip): the accumulation difference to mm/hr, the combination of the
 *        scheme-native accumulators without double counting, the per-column rate
 *        with its rolling snapshot, and the map onto the fire grid.
 */

/// Round-off tolerance of the build precision: 1e-12 in double, 1e-5 in single.
static constexpr double TOL = (sizeof(amrex::Real) == 8) ? 1e-12 : 1e-5;

TEST(FirePrecip, AccumulationDifferenceToMmPerHour)
{
    // 0.5 kg/m2 (= 0.5 mm of water) over 600 s is 3 mm/hr
    EXPECT_NEAR(fire_precip_rate_mm_hr(2.5, 2.0, 600.0), 3.0, TOL);
    // 1 mm over one hour is 1 mm/hr; over a 0.5 s step it is 7200 mm/hr
    EXPECT_NEAR(fire_precip_rate_mm_hr(1.0, 0.0, 3600.0), 1.0, TOL);
    EXPECT_NEAR(fire_precip_rate_mm_hr(1.0, 0.0, 0.5), 7200.0, 1e-9);
    // no change, a reset accumulator, or a non-positive step: no rain
    EXPECT_EQ(fire_precip_rate_mm_hr(2.0, 2.0, 600.0), 0.0);
    EXPECT_EQ(fire_precip_rate_mm_hr(1.0, 2.0, 600.0), 0.0);
    EXPECT_EQ(fire_precip_rate_mm_hr(2.0, 1.0, 0.0), 0.0);
    EXPECT_EQ(fire_precip_rate_mm_hr(2.0, 1.0, -1.0), 0.0);
    // a NaN accumulation gives no rain rather than a NaN moisture
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    EXPECT_EQ(fire_precip_rate_mm_hr(nan, 1.0, 600.0), 0.0);
}

namespace {

/// A 4 x 2 x 3 level of 100 m cells in x and y, one box, and its 2D slab.
struct Slab
{
    amrex::BoxArray ba3d, ba2d;
    amrex::DistributionMapping dm;
    Slab ()
    {
        ba3d = amrex::BoxArray(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 1, 2)));
        dm   = amrex::DistributionMapping(ba3d);
        ba2d = ba3d;
        ba2d.coarsen(amrex::IntVect(1, 1, 3));
    }
};

/// Fill a 3D accumulator with a value linear in the column and the level.
void fill_linear (amrex::MultiFab& mf, amrex::Real a, amrex::Real b, amrex::Real c)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto arr = mf.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i, j, k) = a * i + b * j + c * k;
        });
    }
}

} // namespace

TEST(FirePrecip, SpeciesSumWithoutDoubleCounting)
{
    Slab s;
    amrex::MultiFab rain(s.ba3d, s.dm, 1, 0), snow(s.ba3d, s.dm, 1, 0), total(s.ba3d, s.dm, 1, 0);
    // rain in mm (factor 1), snow in a native unit worth 0.1 kg/m2, and k = 1 differs
    // from k = 0, so a read at the wrong level shows
    fill_linear(rain, 1.0, 10.0, 100.0);
    fill_linear(snow, 2.0, 20.0, 200.0);
    fill_linear(total, 5.0, 50.0, 500.0);
    amrex::MultiFab accum(s.ba2d, s.dm, 1, 0);

    SurfacePrecipAccumulationSources species;
    species.rain = {&rain, amrex::Real(1.0)};
    species.snow = {&snow, amrex::Real(0.1)};
    fire_surface_precip_accum_k0(accum, species, 0);
    for (amrex::MFIter mfi(accum); mfi.isValid(); ++mfi) {
        auto a = accum.array(mfi);
        EXPECT_NEAR(a(3, 1, 0), (3.0 + 10.0) + 0.1 * (6.0 + 20.0), TOL);
        EXPECT_NEAR(a(0, 0, 0), 0.0, TOL);
    }

    // a scheme with a total (and a frozen subset) reports the total only
    SurfacePrecipAccumulationSources with_total = species;
    with_total.total = {&total, amrex::Real(1.0)};
    fire_surface_precip_accum_k0(accum, with_total, 0);
    for (amrex::MFIter mfi(accum); mfi.isValid(); ++mfi) {
        auto a = accum.array(mfi);
        EXPECT_NEAR(a(3, 1, 0), 15.0 + 50.0, TOL);
    }

    // no source at all: zero everywhere
    fire_surface_precip_accum_k0(accum, SurfacePrecipAccumulationSources{}, 0);
    EXPECT_EQ(accum.max(0), 0.0);
    EXPECT_EQ(accum.min(0), 0.0);
}

TEST(FirePrecip, RatePerColumnAndRollingSnapshot)
{
    Slab s;
    amrex::MultiFab prev(s.ba2d, s.dm, 1, 0), now(s.ba2d, s.dm, 1, 0), rate(s.ba2d, s.dm, 1, 0);
    prev.setVal(0.0);
    // column (i, j) received i + 2 j mm over a 0.5 s step
    fill_linear(now, 1.0, 2.0, 0.0);
    fire_precip_rate_from_accum(rate, prev, now, 0.5, true);
    for (amrex::MFIter mfi(rate); mfi.isValid(); ++mfi) {
        auto r = rate.array(mfi);
        auto p = prev.array(mfi);
        EXPECT_NEAR(r(3, 1, 0), 5.0 * 7200.0, 1e-8);
        EXPECT_EQ(r(0, 0, 0), 0.0);
        EXPECT_NEAR(p(3, 1, 0), 5.0, TOL);      // the snapshot rolled forward
    }
    // the same accumulation on the next step: no rain, snapshot unchanged
    fire_precip_rate_from_accum(rate, prev, now, 0.5, true);
    EXPECT_EQ(rate.max(0), 0.0);
    // an invalid snapshot (restart without one): no rain, but the snapshot is taken
    prev.setVal(-7.0);
    fire_precip_rate_from_accum(rate, prev, now, 0.5, false);
    EXPECT_EQ(rate.max(0), 0.0);
    for (amrex::MFIter mfi(prev); mfi.isValid(); ++mfi) {
        auto p = prev.array(mfi);
        EXPECT_NEAR(p(3, 1, 0), 5.0, TOL);
    }
}

TEST(FirePrecip, ColumnRateReachesEveryFireCellOfTheColumn)
{
    // The 4 x 2 level refined by C = 2: fire cell (i_f, j_f) lies in column (i_f / 2, j_f / 2)
    Slab s;
    amrex::RealBox rb({0.0, 0.0, 0.0}, {400.0, 200.0, 300.0});
    amrex::Geometry geom(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 1, 2)), rb,
                         amrex::CoordSys::cartesian, {1, 1, 0});
    FireGrid fg = create_fire_grid(s.ba3d, s.dm, geom, 2);
    amrex::MultiFab rate(s.ba2d, s.dm, 1, 0), fire_rate(fg.ba, fg.dm, 1, 0);
    fill_linear(rate, 1.0, 10.0, 0.0);
    fill_fire_from_atm_k0(fire_rate, rate, fg);
    for (amrex::MFIter mfi(fire_rate); mfi.isValid(); ++mfi) {
        auto f = fire_rate.array(mfi);
        EXPECT_NEAR(f(7, 3, 0), 3.0 + 10.0, TOL);
        EXPECT_NEAR(f(6, 2, 0), 3.0 + 10.0, TOL);
        EXPECT_NEAR(f(1, 0, 0), 0.0, TOL);
        EXPECT_NEAR(f(2, 1, 0), 1.0, TOL);
    }
}
