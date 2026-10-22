#include <gtest/gtest.h>

#include <string>

#include "ERF_MetgridLandUtils.H"

using metgrid_land::interp_in_depth;
using metgrid_land::monthly_interp_weights;

namespace {

// Midpoints of the Noah layers SLM is initialized on, in m
const amrex::Vector<amrex::Real> kNoahMid = {0.05, 0.25, 0.7, 1.5};

constexpr amrex::Real kTol = 1.0e-12;

} // namespace

// Motivation: GFS's soil layers (0-10, 10-40, 40-100, 100-200 cm) are the Noah layers, so
// interpolating them to the Noah midpoints must copy them.
TEST(MetgridLandUtils, GFSLayersAreCopied)
{
    const amrex::Vector<amrex::Real> v = {290.0, 288.0, 285.0, 281.0};
    for (int k = 0; k < 4; ++k) {
        EXPECT_NEAR(interp_in_depth(kNoahMid, v, kNoahMid[k]), v[k], kTol);
    }
}

// Motivation: ERA5's layers (0-7, 7-28, 28-100, 100-289 cm) are not, and must be
// interpolated linearly between their midpoints and held constant beyond the outer ones.
TEST(MetgridLandUtils, ERA5LayersAreInterpolatedInDepth)
{
    const amrex::Vector<amrex::Real> z = {0.035, 0.175, 0.64, 1.945};
    const amrex::Vector<amrex::Real> v = {300.0, 293.0, 286.0, 280.0};

    // 0.05 m lies 0.015 of the 0.14 m between the first two midpoints
    EXPECT_NEAR(interp_in_depth(z, v, 0.05), 300.0 - 7.0 * 0.015 / 0.14, kTol);
    // 0.7 m lies 0.06 of the 1.305 m between the last two
    EXPECT_NEAR(interp_in_depth(z, v, 0.7), 286.0 - 6.0 * 0.06 / 1.305, kTol);
    // Above the shallowest midpoint and below the deepest the end values are kept
    EXPECT_NEAR(interp_in_depth(z, v, 0.01), 300.0, kTol);
    EXPECT_NEAR(interp_in_depth(z, v, 2.5), 280.0, kTol);
    // At a midpoint the value there
    EXPECT_NEAR(interp_in_depth(z, v, 0.64), 286.0, kTol);
}

// Motivation: each monthly value applies on the 15th of its month, as in real.exe, so on
// the 15th at 00 UTC the weight of the next month is zero.
TEST(MetgridLandUtils, MidMonthTakesThatMonth)
{
    int m0 = -1, m1 = -1;
    amrex::Real w1 = -1.0;

    monthly_interp_weights("2023-01-15_00:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 0);
    EXPECT_EQ(m1, 1);
    EXPECT_NEAR(w1, 0.0, kTol);

    monthly_interp_weights("2026-05-15_00:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 4);
    EXPECT_EQ(m1, 5);
    EXPECT_NEAR(w1, 0.0, kTol);

    // Twelve hours later is half a day of the 31 between mid-May and mid-June
    monthly_interp_weights("2026-05-15_12:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 4);
    EXPECT_NEAR(w1, 0.5 / 31.0, kTol);
}

// Motivation: between mid-December and mid-January the interpolation wraps around the end
// of the year, on either side of it.
TEST(MetgridLandUtils, WrapsFromDecemberToJanuary)
{
    int m0 = -1, m1 = -1;
    amrex::Real w1 = -1.0;

    // In the date's own (leap) year mid-December is day 349 and the next mid-January day
    // 14 + 366, 31 days on; January 1 lies 17 days after mid-December
    monthly_interp_weights("2024-01-01_00:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 11);
    EXPECT_EQ(m1, 0);
    EXPECT_NEAR(w1, 17.0 / 31.0, kTol);

    monthly_interp_weights("2023-12-31_12:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 11);
    EXPECT_EQ(m1, 0);
    EXPECT_NEAR(w1, 16.5 / 31.0, kTol);

    monthly_interp_weights("2023-12-15_00:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 11);
    EXPECT_NEAR(w1, 0.0, kTol);
}

// Motivation: February 29 lengthens the interval from mid-February to mid-March.
TEST(MetgridLandUtils, LeapYearFebruary)
{
    int m0 = -1, m1 = -1;
    amrex::Real w1 = -1.0;

    // Mid-February is day 45; mid-March is day 74 in a leap year and 73 otherwise
    monthly_interp_weights("2024-03-01_00:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 1);
    EXPECT_EQ(m1, 2);
    EXPECT_NEAR(w1, 15.0 / 29.0, kTol);

    monthly_interp_weights("2023-03-01_00:00:00", m0, m1, w1);
    EXPECT_EQ(m0, 1);
    EXPECT_EQ(m1, 2);
    EXPECT_NEAR(w1, 14.0 / 28.0, kTol);
}
