#include <gtest/gtest.h>

#include <cmath>
#include <string>

// ERF_EpochTime.H uses the AMReX qualifier macros but does not include them itself,
// so pull them in before it here
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include "ERF_EpochTime.H"

namespace {

const std::string kFormat = "%Y-%m-%d_%H:%M:%S";

} // namespace

// Motivation: the whole seconds and the printed fraction must come from the same
// rounded instant.  Truncating the seconds and then printing the leftover fraction
// to six places emits "1.000000" next to the un-incremented second, which reads as
// a full second earlier than the true time.
TEST(EpochTime, FractionThatRoundsUpCarriesIntoTheSeconds)
{
    // This is what erf.fixed_dt = 0.1 accumulates to after ten steps
    double t = 0.0;
    for (int n = 0; n < 10; ++n) { t += 0.1; }
    ASSERT_LT(t, 1.0);
    ASSERT_GT(t, 0.9999995);

    EXPECT_EQ(getTimestamp(t, kFormat), "1970-01-01_00:00:01.000000");

    // The same one ulp below a whole second at a realistic epoch
    const double e = std::nextafter(1262368800.0 + 1.0, 1262368800.0);
    EXPECT_EQ(getTimestamp(e, kFormat), "2010-01-01_18:00:01.000000");
}

// Motivation: the carry must not perturb ordinary times.
TEST(EpochTime, RepresentableFractionsAreUnchanged)
{
    EXPECT_EQ(getTimestamp(0.5,  kFormat), "1970-01-01_00:00:00.500000");
    EXPECT_EQ(getTimestamp(1.25, kFormat), "1970-01-01_00:00:01.250000");
    EXPECT_EQ(getTimestamp(1262368800.123456, kFormat), "2010-01-01_18:00:00.123456");
    EXPECT_EQ(getTimestamp(0.0,  kFormat), "1970-01-01_00:00:00.000000");
}

// Motivation: plotfile names are built with add_long_frac = false, which reports the
// second the time falls in; that must keep truncating rather than rounding.
TEST(EpochTime, WholeSecondsOnlyStillTruncates)
{
    EXPECT_EQ(getTimestamp(10.6, kFormat, false), "1970-01-01_00:00:10");
    EXPECT_EQ(getTimestamp(10.4, kFormat, false), "1970-01-01_00:00:10");
}

// Motivation: times before the epoch used to report the magnitude of the fraction
// alongside the truncated-toward-zero second, giving a timestamp of the wrong sign.
TEST(EpochTime, TimesBeforeTheEpochCarryDownward)
{
    EXPECT_EQ(getTimestamp(-0.5,  kFormat), "1969-12-31_23:59:59.500000");
    EXPECT_EQ(getTimestamp(-1.25, kFormat), "1969-12-31_23:59:58.750000");
}

// Motivation: getEpochTime and getTimestamp should round-trip a whole second.
TEST(EpochTime, RoundTripsAWholeSecond)
{
    const std::time_t epoch = getEpochTime("2010-01-01_18:00:00", kFormat);
    ASSERT_NE(epoch, static_cast<std::time_t>(-1));
    EXPECT_EQ(getTimestamp(static_cast<double>(epoch), kFormat, false), "2010-01-01_18:00:00");
}
