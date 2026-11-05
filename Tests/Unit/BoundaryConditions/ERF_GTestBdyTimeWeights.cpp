#include <gtest/gtest.h>

#include <cmath>

#include "ERF_Utils.H"

namespace {

// 2024-08-05 00:00:00 UTC, a six-hourly wrfbdy with 21 slices
constexpr double start_bdy_time = 1722816000.0;
constexpr double dT             = 21600.0;
constexpr double data_span      = 20 * dT;

} // namespace

// Motivation: the weight of the next boundary slice must follow the elapsed time to roundoff.
// Forming the time as (elapsed + start_time) - start_bdy_time rounds the elapsed time to the
// ~2.4e-7 s spacing of doubles at epoch times, so the substep times of a refined level, such
// as 10.8333... s, mapped to a weight with only ~7 good digits.  The weight is an amrex::Real,
// so the tolerance is a multiple of the working precision: in double precision that is ~1e-15,
// well below the ~1e-11 error of the old form; in single precision the old form's error is
// below the roundoff of the weight itself, and the test only checks the weight to that roundoff.
TEST(BdyTimeWeights, WeightFollowsElapsedTimeToRoundoff)
{
    const double tol = static_cast<double>(amrex::Real(4) * real_eps);
    const double start_time = start_bdy_time + 3600.0; // run starts an hour into the data
    for (const double elapsed : {10.0, 10.0 + 5.0/6.0, 11.0 + 2.0/3.0, 13.0 + 1.0/3.0, 1.0e-3}) {
        const double tsb = elapsed + 3600.0;
        const BdyTimeWeights w = bdy_time_weights(elapsed + (start_time - start_bdy_time),
                                                  data_span, dT);
        EXPECT_EQ(w.n_lo, 0);
        EXPECT_EQ(w.n_hi, 1);
        EXPECT_NEAR(w.alpha, tsb / dT, tol) << "elapsed = " << elapsed;
    }
}

// Motivation: the same elapsed time must give the same weight whatever the date, as long as
// the run starts at the same point in the boundary data.
TEST(BdyTimeWeights, WeightDoesNotDependOnTheDate)
{
    const double elapsed = 11.0 + 2.0/3.0;
    const BdyTimeWeights a = bdy_time_weights(elapsed + ((start_bdy_time + 60.0) - start_bdy_time),
                                              data_span, dT);
    const double later_bdy = start_bdy_time + 365.0 * 86400.0;
    const BdyTimeWeights b = bdy_time_weights(elapsed + ((later_bdy + 60.0) - later_bdy),
                                              data_span, dT);
    EXPECT_EQ(a.n_lo, b.n_lo);
    EXPECT_EQ(a.alpha, b.alpha);
}

// Motivation: at a slice time the later slice must be picked with zero weight on the next one,
// and past the last slice the last one must be used alone.
TEST(BdyTimeWeights, SliceTimesAndTheEndOfTheData)
{
    const BdyTimeWeights at_slice = bdy_time_weights(3 * dT, data_span, dT);
    EXPECT_EQ(at_slice.n_lo, 3);
    EXPECT_EQ(at_slice.n_hi, 4);
    EXPECT_EQ(at_slice.alpha, 0.0);

    for (const double tsb : {data_span, data_span + 1.0}) {
        const BdyTimeWeights w = bdy_time_weights(tsb, data_span, dT);
        EXPECT_EQ(w.n_lo, 20);
        EXPECT_EQ(w.n_hi, 20);
        EXPECT_EQ(w.alpha, 0.0);
    }
}
