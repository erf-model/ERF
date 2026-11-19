// Unit tests for the CFL limit of a fixed time step (ERF_FixedDtCFLLimit.H): erf.fixed_dt must be
// lowered exactly when its Courant number is above erf.cfl, and then to the Courant number
// min(erf.cfl, 0.8) (erf-model/ERF#5238).

#include <ERF_FixedDtCFLLimit.H>

#include <gtest/gtest.h>

TEST(FixedDtCFLLimit, WithinLimit)
{
    // Courant number 0.5 and 0.8 (at the limit) with erf.cfl = 0.8
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.25, 2.0, 0.8), 1.0);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.4,  2.0, 0.8), 1.0);
    // no fixed dt, no constraining rate, no CFL limit
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(-1.0, 2.0, 0.8), 1.0);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.5,  0.0, 0.8), 1.0);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.5,  2.0, 0.0), 1.0);
}

TEST(FixedDtCFLLimit, LoweredToCappedCFL)
{
    // erf.cfl = 1 is capped at 0.8: Courant number 1.2 -> 0.8
    EXPECT_DOUBLE_EQ(0.6 * FixedDtCFLScale(0.6, 2.0, 1.0) * 2.0, 0.8);
    // the default erf.cfl = 0.8 is used as given: Courant number 1 -> 0.8
    EXPECT_DOUBLE_EQ(0.5 * FixedDtCFLScale(0.5, 2.0, 0.8) * 2.0, 0.8);
    // a smaller erf.cfl is used as given: Courant number 1 -> 0.5
    EXPECT_DOUBLE_EQ(0.5 * FixedDtCFLScale(0.5, 2.0, 0.5) * 2.0, 0.5);
    // a step set to the limit by the scaling is not scaled again
    const double dt = 0.5 * FixedDtCFLScale(0.5, 2.0, 0.5);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(dt, 2.0, 0.5), 1.0);
}
