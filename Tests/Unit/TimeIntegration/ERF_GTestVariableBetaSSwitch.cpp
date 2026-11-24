// Unit tests for the automatic switch to erf.use_variable_beta_s (ERF_VariableBetaSSwitch.H):
// the constant beta_s must be replaced exactly when it is below the large-Courant stability
// limit min(sigma^2, 1), with the limit capped at 1 on very steep terrain (erf-model/ERF#5082).

#include <ERF_VariableBetaSSwitch.H>

#include <gtest/gtest.h>

using amrex::Real;

TEST(VariableBetaSSwitch, BelowSlopeLimit)
{
    // flat or gentle terrain: the default beta_s = 0.1 is stable
    EXPECT_FALSE(VariableBetaSSwitchNeeded(Real(0.1), Real(0.0)));
    EXPECT_FALSE(VariableBetaSSwitchNeeded(Real(0.1), Real(0.05)));
    EXPECT_FALSE(VariableBetaSSwitchNeeded(Real(0.1), Real(0.1)));    // at the limit
    // steeper terrain: the default beta_s is unstable
    EXPECT_TRUE (VariableBetaSSwitchNeeded(Real(0.1), Real(0.18)));
    EXPECT_TRUE (VariableBetaSSwitchNeeded(Real(0.1), Real(0.42)));
    // a larger beta_s covers the steeper slope
    EXPECT_FALSE(VariableBetaSSwitchNeeded(Real(0.5), Real(0.42)));
}

TEST(VariableBetaSSwitch, LimitCappedAtOne)
{
    // sigma^2 > 1: the limit is 1, not sigma^2
    EXPECT_TRUE (VariableBetaSSwitchNeeded(Real(0.95), Real(3.0)));
    EXPECT_FALSE(VariableBetaSSwitchNeeded(Real(1.0),  Real(3.0)));
}
