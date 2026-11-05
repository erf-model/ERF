// The stability coordinate z/L used by the 2 m MOST profile is bounded on the
// stable side, so a very small positive Obukhov length cannot drive the
// diagnostic to an unphysical value.
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "Diagnostics/ERF_NearSurfaceDiagnostics.H"

namespace
{

using amrex::Real;
using near_surface_diagnostics::MostProfile;
using near_surface_diagnostics::MostState;
using near_surface_diagnostics::evaluate_most_factor;
using near_surface_diagnostics::evaluate_most_profile;
using near_surface_diagnostics::most_profile_zeta_max;

MostState state_with (Real roughness, Real obukhov_length)
{
    return MostState{Real(270.0), Real(0.56), Real(0.002), Real(0.0),
                     roughness, obukhov_length};
}

} // namespace

// Motivation: with z0 at its floor and L of a few centimetres, the unbounded
// stable correction gave a 2 m temperature of several hundred K.
TEST(NearSurfaceStability, TinyPositiveObukhovLengthGivesABoundedFactor)
{
    const similarity_funs sfuns{};
    const Real z0 = Real(1.0e-7);
    const Real height = Real(2.0);
    Real factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, Real(0.04)), height, factor, sfuns));
    const Real expected = std::log(height / z0) - sfuns.calc_psi_h(most_profile_zeta_max(sfuns));
    EXPECT_NEAR(factor, expected, Real(1.0e-4));

    MostProfile profile{};
    ASSERT_TRUE(evaluate_most_profile(state_with(z0, Real(0.04)), height, profile, sfuns));
    EXPECT_LT(profile.potential_temperature, Real(400.0));
}

// Motivation: the bound must not change any cell whose z/L is below it.
TEST(NearSurfaceStability, ModerateStabilityIsUnchanged)
{
    const similarity_funs sfuns{};
    const Real z0 = Real(0.05);
    const Real height = Real(2.0);
    const Real obukhov_length = Real(20.0);
    Real factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, obukhov_length), height, factor, sfuns));
    EXPECT_NEAR(factor,
                std::log(height / z0) - sfuns.calc_psi_h(height / obukhov_length),
                Real(1.0e-6));
}

// Motivation: the unstable branch is not bounded here and must stay as it was.
TEST(NearSurfaceStability, UnstableBranchIsUnchanged)
{
    const similarity_funs sfuns{};
    const Real z0 = Real(0.05);
    const Real height = Real(2.0);
    const Real obukhov_length = Real(-5.0);
    Real factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, obukhov_length), height, factor, sfuns));
    EXPECT_NEAR(factor,
                std::log(height / z0) - sfuns.calc_psi_h(height / obukhov_length),
                Real(1.0e-6));
}

// Motivation: the factor must not decrease as stability grows, and must reach
// its bound.
TEST(NearSurfaceStability, FactorIsMonotoneAndSaturatesAtTheBound)
{
    const similarity_funs sfuns{};
    const Real z0 = Real(0.01);
    const Real height = Real(2.0);
    Real previous = -std::numeric_limits<Real>::infinity();
    Real saturated = Real(0.0);
    for (const Real obukhov_length : {Real(1000.0), Real(10.0), Real(2.0), Real(1.0),
                                      Real(0.5), Real(0.04)}) {
        Real factor = Real(0.0);
        ASSERT_TRUE(evaluate_most_factor(state_with(z0, obukhov_length), height, factor, sfuns));
        EXPECT_GE(factor, previous);
        previous = factor;
        saturated = factor;
    }
    Real at_bound = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, Real(2.0)), height, at_bound, sfuns));
    EXPECT_EQ(saturated, at_bound);
}

// Motivation: the QNSE flux kernels run up to z/L = most_zeta_lim (their fold can lie above
// 1), so the profiles must take the same bound and evaluate psi at the z/L of the flux: with
// L = 0.5 m the 2 m factor uses z/L = 4, and only a tiny L (z/L above 100) is held at the
// bound, where ln(1 + b z/L) keeps the factor finite and the 2 m temperature moderate.
TEST(NearSurfaceStability, QnseProfilesTakeTheKernelBound)
{
    similarity_funs qnse{};
    qnse.set_stable_functions(true, Real(2.5), Real(0.2), Real(2.5), Real(0.2));
    const similarity_funs bd{};
    EXPECT_EQ(most_profile_zeta_max(qnse), most_zeta_lim);
    EXPECT_EQ(most_profile_zeta_max(bd), Real(1.0));

    const Real z0 = Real(0.05);
    const Real height = Real(2.0);
    Real factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, Real(0.5)), height, factor, qnse));
    EXPECT_NEAR(factor, std::log(height / z0) - qnse.calc_psi_h(Real(4.0)), Real(1.0e-6));
    Real bd_factor = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(z0, Real(0.5)), height, bd_factor, bd));
    EXPECT_NEAR(bd_factor, std::log(height / z0) - bd.calc_psi_h(Real(1.0)), Real(1.0e-6));

    Real held = Real(0.0);
    ASSERT_TRUE(evaluate_most_factor(state_with(Real(1.0e-7), Real(0.004)), height, held, qnse));
    EXPECT_NEAR(held, std::log(height / Real(1.0e-7)) - qnse.calc_psi_h(most_zeta_lim), Real(1.0e-6));
    MostProfile profile{};
    ASSERT_TRUE(evaluate_most_profile(state_with(Real(1.0e-7), Real(0.004)), height, profile, qnse));
    EXPECT_LT(profile.potential_temperature, Real(400.0));
}
