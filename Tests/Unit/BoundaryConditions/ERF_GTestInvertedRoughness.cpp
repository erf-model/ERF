// The roughness recovered from surface fluxes is bounded above and below.
#include <gtest/gtest.h>

#include "BoundaryConditions/ERF_MOSTUtils.H"

using amrex::Real;

// Motivation: the inversion can return a roughness far above the reference
// height, which makes the 2 m diagnostic fail.
TEST(InvertedRoughness, LargeValuesAreClampedBelowTheReferenceHeight)
{
    EXPECT_EQ(bound_inverted_roughness(Real(1.0e29), Real(10.0)), Real(1.5));
    EXPECT_EQ(bound_inverted_roughness(Real(1.0e29), Real(2.0)), Real(1.0));
    EXPECT_LT(bound_inverted_roughness(Real(5.0), Real(1.0)), Real(1.0));
}

// Motivation: the lower floor is unchanged.
TEST(InvertedRoughness, SmallValuesAreRaisedToTheFloor)
{
    EXPECT_EQ(bound_inverted_roughness(Real(0.0), Real(10.0)), Real(1.0e-7));
    EXPECT_EQ(bound_inverted_roughness(Real(-1.0), Real(10.0)), Real(1.0e-7));
}

// Motivation: values inside the bounds must pass through untouched.
TEST(InvertedRoughness, ValuesInsideTheBoundsAreUnchanged)
{
    EXPECT_EQ(bound_inverted_roughness(Real(0.03), Real(10.0)), Real(0.03));
    EXPECT_EQ(bound_inverted_roughness(Real(1.1), Real(10.0)), Real(1.1));
}
