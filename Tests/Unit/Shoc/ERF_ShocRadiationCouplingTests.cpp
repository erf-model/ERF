#include "ERF_RadiationCloudFraction.H"

#include <gtest/gtest.h>

TEST(ShocRadiationCloudFraction, BinaryClearAndLiquidFallback)
{
    const auto clear = radiation_cloud_fractions(0.0, 0.0, 0.0, false);
    EXPECT_EQ(clear.liquid, 0.0);
    EXPECT_EQ(clear.ice, 0.0);
    EXPECT_EQ(clear.total, 0.0);

    const auto liquid = radiation_cloud_fractions(1.0e-4, 0.0, 0.0, false);
    EXPECT_EQ(liquid.liquid, 1.0);
    EXPECT_EQ(liquid.ice, 0.0);
    EXPECT_EQ(liquid.total, 1.0);
}

TEST(ShocRadiationCloudFraction, SuppliedLiquidFractionIsBounded)
{
    const auto fractional = radiation_cloud_fractions(1.0e-4, 0.0, 0.25, true);
    EXPECT_DOUBLE_EQ(fractional.liquid, 0.25);
    EXPECT_EQ(fractional.ice, 0.0);
    EXPECT_DOUBLE_EQ(fractional.total, 0.25);

    const auto below_zero = radiation_cloud_fractions(1.0e-4, 0.0, -0.5, true);
    EXPECT_DOUBLE_EQ(below_zero.liquid, 1.0e-4);

    const auto above_one = radiation_cloud_fractions(1.0e-4, 0.0, 1.5, true);
    EXPECT_DOUBLE_EQ(above_one.liquid, 1.0);
}

TEST(ShocRadiationCloudFraction, IceRemainsBinaryAndPreservesIceOnlyClouds)
{
    const auto ice_only = radiation_cloud_fractions(0.0, 1.0e-4, 0.0, true);
    EXPECT_EQ(ice_only.liquid, 0.0);
    EXPECT_EQ(ice_only.ice, 1.0);
    EXPECT_EQ(ice_only.total, 1.0);

    const auto mixed = radiation_cloud_fractions(1.0e-4, 1.0e-4, 0.3, true);
    EXPECT_DOUBLE_EQ(mixed.liquid, 0.3);
    EXPECT_EQ(mixed.ice, 1.0);
    EXPECT_EQ(mixed.total, 1.0);

    const auto diagnosed_clear = radiation_cloud_fractions(0.0, 0.0, 0.7, true);
    EXPECT_EQ(diagnosed_clear.liquid, 0.0);
    EXPECT_EQ(diagnosed_clear.ice, 0.0);
    EXPECT_EQ(diagnosed_clear.total, 0.0);
}
