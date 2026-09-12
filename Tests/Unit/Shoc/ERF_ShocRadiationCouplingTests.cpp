#include "ERF_RadiationCloudFraction.H"

#include <gtest/gtest.h>
#include <limits>

namespace {

using amrex::Real;

Real
contract_tolerance (const Real expected)
{
    // The cloud-mass helper performs division, a clamp, and two products.
    // This is deliberately much smaller than the mixed-phase factor-of-four
    // defect that the witness is intended to catch.
    const Real scale = (expected > Real(1.0)) ? expected : Real(1.0);
    return Real(128.0) * std::numeric_limits<Real>::epsilon() * scale;
}

void
expect_mass_contract (const Real mixing_ratio,
                      const Real fraction,
                      const Real expected_grid_mass,
                      const Real expected_in_cloud_mass,
                      const Real rho = Real(1.0),
                      const Real dz = Real(100.0))
{
    const Real in_cloud_mass = radiation_cloud_mass(mixing_ratio, fraction, rho, dz);
    EXPECT_NEAR(in_cloud_mass, expected_in_cloud_mass,
                contract_tolerance(expected_in_cloud_mass));
    EXPECT_NEAR(fraction * in_cloud_mass, expected_grid_mass,
                contract_tolerance(expected_grid_mass));
}

} // namespace

TEST(ShocRadiationCloudFraction, BinaryClearAndLiquidFallback)
{
    const auto clear = radiation_cloud_fractions(Real(0.0), Real(0.0), Real(0.0), false);
    EXPECT_EQ(clear.liquid, Real(0.0));
    EXPECT_EQ(clear.ice, Real(0.0));
    EXPECT_EQ(clear.total, Real(0.0));

    const auto liquid = radiation_cloud_fractions(Real(1.0e-4), Real(0.0), Real(0.0), false);
    EXPECT_EQ(liquid.liquid, Real(1.0));
    EXPECT_EQ(liquid.ice, Real(0.0));
    EXPECT_EQ(liquid.total, Real(1.0));
}

TEST(ShocRadiationCloudFraction, SuppliedLiquidFractionIsBounded)
{
    const auto fractional = radiation_cloud_fractions(Real(1.0e-4), Real(0.0), Real(0.25), true);
    EXPECT_EQ(fractional.liquid, Real(0.25));
    EXPECT_EQ(fractional.ice, Real(0.0));
    EXPECT_EQ(fractional.total, Real(0.25));

    const auto below_zero = radiation_cloud_fractions(Real(1.0e-4), Real(0.0), Real(-0.5), true);
    EXPECT_EQ(below_zero.liquid, Real(1.0e-4));

    const auto above_one = radiation_cloud_fractions(Real(1.0e-4), Real(0.0), Real(1.5), true);
    EXPECT_EQ(above_one.liquid, Real(1.0));
}

TEST(ShocRadiationCloudFraction, IceRemainsBinaryAndPreservesIceOnlyClouds)
{
    const auto ice_only = radiation_cloud_fractions(Real(0.0), Real(1.0e-4), Real(0.0), true);
    EXPECT_EQ(ice_only.liquid, Real(0.0));
    EXPECT_EQ(ice_only.ice, Real(1.0));
    EXPECT_EQ(ice_only.total, Real(1.0));

    const auto mixed = radiation_cloud_fractions(Real(1.0e-4), Real(1.0e-4), Real(0.3), true);
    EXPECT_EQ(mixed.liquid, Real(0.3));
    EXPECT_EQ(mixed.ice, Real(1.0));
    EXPECT_EQ(mixed.total, Real(1.0));

    const auto diagnosed_clear = radiation_cloud_fractions(Real(0.0), Real(0.0), Real(0.7), true);
    EXPECT_EQ(diagnosed_clear.liquid, Real(0.0));
    EXPECT_EQ(diagnosed_clear.ice, Real(0.0));
    EXPECT_EQ(diagnosed_clear.total, Real(0.0));
}

TEST(ShocRadiationCloudFraction, SharedMaskMassContract)
{
    constexpr Real rho = Real(1.0);
    constexpr Real dz = Real(100.0);

    const auto clear = radiation_cloud_fractions(Real(0.0), Real(0.0), Real(0.7), true);
    expect_mass_contract(Real(0.0), clear.total, Real(0.0), Real(0.0), rho, dz);

    const auto liquid = radiation_cloud_fractions(Real(1.0e-4), Real(0.0), Real(0.25), true);
    expect_mass_contract(Real(1.0e-4), liquid.total, Real(0.01), Real(0.04), rho, dz);

    const auto ice = radiation_cloud_fractions(Real(0.0), Real(1.0e-4), Real(0.0), true);
    expect_mass_contract(Real(1.0e-4), ice.total, Real(0.01), Real(0.01), rho, dz);

    const auto mixed = radiation_cloud_fractions(Real(1.0e-4), Real(1.0e-4), Real(0.25), true);
    const Real mixed_lwp = radiation_cloud_mass(Real(1.0e-4), mixed.total, rho, dz);
    const Real mixed_iwp = radiation_cloud_mass(Real(1.0e-4), mixed.total, rho, dz);
    EXPECT_NEAR(mixed.total * mixed_lwp, Real(0.01), contract_tolerance(Real(0.01)));
    EXPECT_NEAR(mixed.total * mixed_iwp, Real(0.01), contract_tolerance(Real(0.01)));

    const auto floor_case = radiation_cloud_fractions(Real(1.0e-8), Real(0.0), Real(0.0), true);
    EXPECT_EQ(floor_case.total, Real(1.0e-4));
    expect_mass_contract(Real(1.0e-8), floor_case.total, Real(1.0e-6), Real(0.01), rho, dz);

    const auto cap_case = radiation_cloud_fractions(Real(1.0e-2), Real(1.0e-4), Real(0.25), true);
    EXPECT_EQ(cap_case.total, Real(1.0));
    expect_mass_contract(Real(1.0e-2), cap_case.total, Real(0.5), Real(0.5), rho, dz);
}

TEST(ShocRadiationCloudFraction, SamplingFlagSelectsProductionDispatch)
{
    EXPECT_EQ(radiation_cloud_sampling_mode(true), RadiationCloudSamplingMode::MCICA);
    EXPECT_EQ(radiation_cloud_sampling_mode(false), RadiationCloudSamplingMode::DeterministicBinary);

    EXPECT_FALSE(radiation_sampling_config_unsupported(true, true, true));
    EXPECT_TRUE(radiation_sampling_config_unsupported(false, true, true));
    EXPECT_FALSE(radiation_sampling_config_unsupported(false, true, false));
    EXPECT_FALSE(radiation_sampling_config_unsupported(false, false, true));

    // Binary fallback ignores a supplied diagnostic value, as the production
    // path does when rad_use_shoc_cldfrac is false.
    const auto binary = radiation_cloud_fractions(Real(1.0e-4), Real(0.0), Real(0.25), false);
    EXPECT_EQ(binary.total, Real(1.0));
}
