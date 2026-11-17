#include <gtest/gtest.h>
#include <limits>

#include <ERF_SuperDropletPC.H>

using namespace amrex;

// Sampled-multiplicity weight of a super-droplet over its independently sampled components.
// Every varying weight carries the factor ndv = numdens*cell_volume = 10 once.

namespace {
constexpr Real sdw_ndv = Real(10);
// a few multiplications and divisions: a few float ulps
constexpr Real sdw_rtol = Real(1.0e-6);
}

TEST(SDMSampledWeight, ConstantsAloneAdd)
{
    SDSampledWeight w;
    w.add(SDDistributionType::mass_constant, Real(0.25), sdw_ndv);
    w.add(SDDistributionType::mass_constant, Real(0.5), sdw_ndv);
    EXPECT_EQ(w.value(), Real(0.75));
}

TEST(SDMSampledWeight, OneVaryingComponentIgnoresConstants)
{
    SDSampledWeight w;
    w.add(SDDistributionType::mass_constant, Real(0.25), sdw_ndv);
    w.add(SDDistributionType::radius_log_normal, Real(4), sdw_ndv);
    w.add(SDDistributionType::mass_constant, Real(0.5), sdw_ndv);
    EXPECT_EQ(w.value(), Real(4));
}

TEST(SDMSampledWeight, VaryingWeightsMultiplyInAnyOrder)
{
    const SDDistributionType types[3] = {SDDistributionType::radius_log_normal,
                                         SDDistributionType::mass_exponential,
                                         SDDistributionType::radius_lognormal_autorange};
    const Real weights[3] = {Real(4), Real(5), Real(2)};
    const Real expected = Real(4) * Real(5) * Real(2) / (sdw_ndv * sdw_ndv);
    const int orders[3][3] = {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}};
    for (const auto& order : orders) {
        SDSampledWeight w;
        for (int c : order) { w.add(types[c], weights[c], sdw_ndv); }
        w.add(SDDistributionType::mass_constant, Real(0.5), sdw_ndv);
        EXPECT_NEAR(w.value(), expected, sdw_rtol * expected);
    }
}

TEST(SDMSampledWeight, NegativeNumberDensityKeepsTheSign)
{
    // an injection reads no number density (-1), so its weights are negative
    SDSampledWeight w;
    w.add(SDDistributionType::radius_log_normal, Real(-4), -sdw_ndv);
    w.add(SDDistributionType::mass_exponential, Real(-5), -sdw_ndv);
    EXPECT_NEAR(w.value(), Real(-2), sdw_rtol * Real(2));
}

TEST(SDMSampledWeight, CellScaleForEitherSign)
{
    Real scale = Real(0);
    EXPECT_TRUE(SD_sampled_scale(Real(4), Real(8), scale));
    EXPECT_EQ(scale, Real(2));
    EXPECT_TRUE(SD_sampled_scale(Real(-4), Real(8), scale));
    EXPECT_EQ(scale, Real(-2));
}

TEST(SDMSampledWeight, CellScaleFallsBackForVanishingWeights)
{
    const Real num = Real(1.0e11);
    Real scale = Real(3);
    EXPECT_FALSE(SD_sampled_scale(Real(0), num, scale));
    EXPECT_FALSE(SD_sampled_scale(Real(0.5) * num / std::numeric_limits<Real>::max(), num, scale));
    EXPECT_EQ(scale, Real(3));
}
