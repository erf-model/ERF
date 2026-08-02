#include "ERF_PrescribedSurfaceFlux.H"

#include <gtest/gtest.h>

#include <limits>
#include <string>
#include <vector>

using prescribed_surface_flux::TimeBracket;

namespace {

TimeBracket bracket_at (const std::vector<double>& time, const double query)
{
    TimeBracket bracket;
    std::string error;
    EXPECT_TRUE(prescribed_surface_flux::find_time_bracket(
        time, query, bracket, error)) << error;
    return bracket;
}

} // namespace

TEST(PrescribedSurfaceFlux, InterpolatesAndAcceptsExactEndpoints)
{
    const std::vector<double> time{0.0, 10.0, 25.0};
    auto bracket = bracket_at(time, 0.0);
    EXPECT_EQ(bracket.lower, 0U);
    EXPECT_EQ(bracket.upper, 0U);

    bracket = bracket_at(time, 5.0);
    EXPECT_EQ(bracket.lower, 0U);
    EXPECT_EQ(bracket.upper, 1U);
    EXPECT_DOUBLE_EQ(bracket.weight, 0.5);
    EXPECT_DOUBLE_EQ(prescribed_surface_flux::interpolate(2.0, 6.0,
                                                          bracket.weight),
                     4.0);

    bracket = bracket_at(time, 25.0);
    EXPECT_EQ(bracket.lower, 2U);
    EXPECT_EQ(bracket.upper, 2U);
}

TEST(PrescribedSurfaceFlux, SingleTimeAllowsOnlyItsEndpoint)
{
    const auto bracket = bracket_at({4.0}, 4.0);
    EXPECT_EQ(bracket.lower, 0U);
    EXPECT_EQ(bracket.upper, 0U);
}

TEST(PrescribedSurfaceFlux, RejectsTimesOutsideCoverage)
{
    const std::vector<double> time{1.0, 2.0};
    TimeBracket bracket;
    std::string error;
    EXPECT_FALSE(prescribed_surface_flux::find_time_bracket(
        time, 0.9, bracket, error));
    EXPECT_FALSE(error.empty());
    EXPECT_FALSE(prescribed_surface_flux::find_time_bracket(
        time, 2.1, bracket, error));
    EXPECT_FALSE(error.empty());
}

TEST(PrescribedSurfaceFlux, ValidatesStrictFiniteTimeAxis)
{
    EXPECT_TRUE(prescribed_surface_flux::validate_time_axis(
        {0.0, 1.0, 2.0}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_time_axis({}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_time_axis(
        {0.0, 0.0}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_time_axis(
        {1.0, 0.0}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_time_axis(
        {0.0, std::numeric_limits<double>::infinity()}).empty());
}

TEST(PrescribedSurfaceFlux, ValidatesShapeMissingNonfiniteAndNegativeUstar)
{
    EXPECT_TRUE(prescribed_surface_flux::validate_field(
        "theta_flux", {1.0, 2.0}, 2, false).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_field(
        "theta_flux", {1.0}, 2, false).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_field(
        "theta_flux", {1.0, 999.0}, 2, false, {999.0}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_field(
        "qv_flux", {0.0, std::numeric_limits<double>::quiet_NaN()},
        2, false).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_field(
        "ustar", {0.1, -0.1}, 2, true).empty());
}

TEST(PrescribedSurfaceFlux, ValidatesDimensionNamesOrderAndShape)
{
    EXPECT_TRUE(prescribed_surface_flux::validate_time_layout({"time"}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_time_layout({"record"}).empty());
    EXPECT_TRUE(prescribed_surface_flux::validate_field_layout(
        "ustar", {"time", "y", "x"}, {3, 4, 5}, {3, 4, 5}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_field_layout(
        "ustar", {"time", "x", "y"}, {3, 5, 4}, {3, 4, 5}).empty());
    EXPECT_FALSE(prescribed_surface_flux::validate_field_layout(
        "ustar", {"time", "y", "x"}, {3, 4, 6}, {3, 4, 5}).empty());
}
