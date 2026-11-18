// Two host-side rules that keep a missing surface value from reaching output as a number:
// the Noah-MP 2 m tile bundle ignores a tile with zero weight, and the station sampler
// leaves missing neighbours out of its bilinear interpolation.
#include <algorithm>
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "Diagnostics/ERF_NearSurfaceDiagnostics.H"
#include "ERF_StationSampler.H"

namespace
{

using amrex::Real;
using near_surface_diagnostics::NativeBundle;
using near_surface_diagnostics::aggregate_native_bundle;
using near_surface_diagnostics::native_bundle_is_valid;
using near_surface_diagnostics::native_bundle_mixing_ratio_is_valid;
using near_surface_diagnostics::native_bundle_temperature_is_valid;

Real tolerance (Real expected)
{
    return Real(128.0) * std::numeric_limits<Real>::epsilon() *
           std::max(std::abs(expected), Real(1.0));
}

Real plain_bilinear (Real v00, Real v10, Real v01, Real v11, Real wx, Real wy)
{
    return (Real(1.0)-wy) * ((Real(1.0)-wx)*v00 + wx*v10)
         +            wy  * ((Real(1.0)-wx)*v01 + wx*v11);
}

} // namespace

// Motivation: Noah-MP writes T2MVXY = 0 where FVEGXY = 0; rejecting the bundle there
// sent those land cells to the MOST fallback, which filled them with -999.
TEST(NearSurfaceDiagnostics, NativeBundleIgnoresTileWithZeroWeight)
{
    // Bare cell: the vegetated tile's placeholder 0 is not checked or used
    const NativeBundle bare{0.0, 268.8, 0.0, 0.0021, 0.0};
    EXPECT_TRUE(native_bundle_temperature_is_valid(bare));
    EXPECT_TRUE(native_bundle_mixing_ratio_is_valid(bare));
    EXPECT_TRUE(native_bundle_is_valid(bare));
    const auto bare_result = aggregate_native_bundle(bare);
    EXPECT_NEAR(bare_result.temperature, 268.8, tolerance(268.8));
    EXPECT_NEAR(bare_result.mixing_ratio, 0.0021, tolerance(0.0021));

    // Fully vegetated cell: a sentinel in the bare tile cannot leak in through 0 * value
    const NativeBundle vegetated{271.5, -9999.0, 0.0030, -9999.0, 1.0};
    EXPECT_TRUE(native_bundle_is_valid(vegetated));
    const auto veg_result = aggregate_native_bundle(vegetated);
    EXPECT_NEAR(veg_result.temperature, 271.5, tolerance(271.5));
    EXPECT_NEAR(veg_result.mixing_ratio, 0.0030, tolerance(0.0030));
    const NativeBundle vegetated_nan{271.5, std::numeric_limits<Real>::quiet_NaN(),
                                     0.0030, 0.0020, 1.0};
    EXPECT_TRUE(native_bundle_temperature_is_valid(vegetated_nan));
    EXPECT_NEAR(aggregate_native_bundle(vegetated_nan).temperature, 271.5, tolerance(271.5));
}

// Motivation: a tile that carries weight must still be valid, so a partial bundle is
// rejected as before and the fallback decides.
TEST(NearSurfaceDiagnostics, NativeBundleStillRejectsWeightedInvalidTile)
{
    EXPECT_FALSE(native_bundle_temperature_is_valid(NativeBundle{0.0, 268.8, 0.003, 0.002, 0.5}));
    EXPECT_FALSE(native_bundle_temperature_is_valid(NativeBundle{0.0, 268.8, 0.003, 0.002, 1.0}));
    EXPECT_FALSE(native_bundle_temperature_is_valid(NativeBundle{271.5, -9999.0, 0.003, 0.002, 0.0}));
    EXPECT_FALSE(native_bundle_mixing_ratio_is_valid(NativeBundle{271.5, 268.8, -0.001, 0.002, 0.25}));
    EXPECT_FALSE(native_bundle_mixing_ratio_is_valid(NativeBundle{271.5, 268.8, 0.003, -9999.0, 0.0}));
    EXPECT_FALSE(native_bundle_is_valid(NativeBundle{271.5, 268.8, 0.003, 0.002, -9999.0}));
    EXPECT_FALSE(native_bundle_is_valid(NativeBundle{271.5, 268.8, 0.003, 0.002, 1.1}));

    const auto mixed = aggregate_native_bundle(NativeBundle{300.0, 280.0, 0.004, 0.001, 0.25});
    EXPECT_NEAR(mixed.temperature, 285.0, tolerance(285.0));
}

// Motivation: a missing (-999) neighbour with a small weight must not shift a 2 m
// station value; missing neighbours must be left out, not averaged in.
TEST(StationSampler, BilinearSkipsMissingNeighbours)
{
    const Real missing = Real(-999.0);
    const Real wx = Real(0.3), wy = Real(0.6);

    // No missing value: identical to the plain bilinear form
    EXPECT_NEAR(station_bilinear_skip_missing(270.0, 271.0, 272.0, 273.0, wx, wy, missing),
                plain_bilinear(270.0, 271.0, 272.0, 273.0, wx, wy),
                tolerance(272.0));

    // One missing neighbour: the others are renormalised
    const Real w10 = wx*(Real(1.0)-wy), w01 = (Real(1.0)-wx)*wy, w11 = wx*wy;
    const Real expected = (w10*271.0 + w01*272.0 + w11*273.0) / (w10 + w01 + w11);
    EXPECT_NEAR(station_bilinear_skip_missing(missing, 271.0, 272.0, 273.0, wx, wy, missing),
                expected, tolerance(expected));

    // A nearly-zero-weight missing neighbour no longer pulls the value toward -999
    const Real near_one = Real(1.0) - Real(1.0e-4);
    const Real value = station_bilinear_skip_missing(270.0, 270.0, 270.0, missing,
                                                     near_one, near_one, missing);
    EXPECT_NEAR(value, 270.0, tolerance(270.0));

    // A non-finite neighbour is treated as missing
    EXPECT_NEAR(station_bilinear_skip_missing(270.0, std::numeric_limits<Real>::quiet_NaN(),
                                              270.0, 270.0, wx, wy, missing),
                270.0, tolerance(270.0));

    // Missing only when every neighbour with weight is missing; zero-weight values do not count
    EXPECT_EQ(station_bilinear_skip_missing(missing, missing, missing, missing, wx, wy, missing),
              missing);
    EXPECT_EQ(station_bilinear_skip_missing(missing, 271.0, 272.0, 273.0,
                                            Real(0.0), Real(0.0), missing),
              missing);
    EXPECT_NEAR(station_bilinear_skip_missing(270.0, missing, missing, missing,
                                              Real(0.0), Real(0.0), missing),
                270.0, tolerance(270.0));
}
