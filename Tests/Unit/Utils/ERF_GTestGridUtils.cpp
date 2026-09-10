#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GridUtils.H"

using erf_grid_utils::UniformGridMetadata;

// Motivation: every gridded forest field is interpolated with LAI metadata,
// so a shifted height grid must fail in constant-Cd mode before interpolation.
TEST(ForestGridMetadata, ConstantCdRejectsShiftedHeightOrigin)
{
    const UniformGridMetadata lai{4, 3, 10.0, 20.0, 100.0, 200.0};
    const UniformGridMetadata height{4, 3, 10.0, 20.0, 101.0, 200.0};

    const auto error = erf_grid_utils::validate_matching_grid(
        lai, height, "forest field 'LAI' in 'lai.nc'",
        "forest field 'height' in 'height.nc'");
    EXPECT_NE(error.find("height.nc"), std::string::npos);
    EXPECT_NE(error.find("x origin"), std::string::npos);
}

// Motivation: file-Cd mode must reject a coefficient field whose dimensions
// happen to match but whose physical spacing differs from LAI.
TEST(ForestGridMetadata, FileCdRejectsDifferentSpacing)
{
    const UniformGridMetadata lai{4, 3, 10.0, 20.0, 100.0, 200.0};
    const UniformGridMetadata cd{4, 3, 10.5, 20.0, 100.0, 200.0};

    const auto error = erf_grid_utils::validate_matching_grid(
        lai, cd, "forest field 'LAI' in 'lai.nc'",
        "forest field 'cd' in 'cd.nc'");
    EXPECT_NE(error.find("cd.nc"), std::string::npos);
    EXPECT_NE(error.find("x spacing"), std::string::npos);
}

TEST(ForestGridMetadata, MatchingConstantAndFileCdGridsAreAccepted)
{
    const UniformGridMetadata lai{4, 3, 10.0, 20.0, 100.0, 200.0};
    const UniformGridMetadata matching{4, 3, 10.0, 20.0, 100.0, 200.0};

    EXPECT_TRUE(erf_grid_utils::validate_matching_grid(
        lai, matching, "LAI", "height").empty());
    EXPECT_TRUE(erf_grid_utils::validate_matching_grid(
        lai, matching, "LAI", "cd").empty());
}

TEST(ForestGridMetadata, CoordinateValidationRejectsUnsupportedAxes)
{
    amrex::Real origin = 0.0;
    amrex::Real spacing = 0.0;

    EXPECT_TRUE(erf_grid_utils::validate_uniform_axis(
        {0.0, 1.0, 2.0}, 3, "x", "field", origin, spacing).empty());
    EXPECT_FALSE(erf_grid_utils::validate_uniform_axis(
        {0.0}, 1, "x", "field", origin, spacing).empty());
    EXPECT_FALSE(erf_grid_utils::validate_uniform_axis(
        {2.0, 1.0}, 2, "x", "field", origin, spacing).empty());
    EXPECT_FALSE(erf_grid_utils::validate_uniform_axis(
        {0.0, 1.0, 2.25}, 3, "x", "field", origin, spacing).empty());
    EXPECT_FALSE(erf_grid_utils::validate_uniform_axis(
        {0.0, 1.0}, 3, "x", "field", origin, spacing).empty());
}

TEST(ForestGridMetadata, AcceptsFloatStoredLongitudeCoordinates)
{
    amrex::Vector<amrex::Real> longitude;
    for (int index = 0; index < 128; ++index) {
        // Simulate a NetCDF float variable before it is converted to Real.
        longitude.push_back(static_cast<amrex::Real>(
            static_cast<float>(-105.0f + static_cast<float>(index) * 0.01f)));
    }

    amrex::Real origin = 0.0;
    amrex::Real spacing = 0.0;
    EXPECT_TRUE(erf_grid_utils::validate_uniform_axis(
        longitude, static_cast<int>(longitude.size()), "longitude", "float-backed field",
        origin, spacing).empty());
}

TEST(ForestGridMetadata, RejectsMateriallyNonuniformFloatBackedCoordinates)
{
    amrex::Vector<amrex::Real> longitude;
    for (int index = 0; index < 32; ++index) {
        longitude.push_back(static_cast<amrex::Real>(
            static_cast<float>(-105.0f + static_cast<float>(index) * 0.01f)));
    }
    longitude[16] += amrex::Real(0.002);

    amrex::Real origin = 0.0;
    amrex::Real spacing = 0.0;
    EXPECT_FALSE(erf_grid_utils::validate_uniform_axis(
        longitude, static_cast<int>(longitude.size()), "longitude", "float-backed field",
        origin, spacing).empty());
}

// Build a float32-quantized axis, as a NetCDF float coordinate variable
// reaches ERF after conversion to Real.
static amrex::Vector<amrex::Real>
float_backed_axis (double axis_origin, double axis_spacing, int point_count)
{
    amrex::Vector<amrex::Real> axis;
    axis.reserve(point_count);
    for (int index = 0; index < point_count; ++index) {
        axis.push_back(static_cast<amrex::Real>(
            static_cast<float>(axis_origin + double(index) * axis_spacing)));
    }
    return axis;
}

// Motivation: a sub-metre canopy raster in a projected coordinate system is
// where float32 quantization is largest relative to the spacing.  At a UTM
// easting near 5e5 the float32 ulp is 0.03125 m, so on a 0.4 m grid
// neighbouring intervals deviate from the nominal spacing by up to 0.025 m --
// more than 5% of the spacing, but far less than the ulp allowance.  The
// spacing tolerance must be the larger of the two allowances, or this
// perfectly uniform grid is rejected with no way for the user to work around
// it.
TEST(ForestGridMetadata, AcceptsSubMetreFloatBackedGridAtLargeProjectedOffset)
{
    const int point_count = 1000;
    const auto easting = float_backed_axis(500000.0, 0.4, point_count);

    amrex::Real origin = 0.0;
    amrex::Real spacing = 0.0;
    EXPECT_TRUE(erf_grid_utils::validate_uniform_axis(
        easting, point_count, "x", "float-backed UTM field",
        origin, spacing).empty());
}

// Motivation: the origin/spacing pair returned here is exactly what
// uniform_interpolation_stencil uses to reconstruct the axis as
// origin + index*spacing.  A spacing taken from the first interval alone
// carries that interval's quantization error into every index, so the
// reconstructed far end drifts.  For a 1.1 m float32 axis at a UTM easting
// near 5e5 the drift reaches 6.25 m -- more than five cells -- on an axis
// that passes validation, silently shifting the interpolation stencil near
// the upper edge.  Pin the reconstructed endpoint to the stored one.
TEST(ForestGridMetadata, ReturnedSpacingReconstructsTheStoredFarEndpoint)
{
    const int point_count = 1000;
    const auto easting = float_backed_axis(500000.0, 1.1, point_count);

    amrex::Real origin = 0.0;
    amrex::Real spacing = 0.0;
    ASSERT_TRUE(erf_grid_utils::validate_uniform_axis(
        easting, point_count, "x", "float-backed UTM field",
        origin, spacing).empty());

    const amrex::Real stored_upper = easting[point_count-1];
    const amrex::Real reconstructed_upper =
        origin + static_cast<amrex::Real>(point_count - 1) * spacing;
    EXPECT_LE(std::abs(reconstructed_upper - stored_upper),
              erf_grid_utils::comparison_tolerance(reconstructed_upper, stored_upper));

    const auto upper_stencil = erf_grid_utils::uniform_interpolation_stencil(
        stored_upper, origin, spacing, point_count);
    EXPECT_TRUE(upper_stencil.inside);
    EXPECT_EQ(upper_stencil.lower, point_count - 2);
}

TEST(ForestGridMetadata, FloatAndDoubleGridMetadataMatchWithinStorageQuantization)
{
    const UniformGridMetadata double_grid{8, 8, 0.01, 0.01, -105.0, 40.0};
    const UniformGridMetadata float_grid{8, 8, 0.010000001, 0.01, -104.999996, 40.000003};

    EXPECT_TRUE(erf_grid_utils::validate_matching_grid(
        double_grid, float_grid, "double grid", "float grid").empty());

    const UniformGridMetadata shifted_grid{8, 8, 0.01, 0.01, -104.99, 40.0};
    EXPECT_FALSE(erf_grid_utils::validate_matching_grid(
        double_grid, shifted_grid, "double grid", "shifted grid").empty());
}

struct EndpointSamples
{
    int lower_endpoint;
    amrex::Real endpoint_weight;
    int endpoint_inside;
    int outside_inside;
};

// Keep the extended device lambda outside GoogleTest's private TestBody so
// nvcc can compile this test with CUDA extended-lambda checks enabled.
EndpointSamples
sample_device_endpoints ()
{
    amrex::Gpu::DeviceVector<int> device_lower(2, -1);
    amrex::Gpu::DeviceVector<amrex::Real> device_weight(2, -1.0);
    amrex::Gpu::DeviceVector<int> device_inside(2, 0);
    int* lower = device_lower.data();
    amrex::Real* weight = device_weight.data();
    int* inside = device_inside.data();

    amrex::ParallelFor(2, [=] AMREX_GPU_DEVICE (int index) noexcept {
        const amrex::Real coordinate = (index == 0) ? amrex::Real(2.0) : amrex::Real(2.25);
        const auto stencil = erf_grid_utils::uniform_interpolation_stencil(
            coordinate, amrex::Real(0.0), amrex::Real(1.0), 3);
        lower[index] = stencil.lower;
        weight[index] = stencil.weight;
        inside[index] = stencil.inside ? 1 : 0;
    });

    amrex::Gpu::HostVector<int> host_lower(2);
    amrex::Gpu::HostVector<amrex::Real> host_weight(2);
    amrex::Gpu::HostVector<int> host_inside(2);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_lower.begin(), device_lower.end(), host_lower.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_weight.begin(), device_weight.end(), host_weight.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_inside.begin(), device_inside.end(), host_inside.begin());
    amrex::Gpu::streamSynchronize();

    return EndpointSamples{host_lower[0], host_weight[0], host_inside[0], host_inside[1]};
}

// Motivation: the endpoint branch executes in terrain and canopy device
// kernels. Verify both upper-endpoint clamping and the no-extrapolation policy
// on the active AMReX backend.
TEST(UniformGridInterpolation, DeviceEndpointIsInsideAndOutsidePointIsRejected)
{
    const auto samples = sample_device_endpoints();
    EXPECT_EQ(samples.endpoint_inside, 1);
    EXPECT_EQ(samples.lower_endpoint, 1);
    EXPECT_EQ(samples.endpoint_weight, amrex::Real(1.0));
    EXPECT_EQ(samples.outside_inside, 0);
}

TEST(UniformGridInterpolation, FloatQuantizedNonExactEndpointsAreClamped)
{
    const amrex::Real origin = amrex::Real(1234.5678);
    const amrex::Real spacing = amrex::Real(0.125);
    const int point_count = 4;
    const amrex::Real upper = origin + amrex::Real(point_count - 1) * spacing;

    const amrex::Real float_lower = static_cast<amrex::Real>(static_cast<float>(origin));
    const amrex::Real float_upper = static_cast<amrex::Real>(static_cast<float>(upper));
    const auto lower = erf_grid_utils::uniform_interpolation_stencil(
        float_lower, origin, spacing, point_count);
    const auto upper_stencil = erf_grid_utils::uniform_interpolation_stencil(
        float_upper, origin, spacing, point_count);
    EXPECT_TRUE(lower.inside);
    EXPECT_EQ(lower.lower, 0);
    EXPECT_EQ(lower.weight, amrex::Real(0.0));
    EXPECT_TRUE(upper_stencil.inside);
    EXPECT_EQ(upper_stencil.lower, point_count - 2);
    EXPECT_EQ(upper_stencil.weight, amrex::Real(1.0));

    const amrex::Real outside = upper +
        amrex::Real(2.0) * erf_grid_utils::comparison_tolerance(upper, upper);
    EXPECT_FALSE(erf_grid_utils::uniform_interpolation_stencil(
        outside, origin, spacing, point_count).inside);
}
