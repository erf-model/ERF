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
