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

// ---------------------------------------------------------------------------
// Geometric vertical stretch
//
// This is the construction ERF uses when the first layer thickness and the
// domain height are the knowns -- an idealized case, or a wrfinput level whose
// grids do not reach the domain top -- so the properties worth pinning are the
// ones a mesh has to have: it closes on z_top, it starts at the requested dz0,
// and it is monotone.
// ---------------------------------------------------------------------------

namespace {

// The same tolerance init_from_wrfinput passes, so these tests exercise the
// solve at the precision the solver is actually used at.
#ifdef AMREX_USE_FLOAT
constexpr amrex::Real stretch_tol = amrex::Real(1.e-4);
#else
constexpr amrex::Real stretch_tol = amrex::Real(1.e-8);
#endif

amrex::Real
roundoff_allowance (amrex::Real scale, amrex::Real layers)
{
    return amrex::Real(16.0) * layers *
           std::numeric_limits<amrex::Real>::epsilon() * std::abs(scale);
}

} // namespace

// Motivation: the whole point of the solve is that the layers add up to the
// domain height, so a factor that does not close the column is a broken mesh.
TEST(GeometricStretch, ClosesOnTheDomainTop)
{
    const amrex::Real z_top = amrex::Real(16500.0);
    const amrex::Real dz0   = amrex::Real(19.57);
    const int         nz    = 176;

    const auto gs = erf_grid_utils::geometric_stretch(dz0, z_top, nz, stretch_tol);

    EXPECT_TRUE(gs.converged);
    EXPECT_FALSE(gs.uniform);
    EXPECT_EQ(gs.dz0, dz0);
    EXPECT_GT(gs.ratio, amrex::Real(1.0));
    EXPECT_LE(std::abs(gs.residual), stretch_tol);
}

// Motivation: make_terrain_fitted_coords drapes this array over the terrain, so
// it has to be strictly increasing, start one dz0 above zero and end exactly on
// z_top -- the caller asserts the last of those.
TEST(GeometricStretch, BuildsAMonotoneArrayThatEndsOnZTop)
{
    const amrex::Real z_top = amrex::Real(16500.0);
    const amrex::Real dz0   = amrex::Real(19.57);
    const int         nz    = 176;

    amrex::Vector<amrex::Real> zlevels(nz+1);
    const auto gs = erf_grid_utils::build_geometric_zlevels(zlevels, dz0, z_top, stretch_tol);

    ASSERT_TRUE(gs.converged);
    EXPECT_EQ(zlevels[0], amrex::Real(0.0));
    EXPECT_EQ(zlevels[nz], z_top);
    EXPECT_EQ(zlevels[1] - zlevels[0], dz0);

    for (int k(1); k<=nz; ++k) {
        EXPECT_GT(zlevels[k], zlevels[k-1]) << "at k = " << k;
    }

    // Every interior layer is the previous one times the factor.  The top layer
    // is excluded because build_geometric_zlevels assigns z_top there exactly
    // rather than by accumulation.
    const amrex::Real allowance = roundoff_allowance(z_top, amrex::Real(nz));
    for (int k(2); k<nz; ++k) {
        const amrex::Real dz_prev = zlevels[k-1] - zlevels[k-2];
        const amrex::Real dz_this = zlevels[k  ] - zlevels[k-1];
        EXPECT_NEAR(dz_this, dz_prev*gs.ratio, allowance) << "at k = " << k;
    }
}

// Motivation: the thicknesses only ever grow, so a dz0 at or above the uniform
// spacing cannot be stretched to fit.  Returning a factor below one would thin
// the layers aloft, which is the opposite of what a stretched mesh is for; the
// uniform grid is the honest answer and the flag says so.
TEST(GeometricStretch, FallsBackToUniformWhenDz0IsTooThick)
{
    const amrex::Real z_top = amrex::Real(16500.0);
    const int         nz    = 176;
    const amrex::Real dz_unif = z_top / static_cast<amrex::Real>(nz);

    for (const amrex::Real dz0 : {dz_unif, amrex::Real(2.0)*dz_unif}) {
        const auto gs = erf_grid_utils::geometric_stretch(dz0, z_top, nz, stretch_tol);
        EXPECT_TRUE(gs.uniform);
        EXPECT_TRUE(gs.converged);
        EXPECT_EQ(gs.ratio, amrex::Real(1.0));
        EXPECT_EQ(gs.dz0, dz_unif);
    }

    amrex::Vector<amrex::Real> zlevels(nz+1);
    erf_grid_utils::build_geometric_zlevels(zlevels, amrex::Real(500.0), z_top, stretch_tol);
    const amrex::Real allowance = roundoff_allowance(z_top, amrex::Real(nz));
    for (int k(1); k<=nz-1; ++k) {
        EXPECT_NEAR(zlevels[k] - zlevels[k-1], dz_unif, allowance) << "at k = " << k;
    }
    EXPECT_EQ(zlevels[nz], z_top);
}

// Motivation: a finer first layer needs a larger factor to cover the same
// domain with the same number of layers, and the solve has to track that rather
// than sitting at its 1.03 starting guess.
TEST(GeometricStretch, AFinerFirstLayerNeedsAStrongerStretch)
{
    const amrex::Real z_top = amrex::Real(16500.0);
    const int         nz    = 176;

    const auto coarse = erf_grid_utils::geometric_stretch(amrex::Real(40.0), z_top, nz, stretch_tol);
    const auto fine   = erf_grid_utils::geometric_stretch(amrex::Real(5.0),  z_top, nz, stretch_tol);

    ASSERT_TRUE(coarse.converged);
    ASSERT_TRUE(fine.converged);
    EXPECT_GT(fine.ratio, coarse.ratio);
}

// Motivation: the caller aborts on !converged, so a degenerate request has to
// report failure rather than hand back the default-constructed factor as if it
// had been solved for.
TEST(GeometricStretch, ReportsFailureOnDegenerateRequests)
{
    const amrex::Real z_top = amrex::Real(16500.0);

    EXPECT_FALSE(erf_grid_utils::geometric_stretch(
        amrex::Real(-1.0), z_top, 176, stretch_tol).converged);
    EXPECT_FALSE(erf_grid_utils::geometric_stretch(
        amrex::Real(0.0), z_top, 176, stretch_tol).converged);
    EXPECT_FALSE(erf_grid_utils::geometric_stretch(
        amrex::Real(19.57), amrex::Real(0.0), 176, stretch_tol).converged);
    EXPECT_FALSE(erf_grid_utils::geometric_stretch(
        amrex::Real(19.57), z_top, 0, stretch_tol).converged);
}

// Motivation: a single-layer column is the degenerate end of the uniform
// branch, and it must not divide by a zero stretch denominator on the way.
TEST(GeometricStretch, HandlesASingleLayer)
{
    const amrex::Real z_top = amrex::Real(1000.0);

    const auto gs = erf_grid_utils::geometric_stretch(amrex::Real(10.0), z_top, 1, stretch_tol);
    EXPECT_TRUE(gs.uniform);
    EXPECT_TRUE(gs.converged);
    EXPECT_EQ(gs.dz0, z_top);

    amrex::Vector<amrex::Real> zlevels(2);
    erf_grid_utils::build_geometric_zlevels(zlevels, amrex::Real(10.0), z_top, stretch_tol);
    EXPECT_EQ(zlevels[0], amrex::Real(0.0));
    EXPECT_EQ(zlevels[1], z_top);
}
