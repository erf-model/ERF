#ifndef ERF_GRID_UTILS_H_
#define ERF_GRID_UTILS_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace erf_grid_utils {

struct UniformGridMetadata
{
    int nx{0};
    int ny{0};
    amrex::Real dx{amrex::Real(0.0)};
    amrex::Real dy{amrex::Real(0.0)};
    amrex::Real xmin{amrex::Real(0.0)};
    amrex::Real ymin{amrex::Real(0.0)};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
comparison_tolerance (amrex::Real lhs, amrex::Real rhs) noexcept
{
    amrex::Real scale = amrex::Real(1.0);
    const amrex::Real lhs_abs = amrex::Math::abs(lhs);
    const amrex::Real rhs_abs = amrex::Math::abs(rhs);
    if (lhs_abs > scale) scale = lhs_abs;
    if (rhs_abs > scale) scale = rhs_abs;
    // NetCDF coordinate variables are commonly stored as float32 even when
    // ERF is built with double precision.  A few float32 ulps at the
    // coordinate magnitude are therefore source-file quantization, not grid
    // nonuniformity.  The Real tolerance is retained for native Real data;
    // the float tolerance is deliberately a small ulp-based allowance rather
    // than an arbitrary absolute threshold.
    const amrex::Real real_ulp =
        amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() * scale;
    const amrex::Real float_ulp =
        amrex::Real(8.0) * static_cast<amrex::Real>(std::numeric_limits<float>::epsilon()) * scale;
    return (real_ulp > float_ulp) ? real_ulp : float_ulp;
}

inline bool
nearly_equal (amrex::Real lhs, amrex::Real rhs)
{
    return std::abs(lhs - rhs) <= comparison_tolerance(lhs, rhs);
}

inline std::string
validate_uniform_axis (const amrex::Vector<amrex::Real>& coordinates,
                       int expected_size,
                       const std::string& axis_name,
                       const std::string& field_description,
                       amrex::Real& origin,
                       amrex::Real& spacing)
{
    if (expected_size < 2) {
        return field_description + " requires at least two " + axis_name +
               " coordinates for bilinear interpolation";
    }
    if (coordinates.size() != static_cast<amrex::Long>(expected_size)) {
        return field_description + " has " + axis_name + " coordinate length " +
               std::to_string(coordinates.size()) + ", expected " +
               std::to_string(expected_size);
    }

    origin = coordinates.front();
    // Derive the spacing from the two endpoints rather than from the first
    // interval alone.  The returned origin/spacing pair is what
    // uniform_interpolation_stencil uses to reconstruct the axis as
    // origin + index*spacing, so a first interval that is rounded slightly
    // short or long would otherwise accumulate that error over the whole axis
    // -- for a float32 axis with many points the reconstructed far end can
    // drift by (point_count-1) ulps, which is large enough to shift the
    // stencil by a cell near the upper edge.  Averaging over the full extent
    // pins both endpoints instead.  A non-finite last coordinate propagates
    // into spacing and is caught by the isfinite check below.
    spacing = (coordinates[expected_size-1] - coordinates.front()) /
              static_cast<amrex::Real>(expected_size - 1);
    if (!std::isfinite(origin) || !std::isfinite(spacing) || spacing <= amrex::Real(0.0)) {
        return field_description + " must have finite, strictly increasing " +
               axis_name + " coordinates with nonzero spacing";
    }

    for (int index = 1; index < expected_size; ++index) {
        const amrex::Real coordinate = coordinates[index];
        if (!std::isfinite(coordinate) || coordinate <= coordinates[index-1]) {
            return field_description + " must have finite, strictly increasing " +
                   axis_name + " coordinates";
        }
        // Compare each interval against the mean spacing instead of
        // accumulating the first spacing over the whole axis.  This avoids
        // rejecting a valid float32 axis because each stored coordinate was
        // rounded independently, while still rejecting materially nonuniform
        // grids.
        //
        // The tolerance is the *larger* of the two allowances, not the
        // smaller: the ulp term covers source-file quantization, which scales
        // with the coordinate magnitude, while the fraction-of-spacing term
        // covers grids whose coordinates are small enough that a few ulps are
        // negligible.  Taking the minimum cancels the ulp allowance in exactly
        // the case it exists for -- a fine grid at a large offset.  Concretely,
        // at a UTM easting near 5e5 the float32 ulp is 0.03125 m, so a 0.4 m
        // raster stores intervals that deviate from the nominal spacing by up
        // to 0.025 m: more than 5% of the spacing, yet far below the ulp
        // allowance.  Integer-metre spacings are unaffected either way, since
        // integers are exactly representable in float32 at UTM magnitudes.
        const amrex::Real local_spacing = coordinate - coordinates[index-1];
        const amrex::Real spacing_tolerance = std::max(
            amrex::Real(0.05) * std::abs(spacing),
            comparison_tolerance(coordinate, coordinates[index-1]));
        if (std::abs(local_spacing - spacing) > spacing_tolerance) {
            return field_description + " has nonuniform " + axis_name +
                   " coordinates, which are not supported";
        }
    }

    return {};
}

inline std::string
validate_matching_grid (const UniformGridMetadata& reference,
                        const UniformGridMetadata& candidate,
                        const std::string& reference_description,
                        const std::string& candidate_description)
{
    if (candidate.nx != reference.nx || candidate.ny != reference.ny) {
        return candidate_description + " dimensions (" + std::to_string(candidate.nx) +
               " x " + std::to_string(candidate.ny) + ") do not match " +
               reference_description + " dimensions (" + std::to_string(reference.nx) +
               " x " + std::to_string(reference.ny) + ")";
    }

    struct MetadataValue {
        const char* name;
        amrex::Real reference_value;
        amrex::Real candidate_value;
    };
    const MetadataValue values[] = {
        {"x spacing", reference.dx, candidate.dx},
        {"y spacing", reference.dy, candidate.dy},
        {"x origin", reference.xmin, candidate.xmin},
        {"y origin", reference.ymin, candidate.ymin}
    };
    for (const auto& value : values) {
        if (!nearly_equal(value.reference_value, value.candidate_value)) {
            std::ostringstream message;
            message << candidate_description << " " << value.name << " ("
                    << value.candidate_value << ") does not match "
                    << reference_description << " (" << value.reference_value << ")";
            return message.str();
        }
    }

    return {};
}

struct InterpolationStencil
{
    int lower{0};
    amrex::Real weight{amrex::Real(0.0)};
    bool inside{false};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
InterpolationStencil
uniform_interpolation_stencil (amrex::Real coordinate,
                               amrex::Real origin,
                               amrex::Real spacing,
                               int point_count) noexcept
{
    InterpolationStencil stencil;
    if (point_count < 2 || spacing <= amrex::Real(0.0)) {
        return stencil;
    }

    const amrex::Real upper = origin + static_cast<amrex::Real>(point_count - 1) * spacing;
    const amrex::Real lower_tolerance = comparison_tolerance(coordinate, origin);
    const amrex::Real upper_tolerance = comparison_tolerance(coordinate, upper);

    if (coordinate < origin - lower_tolerance || coordinate > upper + upper_tolerance) {
        return stencil;
    }

    const amrex::Real bounded_coordinate =
        (coordinate <= origin + lower_tolerance) ? origin :
        ((coordinate >= upper - upper_tolerance) ? upper : coordinate);
    const amrex::Real fractional_index = (bounded_coordinate - origin) / spacing;
    int lower = static_cast<int>(amrex::Math::floor(fractional_index));
    lower = (lower < point_count - 1) ? lower : point_count - 2;

    amrex::Real weight = fractional_index - static_cast<amrex::Real>(lower);
    weight = (weight < amrex::Real(0.0)) ? amrex::Real(0.0) :
             ((weight > amrex::Real(1.0)) ? amrex::Real(1.0) : weight);

    stencil.lower = lower;
    stencil.weight = weight;
    stencil.inside = true;
    return stencil;
}

} // namespace erf_grid_utils

#endif
