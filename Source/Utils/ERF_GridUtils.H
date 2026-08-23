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
    amrex::Real dx{0.0};
    amrex::Real dy{0.0};
    amrex::Real xmin{0.0};
    amrex::Real ymin{0.0};
};

inline amrex::Real
comparison_tolerance (amrex::Real lhs, amrex::Real rhs)
{
    const amrex::Real scale = std::max({amrex::Real(1.0), std::abs(lhs), std::abs(rhs)});
    return amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() * scale;
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
    if (coordinates.size() != static_cast<std::size_t>(expected_size)) {
        return field_description + " has " + axis_name + " coordinate length " +
               std::to_string(coordinates.size()) + ", expected " +
               std::to_string(expected_size);
    }

    origin = coordinates.front();
    spacing = coordinates[1] - coordinates[0];
    if (!std::isfinite(origin) || !std::isfinite(spacing) || spacing <= amrex::Real(0.0)) {
        return field_description + " must have finite, strictly increasing " +
               axis_name + " coordinates with nonzero spacing";
    }

    for (int index = 1; index < expected_size; ++index) {
        const amrex::Real coordinate = coordinates[index];
        const amrex::Real expected = origin + static_cast<amrex::Real>(index) * spacing;
        if (!std::isfinite(coordinate) ||
            (index > 0 && coordinate <= coordinates[index-1])) {
            return field_description + " must have finite, strictly increasing " +
                   axis_name + " coordinates";
        }
        if (!nearly_equal(coordinate, expected)) {
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
    amrex::Real weight{0.0};
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
    const amrex::Real magnitude =
        ((amrex::Real(1.0) > std::abs(origin)) ? amrex::Real(1.0) : std::abs(origin));
    const amrex::Real scale = (magnitude > std::abs(upper)) ? magnitude : std::abs(upper);
    const amrex::Real tolerance =
        amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() * scale;

    if (coordinate < origin - tolerance || coordinate > upper + tolerance) {
        return stencil;
    }

    const amrex::Real bounded_coordinate =
        (coordinate < origin) ? origin : ((coordinate > upper) ? upper : coordinate);
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
