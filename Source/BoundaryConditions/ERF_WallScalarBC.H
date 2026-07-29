#ifndef ERF_WALL_SCALAR_BC_H_
#define ERF_WALL_SCALAR_BC_H_

#include <AMReX_REAL.H>

#include <string>

namespace erf_wall_scalar_bc {

// Host-side intent for a scalar supplied on one physical solid-wall face.
// Presence is represented by this value, never by the numerical value itself.
enum class WallScalarBCIntent : unsigned char {
    Unspecified,
    DirichletConserved,
    DirichletPrimitive,
    Neumann
};

struct ParsedWallScalarBC {
    WallScalarBCIntent intent = WallScalarBCIntent::Unspecified;
    amrex::Real stored_value = amrex::Real(0.0);
};

inline ParsedWallScalarBC
resolve_wall_scalar_value (bool value_specified,
                           amrex::Real primitive_value,
                           bool density_specified_on_same_face,
                           amrex::Real density_value)
{
    if (!value_specified) {
        return {};
    }

    if (density_specified_on_same_face) {
        return {WallScalarBCIntent::DirichletConserved,
                density_value * primitive_value};
    }

    return {WallScalarBCIntent::DirichletPrimitive, primitive_value};
}

inline ParsedWallScalarBC
resolve_wall_scalar_gradient (bool gradient_specified,
                              amrex::Real gradient)
{
    if (!gradient_specified) {
        return {};
    }
    return {WallScalarBCIntent::Neumann, gradient};
}

inline bool
wall_scalar_value_gradient_conflict (bool value_specified,
                                     bool gradient_specified)
{
    return value_specified && gradient_specified;
}

inline std::string
wall_scalar_value_gradient_error (const std::string& face,
                                  const std::string& variable)
{
    return face + "." + variable + " and " + face + "." + variable
        + "_grad are mutually exclusive";
}

inline std::string
anelastic_wall_density_error (const std::string& face,
                              bool density_specified,
                              bool density_gradient_specified)
{
    if (!density_specified && !density_gradient_specified) {
        return {};
    }

    std::string message = face + ": density and density_grad must be omitted "
        "on solid walls in anelastic mode (density is fixed by the base state)";
    if (density_specified && density_gradient_specified) {
        message += "; both keywords were supplied";
    } else if (density_specified) {
        message += "; " + face + ".density was supplied";
    } else {
        message += "; " + face + ".density_grad was supplied";
    }
    return message;
}

} // namespace erf_wall_scalar_bc

#endif
