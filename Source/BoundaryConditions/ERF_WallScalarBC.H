#ifndef ERF_WALL_SCALAR_BC_H_
#define ERF_WALL_SCALAR_BC_H_

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

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

enum class SolidWallKind : unsigned char {
    NoSlip,
    Slip
};

struct ParsedWallScalarBC {
    WallScalarBCIntent intent = WallScalarBCIntent::Unspecified;
    amrex::Real stored_value = amrex::Real(0.0);
};

struct ParsedWallFaceScalars {
    ParsedWallScalarBC density;
    ParsedWallScalarBC theta;
    ParsedWallScalarBC qv;
};

struct WallFaceParseResult {
    ParsedWallFaceScalars scalars;
    std::string error;

    /**
     * Return whether parsing completed without an error.
     */
    [[nodiscard]] bool ok () const noexcept { return error.empty(); }
};

/**
 * Resolve a wall scalar value into primitive or conserved storage.
 *
 * @param[in] value_specified whether the scalar value was specified
 * @param[in] primitive_value primitive scalar value from inputs
 * @param[in] density_specified_on_same_face whether density was specified on the same face
 * @param[in] density_value density value used for conserved conversion
 */
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

/**
 * Resolve a wall scalar gradient into Neumann storage.
 *
 * @param[in] gradient_specified whether the scalar gradient was specified
 * @param[in] gradient gradient value from inputs
 */
inline ParsedWallScalarBC
resolve_wall_scalar_gradient (bool gradient_specified,
                              amrex::Real gradient)
{
    if (!gradient_specified) {
        return {};
    }
    return {WallScalarBCIntent::Neumann, gradient};
}

/**
 * Test whether value and gradient inputs conflict.
 *
 * @param[in] value_specified whether a scalar value was specified
 * @param[in] gradient_specified whether a scalar gradient was specified
 */
inline bool
wall_scalar_value_gradient_conflict (bool value_specified,
                                     bool gradient_specified)
{
    return value_specified && gradient_specified;
}

/**
 * Build the error message for mutually exclusive wall scalar inputs.
 *
 * @param[in] face face name used in the input keys
 * @param[in] variable variable name used in the input keys
 */
inline std::string
wall_scalar_value_gradient_error (const std::string& face,
                                  const std::string& variable)
{
    return face + "." + variable + " and " + face + "." + variable
        + "_grad are mutually exclusive";
}

/**
 * Build the anelastic-mode density error message for a wall face.
 *
 * @param[in] face face name used in the input keys
 * @param[in] density_specified whether density was specified
 * @param[in] density_gradient_specified whether density_grad was specified
 */
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

/**
 * Build the unsupported no-slip density-gradient error message.
 *
 * @param[in] face face name used in the input keys
 */
inline std::string
no_slip_density_gradient_error (const std::string& face)
{
    return face + ".density_grad is not supported for NoSlipWall; use "
        + face + ".density for a fixed compressible wall density or omit "
        "wall density in anelastic mode";
}

/**
 * Parse scalar wall boundary-condition inputs for one face.
 *
 * @param[in] pp ParmParse object scoped to the wall face
 * @param[in] face face name used in errors
 * @param[in] wall_kind no-slip or slip wall type
 * @param[in] anelastic whether anelastic mode is active
 */
inline WallFaceParseResult
parse_wall_face_scalars (amrex::ParmParse& pp,
                         const std::string& face,
                         SolidWallKind wall_kind,
                         bool anelastic)
{
    amrex::Real density = amrex::Real(0.0);
    const bool density_specified = pp.query("density", density);
    amrex::Real density_grad = amrex::Real(0.0);
    const bool density_gradient_specified = pp.query("density_grad", density_grad);
    amrex::Real theta = amrex::Real(0.0);
    const bool theta_specified = pp.query("theta", theta);
    amrex::Real theta_grad = amrex::Real(0.0);
    const bool theta_gradient_specified = pp.query("theta_grad", theta_grad);

    amrex::Real qv = amrex::Real(0.0);
    const bool qv_specified = (wall_kind == SolidWallKind::NoSlip) &&
                              pp.query("qv", qv);

    WallFaceParseResult result;
    if (anelastic && (density_specified || density_gradient_specified)) {
        result.error = anelastic_wall_density_error(
            face, density_specified, density_gradient_specified);
        return result;
    }

    if (wall_kind == SolidWallKind::NoSlip && density_gradient_specified) {
        result.error = no_slip_density_gradient_error(face);
        return result;
    }

    if (wall_kind == SolidWallKind::Slip &&
        wall_scalar_value_gradient_conflict(density_specified,
                                            density_gradient_specified)) {
        result.error = wall_scalar_value_gradient_error(face, "density");
        return result;
    }

    if (wall_scalar_value_gradient_conflict(theta_specified,
                                            theta_gradient_specified)) {
        result.error = wall_scalar_value_gradient_error(face, "theta");
        return result;
    }

    if (density_specified) {
        result.scalars.density = {
            WallScalarBCIntent::DirichletConserved, density};
    }
    result.scalars.theta = resolve_wall_scalar_value(
        theta_specified, theta, density_specified, density);
    if (theta_gradient_specified) {
        result.scalars.theta = resolve_wall_scalar_gradient(true, theta_grad);
    }

    if (wall_kind == SolidWallKind::Slip && density_gradient_specified) {
        result.scalars.density = resolve_wall_scalar_gradient(true, density_grad);
    }

    if (qv_specified) {
        result.scalars.qv = resolve_wall_scalar_value(
            true, qv, density_specified, density);
    }

    return result;
}

} // namespace erf_wall_scalar_bc

#endif
