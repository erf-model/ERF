#ifndef ERF_WALL_GRADIENT_COEFFS_H_
#define ERF_WALL_GRADIENT_COEFFS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include "ERF_NumericalConstants.H"

/**
 * Coefficients of the one-sided z-gradient of a cell-centered velocity at a Dirichlet
 * (no-slip) z wall,
 *
 *   du/dz|wall = (c_wall * u_wall + c_0 * u_0 + c_1 * u_1) / dz0,
 *
 * where u_wall is the wall value, u_0 and u_1 are the first and second cells off the wall,
 * and dz0 and dz1 are their heights.  This is the stencil ComputeStrain_S and
 * ComputeStrain_T use for the wall stress, with the same arithmetic; on a uniform mesh it
 * reduces to (-8/3, 3, -1/3), the stencil of ComputeStrain_N.  It is second order, so it
 * is exact for a quadratic profile.  The coefficients sum to zero.
 *
 * The implicit vertical momentum solve closes its wall row with this same gradient: the
 * solve removes the explicit wall stress (from tau_corr) and puts back its own, so any
 * other wall stencil changes the operator by implicit_fac * (L_solve - L_explicit) at the
 * wall, and moves the steady state of a wall-bounded flow.
 */
struct WallGradCoeffs {
    amrex::Real c_wall;
    amrex::Real c_0;
    amrex::Real c_1;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
WallGradCoeffs
WallGradientCoeffs (const amrex::Real dz0, const amrex::Real dz1) noexcept
{
    const amrex::Real f  = (dz1 / dz0) + two;
    const amrex::Real f2 = f*f;
    const amrex::Real c3 = two / (f - f2);
    const amrex::Real c2 = -f2*c3;
    const amrex::Real c1 = -(one-f2)*c3;
    return {c1, c2, c3};
}

/** The uniform-mesh stencil, with the constants ComputeStrain_N uses. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
WallGradCoeffs
WallGradientCoeffsUniform () noexcept
{
    return {-eightthirds, three, -third};
}

#endif
