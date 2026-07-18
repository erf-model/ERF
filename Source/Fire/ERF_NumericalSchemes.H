#ifndef ERF_NUMERICAL_SCHEMES_H
#define ERF_NUMERICAL_SCHEMES_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <cmath>

/**
 * @file ERF_NumericalSchemes.H
 * @brief Numerical schemes for level-set advection: WENO3 and Godunov norm gradient
 *
 * Ported from wildfire_levelset/src/numerical_schemes.H
 */

namespace fire_levelset {

/**
 * @brief WENO-3 reconstruction on grid point (i,j) from left stencil (i-2,i-1,i,i+1)
 *
 * Reconstructs the third-order accurate flux on the right side of cell (i-1),
 * which is the left side of cell i.
 *
 * Note: v0 (value at i-2) is accepted for interface consistency but is not used
 * in the 2-stencil WENO3 formulation.
 *
 * @param v0  Function value at i-2 (unused in current formulation)
 * @param v1  Function value at i-1
 * @param v2  Function value at i
 * @param v3  Function value at i+1
 * @return Reconstructed value on the right of i-1 (left of i)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real weno3_reconstruct(amrex::Real /*v0*/, amrex::Real v1,
                               amrex::Real v2, amrex::Real v3) noexcept
{
    constexpr amrex::Real eps = 1.0e-6;

    // Smoothness indicators
    amrex::Real IS0 = (v2 - v1) * (v2 - v1);
    amrex::Real IS1 = (v3 - v2) * (v3 - v2);

    // Weights
    amrex::Real w0 = 1.0 / (3.0 * (eps + IS0) * (eps + IS0));
    amrex::Real w1 = 2.0 / (3.0 * (eps + IS1) * (eps + IS1));
    amrex::Real w_sum = w0 + w1;

    // Normalize weights
    w0 /= w_sum;
    w1 /= w_sum;

    // ENO stencils
    amrex::Real u0 = (v1 + v2) / 2.0;
    amrex::Real u1 = (v2 + v3) / 2.0;

    // Weighted combination
    return w0 * u0 + w1 * u1;
}

/**
 * @brief Godunov norm of the gradient, with optional terrain arc-length correction
 *
 * Computes |∇φ| using Godunov's scheme on a staggered grid.
 * With slope_x and slope_y parameters, applies terrain arc-length correction:
 *   |∇φ|_corrected = sqrt(φ_x^2 + φ_y^2 + (φ_x * slope_x + φ_y * slope_y)^2)
 *
 * @param[in] phi         Level-set field (assumed with sufficient ghosts for stencil)
 * @param[in] dx, dy      Cell sizes [m]
 * @param[in] i, j, k     Cell indices
 * @param[in] slope_x     Terrain slope dz/dx (default 0)
 * @param[in] slope_y     Terrain slope dz/dy (default 0)
 * @return Godunov norm of the gradient |∇φ|
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real godunov_norm_grad_phi(const amrex::Array4<const amrex::Real>& phi,
                                   amrex::Real dx, amrex::Real dy,
                                   int i, int j, int k,
                                   amrex::Real slope_x = 0.0,
                                   amrex::Real slope_y = 0.0) noexcept
{
    // One-sided derivatives (Godunov scheme)
    amrex::Real dxm = (phi(i,   j, k) - phi(i-1, j, k)) / dx;
    amrex::Real dxp = (phi(i+1, j, k) - phi(i,   j, k)) / dx;
    amrex::Real dym = (phi(i, j,   k) - phi(i, j-1, k)) / dy;
    amrex::Real dyp = (phi(i, j+1, k) - phi(i, j,   k)) / dy;

    // Take max of upwind/downwind (Godunov)
    amrex::Real phi_x = (dxm * dxm > dxp * dxp) ? dxm : dxp;
    amrex::Real phi_y = (dym * dym > dyp * dyp) ? dym : dyp;

    // Godunov norm
    amrex::Real norm_sq = phi_x * phi_x + phi_y * phi_y;

    // Terrain arc-length correction
    if (slope_x != 0.0 || slope_y != 0.0) {
        amrex::Real slope_term = phi_x * slope_x + phi_y * slope_y;
        norm_sq += slope_term * slope_term;
    }

    return std::sqrt(amrex::max(norm_sq, 1.0e-14));
}

} // namespace fire_levelset

#endif