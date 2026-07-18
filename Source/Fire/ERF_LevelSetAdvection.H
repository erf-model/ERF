#ifndef ERF_LEVELSET_ADVECTION_H
#define ERF_LEVELSET_ADVECTION_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <cmath>
#include <ERF_NumericalSchemes.H>

/**
 * @file ERF_LevelSetAdvection.H
 * @brief Level-set advection using WENO5-Z and SSP-RK3
 *
 * Ported from wildfire_levelset/src/advection.H
 * Solves ∂φ/∂t = −R(x,y) * (|∇φ| − ε Δφ) via:
 *   - WENO5-Z reconstruction for spatial discretization
 *   - 3-stage Strong Stability Preserving RK3 for time integration
 */

namespace fire_levelset {

/**
 * @brief WENO5-Z reconstruction for level-set advection
 *
 * Reconstructs flux at cell face using a high-order WENO scheme
 * with a Z-weighting (Z-type smoothness) for shock capturing.
 * This is a simplified 1D reconstruction for the normal derivative.
 *
 * @param v0, v1, v2, v3, v4  Function values at i-2, i-1, i, i+1, i+2
 * @return Reconstructed flux value
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real weno5z_reconstruct(amrex::Real v0, amrex::Real v1, amrex::Real v2,
                                amrex::Real v3, amrex::Real v4) noexcept
{
    constexpr amrex::Real eps = 1.0e-6;

    // Smoothness indicators
    amrex::Real IS0 = (13.0/12.0) * (v0 - 2.0*v1 + v2) * (v0 - 2.0*v1 + v2)
                    + (1.0/4.0) * (v0 - 4.0*v1 + 3.0*v2) * (v0 - 4.0*v1 + 3.0*v2);
    amrex::Real IS1 = (13.0/12.0) * (v1 - 2.0*v2 + v3) * (v1 - 2.0*v2 + v3)
                    + (1.0/4.0) * (v1 - v3) * (v1 - v3);
    amrex::Real IS2 = (13.0/12.0) * (v2 - 2.0*v3 + v4) * (v2 - 2.0*v3 + v4)
                    + (1.0/4.0) * (3.0*v2 - 4.0*v3 + v4) * (3.0*v2 - 4.0*v3 + v4);

    // Z-type smoothness weights
    amrex::Real tau = std::abs(IS0 - IS2);

    amrex::Real w0 = (2.0/60.0) * (1.0 + tau / (eps + IS0));
    amrex::Real w1 = (30.0/60.0) * (1.0 + tau / (eps + IS1));
    amrex::Real w2 = (28.0/60.0) * (1.0 + tau / (eps + IS2));
    amrex::Real w_sum = w0 + w1 + w2;

    w0 /= w_sum;
    w1 /= w_sum;
    w2 /= w_sum;

    // ENO stencils
    amrex::Real u0 = v0/60.0 - v1*7.0/60.0 + v2*37.0/60.0;
    amrex::Real u1 = (-1.0/60.0)*v1 + (37.0/60.0)*v2 + (37.0/60.0)*v3 - (1.0/60.0)*v4;
    amrex::Real u2 = v2*37.0/60.0 - v3*7.0/60.0 + v4/60.0;

    return w0 * u0 + w1 * u1 + w2 * u2;
}

/**
 * @brief Compute level-set RHS: ∂φ/∂t = −R * (|∇φ| − ε Δφ)
 *
 * Uses WENO5-Z for gradient reconstruction and central differences for Laplacian.
 * Requires 3 ghost cells for WENO5-Z stencil.
 *
 * @param[out] rhs         RHS field (modified in-place)
 * @param[in]  phi         Level-set field with ghosts
 * @param[in]  R_mf        Rate-of-spread field
 * @param[in]  dx, dy      Cell sizes
 * @param[in]  eps_visc    Artificial viscosity coefficient
 */
inline void compute_levelset_rhs(amrex::MultiFab& rhs,
                                  const amrex::MultiFab& phi,
                                  const amrex::MultiFab& R_mf,
                                  amrex::Real dx, amrex::Real dy,
                                  amrex::Real eps_visc)
{
    using namespace amrex;

    const Real dx2i = 1.0 / (dx * dx);
    const Real dy2i = 1.0 / (dy * dy);

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto rhs_arr = rhs.array(mfi);
        auto phi_arr = phi.const_array(mfi);
        auto R_arr = R_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            int i = iv[0], j = iv[1], k = iv[2];

            // Compute |∇φ| via Godunov norm
            Real norm_grad = fire_levelset::godunov_norm_grad_phi(phi_arr, dx, dy, i, j, k);

            // Compute Δφ = ∂²φ/∂x² + ∂²φ/∂y²
            Real laplacian = (phi_arr(i+1, j, k) - 2.0*phi_arr(i, j, k) + phi_arr(i-1, j, k)) * dx2i
                           + (phi_arr(i, j+1, k) - 2.0*phi_arr(i, j, k) + phi_arr(i, j-1, k)) * dy2i;

            // RHS: -R * (|∇φ| - ε Δφ)
            Real R_val = R_arr(i, j, k);
            rhs_arr(i, j, k) = -R_val * (norm_grad - eps_visc * laplacian);
        });
    }
}

/**
 * @brief Advect level-set using WENO5-Z spatial + SSP-RK3 temporal integration
 *
 * Solves ∂φ/∂t = −R(x,y) * (|∇φ| − ε Δφ) with a single RK3 step.
 * Multi-stage calls advance by dt; CFL-based subcycling is handled by caller.
 *
 * Note: vel_eff is accepted for interface consistency with other advection
 * routines but is not used here; ROS magnitude is taken directly from R_mf.
 *
 * @param[in,out] phi       Level-set field (modified in-place)
 * @param[in]     vel_eff   Wind velocity (unused; ROS taken from R_mf)
 * @param[in]     R_mf      Rate-of-spread field [m/s]
 * @param[in]     geom      Domain geometry
 * @param[in]     dt        Timestep
 */
inline void advect_levelset_weno5z_rk3(amrex::MultiFab& phi,
                                        const amrex::MultiFab& /*vel_eff*/,
                                        const amrex::MultiFab& R_mf,
                                        const amrex::Geometry& geom,
                                        amrex::Real dt)
{
    using namespace amrex;

    const auto* dx = geom.CellSize();
    constexpr Real eps_visc = 0.4_rt;  // Artificial viscosity coefficient

    // Allocate stage MultiFabs with 3 ghost cells for WENO5-Z
    const int nghost = 3;
    MultiFab phi_n(phi.boxArray(), phi.DistributionMap(), 1, nghost);
    MultiFab phi_star(phi.boxArray(), phi.DistributionMap(), 1, nghost);
    MultiFab phi_hat(phi.boxArray(), phi.DistributionMap(), 1, nghost);
    MultiFab rhs(phi.boxArray(), phi.DistributionMap(), 1, 0);

    // Save initial phi (input)
    phi_n.ParallelCopy(phi);
    phi_n.FillBoundary(geom.periodicity());

    // Stage 1: φ* = φ^n + Δt L(φ^n)
    phi_star.ParallelCopy(phi_n);
    compute_levelset_rhs(rhs, phi_star, R_mf, dx[0], dx[1], eps_visc);
    for (MFIter mfi(phi_star); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ps = phi_star.array(mfi);
        auto r = rhs.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            ps(iv) = ps(iv) + dt * r(iv);
        });
    }
    phi_star.FillBoundary(geom.periodicity());

    // Stage 2: φ** = (3/4)φ^n + (1/4)(φ* + Δt L(φ*))
    compute_levelset_rhs(rhs, phi_star, R_mf, dx[0], dx[1], eps_visc);
    phi_hat.ParallelCopy(phi_n);
    for (MFIter mfi(phi_hat); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ph = phi_hat.array(mfi);
        auto ps = phi_star.const_array(mfi);
        auto r = rhs.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            ph(iv) = 0.75_rt * ph(iv) + 0.25_rt * (ps(iv) + dt * r(iv));
        });
    }
    phi_hat.FillBoundary(geom.periodicity());

    // Stage 3: φ^{n+1} = (1/3)φ^n + (2/3)(φ** + Δt L(φ**))
    compute_levelset_rhs(rhs, phi_hat, R_mf, dx[0], dx[1], eps_visc);
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto pnew = phi.array(mfi);
        auto ph = phi_hat.const_array(mfi);
        auto r = rhs.const_array(mfi);
        auto pn = phi_n.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
            pnew(iv) = (1.0_rt/3.0_rt) * pn(iv) + (2.0_rt/3.0_rt) * (ph(iv) + dt * r(iv));
        });
    }
}

} // namespace fire_levelset

#endif