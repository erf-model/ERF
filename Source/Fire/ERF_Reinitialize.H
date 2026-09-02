#ifndef ERF_REINITIALIZE_H
#define ERF_REINITIALIZE_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>

/**
 * @file ERF_Reinitialize.H
 * @brief Sussman reinitialization of the level-set function
 *
 * Ported from wildfire_levelset/src/reinitialization.H
 * Resets φ to a normalized signed-distance function via Sussman's method.
 *
 * NOTE ON UNITS. The textbook Sussman update drives |∇φ| → 1, which presumes φ
 * is a signed distance measured in metres. ERF-Hazard's φ is not: it is
 * normalized and clamped to [-1, 1] (initialize_ignition sets φ = (d - r)/r,
 * and the FARSITE path and ember spotting use the same convention). Applying
 * the textbook form to a normalized φ leaves the bracket (1 - |∇φ|) near 1 for
 * every cell, so each iteration added ~dtau and the field diverged by
 * n_iters * dtau per call. The update below instead targets |∇φ| = 1/L for an
 * explicit band length L, so φ reaches ±1 at distance L from the front, which
 * is the convention the rest of the module already uses.
 */

namespace fire_levelset {

/**
 * @brief Reinitialize φ to a normalized signed-distance function
 *
 * Solves the band-normalized reinitialization PDE
 *
 *   ∂φ/∂τ = sgn(φ₀) * (1 - L |∇φ|) / L
 *
 * by explicit Euler in pseudo-time. The fixed point is |∇φ| = 1/L, i.e. φ
 * varies linearly from 0 at the front to ±1 at distance L, matching the
 * normalization initialize_ignition applies. φ is clamped to [-1, 1] after
 * every iteration, which pins the far field (where |∇φ| = 0 and the update
 * would otherwise push φ past ±1 indefinitely).
 *
 * Cells whose original neighbourhood straddles the interface use the
 * Russo & Smereka (2000) subcell update instead, which fixes the interface
 * location from φ₀ rather than letting the pseudo-time iteration move it.
 * Without it each pass erodes the burned area, and because the level-set path
 * never rebuilds φ from fire_arrival_time the loss is permanent.
 *
 * Each iteration moves φ by at most dtau/L, so with the default
 * dtau = 0.5*min(dx,dy) the change per iteration is half a cell measured in
 * band lengths.
 *
 * @param[in,out] phi      Level-set field; modified in-place
 * @param[in]     geom     Domain geometry
 * @param[in]     n_iters  Number of Sussman iterations
 * @param[in]     dtau     Pseudo-timestep [m] for reinitialization
 * @param[in]     band_m   Band half-width L [m]: the distance from the front at
 *                         which φ reaches ±1. Values <= 0 select 3*min(dx,dy).
 */
inline void reinitialize_phi(amrex::MultiFab& phi,
                              const amrex::Geometry& geom,
                              int n_iters,
                              amrex::Real dtau,
                              amrex::Real band_m = -1.0,
                              bool normalized = true)
{
    using namespace amrex;
    
    const auto* dx = geom.CellSize();
    const Real dxi = 1.0 / dx[0];
    const Real dyi = 1.0 / dx[1];

    // Band half-width. The default spreads the [-1, 1] range over three cells,
    // the usual interface thickness for a level-set band.
    // In signed-distance mode phi is already in metres, so the band length is
    // 1 m: the update targets |grad phi| = 1 and nothing is clamped. In
    // normalized mode phi is dimensionless and reaches +/-1 at distance L.
    const Real L    = normalized ? ((band_m > 0.0) ? band_m
                                                  : 3.0 * amrex::min(dx[0], dx[1]))
                                 : 1.0;
    const Real Linv = 1.0 / L;

    // Isotropic spacing for the Russo-Smereka subcell term. The fire grid is
    // square in practice; min() is the conservative choice if it ever is not.
    const Real h    = amrex::min(dx[0], dx[1]);
    const Real hinv = 1.0 / h;
    
    // Make a working copy for pseudo-time iteration
    MultiFab phi_work(phi.boxArray(), phi.DistributionMap(), 1, 3);
    phi_work.ParallelCopy(phi);

    // Jacobi scratch. The update below reads phi_work at i±1 / j±1 while writing
    // the centre cell, so writing back into phi_work in place would be a
    // read-write race between neighbouring GPU threads -- each cell would see
    // either the old or the new neighbour depending on scheduling, making the
    // result non-deterministic and different from the CPU sweep. Write into a
    // separate fab and swap once the pass is complete.
    MultiFab phi_next(phi.boxArray(), phi.DistributionMap(), 1, 3);
    // Freshly allocated fabs hold indeterminate values. Only the valid region is
    // written and only the valid region is copied back, but leaving the ghosts
    // indeterminate makes the run depend on whatever the allocator hands over.
    phi_next.setVal(0.0);

    for (int iter = 0; iter < n_iters; ++iter) {
        // Fill boundary ghosts
        phi_work.FillBoundary(geom.periodicity());
        
        // Update RHS via Sussman reinitialization
        for (MFIter mfi(phi_work); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto phi_arr = phi_work.const_array(mfi);
            auto phi_out = phi_next.array(mfi);
            auto phi_in  = phi.const_array(mfi);
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
                int i = iv[0], j = iv[1], k = iv[2];
                
                // Central differences for gradient
                Real dxm = (phi_arr(i, j, k) - phi_arr(i-1, j, k)) * dxi;
                Real dxp = (phi_arr(i+1, j, k) - phi_arr(i, j, k)) * dxi;
                Real dym = (phi_arr(i, j, k) - phi_arr(i, j-1, k)) * dyi;
                Real dyp = (phi_arr(i, j+1, k) - phi_arr(i, j, k)) * dyi;
                
                // Sign of original level-set
                Real sgn_phi0 = (phi_in(i, j, k) > 0.0) ? 1.0 :
                                (phi_in(i, j, k) < 0.0) ? -1.0 : 0.0;

                // Godunov upwind norm for the reinitialization Hamiltonian.
                // The one-sided difference is selected by the sign of phi_0, not
                // by which has the larger magnitude: taking the larger magnitude
                // unconditionally is not the Godunov scheme and is unstable near
                // extrema, which showed up as the front flooding the domain once
                // enough pseudo-time iterations were taken.
                Real gx2, gy2;
                if (sgn_phi0 > 0.0_rt) {
                    Real a = amrex::max(dxm, 0.0_rt), b = amrex::min(dxp, 0.0_rt);
                    gx2 = amrex::max(a * a, b * b);
                    Real c = amrex::max(dym, 0.0_rt), d = amrex::min(dyp, 0.0_rt);
                    gy2 = amrex::max(c * c, d * d);
                } else {
                    Real a = amrex::min(dxm, 0.0_rt), b = amrex::max(dxp, 0.0_rt);
                    gx2 = amrex::max(a * a, b * b);
                    Real c = amrex::min(dym, 0.0_rt), d = amrex::max(dyp, 0.0_rt);
                    gy2 = amrex::max(c * c, d * d);
                }
                Real norm_grad = std::sqrt(gx2 + gy2 + 1.0e-14);
                
                // Band-normalized Sussman RHS: sgn(φ₀) * (1 - L|∇φ|) / L.
                // The fixed point is |∇φ| = 1/L rather than 1, because φ here is
                // normalized to [-1, 1] rather than measured in metres.
                Real rhs = sgn_phi0 * (1.0_rt - L * norm_grad) * Linv;

                Real phi_new = phi_arr(i, j, k) + dtau * rhs;

                // Russo-Smereka (2000) subcell fix.
                //
                // The plain Sussman iteration does not preserve the zero level
                // set, so each pass walks the front and permanently erodes the
                // burned area -- the level-set path never rebuilds phi from
                // fire_arrival_time the way the FARSITE path does, so nothing
                // puts it back. For cells whose ORIGINAL neighbourhood straddles
                // the interface, replace the Godunov update with one driven by
                // the distance to the interface implied by phi_0 itself, so the
                // interface location is set by phi_0 and not by the iteration.
                const Real p0c = phi_in(i, j, k);
                const bool near_front =
                    (p0c * phi_in(i+1, j, k) < 0.0_rt) ||
                    (p0c * phi_in(i-1, j, k) < 0.0_rt) ||
                    (p0c * phi_in(i, j+1, k) < 0.0_rt) ||
                    (p0c * phi_in(i, j-1, k) < 0.0_rt);

                if (near_front) {
                    // Robust local |grad phi_0|*h estimate: the central
                    // difference, floored by every one-sided difference so a
                    // kink at the interface cannot make it vanish.
                    const Real dcx = 0.5_rt * (phi_in(i+1, j, k) - phi_in(i-1, j, k));
                    const Real dcy = 0.5_rt * (phi_in(i, j+1, k) - phi_in(i, j-1, k));
                    Real dphi0 = std::sqrt(dcx * dcx + dcy * dcy);
                    dphi0 = amrex::max(dphi0, std::abs(phi_in(i+1, j, k) - p0c));
                    dphi0 = amrex::max(dphi0, std::abs(p0c - phi_in(i-1, j, k)));
                    dphi0 = amrex::max(dphi0, std::abs(phi_in(i, j+1, k) - p0c));
                    dphi0 = amrex::max(dphi0, std::abs(p0c - phi_in(i, j-1, k)));
                    dphi0 = amrex::max(dphi0, 1.0e-12_rt);

                    // Signed distance to the interface [m], then expressed in the
                    // normalized units phi uses.
                    const Real D = h * p0c / dphi0;

                    phi_new = phi_arr(i, j, k)
                            - dtau * hinv * (sgn_phi0 * std::abs(phi_arr(i, j, k)) - D * Linv);
                }

                // Clamp only in normalized mode, where it pins the far field:
                // there |grad phi| = 0 and the update would otherwise drive phi
                // past +/-1 without limit. In signed-distance mode phi is meant
                // to keep growing with distance, so clamping is exactly the
                // thing that breaks it.
                phi_out(i, j, k) = normalized
                    ? amrex::max(-1.0_rt, amrex::min(1.0_rt, phi_new))
                    : phi_new;
            });
        }

        // phi_next now holds the whole updated valid region; make it the input
        // for the next sweep. Ghost cells are refilled at the top of the loop.
        MultiFab::Copy(phi_work, phi_next, 0, 0, 1, 0);
    }
    
    // Copy result back to phi
    phi.ParallelCopy(phi_work);
}

} // namespace fire_levelset

#endif
