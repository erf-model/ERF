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
 * Resets φ to be a signed-distance function via Sussman's method.
 */

namespace fire_levelset {

/**
 * @brief Reinitialize φ to be a signed-distance function
 *
 * Solves the reinitialization PDE:
 *   ∂φ/∂τ = sgn(φ₀) * (1 - |∇φ|)
 *
 * using a pseudo-time iterate with explicit Euler or Runge-Kutta stepping.
 *
 * @param[in,out] phi      Level-set field; modified in-place
 * @param[in]     geom     Domain geometry
 * @param[in]     n_iters  Number of Sussman iterations
 * @param[in]     dtau     Pseudo-timestep for reinitialization
 */
inline void reinitialize_phi(amrex::MultiFab& phi,
                              const amrex::Geometry& geom,
                              int n_iters,
                              amrex::Real dtau)
{
    using namespace amrex;
    
    const auto* dx = geom.CellSize();
    const Real dxi = 1.0 / dx[0];
    const Real dyi = 1.0 / dx[1];
    
    // Make a working copy for pseudo-time iteration
    MultiFab phi_work(phi.boxArray(), phi.DistributionMap(), 1, 3);
    phi_work.ParallelCopy(phi);
    
    for (int iter = 0; iter < n_iters; ++iter) {
        // Fill boundary ghosts
        phi_work.FillBoundary(geom.periodicity());
        
        // Update RHS via Sussman reinitialization
        for (MFIter mfi(phi_work); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto phi_arr = phi_work.array(mfi);
            auto phi_in  = phi.const_array(mfi);
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept {
                int i = iv[0], j = iv[1], k = iv[2];
                
                // Central differences for gradient
                Real dxm = (phi_arr(i, j, k) - phi_arr(i-1, j, k)) * dxi;
                Real dxp = (phi_arr(i+1, j, k) - phi_arr(i, j, k)) * dxi;
                Real dym = (phi_arr(i, j, k) - phi_arr(i, j-1, k)) * dyi;
                Real dyp = (phi_arr(i, j+1, k) - phi_arr(i, j, k)) * dyi;
                
                // Use upwind Godunov for norm
                Real phi_x = (dxm * dxm > dxp * dxp) ? dxm : dxp;
                Real phi_y = (dym * dym > dyp * dyp) ? dym : dyp;
                Real norm_grad = std::sqrt(phi_x * phi_x + phi_y * phi_y + 1.0e-14);
                
                // Sign of original level-set
                Real sgn_phi0 = (phi_in(i, j, k) > 0.0) ? 1.0 :
                                (phi_in(i, j, k) < 0.0) ? -1.0 : 0.0;
                
                // Sussman RHS: sgn(φ₀) * (1 - |∇φ|)
                Real rhs = sgn_phi0 * (1.0_rt - norm_grad);
                
                // Forward Euler step
                phi_arr(i, j, k) = phi_arr(i, j, k) + dtau * rhs;
            });
        }
    }
    
    // Copy result back to phi
    phi.ParallelCopy(phi_work);
}

} // namespace fire_levelset

#endif
