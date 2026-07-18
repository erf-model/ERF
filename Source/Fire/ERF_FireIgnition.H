#ifndef ERF_FIRE_IGNITION_H
#define ERF_FIRE_IGNITION_H

/**
 * @file ERF_FireIgnition.H
 * @brief Fire ignition initialization for level-set field
 *
 * Initializes the level-set field phi for a circular ignition point.
 */

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <cmath>

/**
 * @brief Initialize the level-set field phi for a circular ignition point.
 *
 * Sets phi = -(r - dist(x,y, cx, cy)) / r for cells within r metres of (cx, cy),
 * and phi = +1 elsewhere. This gives a smooth signed-distance initialisation
 * within the ignition circle and +1 (unburned) outside.
 *
 * Parameters read from ParmParse("erf.fire"):
 *   erf.fire.ignition_x   [m]
 *   erf.fire.ignition_y   [m]
 *   erf.fire.ignition_r   [m]  (default 20.0)
 *
 * @param[in,out] phi     Level-set field to initialize
 * @param[in]     geom_fire Fire grid geometry
 * @param[in]     cx      Ignition center x [m]
 * @param[in]     cy      Ignition center y [m]
 * @param[in]     radius  Ignition radius [m]
 */
inline
void initialize_ignition(
    amrex::MultiFab& phi,
    const amrex::Geometry& geom_fire,
    amrex::Real cx, amrex::Real cy, amrex::Real radius)
{
    using namespace amrex;
    
    const auto* prob_lo = geom_fire.ProbLo();
    const auto* dx      = geom_fire.CellSize();
    
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
            Real dist = std::sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy));
            
            // Signed distance: negative inside ignition circle, positive outside.
            // Clamped to [-1, 1] to bound the level-set values.
            p(i,j,k) = amrex::max(amrex::min((dist - radius) / radius, 1.0_rt), -1.0_rt);
        });
    }
}

#endif // ERF_FIRE_IGNITION_H
