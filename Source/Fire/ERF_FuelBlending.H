#ifndef ERF_FUEL_BLENDING_H
#define ERF_FUEL_BLENDING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>

/**
 * @file ERF_FuelBlending.H
 * @brief Fuel boundary blending for multi-fuel simulations
 *
 * Ported from wildfire_levelset/src/two_fuel_blending.H
 *
 * Blend rate-of-spread (ROS) at spatial boundaries between different fuel model zones.
 * For cells at fuel boundaries, ROS is smoothly transitioned to reduce discontinuities
 * in fire spread behavior.
 *
 * Reference:
 *   Finney, M.A. (1998). FARSITE: Fire Area Simulator. RMRS-RP-4, Appendix B.
 */

/**
 * @brief Blend ROS at spatial boundaries between different fuel model zones.
 *
 * For each cell that has a different fuel code from at least one of its
 * 4 neighbours (i±1, j) and (i, j±1):
 *   ros(i,j,0) = (1 - f) * ros(i,j,0) + f * mean(ros of differently-coded neighbours)
 *
 * where f = blending_fraction (typically 0.1–0.5).
 *
 * When blending_fraction = 0 this function returns immediately (no-op).
 *
 * Guard all neighbour accesses against mfi.validbox() bounds.
 *
 * Reference:
 *   Finney, M.A. (1998). FARSITE: Fire Area Simulator. RMRS-RP-4, Appendix B.
 *
 * @param[inout] fire_ros          Rate-of-spread field to blend
 * @param[in]    fire_fuel_model   Per-cell fuel model codes (integer stored as Real)
 * @param[in]    blending_fraction Blending weight [0-1]
 */
inline void apply_fuel_boundary_blending(amrex::MultiFab& fire_ros,
                                          const amrex::MultiFab& fire_fuel_model,
                                          amrex::Real blending_fraction)
{
    // No-op if blending is disabled
    if (blending_fraction <= 0.0_rt) { return; }

    // Clamp blending fraction to [0, 1]
    blending_fraction = amrex::min(blending_fraction, 1.0_rt);

    // ParallelFor over ROS field
    for (amrex::MFIter mfi(fire_ros); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::Box& valid_box = bx;

        auto& ros_arr = fire_ros.array(mfi);
        const auto& fuel_arr = fire_fuel_model.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Get fuel code at this cell
            int fuel_c = static_cast<int>(fuel_arr(i, j, k) + 0.5_rt);

            // Check 4 neighbours: (i-1,j), (i+1,j), (i,j-1), (i,j+1)
            bool is_boundary = false;
            int fuel_i_minus = fuel_c;
            int fuel_i_plus = fuel_c;
            int fuel_j_minus = fuel_c;
            int fuel_j_plus = fuel_c;

            if (valid_box.contains(amrex::IntVect(i - 1, j, k))) {
                fuel_i_minus = static_cast<int>(fuel_arr(i - 1, j, k) + 0.5_rt);
                if (fuel_i_minus != fuel_c) { is_boundary = true; }
            }
            if (valid_box.contains(amrex::IntVect(i + 1, j, k))) {
                fuel_i_plus = static_cast<int>(fuel_arr(i + 1, j, k) + 0.5_rt);
                if (fuel_i_plus != fuel_c) { is_boundary = true; }
            }
            if (valid_box.contains(amrex::IntVect(i, j - 1, k))) {
                fuel_j_minus = static_cast<int>(fuel_arr(i, j - 1, k) + 0.5_rt);
                if (fuel_j_minus != fuel_c) { is_boundary = true; }
            }
            if (valid_box.contains(amrex::IntVect(i, j + 1, k))) {
                fuel_j_plus = static_cast<int>(fuel_arr(i, j + 1, k) + 0.5_rt);
                if (fuel_j_plus != fuel_c) { is_boundary = true; }
            }

            // If at a boundary, blend ROS with differently-coded neighbours
            if (is_boundary) {
                amrex::Real ros_blend = 0.0_rt;
                int n_diff = 0;

                if (valid_box.contains(amrex::IntVect(i - 1, j, k)) && fuel_i_minus != fuel_c) {
                    ros_blend += ros_arr(i - 1, j, k);
                    ++n_diff;
                }
                if (valid_box.contains(amrex::IntVect(i + 1, j, k)) && fuel_i_plus != fuel_c) {
                    ros_blend += ros_arr(i + 1, j, k);
                    ++n_diff;
                }
                if (valid_box.contains(amrex::IntVect(i, j - 1, k)) && fuel_j_minus != fuel_c) {
                    ros_blend += ros_arr(i, j - 1, k);
                    ++n_diff;
                }
                if (valid_box.contains(amrex::IntVect(i, j + 1, k)) && fuel_j_plus != fuel_c) {
                    ros_blend += ros_arr(i, j + 1, k);
                    ++n_diff;
                }

                if (n_diff > 0) {
                    ros_blend /= static_cast<amrex::Real>(n_diff);
                    amrex::Real ros_orig = ros_arr(i, j, k);
                    ros_arr(i, j, k) = (1.0_rt - blending_fraction) * ros_orig +
                                       blending_fraction * ros_blend;
                }
            }
        });
    }
}

#endif  // ERF_FUEL_BLENDING_H
