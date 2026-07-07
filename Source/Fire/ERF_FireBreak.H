#ifndef ERF_FIREBREAK_H
#define ERF_FIREBREAK_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <string>
#include <vector>
#include <cmath>

/**
 * @file ERF_FireBreak.H
 * @brief Firebreak barriers for fire simulation
 *
 * Ported from wildfire_levelset/src/barrier_polygons.H
 *
 * Implements rectangular and circular non-burnable zones that are stamped
 * into the level-set phi at initialization. Firebreak cells remain non-burnable
 * throughout the simulation by setting phi to FIREBREAK_PHI_SENTINEL.
 *
 * Reference:
 *   Finney, M.A. (1998). FARSITE: Fire Area Simulator. RMRS-RP-4.
 */

/**
 * @brief Sentinel value for firebreak cells in fire_phi.
 *
 * Must be larger than any phi value produced during FARSITE propagation
 * (typically φ ∈ [-large, +large]) and larger than farsite_phi_threshold
 * (default 0.1).
 */
constexpr amrex::Real FIREBREAK_PHI_SENTINEL = 1.0e6_rt;

/**
 * @brief A firebreak is a rectangular or circular non-burnable zone.
 *
 * Cells within the barrier have phi set to FIREBREAK_PHI_SENTINEL at
 * initialization and remain non-burnable for the duration of the simulation.
 *
 * Reference:
 *   Finney, M.A. (1998). FARSITE: Fire Area Simulator. RMRS-RP-4.
 */
struct FireBreak {
    std::string type   = "rect";  ///< "rect" or "circle"
    amrex::Real x_lo   = 0.0;     ///< rectangle: west edge [m]
    amrex::Real y_lo   = 0.0;     ///< rectangle: south edge [m]
    amrex::Real x_hi   = 0.0;     ///< rectangle: east edge [m]
    amrex::Real y_hi   = 0.0;     ///< rectangle: north edge [m]
    amrex::Real cx     = 0.0;     ///< circle: centre x [m]
    amrex::Real cy     = 0.0;     ///< circle: centre y [m]
    amrex::Real radius = 0.0;     ///< circle: radius [m]
};

/**
 * @brief Apply firebreak barriers to fire_phi using a GPU ParallelFor kernel.
 *
 * For each cell whose centre falls within any barrier:
 *   phi(i,j,0) = FIREBREAK_PHI_SENTINEL
 *
 * Cells outside all barriers are not modified.
 *
 * Called once in FireLayer::initialize() after ignition stamping.
 *
 * Cell centre: x = prob_lo[0] + (i + 0.5) * dx[0],
 *              y = prob_lo[1] + (j + 0.5) * dx[1]
 *
 * Rectangle check: x_lo <= x <= x_hi && y_lo <= y <= y_hi
 * Circle check: sqrt((x-cx)^2 + (y-cy)^2) <= radius
 *
 * @param[inout] fire_phi      Level-set phi to stamp
 * @param[in]    breaks        Vector of firebreak barriers
 * @param[in]    geom          Fire grid geometry
 */
inline void apply_firebreaks(amrex::MultiFab& fire_phi,
                              const std::vector<FireBreak>& breaks,
                              const amrex::Geometry& geom)
{
    if (breaks.empty()) { return; }

    // Copy firebreaks to device vectors for GPU access
    int n_breaks = static_cast<int>(breaks.size());

    // Allocate device vectors for barrier parameters
    amrex::Gpu::DeviceVector<std::string> d_types(n_breaks);
    amrex::Gpu::DeviceVector<amrex::Real> d_x_lo(n_breaks);
    amrex::Gpu::DeviceVector<amrex::Real> d_y_lo(n_breaks);
    amrex::Gpu::DeviceVector<amrex::Real> d_x_hi(n_breaks);
    amrex::Gpu::DeviceVector<amrex::Real> d_y_hi(n_breaks);
    amrex::Gpu::DeviceVector<amrex::Real> d_cx(n_breaks);
    amrex::Gpu::DeviceVector<amrex::Real> d_cy(n_breaks);
    amrex::Gpu::DeviceVector<amrex::Real> d_radius(n_breaks);

    // Host copies for MPI transfer
    std::vector<std::string> h_types(n_breaks);
    std::vector<amrex::Real> h_x_lo(n_breaks);
    std::vector<amrex::Real> h_y_lo(n_breaks);
    std::vector<amrex::Real> h_x_hi(n_breaks);
    std::vector<amrex::Real> h_y_hi(n_breaks);
    std::vector<amrex::Real> h_cx(n_breaks);
    std::vector<amrex::Real> h_cy(n_breaks);
    std::vector<amrex::Real> h_radius(n_breaks);

    // Copy data to host vectors
    for (int nb = 0; nb < n_breaks; ++nb) {
        h_types[nb] = breaks[nb].type;
        h_x_lo[nb] = breaks[nb].x_lo;
        h_y_lo[nb] = breaks[nb].y_lo;
        h_x_hi[nb] = breaks[nb].x_hi;
        h_y_hi[nb] = breaks[nb].y_hi;
        h_cx[nb] = breaks[nb].cx;
        h_cy[nb] = breaks[nb].cy;
        h_radius[nb] = breaks[nb].radius;
    }

    // Copy to device
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_types.begin(), h_types.end(),
                     d_types.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_x_lo.begin(), h_x_lo.end(), d_x_lo.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_y_lo.begin(), h_y_lo.end(), d_y_lo.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_x_hi.begin(), h_x_hi.end(), d_x_hi.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_y_hi.begin(), h_y_hi.end(), d_y_hi.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_cx.begin(), h_cx.end(), d_cx.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_cy.begin(), h_cy.end(), d_cy.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     h_radius.begin(), h_radius.end(), d_radius.begin());

    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* dx = geom.CellSize();

    const auto& d_types_p = d_types.data();
    const auto& d_x_lo_p = d_x_lo.data();
    const auto& d_y_lo_p = d_y_lo.data();
    const auto& d_x_hi_p = d_x_hi.data();
    const auto& d_y_hi_p = d_y_hi.data();
    const auto& d_cx_p = d_cx.data();
    const auto& d_cy_p = d_cy.data();
    const auto& d_radius_p = d_radius.data();

    // ParallelFor over fire grid
    for (amrex::MFIter mfi(fire_phi); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto& phi_arr = fire_phi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            amrex::Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];

            // Check all barriers
            for (int nb = 0; nb < n_breaks; ++nb) {
                bool in_barrier = false;

                // Check barrier type
                if (d_types_p[nb][0] == 'r') {  // "rect"
                    in_barrier = (x >= d_x_lo_p[nb] && x <= d_x_hi_p[nb] &&
                                 y >= d_y_lo_p[nb] && y <= d_y_hi_p[nb]);
                } else if (d_types_p[nb][0] == 'c') {  // "circle"
                    amrex::Real dx_c = x - d_cx_p[nb];
                    amrex::Real dy_c = y - d_cy_p[nb];
                    amrex::Real dist_sq = dx_c * dx_c + dy_c * dy_c;
                    in_barrier = (dist_sq <= d_radius_p[nb] * d_radius_p[nb]);
                }

                // Stamp firebreak sentinel if cell is in barrier
                if (in_barrier) {
                    phi_arr(i, j, k) = FIREBREAK_PHI_SENTINEL;
                }
            }
        });
    }
}

#endif  // ERF_FIREBREAK_H
