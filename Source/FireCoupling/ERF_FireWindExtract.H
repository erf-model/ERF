#ifndef ERF_FIRE_WIND_EXTRACT_H
#define ERF_FIRE_WIND_EXTRACT_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <cmath>
#include <ERF_FireGrid.H>

/**
 * @file ERF_FireWindExtract.H
 * @brief Atmospheric wind sampled onto the fire grid
 *
 * The atmosphere side of the fire coupling: the wind of the atmospheric
 * columns at a height above ground (a uniform height, or one per fire cell
 * that the fire engine supplies), interpolated linearly in height and
 * bilinearly or nearest between the surrounding columns, with an optional
 * two-height log-law resampling; no MOST reconstruction is involved. A fire
 * engine may read its wind through it. What the engine does with the wind
 * afterwards (wind adjustment factors, terrain corrections) is the engine's.
 */


/**
 * @brief Extract wind from direct vertical interpolation on fire grid
 *
 * Maps atmospheric wind to fire grid at reference height by interpolating
 * vertically through the atmospheric column at the fire cell location.
 * Uses cell-center wind values (averaged from face-staggered data) and
 * linearly interpolates to the target height above terrain.
 *
 * When z_ref_cell_mf is non-null, each fire cell is sampled at its own height
 * from that field (for example a per-fuel wind height, as in the WRF-SFIRE
 * fcwh convention) instead of the global z_ref.
 *
 * @param[out] fire_wind_ref   u,v at z_ref on fire grid [m/s]
 * @param[out] fire_wind_extract_z_mf  Wind extraction height z_target [m]
 * @param[in]  xvel_mf         x-velocity (face-staggered at i+1/2, j, k)
 * @param[in]  yvel_mf         y-velocity (face-staggered at i, j+1/2, k)
 * @param[in]  z_phys_cc_mf    Cell-center terrain height z_phys_cc[0] on atm grid
 * @param[in]  fg              Fire grid
 * @param[in]  z_ref           Global fallback wind reference height [m]
 * @param[in]  nz              Number of vertical levels
 * @param[in]  z_ref_cell_mf   Height above ground [m] to sample each fire cell
 *                             at, on the fire grid; nullptr = the global z_ref
 * @param[in]  below_first_log Target below the lowest cell centre: false takes the
 *                             lowest centre's wind as it is; true scales it
 *                             with a neutral log profile of roughness z0_log,
 *                             U(z) = U(z_1) ln(z/z0)/ln(z_1/z0)
 *
 * Algorithm per fire cell (i_f, j_f):
 * 1. Map to atmospheric column: i_a = i_f / C, j_a = j_f / C
 * 2. Get ground elevation: z_surf = fire_surface_z(i_f, j_f), the average of the
 *    four k = 0 nodes of atmospheric cell (i_a, j_a). Using z_phys_cc(i_a, j_a, 0)
 *    instead would sit half an atmospheric cell above the ground.
 * 3. Lookup wind height:
 *    - If z_ref_cell_mf != nullptr: z_ref_cell = z_ref_cell_mf(i_f, j_f)
 *    - Otherwise: z_ref_cell = z_ref (global)
 * 4. Compute target height: z_target = z_surf + z_ref_cell
 * 5. Bisection on the column's z_phys_cc profile for k_lo such that
 *    z_phys_cc(i_a, j_a, k_lo) <= z_target < z_phys_cc(i_a, j_a, k_lo+1),
 *    clamped to [0, nz-2]; a target above the top takes the highest centre,
 *    a target below the first centre the lowest centre's wind scaled by the
 *    log profile (below_first_log) or as it is
 * 6. Compute interpolation weight:
 *    alpha = (z_target - z_phys_cc(i_a, j_a, k_lo)) /
 *            (z_phys_cc(i_a, j_a, k_lo+1) - z_phys_cc(i_a, j_a, k_lo))
 *    Clamp alpha to [0, 1]
 * 7. Average u/v from bracketing face values to get cell centers:
 *    u_cc(i,j,k) = 0.5 * (xvel(i,j,k) + xvel(i+1,j,k))
 *    v_cc(i,j,k) = 0.5 * (yvel(i,j,k) + yvel(i,j+1,k))
 * 8. Interpolate: fw(i_f,j_f,0,0) = u_cc(k_lo) + alpha*(u_cc(k_lo+1) - u_cc(k_lo))
 *                fw(i_f,j_f,0,1) = v_cc(k_lo) + alpha*(v_cc(k_lo+1) - v_cc(k_lo))
 * 9. Store z_target in fire_wind_extract_z_mf
 *
 * With wind_interp = 1 (bilinear, the default) steps 3-8 run for each of the
 * four atmospheric columns surrounding the fire cell, every one sampled at the
 * reference height above *its own* ground, and the results are blended with
 * bilinear weights. With wind_interp = 0 only the column containing the cell is
 * used, which makes the fire-grid wind piecewise constant on atmospheric cells.
 *
 * With col_open and col_roof given (2D atmospheric-column fields, one ghost
 * cell) each bilinear weight is multiplied by the column's open fraction when
 * its roof is above the cell's wind height, and the four weights are
 * renormalised, so a column inside a building does not contribute its
 * relaxed in-building velocity to the fire cells around it. A cell whose four
 * columns are all closed keeps the plain weights; it is inside a footprint
 * and masked anyway.
 *
 * With z_sample > 0 (an engine input) every column is sampled at
 * z_sample above its ground instead of at the target height, and the result
 * is brought down to the target height with a neutral log profile of
 * roughness z0_log: U(z_ref) = U(z_sample) ln(z_ref / z0) / ln(z_sample / z0).
 * This is the second interpolation option of the Community Fire Behavior
 * Model (Jimenez y Munoz et al. 2026), meant to drive the fire with a wind
 * taken above its own plume. fire_wind_extract_z then records the sampled
 * height; the open-column test uses that height as well.
 */
void fill_fire_wind_from_interpolation(
    amrex::MultiFab&       fire_wind_ref,
    amrex::MultiFab&       fire_wind_extract_z_mf,
    const amrex::MultiFab& xvel_mf,
    const amrex::MultiFab& yvel_mf,
    const amrex::MultiFab& z_phys_cc_mf,
    const amrex::MultiFab& fire_surface_z_mf,
    const amrex::MultiFab& fire_col_ground_mf,
    const FireGrid&        fg,
    amrex::Real            z_ref,
    int                    nz,
    const amrex::MultiFab* z_ref_cell_mf  = nullptr,
    int                    wind_interp    = 1,
    const amrex::MultiFab* col_open       = nullptr,
    const amrex::MultiFab* col_roof       = nullptr,
    amrex::Real            z_sample       = 0.0,
    amrex::Real            z0_log         = 0.1,
    bool                   below_first_log = false);

/**
 * @brief The wind at the lowest cell centre of each atmospheric column, and
 *        its height above that column's ground
 *
 * For an engine that drives its fire with the k = 0 wind of the atmospheric
 * columns rather than a wind interpolated to a height (LSFIRE does). The face
 * velocities are averaged to the k = 0 cell centre; the height is that centre
 * minus the column's ground, the mean of its four k = 0 nodes (the ground
 * compute_fire_surface_height uses), or minus the domain floor without nodal
 * heights.
 *
 * @param[out] u_k0, v_k0  Wind [m/s] on the k = 0 slabs of the level (one per box of
 *                         the level, the BoxArray of FireForcing::T_k0), with any
 *                         ghost cells: those are filled with FillBoundary
 * @param[out] h_k0        Height of that wind above the column's ground [m], laid out as u_k0
 * @param[in]  xvel, yvel  Face velocities of the level (FireForcing::xvel, yvel)
 * @param[in]  z_phys_cc   Cell-centre height of the level [m]
 * @param[in]  z_phys_nd   Nodal height of the level [m], or nullptr (flat floor)
 * @param[in]  geom        Geometry of the level
 */
void extract_fire_wind_k0(
    amrex::MultiFab&       u_k0,
    amrex::MultiFab&       v_k0,
    amrex::MultiFab&       h_k0,
    const amrex::MultiFab& xvel,
    const amrex::MultiFab& yvel,
    const amrex::MultiFab& z_phys_cc,
    const amrex::MultiFab* z_phys_nd,
    const amrex::Geometry& geom);

#endif
