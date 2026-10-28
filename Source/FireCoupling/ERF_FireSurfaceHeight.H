#ifndef ERF_FIRE_SURFACE_HEIGHT_H
#define ERF_FIRE_SURFACE_HEIGHT_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <ERF_FireGrid.H>

/**
 * @file ERF_FireSurfaceHeight.H
 * @brief Ground of the atmospheric columns, seen from the fire grid
 *
 * The atmosphere side of the fire coupling: the ground elevation the shared
 * wind sampler (ERF_FireWindExtract.H) measures its sampling height from, for
 * the column holding each fire cell and for the four columns of its bilinear
 * stencil. Terrain slopes, which a fire engine may take from its own
 * fire-resolved terrain, are the engine's.
 */

/**
 * @brief Fill the ground elevation of each fire cell's atmospheric column
 *
 * The fire module extracts wind at a height above ground by walking the
 * atmospheric column's z_phys_cc profile, so the ground it measures from has to
 * be the ground of that same column: the average of the four k = 0 nodes of
 * atmospheric cell (i_f/C, j_f/C). Taking the first cell centre instead puts the
 * extraction height half an atmospheric cell too high, and by an amount that
 * changes with vertical resolution.
 *
 * This is deliberately the atmospheric column's ground, not the fire-resolved
 * terrain, even when a finer terrain file is supplied for the slopes: the wind
 * profile being interpolated belongs to the atmospheric column.
 *
 * @param[out] fire_surface_z  1-component MultiFab on the fire grid [m]
 * @param[in]  z_phys_nd       Nodal terrain height from ERF [m], or nullptr for
 *                             flat terrain, in which case the domain floor is used
 * @param[in]  geom_atm        Geometry of the level the fire grid refines
 * @param[in]  fg              Fire grid structure
 *
 * Called once by the fire engine at initialisation; terrain is static.
 */
void compute_fire_surface_height(
    amrex::MultiFab& fire_surface_z,
    const amrex::MultiFab* z_phys_nd,
    const amrex::Geometry& geom_atm,
    const FireGrid& fg);

/**
 * @brief Fill the ground elevations of a fire cell's bilinear wind stencil
 *
 * Bilinear wind interpolation blends the four atmospheric columns surrounding a
 * fire cell. Each of those columns has its own ground, and near the surface the
 * wind profile is anchored to it, so every column has to be sampled at the same
 * height *above its own ground* before the horizontal blend. Over the ridge of
 * the ROS_Slope_Effects case, neighbouring column grounds differ by up to 28 m,
 * so blending at a single absolute height would sample a neighbour at 35 m above
 * ground while asking for 6 m.
 *
 * Component c holds the ground of stencil column c, ordered
 * (i0,j0), (i0+1,j0), (i0,j0+1), (i0+1,j0+1), with i0 = floor((i_f + 0.5)/C - 0.5).
 *
 * @param[out] fire_col_ground  4-component MultiFab on the fire grid [m]
 * @param[in]  z_phys_nd        Nodal terrain height from ERF [m], or nullptr for
 *                              flat terrain, in which case the domain floor is used
 * @param[in]  geom_atm         Geometry of the level the fire grid refines
 * @param[in]  fg               Fire grid structure
 *
 * Called once by the fire engine at initialisation; terrain is static.
 */
void compute_fire_column_grounds(
    amrex::MultiFab& fire_col_ground,
    const amrex::MultiFab* z_phys_nd,
    const amrex::Geometry& geom_atm,
    const FireGrid& fg);

#endif
