#ifndef ERF_TERRAIN_SLOPE_H
#define ERF_TERRAIN_SLOPE_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <ERF_FireGrid.H>

/**
 * @file ERF_TerrainSlope.H
 * @brief Terrain slope computation on fire grid
 *
 * Computes dz/dx and dz/dy at fire grid cell centers via centered finite
 * differences of nodal terrain height from ERF. Result is stored and
 * accessed throughout the fire simulation.
 *
 * Reference: wildfire_levelset/src/terrain_slope.H
 */

/**
 * @brief Compute terrain slopes on fire grid
 *
 * Computes centered differences of z_phys_nd (nodal terrain height) at
 * fire grid cell centers via bilinear interpolation from coarser atmospheric nodes.
 * If a terrain file is provided, uses fine-grid interpolated terrain instead of
 * atmospheric terrain data.
 *
 * @param[out] fire_slopes    2-component MultiFab (dz/dx, dz/dy) on fire grid [m/m]
 * @param[in]  z_phys_nd      Nodal terrain height from ERF [m], or nullptr for
 *                            flat terrain (TerrainType::None), in which case the
 *                            slopes are set to zero
 * @param[in]  geom_atm       Atmospheric level-0 geometry
 * @param[in]  fg             Fire grid structure
 * @param[in]  terrain_fname  Path to terrain file; empty string uses z_phys_nd (default: "")
 *
 * Algorithm (with terrain file):
 * 1. Allocate z_fire_nd on fire grid with 1 ghost cell
 * 2. Call read_terrain_onto_fire_grid to fill z_fire_nd
 * 3. Fill boundary conditions with periodicity
 * 4. Compute slopes using centered bilinear averages at fire-grid resolution
 *
 * Algorithm (without terrain file, fallback):
 * 1. For each fire cell (i_f, j_f) at cell center
 * 2. Compute indices in nodal array via interpolation (i_n = i_f/C, j_n = j_f/C)
 * 3. Extract z values at four corners of the nodal element
 * 4. Compute slopes using centered bilinear average formula with atmospheric spacing
 *
 * Called once at FireLayer::initialize()
 */
void compute_terrain_slopes(
    amrex::MultiFab& fire_slopes,
    const amrex::MultiFab* z_phys_nd,
    const amrex::Geometry& geom_atm,
    const FireGrid& fg,
    const std::string& terrain_fname = "");

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
 * @param[in]  geom_atm        Atmospheric level-0 geometry
 * @param[in]  fg              Fire grid structure
 *
 * Called once at FireLayer::initialize(); terrain is static.
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
 * @param[in]  geom_atm         Atmospheric level-0 geometry
 * @param[in]  fg               Fire grid structure
 *
 * Called once at FireLayer::initialize(); terrain is static.
 */
void compute_fire_column_grounds(
    amrex::MultiFab& fire_col_ground,
    const amrex::MultiFab* z_phys_nd,
    const amrex::Geometry& geom_atm,
    const FireGrid& fg);

#endif
