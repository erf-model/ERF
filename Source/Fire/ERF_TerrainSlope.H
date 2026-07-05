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
 * @param[in]  z_phys_nd      Nodal terrain height from ERF [m]
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
    const amrex::MultiFab& z_phys_nd,
    const amrex::Geometry& geom_atm,
    const FireGrid& fg,
    const std::string& terrain_fname = "");

#endif
