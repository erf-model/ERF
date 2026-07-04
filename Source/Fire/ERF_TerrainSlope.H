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
 *
 * @param[out] fire_slopes    2-component MultiFab (dz/dx, dz/dy) on fire grid [m/m]
 * @param[in]  z_phys_nd      Nodal terrain height from ERF [m]
 * @param[in]  geom_atm       Atmospheric level-0 geometry
 * @param[in]  fg             Fire grid structure
 *
 * Algorithm:
 * 1. For each fire cell (i_f, j_f) at cell center
 * 2. Compute indices in nodal array via interpolation (i_n = i_f/C, j_n = j_f/C)
 * 3. Extract z values at four corners of the nodal element
 * 4. Compute dz/dx = (z_right - z_left) / (2 * dx_atm)
 * 5. Compute dz/dy = (z_top - z_bottom) / (2 * dy_atm)
 *
 * Called once at FireLayer::initialize()
 */
void compute_terrain_slopes(
    amrex::MultiFab& fire_slopes,
    const amrex::MultiFab& z_phys_nd,
    const amrex::Geometry& geom_atm,
    const FireGrid& fg);

#endif
