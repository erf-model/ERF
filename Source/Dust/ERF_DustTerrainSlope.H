#ifndef ERF_DUST_TERRAIN_SLOPE_H
#define ERF_DUST_TERRAIN_SLOPE_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <ERF_DustGrid.H>
#include <string>

/**
 * @file ERF_DustTerrainSlope.H
 * @brief Terrain slope computation on dust grid
 *
 * Computes dz/dx and dz/dy at dust grid cell centers via centered finite
 * differences of nodal terrain height from ERF. Result is stored and
 * accessed throughout the dust simulation.
 *
 * Reference: wildfire_levelset/src/terrain_slope.H
 */

/**
 * @brief Compute terrain slopes on dust grid
 *
 * Computes centered differences of z_phys_nd (nodal terrain height) at
 * dust grid cell centers via bilinear interpolation from coarser atmospheric nodes.
 * If a terrain file is provided, uses fine-grid interpolated terrain instead of
 * atmospheric terrain data.
 *
 * @param[out] dust_slopes    2-component MultiFab (dz/dx, dz/dy) on dust grid [m/m]
 * @param[in]  z_phys_nd      Nodal terrain height from ERF [m]
 * @param[in]  geom_atm       Atmospheric level-0 geometry
 * @param[in]  dg             Dust grid structure
 * @param[in]  terrain_fname  Path to terrain file; empty string uses z_phys_nd (default: "")
 *
 * Algorithm (with terrain file):
 * 1. Allocate z_dust_nd on dust grid with 1 ghost cell
 * 2. Fill z_dust_nd from the optional terrain file
 * 3. Fill boundary conditions with periodicity
 * 4. Compute slopes using centered bilinear averages at dust-grid resolution
 *
 * Algorithm (without terrain file, fallback):
 * 1. For each dust cell (i_d, j_d) at cell center
 * 2. Compute indices in nodal array via interpolation (i_n = i_d/C, j_n = j_d/C)
 * 3. Extract z values at four corners of the nodal element
 * 4. Compute slopes using centered bilinear average formula with atmospheric spacing
 *
 * Called once at DustLayer::initialize()
 */
void compute_dust_terrain_slopes(
    amrex::MultiFab& dust_slopes,
    const amrex::MultiFab& z_phys_nd,
    const amrex::Geometry& geom_atm,
    const DustGrid& dg,
    const std::string& terrain_fname = "");

#endif
