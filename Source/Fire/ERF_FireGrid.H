#ifndef ERF_FIRE_GRID_H
#define ERF_FIRE_GRID_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>

/**
 * @file ERF_FireGrid.H
 * @brief Fire grid data structure and construction
 *
 * Defines a 2D grid for fire simulation on the horizontal plane, refined by
 * factor C relative to the atmospheric level-0 grid. The fire grid maintains
 * x-y domain decomposition aligned with the atmospheric MPI layout.
 */

/**
 * @struct FireGrid
 * @brief 2D horizontal fire computational grid
 *
 * Represents a 2D slab (kmax=1) with dimensions C*nx × C*ny × 1,
 * where C is the refinement factor and (nx, ny) is the atmospheric grid size.
 */
struct FireGrid
{
    amrex::BoxArray              ba;   ///< 2D BoxArray (C*nx × C*ny × 1)
    amrex::DistributionMapping   dm;   ///< Distribution map (x-y decomposition aligned with atm)
    amrex::Geometry              geom; ///< 2D Geometry in physical coordinates [m]
    int                          C;    ///< Refinement factor
};

/**
 * @brief Construct a refined 2D fire grid from atmospheric level-0
 *
 * @param[in] ba_atm       Atmospheric level-0 BoxArray
 * @param[in] dm_atm       Atmospheric level-0 DistributionMapping
 * @param[in] geom_atm     Atmospheric level-0 Geometry
 * @param[in] C            Refinement factor (grid_ratio from FireParams)
 *
 * @return FireGrid with refined 2D domain
 *
 * Algorithm:
 * 1. Extract k=0 2D slice from each box in ba_atm using makeSlab()
 * 2. Refine the 2D BoxArray by {C, C, 1} using amrex::refine()
 * 3. Refine DistributionMapping by {C, C, 1} to align ownership
 * 4. Create 2D Geometry from geom_atm (drop z extent, keep x-y physical domain)
 */
FireGrid create_fire_grid(const amrex::BoxArray& ba_atm,
                          const amrex::DistributionMapping& dm_atm,
                          const amrex::Geometry& geom_atm,
                          int C);

#endif
