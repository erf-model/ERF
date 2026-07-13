/**
 * @file ERF_DustGrid.H
 * @brief Definition of the 2D dust emission grid and factory function.
 *
 * The dust grid is a horizontal slab covering the same physical domain as the
 * ERF level-0 atmospheric grid, optionally refined by grid_ratio in x and y.
 * The MPI decomposition follows the atmospheric DistributionMapping to avoid
 * inter-rank communication during coupling phases.
 *
 * References:
 *   AMReX BoxArray refinement:
 *     https://amrex-codes.github.io/amrex/docs_html/BoxArray.html
 */

#ifndef ERF_DUST_GRID_H
#define ERF_DUST_GRID_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>

/**
 * @struct DustGrid
 * @brief Container for the 2D dust emission grid infrastructure.
 *
 * Holds the BoxArray, DistributionMapping, Geometry, and refinement factor
 * of the 2D dust grid.
 */
struct DustGrid
{
    amrex::BoxArray           ba;          ///< 2D BoxArray (grid_ratio·nx × grid_ratio·ny × 1)
    amrex::DistributionMapping dm;         ///< Distribution map aligned with atmospheric layout
    amrex::Geometry           geom;        ///< 2D geometry in physical coordinates [m]
    int                       grid_ratio;  ///< Refinement factor
};

/**
 * @brief Construct a 2D dust emission grid from the atmospheric level-0 grid.
 *
 * The dust grid is a horizontal slab (k-extent = 1) covering the same
 * physical domain as the atmospheric level-0 grid, optionally refined by
 * grid_ratio in x and y. The MPI decomposition follows the atmospheric
 * DistributionMapping so that each rank owns the same horizontal tiles in
 * both grids.
 *
 * Algorithm:
 *   1. Extract k=0 slab from each box in ba_atm using setSmall/setBig on z.
 *   2. Refine the 2D BoxArray by IntVect{grid_ratio, grid_ratio, 1}.
 *   3. Reuse dm_atm directly (refinement does not change the number of boxes).
 *   4. Construct a refined 2D Geometry whose index-space domain spans
 *      [0, grid_ratio*nx - 1] x [0, grid_ratio*ny - 1] x [0, 0], with the
 *      same physical (RealBox) extent as geom_atm. The z physical extent is
 *      set to 1 m as a dummy value (the dust grid is 2D).
 *
 * References:
 *   AMReX BoxArray refinement:
 *     https://amrex-codes.github.io/amrex/docs_html/BoxArray.html
 *
 * @param[in] ba_atm       Atmospheric level-0 BoxArray.
 * @param[in] dm_atm       Atmospheric level-0 DistributionMapping.
 * @param[in] geom_atm     Atmospheric level-0 Geometry.
 * @param[in] grid_ratio   Horizontal refinement factor (>= 1).
 * @return DustGrid with refined 2D domain.
 */
DustGrid create_dust_grid(const amrex::BoxArray&          ba_atm,
                          const amrex::DistributionMapping& dm_atm,
                          const amrex::Geometry&          geom_atm,
                          int                             grid_ratio);

#endif // ERF_DUST_GRID_H
