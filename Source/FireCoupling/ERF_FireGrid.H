#ifndef ERF_FIRE_GRID_H
#define ERF_FIRE_GRID_H

#include <ERF_SurfaceGrid.H>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

/**
 * @file ERF_FireGrid.H
 * @brief Fire grid data structure and construction
 *
 * Defines a 2D grid for fire simulation on the horizontal plane, refined by
 * factor C relative to one atmospheric AMR level (erf.fire.anchor_level, the
 * finest level by default). The fire grid covers that level's boxes, keeps their
 * x-y decomposition and their MPI owners, and its index space starts at zero at
 * the lower corner of the region they cover. It is the surface grid of the
 * surface-coupling layer (ERF_SurfaceGrid.H) under the fire coupling's name;
 * the functions below call that layer with the fire coupling's names.
 */

/// The fire coupling's names: erf.fire.<key>, "[FIRE]", "fire grid"
inline const SurfaceCouplingNames& fire_coupling_names ()
{
    static const SurfaceCouplingNames names{"fire", "Fire", "FIRE"};
    return names;
}

/**
 * @struct FireGrid
 * @brief 2D horizontal fire computational grid
 *
 * Represents a 2D slab (kmax=1) with dimensions C*nx × C*ny × 1, where C is the
 * refinement factor (erf.fire.grid_ratio) and (nx, ny) is the size of the
 * region the atmospheric level's boxes cover: the whole domain on level 0, the
 * refined region on a finer level. The members (ba, dm, geom, C, atm_lo, lev)
 * are those of SurfaceGrid.
 */
struct FireGrid : SurfaceGrid
{
};

/**
 * @brief Construct a refined 2D fire grid over one atmospheric level's boxes
 *
 * @param[in] ba_atm       BoxArray of the atmospheric level
 * @param[in] dm_atm       DistributionMapping of the atmospheric level
 * @param[in] geom_atm     Geometry of the atmospheric level
 * @param[in] C            Refinement factor (erf.fire.grid_ratio)
 *
 * @return FireGrid over the region the boxes cover; its lev is left at 0 for
 *         the caller to set. See create_surface_grid() for the algorithm.
 */
FireGrid create_fire_grid(const amrex::BoxArray& ba_atm,
                          const amrex::DistributionMapping& dm_atm,
                          const amrex::Geometry& geom_atm,
                          int C);

/**
 * Exchange the ghost cells of a fire-grid MultiFab and fill the ones outside a
 * non-periodic face of the domain by copying the nearest interior value
 * (surface_fill_boundary).
 *
 * FillBoundary alone leaves those cells untouched, so with an inflow-outflow
 * atmosphere (the fire grid inherits its periodicity) the gradient stencil of
 * the level set read whatever the allocator had left there. Whether that was
 * benign depended on the box layout: one box per rank happened to run, two
 * boxes on a rank blew the level set up at the inflow face within three steps.
 * Every fire-grid exchange goes through here.
 */
inline void fire_fill_boundary (amrex::MultiFab& mf, const amrex::Geometry& geom)
{
    surface_fill_boundary(mf, geom, fire_coupling_names());
}

/**
 * Copy a k = 0 atmospheric field onto the fire grid: every fire cell takes the
 * value of the atmospheric column it lies in, (i / C + atm_lo[0], j / C + atm_lo[1])
 * (surface_fill_from_atm_k0). Box n of the fire grid refines box n of the
 * atmospheric level and has the same owner, so one MFIter serves both MultiFabs.
 */
inline void fill_fire_from_atm_k0 (amrex::MultiFab& fire_mf, const amrex::MultiFab& atm_k0,
                                   const FireGrid& fg)
{
    surface_fill_from_atm_k0(fire_mf, atm_k0, fg, fire_coupling_names());
}

#endif
