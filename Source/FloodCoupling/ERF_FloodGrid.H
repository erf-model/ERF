#ifndef ERF_FLOOD_GRID_H
#define ERF_FLOOD_GRID_H

#include <ERF_SurfaceGrid.H>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

/**
 * @file ERF_FloodGrid.H
 * @brief Flood grid data structure, construction and the maps between the
 *        flood grid and the atmospheric columns
 *
 * A 2D grid on the horizontal plane, refined by grid_ratio relative to one
 * atmospheric AMR level (erf.flood.anchor_level, the finest level by default).
 * The flood grid covers that level's boxes, keeps their x-y decomposition and
 * their MPI owners, and its index space starts at zero at the lower corner of
 * the region they cover. It is the surface grid of the surface-coupling layer
 * (ERF_SurfaceGrid.H) under the flood coupling's name, as the fire grid is;
 * the functions below call that layer with the flood coupling's names.
 */

/// The flood coupling's names: erf.flood.<key>, "[FLOOD]", "flood grid"
inline const SurfaceCouplingNames& flood_coupling_names ()
{
    static const SurfaceCouplingNames names{"flood", "Flood", "FLOOD"};
    return names;
}

/**
 * @struct FloodGrid
 * @brief 2D horizontal flood grid
 *
 * A slab one cell thick with C*nx x C*ny cells, where C is erf.flood.grid_ratio
 * and (nx, ny) is the size of the region the atmospheric level's boxes cover:
 * the whole domain on level 0, the refined region on a finer level. The
 * members (ba, dm, geom, C, atm_lo, lev) are those of SurfaceGrid.
 */
struct FloodGrid : SurfaceGrid
{
};

/**
 * @brief Construct a refined 2D flood grid over one atmospheric level's boxes
 *
 * @param[in] ba_atm    BoxArray of the atmospheric level
 * @param[in] dm_atm    DistributionMapping of the atmospheric level
 * @param[in] geom_atm  Geometry of the atmospheric level
 * @param[in] C         Refinement factor (erf.flood.grid_ratio, >= 1)
 *
 * @return FloodGrid over the region the boxes cover; its lev is left at 0 for
 *         the caller to set. See create_surface_grid() for the algorithm.
 */
FloodGrid create_flood_grid (const amrex::BoxArray& ba_atm,
                             const amrex::DistributionMapping& dm_atm,
                             const amrex::Geometry& geom_atm,
                             int C);

/// surface_same_columns: box n of a and box n of b cover the same columns.
inline bool flood_same_columns (const amrex::BoxArray& a, const amrex::BoxArray& b)
{
    return surface_same_columns(a, b);
}

/// surface_pairs_with_atm: box n of the flood grid refines the columns of
/// atmospheric box n for every n, with the same owners.
inline bool flood_pairs_with_atm (const amrex::MultiFab& flood_mf, const amrex::MultiFab& atm_mf,
                                  const FloodGrid& fg)
{
    return surface_pairs_with_atm(flood_mf, atm_mf, fg);
}

/// surface_fill_boundary: exchange the ghost cells of a cell-centred
/// flood-grid MultiFab and fill the ones outside a non-periodic face.
inline void flood_fill_boundary (amrex::MultiFab& mf, const amrex::Geometry& geom)
{
    surface_fill_boundary(mf, geom, flood_coupling_names());
}

/**
 * Copy component acomp of a k = 0 atmospheric field into component fcomp of a
 * flood-grid field (surface_fill_from_atm_k0): every flood cell takes the
 * value of the atmospheric column it lies in. A per-area quantity
 * (precipitation in kg/m2) keeps its column total.
 */
inline void fill_flood_from_atm_k0 (amrex::MultiFab& flood_mf, const amrex::MultiFab& atm_k0,
                                    const FloodGrid& fg, int fcomp = 0, int acomp = 0)
{
    surface_fill_from_atm_k0(flood_mf, atm_k0, fg, flood_coupling_names(), fcomp, acomp);
}

/**
 * Average a flood-grid field over each atmospheric column (the C x C flood
 * cells under it) into the k = 0 slabs of the flood grid's level
 * (surface_coarsen_to_atm). Box n of the output is the k = 0 slab of
 * atmospheric box n.
 */
inline void coarsen_flood_to_atm (amrex::MultiFab& atm_out, const amrex::MultiFab& flood,
                                  const FloodGrid& fg)
{
    surface_coarsen_to_atm(atm_out, flood, fg, flood_coupling_names());
}

#endif
