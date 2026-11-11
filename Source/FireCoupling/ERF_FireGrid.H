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
 * factor C relative to one atmospheric AMR level (erf.fire.anchor_level, the
 * finest level by default). The fire grid covers that level's boxes, keeps their
 * x-y decomposition and their MPI owners, and its index space starts at zero at
 * the lower corner of the region they cover.
 */

/**
 * @struct FireGrid
 * @brief 2D horizontal fire computational grid
 *
 * Represents a 2D slab (kmax=1) with dimensions C*nx × C*ny × 1, where C is the
 * refinement factor and (nx, ny) is the size of the region the atmospheric
 * level's boxes cover: the whole domain on level 0, the refined region on a
 * finer level.
 */
struct FireGrid
{
    amrex::BoxArray              ba;   ///< 2D BoxArray (C*nx × C*ny × 1)
    amrex::DistributionMapping   dm;   ///< Distribution map (x-y decomposition aligned with atm)
    amrex::Geometry              geom; ///< 2D Geometry in physical coordinates [m]
    int                          C;    ///< Refinement factor
    /// Index, on the atmospheric level the grid refines, of the lower corner of
    /// the region the fire grid covers. The fire grid's own index space starts at
    /// zero at that corner, as every AMReX Geometry does, so fire cell (i, j) lies
    /// in atmospheric cell (i / C + atm_lo[0], j / C + atm_lo[1]). Zero on level 0.
    amrex::IntVect               atm_lo {0, 0, 0};
    int                          lev = 0; ///< The atmospheric level (erf.fire.anchor_level, resolved)
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
 *         the caller to set.
 *
 * Algorithm:
 * 1. Take the k = 0 slab of each box in ba_atm, keeping the box order.
 * 2. The region is their bounding rectangle, which must be covered exactly
 *    (aborts otherwise); atm_lo is its lower corner.
 * 3. Shift the slabs by -atm_lo and refine by {C, C, 1}; box n keeps the owner
 *    of atmospheric box n, so dm_atm is reused.
 * 4. The Geometry spans the region: index space 0 to C * size - 1, the
 *    physical extent of the region (the domain's own bounds where the region
 *    spans the domain), periodic only in a direction the region spans.
 */
FireGrid create_fire_grid(const amrex::BoxArray& ba_atm,
                          const amrex::DistributionMapping& dm_atm,
                          const amrex::Geometry& geom_atm,
                          int C);

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>

/**
 * Exchange the ghost cells of a fire-grid MultiFab and fill the ones outside a
 * non-periodic face of the domain by copying the nearest interior value.
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
    mf.FillBoundary(geom.periodicity());
    const auto ng = mf.nGrowVect();
    const bool px = geom.isPeriodic(0) || ng[0] == 0;
    const bool py = geom.isPeriodic(1) || ng[1] == 0;
    const bool pz = geom.isPeriodic(2) || ng[2] == 0;
    if (px && py && pz) { return; }   // nothing outside a non-periodic face to fill
    const amrex::Box& dom = geom.Domain();
    const auto dlo = dom.smallEnd();
    const auto dhi = dom.bigEnd();
    const int ncomp = mf.nComp();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& gbx = mfi.fabbox();
        if (dom.contains(gbx)) { continue; }
        auto a = mf.array(mfi);
        amrex::ParallelFor(gbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            int ic = i, jc = j, kc = k;
            if (!px) { ic = amrex::min(amrex::max(i, dlo[0]), dhi[0]); }
            if (!py) { jc = amrex::min(amrex::max(j, dlo[1]), dhi[1]); }
            if (!pz) { kc = amrex::min(amrex::max(k, dlo[2]), dhi[2]); }
            if (ic != i || jc != j || kc != k) { a(i, j, k, n) = a(ic, jc, kc, n); }
        });
    }
}

/**
 * Copy a k = 0 atmospheric field onto the fire grid: every fire cell takes the
 * value of the atmospheric column it lies in, (i / C + atm_lo[0], j / C + atm_lo[1]).
 * Box n of the fire grid refines box n of the atmospheric level and has the
 * same owner, so one MFIter serves both MultiFabs.
 */
inline void fill_fire_from_atm_k0 (amrex::MultiFab& fire_mf, const amrex::MultiFab& atm_k0,
                                   const FireGrid& fg)
{
    const int C    = fg.C;
    const int lo_x = fg.atm_lo[0];
    const int lo_y = fg.atm_lo[1];
    for (amrex::MFIter mfi(fire_mf, false); mfi.isValid(); ++mfi) {
        auto f = fire_mf.array(mfi);
        auto a = atm_k0.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) noexcept
        {
            f(iv[0], iv[1], 0) = a(iv[0] / C + lo_x, iv[1] / C + lo_y, 0);
        });
    }
}

#endif
