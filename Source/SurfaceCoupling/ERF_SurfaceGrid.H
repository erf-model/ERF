#ifndef ERF_SURFACE_GRID_H
#define ERF_SURFACE_GRID_H

#include <ERF_SurfaceCouplingNames.H>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>

/**
 * @file ERF_SurfaceGrid.H
 * @brief The 2D surface grid of the fire and flood couplings, and the maps
 *        between it and the atmospheric columns
 *
 * A 2D grid on the horizontal plane, refined by C relative to one atmospheric
 * AMR level (erf.fire.anchor_level or erf.flood.anchor_level, the finest level
 * by default). The grid covers that level's boxes, keeps their x-y
 * decomposition and their MPI owners, and its index space starts at zero at
 * the lower corner of the region they cover. FireGrid (ERF_FireGrid.H) and
 * FloodGrid (ERF_FloodGrid.H) are this grid under each coupling's name.
 */

/**
 * @struct SurfaceGrid
 * @brief 2D horizontal surface grid
 *
 * A slab one cell thick with C*nx x C*ny cells, where (nx, ny) is the size of
 * the region the atmospheric level's boxes cover: the whole domain on level 0,
 * the refined region on a finer level.
 */
struct SurfaceGrid
{
    amrex::BoxArray              ba;   ///< 2D BoxArray (C*nx x C*ny x 1), box n refines atmospheric box n
    amrex::DistributionMapping   dm;   ///< The atmospheric level's owners
    amrex::Geometry              geom; ///< 2D Geometry in physical coordinates [m]
    int                          C = 1; ///< Refinement factor (erf.fire.grid_ratio, erf.flood.grid_ratio)
    /// Index, on the atmospheric level the grid refines, of the lower corner of
    /// the region the grid covers. The grid's own index space starts at zero at
    /// that corner, as every AMReX Geometry does, so surface cell (i, j) lies in
    /// atmospheric cell (i / C + atm_lo[0], j / C + atm_lo[1]). Zero on level 0.
    amrex::IntVect               atm_lo {0, 0, 0};
    int                          lev = 0; ///< The atmospheric level (anchor_level, resolved)
};

/**
 * @brief Construct a refined 2D surface grid over one atmospheric level's boxes
 *
 * @param[in] ba_atm    BoxArray of the atmospheric level
 * @param[in] dm_atm    DistributionMapping of the atmospheric level
 * @param[in] geom_atm  Geometry of the atmospheric level
 * @param[in] C         Refinement factor (>= 1, aborts otherwise)
 * @param[in] names     The coupling's names, for the messages
 *
 * @return SurfaceGrid over the region the boxes cover; its lev is left at 0
 *         for the caller to set.
 *
 * 1. Take the k = 0 slab of each box in ba_atm, keeping the box order.
 * 2. The region is their bounding rectangle, which must be covered exactly
 *    (aborts otherwise); atm_lo is its lower corner.
 * 3. Shift the slabs by -atm_lo and refine by {C, C, 1}; box n keeps the owner
 *    of atmospheric box n, so dm_atm is reused.
 * 4. The Geometry spans the region: index space 0 to C * size - 1, the
 *    physical extent of the region (the domain's own bounds where the region
 *    spans the domain), periodic only in a direction the region spans and the
 *    atmosphere is periodic in.
 */
SurfaceGrid create_surface_grid (const amrex::BoxArray& ba_atm,
                                 const amrex::DistributionMapping& dm_atm,
                                 const amrex::Geometry& geom_atm,
                                 int C,
                                 const SurfaceCouplingNames& names);

/**
 * True when box n of a and box n of b cover the same columns (the same x-y
 * extent) for every n; their z extents may differ.
 */
inline bool surface_same_columns (const amrex::BoxArray& a, const amrex::BoxArray& b)
{
    if (a.size() != b.size()) { return false; }
    for (int n = 0; n < a.size(); ++n) {
        amrex::Box ba = a[n];
        amrex::Box bb = b[n];
        ba.setSmall(2, 0); ba.setBig(2, 0);
        bb.setSmall(2, 0); bb.setBig(2, 0);
        if (ba != bb) { return false; }
    }
    return true;
}

/**
 * True when box n of the surface grid refines the columns of box n of the
 * atmospheric level (shifted by -atm_lo, refined by {C, C, 1}) for every n,
 * the two have the same owners, and both are cell-centred (the maps loop over
 * the valid cells, so a nodal box would reach one column past its pair).
 */
inline bool surface_pairs_with_atm (const amrex::MultiFab& surf_mf, const amrex::MultiFab& atm_mf,
                                    const SurfaceGrid& sg)
{
    const amrex::BoxArray& sba = surf_mf.boxArray();
    const amrex::BoxArray& aba = atm_mf.boxArray();
    if (!sba.ixType().cellCentered() || !aba.ixType().cellCentered()) {
        return false;
    }
    if (sba.size() != aba.size() || surf_mf.DistributionMap() != atm_mf.DistributionMap()) {
        return false;
    }
    for (int n = 0; n < aba.size(); ++n) {
        amrex::Box ab = aba[n];
        amrex::Box sb = sba[n];
        ab.setSmall(2, 0); ab.setBig(2, 0);
        sb.setSmall(2, 0); sb.setBig(2, 0);
        ab.shift(-sg.atm_lo);
        ab.refine(amrex::IntVect(sg.C, sg.C, 1));
        if (ab != sb) { return false; }
    }
    return true;
}

/**
 * Exchange the ghost cells of a cell-centred MultiFab on a surface grid (or
 * on the atmospheric columns it covers) and fill the ones outside a
 * non-periodic face of the domain by copying the nearest interior value.
 *
 * FillBoundary alone leaves those cells untouched, so with an inflow-outflow
 * atmosphere (the surface grid inherits its periodicity) a stencil that reads
 * past the region's edge would read whatever the allocator had left there.
 */
inline void surface_fill_boundary (amrex::MultiFab& mf, const amrex::Geometry& geom,
                                   const SurfaceCouplingNames& names)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf.ixType().cellCentered(),
        (names.tag() + " Filling the ghost cells of the " + names.grid()
         + ": the MultiFab must be cell-centred").c_str());
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
 * Copy component acomp of a k = 0 atmospheric field into component scomp of a
 * surface-grid field: every surface cell takes the value of the atmospheric
 * column it lies in, (i / C + atm_lo[0], j / C + atm_lo[1]). Box n of the
 * surface grid refines box n of the atmospheric level and has the same owner,
 * so one MFIter serves both MultiFabs. A per-area quantity (precipitation in
 * kg/m2) keeps its column total: each of the C^2 cells under a column receives
 * the column's value per unit area.
 */
inline void surface_fill_from_atm_k0 (amrex::MultiFab& surf_mf, const amrex::MultiFab& atm_k0,
                                      const SurfaceGrid& sg, const SurfaceCouplingNames& names,
                                      int scomp = 0, int acomp = 0)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_pairs_with_atm(surf_mf, atm_k0, sg),
        (names.tag() + " Copying the atmospheric columns to the " + names.grid()
         + ": the two must pair box for box").c_str());
    const int C    = sg.C;
    const int lo_x = sg.atm_lo[0];
    const int lo_y = sg.atm_lo[1];
    for (amrex::MFIter mfi(surf_mf, false); mfi.isValid(); ++mfi) {
        auto f = surf_mf.array(mfi);
        auto a = atm_k0.const_array(mfi);
        const int ka = atm_k0.boxArray()[mfi.index()].smallEnd(2);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            f(i, j, k, scomp) = a(i / C + lo_x, j / C + lo_y, ka, acomp);
        });
    }
}

/**
 * Average a surface-grid field over each atmospheric column (the C x C surface
 * cells under it) into the k = 0 slabs of the grid's level: surface cells
 * (i - lo_x) C + r, (j - lo_y) C + q for r, q in [0, C) average into
 * atmospheric cell (i, j). Box n of the output is the k = 0 slab of
 * atmospheric box n; an output box taller than one cell gets the average at
 * every k of its valid region. On level 0, with k = 0 slabs, this is
 * amrex::average_down() with ratio (C, C, 1), summed in the same order.
 */
inline void surface_coarsen_to_atm (amrex::MultiFab& atm_out, const amrex::MultiFab& surf,
                                    const SurfaceGrid& sg, const SurfaceCouplingNames& names)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_pairs_with_atm(surf, atm_out, sg),
        (names.tag() + " Averaging the " + names.grid()
         + " onto the atmospheric columns: the two must pair box for box").c_str());
    AMREX_ALWAYS_ASSERT(atm_out.nComp() >= surf.nComp());
    const int C     = sg.C;
    const int lo_x  = sg.atm_lo[0];
    const int lo_y  = sg.atm_lo[1];
    const int ncomp = surf.nComp();
    const amrex::Real inv_cells = amrex::Real(1.0) / amrex::Real(C * C);
    for (amrex::MFIter mfi(atm_out); mfi.isValid(); ++mfi) {
        auto crse = atm_out.array(mfi);
        auto fine = surf.const_array(mfi);
        const int kf = surf.boxArray()[mfi.index()].smallEnd(2);
        amrex::ParallelFor(mfi.validbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ii = (i - lo_x) * C;
            const int jj = (j - lo_y) * C;
            amrex::Real c = amrex::Real(0.0);
            for (int q = 0; q < C; ++q) {
                for (int r = 0; r < C; ++r) {
                    c += fine(ii + r, jj + q, kf, n);
                }
            }
            crse(i, j, k, n) = inv_cells * c;
        });
    }
}

#endif
