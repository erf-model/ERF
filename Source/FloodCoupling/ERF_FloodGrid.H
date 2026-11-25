#ifndef ERF_FLOOD_GRID_H
#define ERF_FLOOD_GRID_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>

/**
 * @file ERF_FloodGrid.H
 * @brief Flood grid data structure, construction and the maps between the
 *        flood grid and the atmospheric columns
 *
 * A 2D grid on the horizontal plane, refined by grid_ratio relative to one
 * atmospheric AMR level (erf.flood.anchor_level, the finest level by default).
 * The flood grid covers that level's boxes, keeps their x-y decomposition and
 * their MPI owners, and its index space starts at zero at the lower corner of
 * the region they cover. The construction is the one of the fire grid
 * (ERF_FireGrid.H), kept separate so that neither build needs the other.
 */

/**
 * @struct FloodGrid
 * @brief 2D horizontal flood grid
 *
 * A slab one cell thick with C*nx x C*ny cells, where C is erf.flood.grid_ratio
 * and (nx, ny) is the size of the region the atmospheric level's boxes cover:
 * the whole domain on level 0, the refined region on a finer level.
 */
struct FloodGrid
{
    amrex::BoxArray              ba;   ///< 2D BoxArray (C*nx x C*ny x 1), box n refines atmospheric box n
    amrex::DistributionMapping   dm;   ///< The atmospheric level's owners
    amrex::Geometry              geom; ///< 2D Geometry in physical coordinates [m]
    int                          C = 1; ///< erf.flood.grid_ratio
    /// Index, on the atmospheric level the grid refines, of the lower corner of
    /// the region the flood grid covers: flood cell (i, j) lies in atmospheric
    /// cell (i / C + atm_lo[0], j / C + atm_lo[1]). Zero on level 0.
    amrex::IntVect               atm_lo {0, 0, 0};
    int                          lev = 0; ///< The atmospheric level (erf.flood.anchor_level, resolved)
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
 *         the caller to set.
 *
 * 1. Take the k = 0 slab of each box in ba_atm, keeping the box order.
 * 2. The region is their bounding rectangle, which must be covered exactly
 *    (aborts otherwise); atm_lo is its lower corner.
 * 3. Shift the slabs by -atm_lo and refine by {C, C, 1}; box n keeps the owner
 *    of atmospheric box n, so dm_atm is reused.
 * 4. The Geometry spans the region: index space 0 to C * size - 1, the
 *    physical extent of the region (the domain's own bounds where the region
 *    spans the domain), periodic only in a direction the region spans.
 */
FloodGrid create_flood_grid (const amrex::BoxArray& ba_atm,
                             const amrex::DistributionMapping& dm_atm,
                             const amrex::Geometry& geom_atm,
                             int C);

/**
 * True when box n of a and box n of b cover the same columns (the same x-y
 * extent) for every n; their z extents may differ.
 */
inline bool flood_same_columns (const amrex::BoxArray& a, const amrex::BoxArray& b)
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
 * True when box n of the flood grid refines the columns of box n of the
 * atmospheric level (shifted by -atm_lo, refined by {C, C, 1}) for every n,
 * and the two have the same owners.
 */
inline bool flood_pairs_with_atm (const amrex::MultiFab& flood_mf, const amrex::MultiFab& atm_mf,
                                  const FloodGrid& fg)
{
    const amrex::BoxArray& fba = flood_mf.boxArray();
    const amrex::BoxArray& aba = atm_mf.boxArray();
    if (fba.size() != aba.size() || flood_mf.DistributionMap() != atm_mf.DistributionMap()) {
        return false;
    }
    for (int n = 0; n < aba.size(); ++n) {
        amrex::Box ab = aba[n];
        amrex::Box fb = fba[n];
        ab.enclosedCells();
        fb.enclosedCells();
        ab.setSmall(2, 0); ab.setBig(2, 0);
        fb.setSmall(2, 0); fb.setBig(2, 0);
        ab.shift(-fg.atm_lo);
        ab.refine(amrex::IntVect(fg.C, fg.C, 1));
        if (ab != fb) { return false; }
    }
    return true;
}

/**
 * Exchange the ghost cells of a flood-grid MultiFab and fill the ones outside
 * a non-periodic face of the flood region by copying the nearest interior value.
 * For cell-centred MultiFabs only.
 */
inline void flood_fill_boundary (amrex::MultiFab& mf, const amrex::Geometry& geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf.ixType().cellCentered(),
        "[FLOOD] flood_fill_boundary: the flood-grid MultiFab must be cell-centred");
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
 * Copy component acomp of a k = 0 atmospheric field into component fcomp of a
 * flood-grid field: every flood cell takes the value of the atmospheric column
 * it lies in, (i / C + atm_lo[0], j / C + atm_lo[1]). Box n of the flood grid
 * refines box n of the atmospheric level and has the same owner, so one MFIter
 * serves both MultiFabs. A per-area quantity (precipitation in kg/m2) keeps its
 * column total: each of the C^2 cells under a column receives the column's
 * value per unit area.
 */
inline void fill_flood_from_atm_k0 (amrex::MultiFab& flood_mf, const amrex::MultiFab& atm_k0,
                                    const FloodGrid& fg, int fcomp = 0, int acomp = 0)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        flood_pairs_with_atm(flood_mf, atm_k0, fg),
        "[FLOOD] fill_flood_from_atm_k0: the flood grid and the atmospheric slabs must pair box for box");
    const int C    = fg.C;
    const int lo_x = fg.atm_lo[0];
    const int lo_y = fg.atm_lo[1];
    for (amrex::MFIter mfi(flood_mf, false); mfi.isValid(); ++mfi) {
        auto f = flood_mf.array(mfi);
        auto a = atm_k0.const_array(mfi);
        const int ka = atm_k0.boxArray()[mfi.index()].smallEnd(2);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            f(i, j, k, fcomp) = a(i / C + lo_x, j / C + lo_y, ka, acomp);
        });
    }
}

/**
 * Average a flood-grid field over each atmospheric column (the C x C flood
 * cells under it) into the k = 0 slabs of the flood grid's level. Box n of the
 * output is the k = 0 slab of atmospheric box n.
 */
inline void coarsen_flood_to_atm (amrex::MultiFab& atm_out, const amrex::MultiFab& flood,
                                  const FloodGrid& fg)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        flood_pairs_with_atm(flood, atm_out, fg),
        "[FLOOD] coarsen_flood_to_atm: the flood grid and the atmospheric slabs must pair box for box");
    AMREX_ALWAYS_ASSERT(atm_out.nComp() >= flood.nComp());
    const int C     = fg.C;
    const int lo_x  = fg.atm_lo[0];
    const int lo_y  = fg.atm_lo[1];
    const int ncomp = flood.nComp();
    const amrex::Real inv_cells = amrex::Real(1.0) / amrex::Real(C * C);
    for (amrex::MFIter mfi(atm_out); mfi.isValid(); ++mfi) {
        auto crse = atm_out.array(mfi);
        auto fine = flood.const_array(mfi);
        const int kf = flood.boxArray()[mfi.index()].smallEnd(2);
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
