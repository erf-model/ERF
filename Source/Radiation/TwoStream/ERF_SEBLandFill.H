#ifndef ERF_SEB_LAND_FILL_H_
#define ERF_SEB_LAND_FILL_H_

#include <AMReX_BCRec.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>
#include <AMReX_Vector.H>
#include <AMReX_iMultiFab.H>

#include "ERF_SEBLandMask.H"

/**
 * @file ERF_SEBLandFill.H
 * @brief The fill of a new fine level's force-restore surface fields from its parent
 * (ERF::fill_seb_from_coarse).
 */

/**
 * @brief Fill the halo of a coarse 2D surface field before it is interpolated.
 *
 * FillBoundary sets the interior and periodic ghost columns. Across a non-periodic face a
 * ghost column takes the nearest column inside the domain: the slope stencil of the
 * interpolater reads it for every boundary type. A periodic ghost column is left as
 * FillBoundary set it. InterpFromCoarseLevel reads it twice, from this halo and from the
 * periodic image of the valid column, and the two must agree.
 */
inline void
seb_extend_halo (amrex::MultiFab& mf, const amrex::Geometry& geom)
{
    using amrex::Real;
    mf.FillBoundary(geom.periodicity());
    const bool clamp_x = !geom.isPeriodic(0);
    const bool clamp_y = !geom.isPeriodic(1);
    const auto dlo = amrex::lbound(geom.Domain());
    const auto dhi = amrex::ubound(geom.Domain());
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& gbx = mfi.growntilebox();
        const amrex::Array4<Real>& a = mf.array(mfi);
        amrex::ParallelFor(gbx, mf.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ic = clamp_x ? amrex::max(dlo.x, amrex::min(dhi.x, i)) : i;
            const int jc = clamp_y ? amrex::max(dlo.y, amrex::min(dhi.y, j)) : j;
            if (ic != i || jc != j) { a(i,j,k,n) = a(ic,jc,k,n); }
        });
    }
}

/**
 * @brief Fill @p fine from @p crse when the coarse land mask holds water.
 *
 * A water column holds a value the balance never updated. The fill keeps it out of the
 * land columns:
 *
 * - A fine column whose parent and the parent's eight neighbors are all land takes the
 *   conservative linear interpolation, as on a level with no water.
 * - Elsewhere no slope is used: the slope of a coarse column next to water is formed with
 *   the held value. A fine land column over a coarse land column takes the parent's value.
 *   A fine land column over a coarse water column takes the mean of the land columns among
 *   the parent's eight neighbors.
 * - A fine land column whose parent and all eight neighbors are water has no land source.
 *   It takes the parent's held value and is counted.
 * - A fine water column in that region takes its parent's value.
 *
 * Neighbors beyond a non-periodic domain face are left out. RemakeLevel copies the old
 * fine values on top afterwards, so only columns new to the level keep this fill.
 *
 * @param fine      the fine field; valid and ghost columns inside the domain are set
 * @param crse      the coarse field, after seb_extend_halo()
 * @param crse_mask the coarse land mask, ghost columns filled, at least one ghost column
 * @param fine_mask the fine land mask with the ghost columns of @p fine, or nullptr (all land)
 * @param bcs       one boundary record per component of @p fine
 * @return the number of fine land columns with no land source on this rank (not reduced)
 */
inline amrex::Long
seb_fill_from_coarse_over_land (amrex::MultiFab& fine, const amrex::MultiFab& crse,
                                const amrex::iMultiFab& crse_mask, const amrex::iMultiFab* fine_mask,
                                const amrex::Geometry& cgeom, const amrex::Geometry& fgeom,
                                const amrex::IntVect& ratio, const amrex::Vector<amrex::BCRec>& bcs)
{
    using amrex::Real;
    const int nc = fine.nComp();
    AMREX_ALWAYS_ASSERT(crse.nComp() == nc && bcs.size() == nc);
    seb_assert_land_mask_layout(&crse_mask, crse, "seb_fill_from_coarse_over_land (coarse)");
    seb_assert_land_mask_layout(fine_mask, fine, "seb_fill_from_coarse_over_land (fine)");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(crse_mask.nGrow(0) >= 1 && crse_mask.nGrow(1) >= 1 &&
                                     crse.nGrow(0) >= 1 && crse.nGrow(1) >= 1,
        "seb_fill_from_coarse_over_land: the coarse field and mask need one ghost column");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine_mask == nullptr || fine_mask->nGrowVect().allGE(fine.nGrowVect()),
        "seb_fill_from_coarse_over_land: the fine mask has fewer ghost columns than the field");

    // The interpolation with slopes, as on a level with no water.
    amrex::InterpFromCoarseLevel(fine, fine.nGrowVect(), amrex::IntVect(0,0,0),
                                 crse, 0, 0, nc, cgeom, fgeom, ratio, &amrex::cell_cons_interp, bcs, 0);

    // Per coarse column: [0, nc) the value a fine land column takes, [nc, 2 nc) the column's
    // own value, 2 nc whether the column or one of its eight neighbors is water, 2 nc + 1
    // whether it is water with no land neighbor.
    const int i_near = 2 * nc;
    const int i_none = 2 * nc + 1;
    amrex::MultiFab src(crse.boxArray(), crse.DistributionMap(), 2 * nc + 2, crse.nGrowVect());
    // The ghost columns at a coarse-fine edge of the coarse level are not set below, and
    // those of crse and crse_mask there are not refreshed from the level under it. With
    // the proper nesting ERF runs with (amr.n_proper = 2) no fine column has its parent
    // there. The zero flag is a stopgap for looser nesting: it leaves such a column with
    // the interpolation above, which still reads the coarse ghost columns as they are.
    // Filling them from the level under the coarse one first would be the full fix.
    src.setBndry(Real(0.0));
    const bool per_x = cgeom.isPeriodic(0);
    const bool per_y = cgeom.isPeriodic(1);
    const auto dlo = amrex::lbound(cgeom.Domain());
    const auto dhi = amrex::ubound(cgeom.Domain());
    for (amrex::MFIter mfi(src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto s = src.array(mfi);
        const auto c = crse.const_array(mfi);
        const SEBLandMaskView cm = seb_land_mask_view(&crse_mask, mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int n_water = 0;
            int n_land_nb = 0;
            for (int jj = j - 1; jj <= j + 1; ++jj) {
                for (int ii = i - 1; ii <= i + 1; ++ii) {
                    const bool inside = (per_x || (ii >= dlo.x && ii <= dhi.x)) &&
                                        (per_y || (jj >= dlo.y && jj <= dhi.y));
                    if (!inside) { continue; }
                    if (cm.water(ii, jj)) { ++n_water; }
                    else if (ii != i || jj != j) { ++n_land_nb; }
                }
            }
            s(i, j, k, i_near) = (n_water > 0) ? Real(1.0) : Real(0.0);
            s(i, j, k, i_none) = (cm.water(i, j) && n_land_nb == 0) ? Real(1.0) : Real(0.0);
        });
        amrex::ParallelFor(bx, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const Real own = c(i, j, k, n);
            s(i, j, k, nc + n) = own;
            if (!cm.water(i, j)) { s(i, j, k, n) = own; return; }
            Real sum = Real(0.0);
            int n_land_nb = 0;
            for (int jj = j - 1; jj <= j + 1; ++jj) {
                for (int ii = i - 1; ii <= i + 1; ++ii) {
                    const bool inside = (per_x || (ii >= dlo.x && ii <= dhi.x)) &&
                                        (per_y || (jj >= dlo.y && jj <= dhi.y));
                    if (!inside || cm.water(ii, jj)) { continue; }
                    sum += c(ii, jj, k, n);
                    ++n_land_nb;
                }
            }
            // The divisor is chosen first: no division by zero is formed.
            const Real divisor = (n_land_nb > 0) ? static_cast<Real>(n_land_nb) : Real(1.0);
            s(i, j, k, n) = (n_land_nb > 0) ? sum / divisor : own;
        });
    }
    seb_extend_halo(src, cgeom);

    // Each fine column takes its parent's entries: no arithmetic, no slope.
    amrex::MultiFab tmp(fine.boxArray(), fine.DistributionMap(), 2 * nc + 2, fine.nGrowVect());
    tmp.setVal(Real(0.0));
    const amrex::Vector<amrex::BCRec> bcs_src(2 * nc + 2, bcs[0]);
    amrex::InterpFromCoarseLevel(tmp, tmp.nGrowVect(), amrex::IntVect(0,0,0),
                                 src, 0, 0, 2 * nc + 2, cgeom, fgeom, ratio, &amrex::pc_interp, bcs_src, 0);

    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (amrex::MFIter mfi(fine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto f = fine.array(mfi);
        const auto t = tmp.const_array(mfi);
        const SEBLandMaskView fm = seb_land_mask_view(fine_mask, mfi);
        // Ghost columns outside the domain were not filled: their flag is zero.
        amrex::ParallelFor(mfi.growntilebox(), nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (t(i, j, k, i_near) > Real(0.5)) {
                f(i, j, k, n) = fm.water(i, j) ? t(i, j, k, nc + n) : t(i, j, k, n);
            }
        });
        reduce_op.eval(mfi.tilebox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            return { (!fm.water(i, j) && t(i, j, k, i_none) > Real(0.5)) ? amrex::Long(1) : amrex::Long(0) };
        });
    }
    return amrex::get<0>(reduce_data.value(reduce_op));
}

#endif
