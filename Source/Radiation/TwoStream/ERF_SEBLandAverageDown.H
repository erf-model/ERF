#ifndef ERF_SEB_LAND_AVERAGE_DOWN_H_
#define ERF_SEB_LAND_AVERAGE_DOWN_H_

#include <AMReX_BoxArray.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_iMultiFab.H>

#include "ERF_SEBLandMask.H"

/**
 * @file ERF_SEBLandAverageDown.H
 * @brief The average-down of the two-stream force-restore balance's surface fields between
 * two levels whose land masks hold water.
 *
 * The balance runs on land columns only (TwoStreamRadiation::advance), and a water column
 * keeps the values it had at the start, at the restart, or when its level was built (the
 * snow layer's water is set to zero there). An average over every fine column would then
 * mix those held values into a coarse land column on a shore, and would overwrite a coarse
 * water column that lies under fine land. Here a coarse land column takes the mean of the fine
 * land columns under it; a coarse water column, and a coarse land column with no fine land
 * column under it, keep their values.
 *
 * With no water column on either level this is not called: ERF then uses
 * amrex::average_down, as before the mask.
 *
 * All fields go in one call: the fine land columns are counted once and one ParallelCopy
 * carries every field. A pair with a null pointer on either side is left out.
 *
 * @param fine      the 2D surface fields of the fine level (any number of components each)
 * @param crse      the same fields of the coarse level; updated under the fine level
 * @param fine_mask the fine level's land mask (0 water), on the grids of @p fine
 * @param crse_mask the coarse level's land mask (0 water), on the grids of @p crse
 * @param ratio     the refinement ratio, 1 in z
 */
inline void
seb_average_down_land (const amrex::Vector<const amrex::MultiFab*>& fine,
                       const amrex::Vector<amrex::MultiFab*>& crse,
                       const amrex::iMultiFab& fine_mask, const amrex::iMultiFab& crse_mask,
                       const amrex::IntVect& ratio)
{
    using amrex::Real;
    AMREX_ALWAYS_ASSERT(fine.size() == crse.size() && ratio[2] == 1);

    // The fields in use, and the first component of each in the packed array.
    amrex::Vector<int> used;
    amrex::Vector<int> first;
    int nc = 0;
    for (int f = 0; f < fine.size(); ++f) {
        if (fine[f] == nullptr || crse[f] == nullptr) { continue; }
        AMREX_ALWAYS_ASSERT(fine[f]->nComp() == crse[f]->nComp());
        seb_assert_land_mask_layout(&fine_mask, *fine[f], "seb_average_down_land (fine)");
        seb_assert_land_mask_layout(&crse_mask, *crse[f], "seb_average_down_land (coarse)");
        used.push_back(f);
        first.push_back(nc);
        nc += fine[f]->nComp();
    }
    if (used.empty()) { return; }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine_mask.boxArray().coarsenable(ratio),
        "seb_average_down_land: the fine grids are not aligned with the refinement ratio");

    const int rx = ratio[0];
    const int ry = ratio[1];

    // The mean over the fine land columns, on the coarsened fine grids; the last component
    // is the number of land columns it was taken over, counted once for all fields.
    amrex::BoxArray cba = fine_mask.boxArray();
    cba.coarsen(ratio);
    amrex::MultiFab cfine(cba, fine_mask.DistributionMap(), nc + 1, 0);
    for (amrex::MFIter mfi(cfine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto c = cfine.array(mfi);
        const SEBLandMaskView fm = seb_land_mask_view(&fine_mask, mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int n_land = 0;
            for (int jj = ry * j; jj < ry * (j + 1); ++jj) {
                for (int ii = rx * i; ii < rx * (i + 1); ++ii) {
                    if (!fm.water(ii, jj)) { ++n_land; }
                }
            }
            c(i, j, k, nc) = static_cast<Real>(n_land);
        });
        for (int u = 0; u < used.size(); ++u) {
            const auto f = fine[used[u]]->const_array(mfi);
            const int nf = fine[used[u]]->nComp();
            const int c0 = first[u];
            amrex::ParallelFor(bx, nf, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                // The divisor is chosen first: no division by zero is formed for a coarse
                // column with no fine land column.
                const Real n_land = c(i, j, k, nc);
                const Real divisor = (n_land > Real(0.5)) ? n_land : Real(1.0);
                Real sum = Real(0.0);
                for (int jj = ry * j; jj < ry * (j + 1); ++jj) {
                    for (int ii = rx * i; ii < rx * (i + 1); ++ii) {
                        if (!fm.water(ii, jj)) { sum += f(ii, jj, k, n); }
                    }
                }
                c(i, j, k, c0 + n) = sum / divisor;
            });
        }
    }

    // On the coarse grids: zero land columns where no fine level lies above.
    const amrex::MultiFab& crse0 = *crse[used[0]];
    amrex::MultiFab cover(crse0.boxArray(), crse0.DistributionMap(), nc + 1, 0);
    cover.setVal(Real(0.0));
    cover.ParallelCopy(cfine, 0, 0, nc + 1);
    for (amrex::MFIter mfi(cover, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto v = cover.const_array(mfi);
        const SEBLandMaskView cm = seb_land_mask_view(&crse_mask, mfi);
        for (int u = 0; u < used.size(); ++u) {
            const auto c = crse[used[u]]->array(mfi);
            const int nf = crse[used[u]]->nComp();
            const int c0 = first[u];
            amrex::ParallelFor(bx, nf, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (cm.water(i, j) || v(i, j, k, nc) < Real(0.5)) { return; }
                c(i, j, k, n) = v(i, j, k, c0 + n);
            });
        }
    }
}

/** @brief One field: the form above with a single pair. */
inline void
seb_average_down_land (const amrex::MultiFab& fine, amrex::MultiFab& crse,
                       const amrex::iMultiFab& fine_mask, const amrex::iMultiFab& crse_mask,
                       const amrex::IntVect& ratio)
{
    seb_average_down_land(amrex::Vector<const amrex::MultiFab*>{&fine},
                          amrex::Vector<amrex::MultiFab*>{&crse}, fine_mask, crse_mask, ratio);
}

#endif
