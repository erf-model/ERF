#ifndef ERF_SEB_LAND_AVERAGE_DOWN_H_
#define ERF_SEB_LAND_AVERAGE_DOWN_H_

#include <AMReX_BoxArray.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

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
 * @param fine      a 2D surface field of the fine level (any number of components)
 * @param crse      the same field of the coarse level; updated under the fine level
 * @param fine_mask the fine level's land mask (0 water), on the grids of @p fine
 * @param crse_mask the coarse level's land mask (0 water), on the grids of @p crse
 * @param ratio     the refinement ratio, 1 in z
 */
inline void
seb_average_down_land (const amrex::MultiFab& fine, amrex::MultiFab& crse,
                       const amrex::iMultiFab& fine_mask, const amrex::iMultiFab& crse_mask,
                       const amrex::IntVect& ratio)
{
    using amrex::Real;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_mask.boxArray() == fine.boxArray() && fine_mask.DistributionMap() == fine.DistributionMap() &&
        crse_mask.boxArray() == crse.boxArray() && crse_mask.DistributionMap() == crse.DistributionMap(),
        "seb_average_down_land: a land mask is not on its field's grids");
    AMREX_ALWAYS_ASSERT(fine.nComp() == crse.nComp() && ratio[2] == 1);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine.boxArray().coarsenable(ratio),
        "seb_average_down_land: the fine grids are not aligned with the refinement ratio");

    const int nc = fine.nComp();
    const int rx = ratio[0];
    const int ry = ratio[1];

    // The mean over the fine land columns, on the coarsened fine grids; the last component
    // is the number of land columns it was taken over.
    amrex::BoxArray cba = fine.boxArray();
    cba.coarsen(ratio);
    amrex::MultiFab cfine(cba, fine.DistributionMap(), nc + 1, 0);
    for (amrex::MFIter mfi(cfine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto c = cfine.array(mfi);
        const auto f = fine.const_array(mfi);
        const auto fm = fine_mask.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int n_land = 0;
            for (int jj = ry * j; jj < ry * (j + 1); ++jj) {
                for (int ii = rx * i; ii < rx * (i + 1); ++ii) {
                    if (fm(ii, jj, 0) != 0) { ++n_land; }
                }
            }
            c(i, j, k, nc) = static_cast<Real>(n_land);
            // The divisor is chosen first: no division by zero is formed for a coarse
            // column with no fine land column.
            const Real divisor = static_cast<Real>(n_land > 0 ? n_land : 1);
            for (int n = 0; n < nc; ++n) {
                Real sum = Real(0.0);
                for (int jj = ry * j; jj < ry * (j + 1); ++jj) {
                    for (int ii = rx * i; ii < rx * (i + 1); ++ii) {
                        if (fm(ii, jj, 0) != 0) { sum += f(ii, jj, k, n); }
                    }
                }
                c(i, j, k, n) = sum / divisor;
            }
        });
    }

    // On the coarse grids: zero land columns where no fine level lies above.
    amrex::MultiFab cover(crse.boxArray(), crse.DistributionMap(), nc + 1, 0);
    cover.setVal(Real(0.0));
    cover.ParallelCopy(cfine, 0, 0, nc + 1);
    for (amrex::MFIter mfi(crse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto c = crse.array(mfi);
        const auto v = cover.const_array(mfi);
        const auto cm = crse_mask.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (cm(i, j, 0) == 0 || v(i, j, k, nc) < Real(0.5)) { return; }
            for (int n = 0; n < nc; ++n) { c(i, j, k, n) = v(i, j, k, n); }
        });
    }
}

#endif
