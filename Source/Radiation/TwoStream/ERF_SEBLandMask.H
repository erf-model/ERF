#ifndef ERF_SEB_LAND_MASK_H_
#define ERF_SEB_LAND_MASK_H_

#include <string>

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_iMultiFab.H>

/**
 * @file ERF_SEBLandMask.H
 * @brief The one place that says how the two-stream force-restore balance reads the land
 * mask (ERF::lmask_lev: 0 water, 1 land, 2 a building's column).
 *
 * The balance is a model of soil. A column is left out of it when the mask is 0 there.
 * The balance, its snow and hydrology steps, the albedo with snow, the average-down and
 * the fill of a new level all ask through this file, so a change of the rule is one edit.
 */

/**
 * @brief The land mask of one box, for a device kernel. Captured by value.
 *
 * With no mask (on == false) no column is water.
 */
struct SEBLandMaskView
{
    bool on = false;
    amrex::Array4<const int> arr;

    /** @brief True when the balance does not run on column (i, j). */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool water (int i, int j) const noexcept { return on && arr(i, j, 0) == 0; }
};

/**
 * @brief Stop when @p mask is given and is not on the grids of @p field.
 * @param who the caller, named in the message
 */
inline void
seb_assert_land_mask_layout (const amrex::iMultiFab* mask, const amrex::FabArrayBase& field,
                             const char* who)
{
    if (mask == nullptr) { return; }
    if (mask->boxArray() != field.boxArray() || mask->DistributionMap() != field.DistributionMap()) {
        amrex::Abort(std::string(who) + ": the land mask is not on the level's 2D grids");
    }
}

/** @brief The view of @p mask on the box of @p mfi; off when @p mask is nullptr. */
[[nodiscard]] inline SEBLandMaskView
seb_land_mask_view (const amrex::iMultiFab* mask, const amrex::MFIter& mfi)
{
    SEBLandMaskView v;
    v.on = (mask != nullptr);
    if (v.on) { v.arr = mask->const_array(mfi); }
    return v;
}

/**
 * @brief True when @p mask holds a water column on any rank; false for nullptr.
 *
 * Read from the mask each time it is asked: the mask can be rewritten during a run (the
 * WW3 coupling does), and a cached answer would then be wrong. One reduction and one MPI
 * collective per call.
 */
[[nodiscard]] inline bool
seb_land_mask_has_water (const amrex::iMultiFab* mask)
{
    return mask != nullptr && mask->min(0) == 0;
}

#endif
