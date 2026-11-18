#ifndef ERF_CHECKPOINT_LANDMASK_H_
#define ERF_CHECKPOINT_LANDMASK_H_

#include <AMReX_iMultiFab.H>
#include <AMReX_GpuLaunch.H>

namespace erf_checkpoint_landmask {

// Older checkpoints can mark building footprints as immersed on levels where
// the current resolution cutoff disables building forcing. Those footprints
// are land; restore their ordinary land mask before MOST and the LSM read it.
inline void restore_coarse_building_landmask (amrex::iMultiFab& mask)
{
    for (amrex::MFIter mfi(mask); mfi.isValid(); ++mfi) {
        const auto arr = mask.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (arr(i,j,k) == 2) { arr(i,j,k) = 1; }
        });
    }
}

} // namespace erf_checkpoint_landmask

#endif
