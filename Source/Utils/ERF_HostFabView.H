#ifndef ERF_HOST_FAB_VIEW_H
#define ERF_HOST_FAB_VIEW_H

#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

/**
 * @file ERF_HostFabView.H
 * @brief Host-readable view of one FAB of a MultiFab.
 *
 * AMReX leaves `the_arena_is_managed` at false by default, so under CUDA a
 * MultiFab's data lives in device memory and an `Array4` obtained from it
 * cannot be dereferenced on the host at all -- a plain `for` loop or
 * `amrex::LoopOnCpu` over it faults, or silently reads garbage if the address
 * happens to be mapped. `Gpu::streamSynchronize()` does not help: it orders
 * kernels, it does not migrate memory.
 *
 * Most of the time the answer is to do the work in a kernel. This helper is for
 * the cases where it genuinely cannot be: serial trajectory integration with a
 * host RNG, scatter-accumulation into named bins, and gathering a field into a
 * flat host vector ahead of an MPI reduce or a file write.
 *
 * On CPU builds the view is a plain alias for the device FAB, so host runs pay
 * nothing for it.
 *
 * Typical use:
 * @code
 *   for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
 *       const ERFHostFabView view(mf[mfi]);
 *       auto arr = view.array();
 *       amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
 *           total += arr(i,j,k);
 *       });
 *   }
 * @endcode
 */
struct ERFHostFabView
{
#ifdef AMREX_USE_GPU
    amrex::FArrayBox fab;

    explicit ERFHostFabView (const amrex::FArrayBox& dev_fab)
        : fab(dev_fab.box(), dev_fab.nComp(), amrex::The_Pinned_Arena())
    {
        // The staging FAB carries the device FAB's own box, so an Array4 taken
        // from it is indexed with the same global (i,j,k).
        amrex::Gpu::dtoh_memcpy_async(fab.dataPtr(), dev_fab.dataPtr(),
                                      dev_fab.size() * sizeof(amrex::Real));
        amrex::Gpu::streamSynchronize();
    }

    [[nodiscard]] amrex::Array4<const amrex::Real> array () const
    {
        return fab.const_array();
    }
#else
    const amrex::FArrayBox* dev_fab;

    explicit ERFHostFabView (const amrex::FArrayBox& d) : dev_fab(&d) {}

    [[nodiscard]] amrex::Array4<const amrex::Real> array () const
    {
        return dev_fab->const_array();
    }
#endif
};

#endif /* ERF_HOST_FAB_VIEW_H */
