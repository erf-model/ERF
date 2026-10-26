#ifndef ERF_HOST_FAB_OPS_H_
#define ERF_HOST_FAB_OPS_H_

#include <AMReX_BaseFab.H>
#include <AMReX_Box.H>
#include <AMReX_Loop.H>

namespace erf_host_fab_ops {

/**
 * \brief Host-side whole-fab fill and copy that are safe on every backend.
 *
 * BaseFab::setVal<RunOn::Host> and BaseFab::copy<RunOn::Host> expand (via the
 * FLAG macros in AMReX_GpuLaunch.H) to a device launch when run_on is Device,
 * and to LoopConcurrentOnCpu on CUDA/HIP in the else -- but amrex::Abort on
 * SYCL builds unconditionally.  These helpers run the same
 * LoopConcurrentOnCpu on every backend, so they are bit-identical to the
 * RunOn::Host forms on CUDA/HIP and to setVal/copy on CPU, and they work on
 * SYCL where the originals abort.
 *
 * The fab data must be host-accessible (managed, pinned, or CPU arena).
 */

//! Fill every component of \p fab over its whole box with \p value.
inline void
host_setval (amrex::BaseFab<amrex::Real>& fab, amrex::Real value) noexcept
{
    auto const f_arr = fab.array();
    amrex::LoopConcurrentOnCpu(fab.box(), f_arr.nComp(),
        [=] (int i, int j, int k, int n) noexcept
    {
        f_arr(i,j,k,n) = value;
    });
}

//! Copy \p ncomp components from \p src to \p dst over aligned whole-fab boxes.
inline void
host_copy (amrex::BaseFab<amrex::Real>& dst, amrex::BaseFab<amrex::Real> const& src,
           int ncomp) noexcept
{
    auto const dst_arr = dst.array();
    auto const src_arr = src.const_array();
    amrex::LoopConcurrentOnCpu(src.box(), ncomp,
        [=] (int i, int j, int k, int n) noexcept
    {
        dst_arr(i,j,k,n) = src_arr(i,j,k,n);
    });
}

} // namespace erf_host_fab_ops

#endif // ERF_HOST_FAB_OPS_H_
