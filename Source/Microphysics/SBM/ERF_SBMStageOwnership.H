#ifndef ERF_SBM_STAGE_OWNERSHIP_H_
#define ERF_SBM_STAGE_OWNERSHIP_H_

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>

namespace erf_sbm {

/** Perform ERF's predictor-to-current copy for one slow state component. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
copy_host_slow_component (const amrex::Array4<amrex::Real>& current,
                          const amrex::Array4<const amrex::Real>& predictor,
                          const int i,
                          const int j,
                          const int k,
                          const int component,
                          const bool sbm_active,
                          const int qc_component,
                          const int qr_component) noexcept
{
    const bool sbm_projected_liquid =
        sbm_active &&
        (component == qc_component || component == qr_component);
    if (!sbm_projected_liquid) {
        current(i, j, k, component) = predictor(i, j, k, component);
    }
}

} // namespace erf_sbm

#endif
