#ifndef ERF_SBM_STAGE_OWNERSHIP_H_
#define ERF_SBM_STAGE_OWNERSHIP_H_

#include <AMReX_GpuQualifiers.H>

namespace erf_sbm {

/** Whether ERF's predictor-to-current slow-state copy owns this component. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
copy_host_slow_component (const int component,
                          const bool sbm_active,
                          const int qc_component,
                          const int qr_component) noexcept
{
    return !sbm_active ||
           (component != qc_component && component != qr_component);
}

} // namespace erf_sbm

#endif
