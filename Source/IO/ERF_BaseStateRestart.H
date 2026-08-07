#ifndef ERF_BASE_STATE_RESTART_H_
#define ERF_BASE_STATE_RESTART_H_

#include <AMReX_Array4.H>

#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"

namespace erf_checkpoint {

// Legacy checkpoints contain p0 but may omit pi0. Current checkpoints carry
// pi0 and must retain their stored value.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void reconstruct_missing_pi0 (const amrex::Array4<amrex::Real>& base,
                              const int ncomp_base_to_read,
                              const amrex::Real rdOcp,
                              const int i, const int j, const int k) noexcept
{
    if (ncomp_base_to_read <= BaseState::pi0_comp) {
        base(i,j,k,BaseState::pi0_comp) =
            getExnergivenP(base(i,j,k,BaseState::p0_comp), rdOcp);
    }
}

} // namespace erf_checkpoint

#endif
