#ifndef ERF_BASE_STATE_HSE_H_
#define ERF_BASE_STATE_HSE_H_

#include <AMReX_Array4.H>

#include "ERF_EOS.H"

namespace erf_base_state {

// Store an HSE pressure and its Exner value as one production operation.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_pressure_and_exner (const amrex::Array4<amrex::Real>& pressure,
                             const amrex::Array4<amrex::Real>& exner,
                             const int i, const int j, const int k,
                             const amrex::Real pressure_value,
                             const amrex::Real rdOcp) noexcept
{
    pressure(i,j,k) = pressure_value;
    exner(i,j,k) = getExnergivenP(pressure_value, rdOcp);
}

} // namespace erf_base_state

#endif
