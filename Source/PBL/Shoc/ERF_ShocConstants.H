#ifndef ERF_SHOC_CONSTANTS_H_
#define ERF_SHOC_CONSTANTS_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace shoc::constants {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real min_tke () noexcept { return amrex::Real(4.0e-4); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real latent_ice () noexcept { return amrex::Real(3.34e5); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real min_temp () noexcept { return amrex::Real(180.0); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real freezing_temp () noexcept { return amrex::Real(273.15); }

} // namespace shoc::constants

#endif