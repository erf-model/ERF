#ifndef ERF_SHOC_CONSTANTS_H_
#define ERF_SHOC_CONSTANTS_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace shoc::constants {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real min_tke () noexcept { return amrex::Real(4.0e-4); }

// Relative margin above min_tke() used to decide whether a level is still
// turbulently active. Above the boundary layer both the explicit TKE update
// and the implicit TKE solve clamp TKE to exactly min_tke(), so a bare
// "tke > min_tke()" test is settled by the last bit of the tridiagonal
// solution. Host and device builds differ there by one ulp, because the CUDA
// build contracts to FMA and the host build does not, and that flips the
// diagnosed active top from step to step. The margin keeps the comparison off
// that measure-zero boundary.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real active_tke_rel_margin () noexcept
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(1.0e-4);
#else
    return amrex::Real(1.0e-6);
#endif
}

// Smallest TKE that counts as turbulently active. See active_tke_rel_margin().
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real active_tke_threshold () noexcept
{
    return min_tke() * (amrex::Real(1.0) + active_tke_rel_margin());
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real latent_ice () noexcept { return amrex::Real(3.34e5); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real min_temp () noexcept { return amrex::Real(180.0); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real freezing_temp () noexcept { return amrex::Real(273.15); }

} // namespace shoc::constants

#endif