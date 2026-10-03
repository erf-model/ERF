#ifndef ERF_IMMERSED_WALL_STABILITY_H_
#define ERF_IMMERSED_WALL_STABILITY_H_

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

/**
 * Bounds on the stability estimate of the immersed-forcing wall law. The wind floor and the zeta
 * bound are the flat-ground surface layer's (ERF_MOSTStress.H: the wind speed floored at
 * WSMIN = 0.1 m/s, and zeta = z/L held within +-100); without them a zero tangential wind one cell
 * from the wall gives u* = 0, an Obukhov length of zero and zeta = z/0 (erf-model/ERF#4016). The
 * psi_h cap is WRF's.
 */
namespace ib_stability
{
    /** Wind speed [m/s] below which the stability estimate does not go. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real wind_floor () noexcept { return amrex::Real(0.1); }

    /** Largest |zeta| = |z / L| the wall law uses. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real zeta_max () noexcept { return amrex::Real(100.0); }

    /** The wind speed the stability estimate uses: the speed, floored at wind_floor(). */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real floored_wind (amrex::Real speed) noexcept
    {
        return amrex::max(speed, wind_floor());
    }

    /**
     * The Obukhov length with |L| no smaller than z_ref / zeta_max(), so that z / L lies within
     * +-zeta_max() for every z <= z_ref. The sign is kept: L = -0 is the unstable bound and
     * L = +0 the stable one. Unchanged, bit for bit, when |L| >= z_ref / zeta_max().
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real bounded_obukhov_length (amrex::Real Olen, amrex::Real z_ref) noexcept
    {
        return std::copysign(amrex::max(std::abs(Olen), z_ref / zeta_max()), Olen);
    }

    /**
     * The largest psi_h the wall law uses at height z: 0.9 ln(z / z0), the cap of WRF's revised
     * surface layer (sf_sfclayrev). It keeps ln(z / z0) - psi_h at no less than 0.1 ln(z / z0), so
     * strongly unstable air still carries heat between the wall and the cell above; a cap of
     * ln(z / z0) itself flattens the profile and the wall heat transfer stops.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real psi_h_cap (amrex::Real z, amrex::Real z0) noexcept
    {
        return amrex::Real(0.9) * std::log(z / z0);
    }
}

#endif
