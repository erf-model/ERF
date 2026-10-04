#ifndef ERF_IMMERSED_WALL_STABILITY_H_
#define ERF_IMMERSED_WALL_STABILITY_H_

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

/**
 * Bounds on the stability estimate of the immersed-forcing wall law. Without them a zero
 * tangential wind one cell from the wall gives u* = 0, an Obukhov length of zero and zeta = z/0
 * (erf-model/ERF#4016).
 *
 * - The wind behind the friction velocity of the stability estimate is floored at
 *   erf.if_stability_wind_floor (default 0.1 m/s, the flat-ground surface layer's WSMIN in
 *   ERF_MOSTStress.H and WRF's floor); 0 leaves the wind as it is.
 * - The Obukhov length is held at |L| >= z_ref / 100, so |zeta| <= 100 as on flat ground
 *   (ERF_MOSTStress.H).
 * - psi_m and psi_h are capped at erf.if_psi_cap_factor ln(z / z0) (default 0.9, WRF's cap in
 *   sf_sfclayrev), in the targets and, for psi_m, in the friction velocity as WRF does; 1 is a cap
 *   of ln(z / z0) in the targets and none in the friction velocity.
 */
namespace ib_stability
{
    /** Largest |zeta| = |z / L| the wall law uses. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real zeta_max () noexcept { return amrex::Real(100.0); }

    /** The wind speed the stability estimate uses: the speed, floored at wind_floor (>= 0). */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real floored_wind (amrex::Real speed, amrex::Real wind_floor) noexcept
    {
        return amrex::max(speed, wind_floor);
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
     * The largest psi (momentum or heat) the wall law uses at height z: factor ln(z / z0) for
     * 0 < factor <= 1, and never more than ln(z / z0). For z > z0 a factor below 1 keeps
     * ln(z / z0) - psi at no less than (1 - factor) ln(z / z0) > 0, so strongly unstable air
     * still carries momentum and heat between the wall and the cell above; a cap of ln(z / z0)
     * itself makes the target zero (momentum) or the value above (heat). For z <= z0, where the
     * log law does not hold, the cap stays ln(z / z0).
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real psi_cap (amrex::Real z, amrex::Real z0, amrex::Real factor) noexcept
    {
        const amrex::Real ln = std::log(z / z0);
        return amrex::min(factor * ln, ln);
    }

    /** psi held at or below psi_cap(z, z0, factor). */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real capped_psi (amrex::Real psi, amrex::Real z, amrex::Real z0, amrex::Real factor) noexcept
    {
        return amrex::min(psi, psi_cap(z, z0, factor));
    }

    /**
     * psi_m in the friction velocity u* = U kappa / (ln(z / z0) - psi_m) at height z: held at or
     * below factor ln(z / z0), as WRF caps PSIM before it forms UST, so the denominator stays at
     * least (1 - factor) ln(z / z0) > 0 and strongly unstable air keeps a positive u*. Uncapped for
     * factor = 1 (a cap of ln(z / z0) would zero the denominator) and for z <= z0, where the log law
     * does not hold; there the wall law's u* clamp handles a denominator of either sign.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real psi_m_for_ustar (amrex::Real psi_m, amrex::Real z, amrex::Real z0, amrex::Real factor) noexcept
    {
        const amrex::Real ln = std::log(z / z0);
        return (factor < amrex::Real(1.0) && ln > amrex::Real(0.0)) ? amrex::min(psi_m, factor * ln) : psi_m;
    }

    /**
     * The friction velocity held within [0, 2] m/s, a NaN going to 0: the wall law's guard
     * against a log-law denominator that cancels or changes sign.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real clamped_ustar (amrex::Real ustar) noexcept
    {
        if (!(ustar > amrex::Real(0.0) && !std::isnan(ustar))) { ustar = amrex::Real(0.0); }
        if (!(ustar < amrex::Real(2.0) && !std::isnan(ustar))) { ustar = amrex::Real(2.0); }
        return ustar;
    }
}

#endif
