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
 * Always on:
 * - An Obukhov length derived from u* and the heat flux is held at |L| >= z_ref / 100, and every
 *   zeta = z / L at |zeta| <= 100, as on flat ground (ERF_MOSTStress.H). A prescribed
 *   erf.if_Olen is used as given: only its zeta is held, so one input means one L on every face
 *   and at every height.
 *
 * Opt-in (WRF's revised surface layer, sf_sfclayrev):
 * - erf.if_stability_wind_floor (default 0, WRF and the flat-ground WSMIN use 0.1 m/s) floors the
 *   wind behind the friction velocity of the stability estimate.
 * - erf.if_psi_cap_factor (default 1, WRF uses 0.9) caps psi_m at factor ln(z / z0) in the
 *   velocity target and before it forms the friction velocity of the momentum law. The
 *   temperature forcing keeps caps of ln(z / z0): its surface-flux branch forms theta* from an
 *   Obukhov length that is not consistent with u*, which those caps keep from reaching the target
 *   in strongly unstable air.
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
     * An Obukhov length derived from u* and the heat flux, with |L| no smaller than
     * z_ref / zeta_max(), so that z / L lies within +-zeta_max() for every z <= z_ref and the
     * zero L of a zero u* is never divided by. The sign is kept: L = -0 is the unstable bound and
     * L = +0 the stable one. Unchanged, bit for bit, when |L| >= z_ref / zeta_max(). Not for a
     * prescribed L, whose bound would then depend on the local cell size; see bounded_zeta.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real bounded_obukhov_length (amrex::Real Olen, amrex::Real z_ref) noexcept
    {
        return std::copysign(amrex::max(std::abs(Olen), z_ref / zeta_max()), Olen);
    }

    /**
     * zeta = z / L held within +-zeta_max(), for an L that is not zero (a prescribed erf.if_Olen,
     * which start-up refuses at 0, or one from bounded_obukhov_length). Unchanged, bit for bit,
     * when |z / L| <= zeta_max(); a NaN stays NaN.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real bounded_zeta (amrex::Real z, amrex::Real Olen) noexcept
    {
        // z / L first in both, so a NaN propagates instead of becoming a bound
        return amrex::min(amrex::max(z / Olen, -zeta_max()), zeta_max());
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
