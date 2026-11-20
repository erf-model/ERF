#ifndef ERF_IBSEB_WALLFUNCTION_H
#define ERF_IBSEB_WALLFUNCTION_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <cmath>

/**
 * \file ERF_IBSEBWallFunction.H
 * \brief Explicit stability factors of the face wall function.
 */
namespace ibseb {

/**
 * Louis (1979) stability factors on the bulk Richardson number.
 *
 * The exchange coefficients are the neutral ones times a factor:
 * ``u*^2 = a^2 U^2 F_m`` and ``u* theta* = a_h^2 U dtheta F_h``. Unstable
 * (``Ri_b < 0``):
 * ``F = 1 - b Ri_b / (1 + C* a^2 b sqrt(|Ri_b| z / z0))`` with ``b = 9.4``,
 * ``C* = 7.4`` for momentum and ``5.3`` for heat. Stable:
 * ``F_m = F_h = 1 / (1 + b/2 Ri_b)^2``. Both are 1 at ``Ri_b = 0`` and
 * continuous there. No iteration: the factors come straight from the
 * resolved state.
 *
 * @param[in]  Rib        Bulk Richardson number between the surface and the reference height.
 * @param[in]  z_over_z0  Reference height over the momentum roughness length.
 * @param[in]  a2         Neutral momentum drag coefficient ``(kappa / ln(z / z0))^2``.
 * @param[out] Fm         Momentum factor.
 * @param[out] Fh         Heat factor.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void louis_factors (amrex::Real Rib, amrex::Real z_over_z0, amrex::Real a2,
                    amrex::Real& Fm, amrex::Real& Fh) noexcept
{
    const amrex::Real one  = amrex::Real(1.0);
    const amrex::Real b    = amrex::Real(9.4);
    const amrex::Real cs_m = amrex::Real(7.4);
    const amrex::Real cs_h = amrex::Real(5.3);
    // The square root is formed from a clipped argument so a speculatively
    // evaluated unstable arm never takes the root of a negative number.
    const amrex::Real s = std::sqrt(amrex::max(-Rib, amrex::Real(0.0)) * z_over_z0);
    if (Rib < amrex::Real(0.0)) {
        Fm = one - b * Rib / (one + cs_m * a2 * b * s);
        Fh = one - b * Rib / (one + cs_h * a2 * b * s);
    } else {
        const amrex::Real q = one + amrex::Real(0.5) * b * Rib;
        Fm = one / (q * q);
        Fh = Fm;
    }
}

} // namespace ibseb

#endif
