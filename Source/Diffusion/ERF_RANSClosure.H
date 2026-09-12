/** \file ERF_RANSClosure.H
 *
 *  Closure relations of the one-equation k RANS model of
 *  Axell & Liungman (2001), Environ. Fluid Mech. 1, 71-106 (AL01), as
 *  device-inline functions so that ComputeTurbulentViscosityRANS and the
 *  unit tests evaluate exactly the same expressions.
 */
#ifndef ERF_RANS_CLOSURE_H_
#define ERF_RANS_CLOSURE_H_

#include <cmath>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include "ERF_Constants.H"

namespace AL01 {

/** Geometric length scale with the harmonic cap (AL01 Eq. 22 plus limiter).
 *  @param l_g_raw kappa times (wall distance + z0)
 *  @param l_g_max upper bound, about kappa times 0.1 zi
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real geom_length (amrex::Real l_g_raw, amrex::Real l_g_max)
{
    return l_g_max * l_g_raw / (l_g_max + l_g_raw);
}

/** Burchard & Petersen smoothing of Rt below Rt_crit, which maps
 *  (-inf, Rt_crit) onto (Rt_min, Rt_crit) so the stability functions stay
 *  finite in strongly unstable air:
 *      Rt_s = Rt - (Rt - Rt_crit)^2 / (Rt + Rt_min - 2 Rt_crit).
 *  Written with x = Rt - Rt_crit and a = Rt_min - Rt_crit as
 *      Rt_s = Rt_crit + a x / (x + a),
 *  which is the same function without the cancellation of two O(|Rt|)
 *  terms that returned garbage for |Rt| beyond about 1e15 (tiny k under a
 *  strong unstable N^2). Limits: Rt_crit as x -> 0, Rt_min as x -> -inf. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real smooth_Rt (amrex::Real Rt, amrex::Real Rt_crit, amrex::Real Rt_min)
{
    if (Rt >= Rt_crit) { return Rt; }
    const amrex::Real x = Rt - Rt_crit;      // < 0
    const amrex::Real a = Rt_min - Rt_crit;  // < 0 by input validation
    return Rt_crit + a * x / (x + a);
}

/** Turbulent length scale: neutral (l_g), stable (Eq. 26) or unstable
 *  (Eq. 28).
 *
 *  Unstable: Eq. 28 is Eq. 26 rewritten with the dissipation of Eq. 19 on
 *  the right-hand side. Iterating it (dissipation from the new length,
 *  length from the new Rt) is the fixed-point map of Eq. 26 with N^2 < 0,
 *  which has no fixed point once l_g^2 |N^2| exceeds Cb^2 k (AL01, p. 78),
 *  so the length grows without bound in strongly convective, weakly
 *  turbulent air. Here Eq. 28 is evaluated once with the dissipation of the
 *  geometric length and the Burchard & Petersen smoothed Rt, as in the
 *  Kynema KLAxell implementation, which bounds the unstable length by
 *  l_g sqrt(1 + Cmu0^6 |Rt_min| / Cb^2), about 1.31 l_g for the defaults.
 *
 *  @param l_g capped geometric length
 *  @param N2 Brunt-Vaisala frequency squared
 *  @param tke turbulent kinetic energy (floored, > 0)
 *  @param Cmu0_pow3 Cmu0 cubed
 *  @param inv_Cb_sq one over Cb squared
 *  @param Rt_crit, Rt_min smoothing parameters
 *  @param eps threshold below which N2 counts as neutral
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real turb_length (amrex::Real l_g, amrex::Real N2, amrex::Real tke,
                         amrex::Real Cmu0_pow3, amrex::Real inv_Cb_sq,
                         amrex::Real Rt_crit, amrex::Real Rt_min,
                         amrex::Real eps)
{
    amrex::Real length;
    if (std::abs(N2) <= eps) {
        length = l_g;
    } else if (N2 > eps) {
        // Stable (AL01, Eqn. 26)
        length = std::sqrt(one /
                (one / (l_g * l_g) + inv_Cb_sq * N2 / tke));
    } else {
        // Unstable (AL01, Eqn. 28) with Rt from the geometric length
        // (Eq. 29 with Eq. 19: Rt = l_g^2 N^2 / (k Cmu0^6)), smoothed
        const amrex::Real Rt = smooth_Rt(l_g*l_g * N2 / (tke * Cmu0_pow3 * Cmu0_pow3),
                                         Rt_crit, Rt_min);
        length = l_g * std::sqrt(one - Cmu0_pow3*Cmu0_pow3 * inv_Cb_sq * Rt);
    }
    return length;
}

/** Upper bound of the unstable length for given smoothing parameters. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real unstable_length_bound (amrex::Real l_g, amrex::Real Cmu0_pow3,
                                   amrex::Real inv_Cb_sq, amrex::Real Rt_min)
{
    return l_g * std::sqrt(one - Cmu0_pow3*Cmu0_pow3 * inv_Cb_sq * Rt_min);
}

/** Dissipation rate per unit volume, rho Cmu0^3 k^1.5 / L (AL01 Eq. 19). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dissipation (amrex::Real rho, amrex::Real Cmu0_pow3,
                         amrex::Real tke, amrex::Real length)
{
    return rho * Cmu0_pow3 * std::pow(tke,amrex::Real(1.5)) / length;
}

/** Turbulent Richardson number, Eq. 29 combined with Eq. 19:
 *  Rt = k^2 N^2 / eps^2 = L^2 N^2 / (k Cmu0^6). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real richardson (amrex::Real length, amrex::Real N2, amrex::Real tke,
                        amrex::Real Cmu0_pow3)
{
    return length*length * N2 / (tke * Cmu0_pow3 * Cmu0_pow3);
}


/** Momentum stability function (AL01 Eq. 31). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cmu (amrex::Real Rt, amrex::Real Cmu0)
{
    return (Cmu0 + amrex::Real(0.108)*Rt)
         / (one + amrex::Real(0.308)*Rt + amrex::Real(0.00837)*Rt*Rt);
}

/** Scalar stability function (AL01 Eq. 32). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cmu_prime (amrex::Real Rt, amrex::Real Cmu0)
{
    return Cmu0 / (1 + amrex::Real(0.277)*Rt);
}

/** Lowest admissible Rt_min: the poles of cmu' (Rt = -1/0.277) and of cmu
 *  (larger root of 1 + 0.308 Rt + 0.00837 Rt^2) both sit near -3.6. */
constexpr amrex::Real Rt_min_lower_bound = amrex::Real(-3.6);

} // namespace AL01

#endif
