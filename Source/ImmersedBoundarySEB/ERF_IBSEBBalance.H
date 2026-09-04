#ifndef ERF_IBSEB_BALANCE_H
#define ERF_IBSEB_BALANCE_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <cmath>

/**
 * \file ERF_IBSEBBalance.H
 * \brief Newton solve of the surface energy balance on one building face.
 *
 * Adapted from the ERF-SLUCM branch's ERF_UCMSEBSolver.H (``solve_facet_seb``).
 * The face skin is massless, so at the end of every atmospheric step its
 * temperature :math:`T_s` satisfies
 *
 * \f[
 *   F(T_s) = SW_{abs} + \epsilon\,Q_{ext} + \epsilon\,[LW_{ext} - (1 - f_b)\,\sigma T_s^4]
 *            - C_H\,(T_s/\Pi - \theta_a) - LE - (a\,T_s - b) = 0 ,
 * \f]
 *
 * where every term keeps the sign convention of the module (radiation into
 * the face positive, H, LE and G out of it positive):
 *  - ``SW_abs`` is the absorbed shortwave of the step (phase 2);
 *  - ``Q_ext`` is an external incident flux (a fire's radiation), absorbed
 *    with the face's longwave emissivity by Kirchhoff's law since such
 *    sources are thermal; zero unless something sets it;
 *  - ``LW_ext`` is the incoming longwave from the sky and the ground
 *    (phase 3); the part reflected back by the surrounding walls is
 *    ``f_b sigma T_s^4`` (isothermal surroundings) and folds into the
 *    emission term as the factor ``(1 - f_b)``;
 *  - ``C_H (T_s / Pi - theta_a)`` is the wall function's sensible flux
 *    (phase 4) with the coefficient ``C_H = rho c_p kappa u* / ln(delta/z0h)``
 *    frozen at its value from the wind of the step, ``Pi`` the Exner
 *    function of the fluid cell so the skin temperature enters as a
 *    potential temperature;
 *  - ``LE`` is the latent flux, held at its stored value (zero so far);
 *  - ``a T_s - b`` is the conduction into the slab that the implicit slab
 *    step of phase 5 will produce for a skin temperature ``T_s``
 *    (ibseb::slab_skin_response()), so the balance and the slab agree
 *    exactly instead of lagging the slab's top layer by one step.
 *
 * Newton's method with the analytic Jacobian
 * ``F' = -4 eps (1 - f_b) sigma T_s^3 - C_H / Pi - a`` (every term negative,
 * so the iteration is monotone), a step cap and bounds on the temperature.
 * The bounds are inputs (``erf.ibseb.T_skin_min`` / ``T_skin_max``) because a
 * face under a fire's radiation legitimately exceeds the 380 K of the urban
 * canopy model; a face pinned at a bound leaves a non-zero residual, which
 * the caller stores and the summary counts.
 */
namespace ibseb {

/// Stefan-Boltzmann constant [W/m2/K4].
constexpr amrex::Real SIGMA_SB = 5.670374419e-8;

/**
 * Solve the face balance for the skin temperature.
 *
 * @param[in]  T_start   Starting guess, the skin temperature of the previous step [K].
 * @param[in]  SW_abs    Absorbed shortwave [W/m2].
 * @param[in]  Q_abs     Absorbed external flux ``eps Q_ext`` [W/m2].
 * @param[in]  LW_ext    Incoming longwave from sky and ground [W/m2].
 * @param[in]  emis      Face emissivity.
 * @param[in]  f_bldg    Building view fraction of the face.
 * @param[in]  C_H       Sensible coefficient ``rho c_p kappa u* / ln(delta/z0h)`` [W/m2/K].
 * @param[in]  exner     Exner function of the fluid cell, ``T_a / theta_a``.
 * @param[in]  theta_air Potential temperature of the fluid cell [K].
 * @param[in]  LE        Latent flux, frozen [W/m2].
 * @param[in]  slab_a, slab_b  Slab response ``G = a T_s - b`` [W/m2/K], [W/m2].
 * @param[in]  T_min, T_max    Bounds on the skin temperature [K].
 * @param[in]  max_step  Largest Newton step [K].
 * @param[in]  tol       Convergence: stop when a step changes T_s by less than this [K].
 * @param[in]  max_iter  Iteration cap.
 * @param[out] n_iter    Iterations taken.
 * @param[out] resid     ``|F|`` at the returned temperature [W/m2]; non-zero
 *                       when a bound was hit or the cap was reached.
 * @return The skin temperature [K], within the bounds.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solve_skin_balance (amrex::Real T_start,
                                amrex::Real SW_abs, amrex::Real Q_abs, amrex::Real LW_ext,
                                amrex::Real emis, amrex::Real f_bldg,
                                amrex::Real C_H, amrex::Real exner, amrex::Real theta_air,
                                amrex::Real LE, amrex::Real slab_a, amrex::Real slab_b,
                                amrex::Real T_min, amrex::Real T_max,
                                amrex::Real max_step, amrex::Real tol, int max_iter,
                                int& n_iter, amrex::Real& resid)
{
    using amrex::Real;
    const Real e_eff = emis * (Real(1.0) - f_bldg) * SIGMA_SB;   // effective emission coefficient
    auto residual = [&](Real T) {
        const Real T4 = T * T * T * T;
        return SW_abs + Q_abs + emis * LW_ext - e_eff * T4
             - C_H * (T / exner - theta_air) - LE - (slab_a * T - slab_b);
    };
    Real T = amrex::max(amrex::min(T_start, T_max), T_min);
    n_iter = 0;
    Real F = residual(T);
    for (int it = 0; it < max_iter; ++it) {
        ++n_iter;
        const Real Fp = -Real(4.0) * e_eff * T * T * T - C_H / exner - slab_a;
        Real step = -F / Fp;                                     // Fp < 0 always
        step = amrex::max(amrex::min(step, max_step), -max_step);
        const Real T_new = amrex::max(amrex::min(T + step, T_max), T_min);
        const Real dT = std::abs(T_new - T);
        T = T_new;
        F = residual(T);
        if (dT < tol) break;
    }
    resid = std::abs(F);
    return T;
}

} // namespace ibseb

#endif
