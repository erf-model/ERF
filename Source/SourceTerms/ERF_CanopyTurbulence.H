#ifndef ERF_CANOPY_TURBULENCE_H_
#define ERF_CANOPY_TURBULENCE_H_

#include <cmath>
#include <limits>
#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <ERF_CanopyTurbStruct.H>
#include <ERF_IndexDefines.H>

/**
 * Canopy terms of the one-equation k RANS model (erf.rans_type = kEqn).
 *
 * The forest drag removes mean kinetic energy from the flow. Part of it becomes
 * turbulence in the wakes of leaves and branches, at scales far below the mesh,
 * and dissipates quickly (the spectral short cut). With
 * erf.forest_canopy_tke = true, two terms act where the forest drag field
 * f = C_d LAD of ForestDrag is positive:
 *
 *  - a TKE source, S = rho f (beta_p |U|^3 - beta_d |U| k)
 *    (Green 1992, Liu et al. 1996, Sanz 2003). Production and sink are
 *    integrated together and exactly over each RK stage, from the k at the
 *    start of the step (see tke_source), so the steady k is the exact balance
 *    at any time step and the sink never removes more than that k;
 *  - a cap on the turbulent length scale, L <= alpha / f (Harman and Finnigan
 *    2007, alpha = 2 beta^3 with beta = u_* / U_h), applied before the
 *    dissipation, the Richardson number and the stability functions are
 *    evaluated from L.
 *
 * Input names, defaults and formulas follow the ForestDrag.canopy_* inputs of
 * the KLAxell model in kynema-sgf (branch f/winsent of rthedin/amr-wind). In
 * neutral air both codes run the same canopy model; Forest.rst lists where the
 * time integration and the stratified closure differ.
 */
namespace canopy_turbulence {

/**
 * Turbulent length scale capped in the canopy, min(length, alpha / f).
 * Cells with f <= 0 are unchanged, and alpha = 0 disables the cap. The test
 * alpha < length f avoids comparing with a large alpha / f. The cap is
 * computed with a denominator that is never zero, because an optimised build
 * may evaluate it before the test (as for 1/theta_ref in ERF_TurbKESources.H);
 * adding the smallest normal Real leaves any f above min / epsilon unchanged.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real limited_length (amrex::Real length, amrex::Real f_drag, amrex::Real alpha)
{
    const amrex::Real cap = alpha / (amrex::max(f_drag, amrex::Real(0.0))
                                     + std::numeric_limits<amrex::Real>::min());
    const bool capped = (alpha > amrex::Real(0.0)) && (f_drag > amrex::Real(0.0))
                        && (alpha < length * f_drag);
    return capped ? cap : length;
}

/**
 * Canopy source of rho k per unit volume and time for an RK stage that updates
 * rho k from its value at the start of the step, (rho k)^n, over dt:
 *   rho beta_p f |U|^3 phi1(x) - (rho k)^n (1 - exp(-x)) / dt,
 * with x = beta_d f |U| dt and phi1(x) = (1 - exp(-x)) / x. The stage update
 * (rho k)^n + dt S is then the exact solution of d(rho k)/dt = rho P - c rho k
 * with P = beta_p f |U|^3 and c = beta_d f |U| frozen over the stage,
 *   (rho k)^n exp(-x) + rho (P / c) (1 - exp(-x)),
 * so the steady state is the exact balance rho P / c at any dt, and the sink
 * never removes more than (rho k)^n. The RK3 stages (dt/3, dt/2, dt from the
 * start of the step) and the trapezoidal anelastic RK2 (both stages over dt)
 * both update from (rho k)^n.
 * @param rho density of the stage
 * @param rhok_old rho k at the start of the step
 * @param f_drag forest drag field f = C_d LAD (1/m)
 * @param wind cell-centred wind speed |U| of the stage
 * @param beta_p wake production coefficient
 * @param beta_d short-circuit dissipation coefficient
 * @param dt time step of the update the source enters (> 0)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real tke_source (amrex::Real rho, amrex::Real rhok_old, amrex::Real f_drag,
                        amrex::Real wind, amrex::Real beta_p, amrex::Real beta_d,
                        amrex::Real dt)
{
    if (f_drag <= amrex::Real(0.0)) {
        return amrex::Real(0.0);
    }
    const amrex::Real x = beta_d * f_drag * wind * dt;
    // 1 - exp(-x), in [0, 1) for x >= 0
    const amrex::Real decay = -std::expm1(-x);
    // phi1 with a divisor that is never zero (an optimised build may evaluate it
    // for x = 0); at the smallest normal Real, phi1 rounds to its limit 1
    const amrex::Real x_safe = amrex::max(x, std::numeric_limits<amrex::Real>::min());
    const amrex::Real phi1 = -std::expm1(-x_safe) / x_safe;
    return rho * beta_p * f_drag * wind * wind * wind * phi1 - rhok_old * decay / dt;
}

} // namespace canopy_turbulence

/**
 * Add the canopy TKE source to the RhoKE component of the staged cell source.
 * @param[in,out] cell_source staged source terms for the conserved state
 * @param[in] S_data conserved state of the stage (Rho_comp is read)
 * @param[in] S_old conserved state at the start of the step (RhoKE_comp is read)
 * @param[in] xvel x-face velocity
 * @param[in] yvel y-face velocity
 * @param[in] zvel z-face velocity
 * @param[in] forest_drag cell-centred forest drag field f = C_d LAD
 * @param[in] canopy canopy turbulence inputs
 * @param[in] dt time step of the stage the source enters
 */
AMREX_INLINE void
AddCanopyTKESources (amrex::MultiFab&        cell_source,
                     const amrex::MultiFab&  S_data,
                     const amrex::MultiFab&  S_old,
                     const amrex::MultiFab&  xvel,
                     const amrex::MultiFab&  yvel,
                     const amrex::MultiFab&  zvel,
                     const amrex::MultiFab&  forest_drag,
                     const CanopyTurbulence& canopy,
                     amrex::Real             dt)
{
    BL_PROFILE("AddCanopyTKESources");

    const amrex::Real beta_p = canopy.beta_p;
    const amrex::Real beta_d = canopy.beta_d;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(cell_source, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        const auto& source = cell_source.array(mfi);
        const auto& state  = S_data.const_array(mfi);
        const auto& state_old = S_old.const_array(mfi);
        const auto& u      = xvel.const_array(mfi);
        const auto& v      = yvel.const_array(mfi);
        const auto& w      = zvel.const_array(mfi);
        const auto& f_arr  = forest_drag.const_array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real f_drag = f_arr(i, j, k);
            if (f_drag <= amrex::Real(0.0)) {
                return;
            }
            const amrex::Real rho  = state(i, j, k, Rho_comp);
            const amrex::Real rhok_old = amrex::max(state_old(i, j, k, RhoKE_comp), amrex::Real(0.0));
            const amrex::Real u_cc = amrex::Real(0.5) * (u(i, j, k) + u(i+1, j, k));
            const amrex::Real v_cc = amrex::Real(0.5) * (v(i, j, k) + v(i, j+1, k));
            const amrex::Real w_cc = amrex::Real(0.5) * (w(i, j, k) + w(i, j, k+1));
            const amrex::Real wind = std::sqrt(u_cc*u_cc + v_cc*v_cc + w_cc*w_cc);
            source(i, j, k, RhoKE_comp) +=
                canopy_turbulence::tke_source(rho, rhok_old, f_drag, wind, beta_p, beta_d, dt);
        });
    }
}

#endif
