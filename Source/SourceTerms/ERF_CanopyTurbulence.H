#ifndef ERF_CANOPY_TURBULENCE_H_
#define ERF_CANOPY_TURBULENCE_H_

#include <cmath>
#include <limits>
#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

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
 *    (Green 1992, Liu et al. 1996, Sanz 2003). The sink is integrated exactly
 *    over the stage time step, -rho k (1 - exp(-beta_d f |U| dt)) / dt, so it
 *    never removes more than the local k;
 *  - a cap on the turbulent length scale, L <= alpha / f (Harman and Finnigan
 *    2007, alpha = 2 beta^3 with beta = u_* / U_h), applied before the
 *    dissipation, the Richardson number and the stability functions are
 *    evaluated from L.
 *
 * Input names, defaults and formulas follow the ForestDrag.canopy_* inputs of
 * the KLAxell model in kynema-sgf (branch f/winsent of rthedin/amr-wind), so
 * that both codes run the same canopy model.
 */
struct CanopyTurbulence
{
    //! erf.forest_canopy_tke: add the canopy terms to the k-eqn RANS model
    bool enabled = false;
    //! erf.forest_canopy_beta_p: fraction of the drag work that becomes TKE
    amrex::Real beta_p = amrex::Real(1.0);
    //! erf.forest_canopy_beta_d: short-circuit dissipation coefficient
    amrex::Real beta_d = amrex::Real(4.0);
    //! erf.forest_canopy_length_alpha: length scale cap L <= alpha / f;
    //! zero disables the cap and keeps the TKE source
    amrex::Real length_alpha = amrex::Real(0.04);

    /**
     * Read the erf.forest_canopy_* inputs and check their ranges. The checks
     * that need the rest of the configuration (a forest and the k-eqn model)
     * are in check_config.
     */
    void init_params (amrex::ParmParse& pp)
    {
        pp.query("forest_canopy_tke", enabled);
        const bool has_beta_p = pp.query("forest_canopy_beta_p", beta_p);
        const bool has_beta_d = pp.query("forest_canopy_beta_d", beta_d);
        const bool has_alpha  = pp.query("forest_canopy_length_alpha", length_alpha);

        if (!enabled && (has_beta_p || has_beta_d || has_alpha)) {
            amrex::Abort("erf.forest_canopy_beta_p, erf.forest_canopy_beta_d and "
                         "erf.forest_canopy_length_alpha need erf.forest_canopy_tke = true");
        }
        if (!(beta_p >= amrex::Real(0.0) && beta_p <= amrex::Real(1.0))) {
            amrex::Abort("erf.forest_canopy_beta_p must be in [0, 1], the fraction of the "
                         "drag work that becomes TKE; got " + std::to_string(beta_p));
        }
        if (!(beta_d >= amrex::Real(0.0) && std::isfinite(beta_d))) {
            amrex::Abort("erf.forest_canopy_beta_d must be finite and not negative; got " +
                         std::to_string(beta_d));
        }
        if (!(length_alpha >= amrex::Real(0.0) && length_alpha <= amrex::Real(1.0))) {
            amrex::Abort("erf.forest_canopy_length_alpha must be in [0, 1] (0 disables the "
                         "length scale cap); got " + std::to_string(length_alpha));
        }
    }

    /**
     * The canopy terms act on the forest drag field and on the k-eqn RANS
     * model only: abort when either is missing, rather than ignore the input.
     * @param do_forest_drag a forest is configured (erf.forest_file or gridded)
     * @param keqn_on_every_level erf.rans_type = kEqn on every level
     */
    void check_config (bool do_forest_drag, bool keqn_on_every_level) const
    {
        if (!enabled) { return; }
        if (!do_forest_drag) {
            amrex::Abort("erf.forest_canopy_tke = true needs a forest: set erf.forest_file, "
                         "or the gridded erf.forest_lai_file/forest_height_file inputs");
        }
        if (!keqn_on_every_level) {
            amrex::Abort("erf.forest_canopy_tke = true needs erf.rans_type = kEqn on every "
                         "level");
        }
    }
};

namespace canopy_turbulence {

/**
 * Turbulent length scale capped in the canopy, min(length, alpha / f).
 * Cells with f <= 0 are unchanged, and alpha = 0 disables the cap. The test
 * alpha < length f avoids comparing with a large alpha / f. The cap is
 * computed with a denominator that is never zero, because an optimised build
 * may evaluate it before the test (as for 1/theta_ref in ERF_TurbKESources.H);
 * adding the smallest normal Real leaves any f above about 1e-30 unchanged.
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
 * Canopy source of rho k per unit volume and time,
 *   rho f beta_p |U|^3 + rho k expm1(-beta_d f |U| dt) / dt.
 * expm1 is in (-1, 0], so the sink never exceeds rho k / dt.
 * @param rho density
 * @param tke turbulent kinetic energy k
 * @param f_drag forest drag field f = C_d LAD (1/m)
 * @param wind cell-centred wind speed |U|
 * @param beta_p wake production coefficient
 * @param beta_d short-circuit dissipation coefficient
 * @param dt time step of the update the source enters (> 0)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real tke_source (amrex::Real rho, amrex::Real tke, amrex::Real f_drag,
                        amrex::Real wind, amrex::Real beta_p, amrex::Real beta_d,
                        amrex::Real dt)
{
    if (f_drag <= amrex::Real(0.0)) {
        return amrex::Real(0.0);
    }
    return rho * (beta_p * f_drag * wind * wind * wind
                  + tke * std::expm1(-beta_d * f_drag * wind * dt) / dt);
}

} // namespace canopy_turbulence

/**
 * Add the canopy TKE source to the RhoKE component of the staged cell source.
 * @param[in,out] cell_source staged source terms for the conserved state
 * @param[in] S_data conserved state (Rho_comp and RhoKE_comp are read)
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
            const amrex::Real tke  = amrex::max(state(i, j, k, RhoKE_comp) / rho, amrex::Real(0.0));
            const amrex::Real u_cc = amrex::Real(0.5) * (u(i, j, k) + u(i+1, j, k));
            const amrex::Real v_cc = amrex::Real(0.5) * (v(i, j, k) + v(i, j+1, k));
            const amrex::Real w_cc = amrex::Real(0.5) * (w(i, j, k) + w(i, j, k+1));
            const amrex::Real wind = std::sqrt(u_cc*u_cc + v_cc*v_cc + w_cc*w_cc);
            source(i, j, k, RhoKE_comp) +=
                canopy_turbulence::tke_source(rho, tke, f_drag, wind, beta_p, beta_d, dt);
        });
    }
}

#endif
