#ifndef ERF_CANOPY_TURB_STRUCT_H_
#define ERF_CANOPY_TURB_STRUCT_H_

#include <cmath>
#include <string>

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

/**
 * Inputs of the canopy terms of the k-eqn RANS model (erf.forest_canopy_*).
 * The terms themselves are in Source/SourceTerms/ERF_CanopyTurbulence.H.
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

        // Unused coefficients are left in place (with a warning) so that a deck
        // can switch the canopy terms on and off with erf.forest_canopy_tke alone;
        // their ranges are still checked below
        if (!enabled && (has_beta_p || has_beta_d || has_alpha)) {
            amrex::Warning("erf.forest_canopy_beta_p, erf.forest_canopy_beta_d and "
                           "erf.forest_canopy_length_alpha are not used without "
                           "erf.forest_canopy_tke = true");
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
        if (enabled && beta_p == amrex::Real(0.0) && beta_d == amrex::Real(0.0) &&
            length_alpha == amrex::Real(0.0)) {
            amrex::Abort("erf.forest_canopy_tke = true with erf.forest_canopy_beta_p, "
                         "erf.forest_canopy_beta_d and erf.forest_canopy_length_alpha all "
                         "zero would change nothing; set erf.forest_canopy_tke = false");
        }
    }

    /**
     * The canopy terms act on the forest drag field and on the k-eqn RANS
     * model only: abort when either is missing, rather than ignore the input,
     * and when the length cap would meet the explicit TKE dissipation.
     * @param do_forest_drag a forest is configured (erf.forest_file or gridded)
     * @param keqn_on_every_level erf.rans_type = kEqn on every level
     * @param implicit_diss_on_every_level erf.implicit_tke_dissipation on every level
     */
    void check_config (bool do_forest_drag, bool keqn_on_every_level,
                       bool implicit_diss_on_every_level) const
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
        // The cap makes the explicit dissipation rho Cmu0^3 k^1.5 / L stiff: with
        // L = alpha / f its rate Cmu0^3 sqrt(k) f / alpha exceeds 1/dt, and k can go
        // negative in one step, at time steps of a fraction of a second in a dense
        // canopy (f = 0.6 1/m and k = 5 m^2/s^2 give 5.8 1/s). The implicit form has
        // no such limit, so the cap requires it rather than warning about it.
        if (length_alpha > amrex::Real(0.0) && !implicit_diss_on_every_level) {
            amrex::Abort("erf.forest_canopy_tke = true with erf.forest_canopy_length_alpha > 0 "
                         "needs erf.implicit_tke_dissipation = true on every level: the canopy "
                         "length cap makes the explicit TKE dissipation stiff (see the Forest "
                         "model section of the docs)");
        }
    }
};

#endif
