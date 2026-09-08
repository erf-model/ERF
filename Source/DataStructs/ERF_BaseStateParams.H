/**
 * \file ERF_BaseStateParams.H
 */
#ifndef ERF_BASE_STATE_PARAMS_H_
#define ERF_BASE_STATE_PARAMS_H_

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Print.H>

#include <ERF_Constants.H>

/**
 * Parameters defining the WRF reference (base) state, along with the layer
 * interface heights derived from them.
 *
 * The base state is piecewise in log-pressure:
 *   (1) troposphere      T = T00  + TLP       * ln(p/P00)         p >  P_iso
 *   (2) isothermal layer T = TISO                       P_STRAT < p <= P_iso
 *   (3) stratosphere     T = TISO + TLP_STRAT * ln(p/P_STRAT)     p <= P_STRAT
 * Each piece inverts to p(z) in closed form under dp/dz = -rho g, so the two
 * interface heights are precomputed once and the inversion branches on z.
 *
 * Because the resulting profile is a pure function of height, it can be
 * evaluated on any subset of a column -- it does not require a box that spans
 * the domain in the vertical.
 *
 * NOTE: The defaults are the WRF namelist defaults for base_temp, base_pres,
 *       base_lapse, iso_temp, base_lapse_strat and base_pres_strat.  metgrid
 *       files do not carry these parameters at all; wrfinput files usually do,
 *       and are read by read_base_state_params_from_wrfinput.
 */
struct BaseStateParams
{
    amrex::Real T00       = amrex::Real(290.0);
    amrex::Real P00       = p_0;
    amrex::Real TLP       = amrex::Real(50.0);
    amrex::Real TISO      = amrex::Real(200.0);
    amrex::Real TLP_STRAT = amrex::Real(-11.0);
    amrex::Real P_STRAT   = zero;

    // Derived by set_layer_interfaces()
    amrex::Real P_iso     = zero;
    amrex::Real z_iso     = zero;
    amrex::Real z_strat   = zero;
    bool        use_strat = false;

    // Whether set_layer_interfaces() has been called, so that consumers can
    // catch an un-initialized profile rather than silently using z_iso = 0.
    bool        is_set    = false;

    //
    // Compute the layer interfaces and report them.
    //
    // NOTE: This must be called from outside of any MFIter (and outside of any
    //       OpenMP parallel region) so the diagnostics are emitted exactly once
    //       rather than once per box per thread.
    //
    void set_layer_interfaces ()
    {
        // The troposphere inverts p(z) as sqrt((T00/TLP)^2 - 2*g*z/(TLP*R_d)), which
        // folds 1/TLP inside the square root and so discards the sign of TLP. A
        // non-positive lapse parameter, or an isothermal temperature at or above the
        // surface temperature, would silently select the wrong root instead of failing.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(TLP > zero,
                                         "Base state TLP must be positive");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(TISO < T00,
                                         "Base state TISO must be less than T00");

        const amrex::Real x_iso = (TISO - T00) / TLP;
        P_iso = P00 * std::exp(x_iso);
        z_iso = -(R_d/CONST_GRAV) * (T00*x_iso + myhalf*TLP*x_iso*x_iso);

        // The upper stratospheric layer is optional (P_STRAT == 0 or TLP_STRAT == 0
        // disables it) and is only meaningful if it begins above the isothermal layer,
        // i.e. if P_STRAT is below the pressure at which the isothermal layer starts.
        const bool want_strat = ((P_STRAT > zero) && (TLP_STRAT != zero));
        use_strat = (want_strat && (P_STRAT < P_iso));
        z_strat   = (use_strat) ? z_iso + (R_d*TISO/CONST_GRAV)*std::log(P_iso/P_STRAT)
                                : z_iso;

        // A configured stratospheric layer that lies at or below the isothermal
        // transition cannot be represented, so say so rather than dropping it quietly.
        if (want_strat && !use_strat) {
            amrex::Print() << "WARNING: the base stratospheric layer is being ignored: P_STRAT = "
                           << P_STRAT << " Pa is not below the pressure at the base of the "
                           << "isothermal layer, P_iso = " << P_iso << " Pa.\n";
            amrex::Print() << "         TLP_STRAT = " << TLP_STRAT << " will have no effect and "
                           << "the atmosphere above z_iso will be isothermal at TISO = "
                           << TISO << " K.\n";
        }

        amrex::Print() << "Base state layer interfaces: z_iso = " << z_iso << " m";
        if (use_strat) amrex::Print() << ", z_strat = " << z_strat << " m";
        amrex::Print() << "\n";

        is_set = true;
    }

    //
    // True if this holds the same six defining parameters as bsp.  Used to warn when
    // a nested wrfinput file disagrees with its parent: the levels of an AMR hierarchy
    // must share one reference profile, or their base states cannot be consistent
    // across the coarse/fine interface.
    //
    bool same_params_as (const BaseStateParams& bsp) const
    {
        return (T00       == bsp.T00 ) && (P00     == bsp.P00 ) &&
               (TLP       == bsp.TLP ) && (TISO    == bsp.TISO) &&
               (TLP_STRAT == bsp.TLP_STRAT) && (P_STRAT == bsp.P_STRAT);
    }
};

#endif
