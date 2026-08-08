#include <ERF.H>

using namespace amrex;

/**
 * @brief Advance radiation diagnostics and heating rates for one time step.
 *
 * **Temporal Semantics (Phase 6):**
 * 
 * This function is called exactly once per ERF::Advance() invocation, after
 * SurfaceLayer/LSM updates but BEFORE the dycore (slow+fast substeps) are
 * executed. It operates on the "old" state (t^n) at the beginning of the
 * slow step.
 *
 * - RRTMGP path (solverChoice.rad_type != None): Full spectral model with its
 *   own time-centering and source-term semantics; produces qheating_rates[lev].
 *
 * - TwoStream path (solverChoice.radChoice.rad_type == TwoStream): Two-stream
 *   SW/LW model that computes heating rates from the old-state atmosphere
 *   (t^n) using clear-sky/cloudy column algorithms. The heating rates are
 *   written to qheating_rates[lev] (2-component MultiFab: [SW, LW]).
 *
 * **Source-Term Application (Phase 6):**
 *
 * The computed qheating_rates are later injected into the RhoTheta source
 * term in ERF_MakeSources.cpp, ONLY during slow-RHS construction (is_slow_step
 * = true). This ensures:
 * 1. Radiation tendencies are applied once per slow step, not per substep.
 * 2. The tendencies represent the old-state atmosphere throughout all fast
 *    substeps of the current slow step.
 * 3. No temporal aliasing from multiple calls to advance_radiation() within a
 *    single slow step (there is only one call per slow step).
 *
 * **Key Contracts (Phase 6):**
 *
 * - R13: Radiation Heating as Old-State Forcing
 *   The qheating_rates computed in this function represent the radiative
 *   heating based on the old-state atmosphere (t^n). These rates are applied
 *   as a source term during the slow RHS construction, providing a single
 *   radiative "kick" per slow step that is consistent with the old state
 *   throughout all fast substeps. No state-dependent (adaptive) radiation
 *   updates occur within a slow-step's fast substeps; radiation is fixed
 *   across the slow step.
 *
 * - R14: Mutually Exclusive Radiation Paths
 *   RRTMGP and TwoStream paths NEVER coexist in the same simulation; they are
 *   mutually exclusive via the if/else_if structure below. Both produce
 *   qheating_rates in the same 2-component format (SW, LW). The source-term
 *   injection gate in ERF_MakeSources.cpp checks both paths; only one gate
 *   will ever match in any given simulation.
 *
 * @param[in] lev Level of refinement (coarsest level is 0)
 * @param[in,out] cons Conservative quantities (Rho, RhoTheta, RhoQ*, RhoRE)
 * @param[in] dt_advance Time step for this slow-step stage [seconds]
 */
void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const double& dt_advance)
{
    if (solverChoice.rad_type != RadiationType::None) {
#ifdef ERF_USE_NETCDF
        MultiFab *lat_ptr = lat_m[lev].get();
        MultiFab *lon_ptr = lon_m[lev].get();
#else
        MultiFab *lat_ptr = nullptr;
        MultiFab *lon_ptr = nullptr;
#endif
        // T surf from SurfaceLayer if we have it
        MultiFab* t_surf = (m_SurfaceLayer) ? m_SurfaceLayer->get_t_surf(lev) : nullptr;

        // RRTMGP inputs names and pointers
        Vector<std::string> lsm_input_names = rad[lev]->get_lsm_input_varnames();
        Vector<MultiFab*> lsm_input_ptrs(lsm_input_names.size(),nullptr);
        for (int i(0); i<lsm_input_ptrs.size(); ++i) {
            int varIdx = lsm.Get_DataIdx(lev,lsm_input_names[i]);
            if (varIdx >= 0) { lsm_input_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx); }
        }

        // RRTMGP output names and pointers
        Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();
        Vector<MultiFab*> lsm_output_ptrs(lsm_output_names.size(),nullptr);
        for (int i(0); i<lsm_output_ptrs.size(); ++i) {
            int varIdx = lsm.Get_DataIdx(lev,lsm_output_names[i]);
            if (varIdx >= 0) { lsm_output_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx); }
        }

        // Force radiation update to sync with lsm?
        bool lsm_updated = (lev==0) ? lsm.Get_LSM_Update_Status(lev) : false;

        // Enter radiation class driver
        double time_for_rad = t_old[lev] + start_time;
        rad[lev]->Run(lev, istep[lev], time_for_rad, dt_advance,
                      cons.boxArray(), geom[lev], &(cons),
                      lmask_lev[lev][0].get(), t_surf,
                      lsm_input_ptrs, lsm_output_ptrs,
                      qheating_rates[lev].get(), rad_fluxes[lev].get(),
                      z_phys_nd[lev].get()     , lat_ptr, lon_ptr,
                      lsm_updated);
    }
    // Phase 6 (Audit + Timing Docs): Two-Stream TwoStream Radiation Driver.
    //
    // Phase 5 (Step 3) wired the Phase 1-5 two-stream radiation driver into
    // the time loop for the first time. Prior to Phase 5,
    // compute_twostream_radiation_diagnostics() was never called anywhere
    // in the codebase -- Phase 1-4 code was fully correct but dead code
    // from the simulation's perspective. This is a separate, independent
    // path from the RRTMGP branch above (mutually exclusive: RRTMGP uses
    // solverChoice.rad_type; TwoStream uses solverChoice.radChoice.rad_type).
    //
    // Phase 6 audit confirms:
    // - This call happens exactly once per slow step (in ERF::Advance)
    // - The heating rates computed here are old-state based (t^n)
    // - The heating rates are injected into RhoTheta source only on is_slow_step
    //   (see ERF_MakeSources.cpp), ensuring no unintended duplicate forcing
    // - istep[lev] is used as both the "current step number" and (Phase 1-4
    //   convention) as the CSV diagnostics row index; dt_advance provides the
    //   time_step diagnostic value logged to the CSV/console output.
    else if (solverChoice.radChoice.rad_type == RadType::TwoStream) {
        compute_twostream_radiation_diagnostics(lev, istep[lev], t_old[lev], "pre_dycore");
    }
}
