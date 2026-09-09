#include <ERF.H>

using namespace amrex;

/**
 * @brief Advance radiation diagnostics and heating rates for one time step.
 *
 * **Temporal semantics**
 *
 * This function is called exactly once per ERF::Advance() invocation, after
 * the SurfaceLayer and LSM updates and before the dycore slow and fast
 * substeps. It operates on the old state (t^n) at the beginning of the slow
 * step.
 *
 * - RRTMGP path (erf.radiation_model, SolverChoice::rad_type != None): a full
 *   spectral model with its own time-centering and source-term semantics.
 *   Produces qheating_rates[lev].
 *
 * - Two-stream path (erf.radiation_type, RadChoice::rad_type == TwoStream): a
 *   shortwave and longwave model that computes heating rates from the
 *   old-state atmosphere (t^n) with clear-sky and cloudy column algorithms.
 *   The heating rates go into qheating_rates[lev], a 2-component MultiFab
 *   holding shortwave and longwave.
 *
 * **Source-term application**
 *
 * The computed qheating_rates are injected into the RhoTheta source term in
 * ERF_MakeSources.cpp only while the slow RHS is being built (is_slow_step is
 * true), which ensures:
 * 1. Radiation tendencies are applied once per slow step, not per substep.
 * 2. The tendencies represent the old-state atmosphere throughout all fast
 *    substeps of the current slow step.
 * 3. There is no temporal aliasing from repeated calls to advance_radiation()
 *    within a slow step, since there is only one call per slow step.
 *
 * **Key contracts**
 *
 * - Radiation heating is an old-state forcing. The qheating_rates computed
 *   here are the radiative heating of the old-state atmosphere (t^n), applied
 *   as a source term while the slow RHS is built. That gives one radiative
 *   increment per slow step, consistent with the old state across every fast
 *   substep. Radiation does not adapt to the state within a slow step.
 *
 * - The two radiation paths are mutually exclusive. RRTMGP and two-stream
 *   never both run in one simulation; the if/else below selects one. Both
 *   produce qheating_rates in the same 2-component (SW, LW) format, and the
 *   source-term gate in ERF_MakeSources.cpp tests both, so exactly one
 *   matches in any given simulation.
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
        bool lsm_updated = (lev==0 && max_level>0) ? lsm.Get_LSM_Update_Status(lev) : false;

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
    // Two-stream radiation driver. This is a separate, mutually exclusive
    // path from the RRTMGP branch above: RRTMGP is selected by
    // erf.radiation_model (SolverChoice::rad_type), two-stream by
    // erf.radiation_type (RadChoice::rad_type).
    //
    // - The call happens exactly once per slow step (from ERF::Advance).
    // - The heating rates computed here are old-state based (t^n).
    // - They are injected into the RhoTheta source only on is_slow_step
    //   (see ERF_MakeSources.cpp), so there is no duplicate forcing.
    // - istep[lev] serves as both the current step number and the CSV
    //   diagnostics row index; dt_advance supplies the time_step value
    //   logged to the CSV and console output.
    else if (solverChoice.radChoice.rad_type == RadType::TwoStream) {
        compute_twostream_radiation_diagnostics(lev, istep[lev], t_old[lev], "pre_dycore");
    }
}
