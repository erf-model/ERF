#include <ERF.H>

using namespace amrex;

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
    // Phase 5 (Step 3): Wire the Phase 1-5 two-stream radiation driver into
    // the time loop for the first time. Prior to this,
    // compute_twostream_radiation_diagnostics() was never called anywhere
    // in the codebase -- Phase 1-4 code was fully correct but dead code
    // from the simulation's perspective. This is a separate, independent
    // path from the RRTMGP branch above (mutually exclusive: RRTMGP uses
    // solverChoice.rad_type; TwoStream uses solverChoice.radChoice.rad_type).
    //
    // istep[lev] is used as both the "current step number" and (Phase 1-4
    // convention) as the CSV diagnostics row index; dt_advance provides the
    // time_step diagnostic value logged to the CSV/console output.
    else if (solverChoice.radChoice.rad_type == RadType::TwoStream) {
        compute_twostream_radiation_diagnostics(lev, istep[lev], dt_advance);
    }
}
