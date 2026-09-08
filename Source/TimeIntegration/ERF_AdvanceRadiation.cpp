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
        Vector<const MultiFab*> lsm_input_ptrs;
        Vector<MultiFab*> lsm_output_ptrs;

        if (m_SurfaceLayer && m_SurfaceLayer->use_sfc_fluxes()) {
            lsm_input_ptrs = {t_surf};
        }

        if (solverChoice.lsm_type == LandSurfaceType::SLM) {
            //rad[lev]->set_lsm_inputs(lsm.get_model_lev<SLM>(lev));

            // TODO: rename SLM variables to match those expected by RRMTGP
            //       for now, pointers are directly returned without name lookup
            lsm_input_ptrs = lsm.get_model_lev<SLM>(lev)->export_to_RRTMGP();

            auto slm_output_name_map = lsm.get_model_lev<SLM>(lev)->get_rad_output_map();
            Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();
            lsm_output_ptrs = Vector<MultiFab*>(lsm_output_names.size(), nullptr);
            for (int i(0); i<lsm_output_ptrs.size(); ++i) {
                std::string slm_name = slm_output_name_map[lsm_output_names[i]];
                if (slm_name == "") continue;

                int varIdx = lsm.Get_DataIdx(lev, slm_name);
                if (varIdx >= 0) { lsm_output_ptrs[i] = lsm.Get_Data_Ptr(lev, varIdx); }
            }
        } else {
            Vector<std::string> lsm_input_names = rad[lev]->get_lsm_input_varnames();
            lsm_input_ptrs.resize(lsm_input_names.size(), nullptr);
            for (int i(0); i<lsm_input_ptrs.size(); ++i) {
                int varIdx = lsm.Get_DataIdx(lev, lsm_input_names[i]);
                if (varIdx >= 0) { lsm_input_ptrs[i] = lsm.Get_Data_Ptr(lev, varIdx); }
            }

            // RRTMGP output names and pointers
            Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();
            lsm_output_ptrs = Vector<MultiFab*>(lsm_output_names.size(), nullptr);
            for (int i(0); i<lsm_output_ptrs.size(); ++i) {
                int varIdx = lsm.Get_DataIdx(lev, lsm_output_names[i]);
                if (varIdx >= 0) { lsm_output_ptrs[i] = lsm.Get_Data_Ptr(lev, varIdx); }
            }
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
}
