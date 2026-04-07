#include <ERF.H>

using namespace amrex;

void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const Real& dt_advance)
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

            Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();
            lsm_output_ptrs = Vector<MultiFab*>(lsm_output_names.size(),nullptr);
            for (int i(0); i<lsm_output_ptrs.size(); ++i) {
                int varIdx = lsm.Get_DataIdx(lev,lsm_output_names[i]);
                if (varIdx >= 0) { lsm_output_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx); }
            }

        } else {
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
        }

        // Enter radiation class driver
        amrex::Real time_for_rad = t_old[lev] + start_time;
        amrex::Print() << " AdvanceRadiation - start_time = " << start_time << " t_old = " << t_old[lev] << " time_for_rad = " << time_for_rad << std::endl;
        rad[lev]->Run(lev, istep[lev], time_for_rad, dt_advance,
                      cons.boxArray(), geom[lev], &(cons),
                      lmask_lev[lev][0].get(), t_surf,
                      lsm_input_ptrs, lsm_output_ptrs,
                      qheating_rates[lev].get(), rad_fluxes[lev].get(),
                      z_phys_nd[lev].get()     , lat_ptr, lon_ptr);
    }
}
