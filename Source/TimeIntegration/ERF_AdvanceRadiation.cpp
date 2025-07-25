#include <ERF.H>

using namespace amrex;

void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const Real& dt_advance)
{
    if (solverChoice.rad_type != RadiationType::None) {
        // TODO: Address issue with lev>0 not spanning all z?
        if (lev == 0) {
#ifdef ERF_USE_NETCDF
            MultiFab *lat_ptr = lat_m[lev].get();
            MultiFab *lon_ptr = lon_m[lev].get();
#else
            MultiFab *lat_ptr = nullptr;
            MultiFab *lon_ptr = nullptr;
#endif
            // RRTMGP inputs names and pointers
            Vector<std::string> lsm_input_names = rad[lev]->get_lsm_input_varnames();
            Vector<MultiFab*> lsm_input_ptrs(lsm_input_names.size(),nullptr);
            for (int i(0); i<lsm_input_ptrs.size(); ++i) {
                int varIdx = lsm.Get_VarIdx(lev,lsm_input_names[i]);
                lsm_input_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx);
            }

            // RRTMGP output names and pointers
            Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();
            Vector<MultiFab*> lsm_output_ptrs(lsm_output_names.size(),nullptr);
            for (int i(0); i<lsm_output_ptrs.size(); ++i) {
                int varIdx = lsm.Get_VarIdx(lev,lsm_output_names[i]);
                lsm_output_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx);
            }

            // Enter radiation class driver
            rad[lev]->Run(lev, istep[lev]  , t_new[lev], dt_advance,
                          cons.boxArray()  , geom[lev], &(cons),
                          sw_lw_fluxes[lev].get(), solar_zenith[lev].get(),
                          lsm_input_ptrs   , lsm_output_ptrs,
                          qheating_rates[lev].get(), z_phys_nd[lev].get()   ,
                          lat_ptr, lon_ptr);
            /*
            // TODO: fix this - can't use set_lsm_inputs with IRadiation::Run()
            if (solverChoice.lsm_type == LandSurfaceType::SLM) {
                rad[lev]->set_lsm_inputs(lsm.get_model_lev<SLM>(lev));
            }
            */
        }
    }
}
