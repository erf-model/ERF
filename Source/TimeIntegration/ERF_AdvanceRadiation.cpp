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
            rad[lev]->Run(lev, istep[lev], t_new[lev], dt_advance,
                          cons.boxArray(), geom[lev], &(cons),
                          sw_lw_fluxes[lev].get()  , solar_zenith[lev].get(),
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
