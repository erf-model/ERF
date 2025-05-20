#include <ERF.H>

using namespace amrex;

void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const Real& dt_advance)
{
    if (solverChoice.rad_type != RadiationType::None) {
        // TODO: Address issue with lev>0 not spanning all z?
        if (lev == 0) {
            rad[lev]->Run(lev, istep[lev], t_new[lev], dt_advance,
                          cons.boxArray(), geom[lev], &(cons),
                          sw_lw_fluxes[lev].get()  , solar_zenith[lev].get(),
                          qheating_rates[lev].get(), z_phys_nd[lev].get()   ,
                          lat_m[lev].get(), lon_m[lev].get());
            /*
            // TODO: fix this - can't use set_lsm_inputs with IRadiation::Run()
            if (solverChoice.lsm_type == LandSurfaceType::SLM) {
                rad[lev]->set_lsm_inputs(lsm.get_model_lev<SLM>(lev));
            }
            */
        }
    }
}
