#include <ERF.H>

using namespace amrex;

#if defined(ERF_USE_RRTMGP)
void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const Real& dt_advance)
{
    // TODO: Address issue with lev>0 not spanning all z?
    if (lev == 0) {
        rad[lev]->set_grids(lev, istep[lev], t_new[lev], dt_advance,
                            cons.boxArray(), geom[lev], &(cons),
                            sw_lw_fluxes[lev].get()  , solar_zenith[lev].get(),
                            qheating_rates[lev].get(), z_phys_nd[lev].get()   ,
                            lat_m[lev].get(), lon_m[lev].get());
        rad[lev]->rad_run_impl();
    }
}
#endif
