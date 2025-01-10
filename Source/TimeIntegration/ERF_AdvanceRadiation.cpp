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
                            ls_lw_fluxes[lev], solar_zenith[lev],
                            qheating_rates[lev], z_phys_nd[lev],
                            lat_m[lev], lon_m[lev]);
        rad[lev]->rad_run_impl();
    }
}
#endif
