#include <ERF.H>

using namespace amrex;

/**
 * Advances the urban-model interface on one AMR level.
 */
void
ERF::advance_urban (int lev,
                    MultiFab& cons_in,
                    MultiFab& u_in,
                    MultiFab& v_in,
                    const Real& dt_advance,
                    const Geometry& geom_in,
                    const MultiFab* z_phys_nd_in,
                    MultiFab& eddyDiffs_in)
{
    if (solverChoice.urban_type != UrbanType::None &&
        solverChoice.urban_enabled_lev[lev] == 1) {
        u_in.FillBoundary(geom_in.periodicity());
        v_in.FillBoundary(geom_in.periodicity());

        urban.Update_Urban_Vars_Lev(lev, cons_in, u_in, v_in);

#ifdef ERF_USE_NETCDF
        MultiFab* lat_ptr = lat_m[lev].get();
        MultiFab* lon_ptr = lon_m[lev].get();
#else
        MultiFab* lat_ptr = nullptr;
        MultiFab* lon_ptr = nullptr;
#endif
        urban.set_urban_terrain_inputs(lev, tsk_lev, lmask_lev, land_type_lev,
                                       urb_frac_lev, lat_ptr, lon_ptr);

        if (solverChoice.rad_type != RadiationType::None) {
            urban.set_urban_flux_inputs(lev, solar_declin);
        }

        urban.Advance(lev, dt_advance, t_new[lev], start_time, geom_in,
                      z_phys_nd_in, eddyDiffs_in, calday);
        urban.Update_State_Vars_Lev(lev, cons_in);
    }
}
