/**
 * \file ERF_init_windfarm.cpp
 */
#include <ERF.H>

using namespace amrex;

/**
 * Read in the turbine locations in latitude-longitude from windturbines.txt
 * and convert it into x and y coordinates in metres
 *
 * @param lev Integer specifying the current level
 */

// Explicit instantiation

void
ERF::init_windfarm (int lev)
{
    if(solverChoice.windfarm_loc_type == WindFarmLocType::lat_lon) {
        windfarm->read_tables(solverChoice.windfarm_loc_table,
                             solverChoice.windfarm_spec_table,
                             false, true,
                             solverChoice.latitude_lo, solverChoice.longitude_lo);
    } else if(solverChoice.windfarm_loc_type == WindFarmLocType::x_y) {
        windfarm->read_tables(solverChoice.windfarm_loc_table,
                             solverChoice.windfarm_spec_table,
                             true, false);
    }

    windfarm->fill_Nturb_multifab(geom[lev], Nturb[lev]);

    windfarm->write_turbine_locations_vtk();

    if(solverChoice.windfarm_type == WindFarmType::SimpleAD) {
        windfarm->write_actuator_disks_vtk();
    }
}

void
ERF::advance_windfarm (const Geometry& a_geom,
                       const Real& dt_advance,
                       MultiFab& cons_in,
                       MultiFab& U_old,
                       MultiFab& V_old,
                       MultiFab& W_old,
                       MultiFab& mf_vars_windfarm,
                       const MultiFab& mf_Nturb)
{
        windfarm->advance(a_geom, dt_advance, cons_in, mf_vars_windfarm,
                          U_old, V_old, W_old, mf_Nturb);
}
