#include <ERF.H>

using namespace amrex;

void ERF::advance_lsm (int lev,
                       MultiFab& cons_in,
                       MultiFab& xvel_in,
                       MultiFab& yvel_in,
                       const Real& dt_advance)
{
    if (solverChoice.lsm_type != LandSurfaceType::None) {
        if (solverChoice.lsm_type == LandSurfaceType::NOAHMP) {
            lsm.Advance(lev, cons_in, xvel_in, yvel_in, SFS_hfx3_lev[lev].get(), SFS_q1fx3_lev[lev].get(), dt_advance, istep[0]);
        } else {
            lsm.Advance(lev, dt_advance);
        }
    }
}
