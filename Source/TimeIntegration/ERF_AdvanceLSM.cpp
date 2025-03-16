#include <ERF.H>

using namespace amrex;

void ERF::advance_lsm (int lev,
                       MultiFab& cons_in,
                       MultiFab& xvel_in,
                       MultiFab& yvel_in,
                       const Real& dt_advance)
{
    if (solverChoice.lsm_type != LandSurfaceType::None) {
        if (solverChoice.lsm_type == LandSurfaceType::NOAH) {
            std::unique_ptr<MultiFab>& hfx3_in = SFS_hfx3_lev[lev];
            std::unique_ptr<MultiFab>& qfx3_in = SFS_q1fx3_lev[lev];
            lsm.Advance(lev, cons_in, xvel_in, yvel_in, hfx3_in, qfx3_in, dt_advance);
        } else {
            lsm.Advance(lev, dt_advance);
        }
    }
}
