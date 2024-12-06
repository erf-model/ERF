#include <ERF.H>

using namespace amrex;

void ERF::advance_lsm (int lev,
                       MultiFab& /*cons*/,
                       const Real& dt_advance)
{
    if (solverChoice.lsm_type != LandSurfaceType::None) {
        lsm.Advance(lev, dt_advance);
    }
}
