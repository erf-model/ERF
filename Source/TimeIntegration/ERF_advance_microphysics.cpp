#include <ERF.H>

using namespace amrex;

void ERF::advance_microphysics (int lev,
                                MultiFab& cons,
                                const Real& dt_advance)
{
    if (solverChoice.moisture_type != MoistureType::None) {
        micro.Init(cons, *(qmoist[lev]),
                   grids[lev],
                   Geom(lev),
                   dt_advance);
        micro.Advance();
        micro.Update(cons, *(qmoist[lev]));
    }
}
