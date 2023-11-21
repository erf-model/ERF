#include <ERF.H>

using namespace amrex;

#if defined(ERF_USE_MOISTURE)
void ERF::advance_microphysics (int lev,
                                MultiFab& cons,
                                const Real& dt_advance)
{
    micro.Init(cons, qmoist[lev],
               grids[lev],
               Geom(lev),
               dt_advance);
    micro.Advance();
	micro.Diagnose();
    micro.Update(cons, qmoist[lev]);
}
#endif
