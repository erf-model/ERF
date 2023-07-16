#include <ERF.H>

using namespace amrex;

#if defined(ERF_USE_MOISTURE)
void ERF::advance_microphysics(int lev,
                               MultiFab& cons,
                               const Real& dt_advance)
{
    micro.Init(cons, qmoist[lev],
               grids_to_evolve[lev],
               Geom(lev),
               dt_advance);

    micro.Cloud();
    micro.Diagnose();
    micro.IceFall();
    micro.Precip();
    micro.MicroPrecipFall();

    micro.Update(cons, qmoist[lev]);
}
#endif
