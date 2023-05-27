#include <ERF.H>

using namespace amrex;

void ERF::advance_microphysics(int lev,
                               MultiFab& cons,
                               const Real& dt_advance)
{
#if defined(ERF_USE_MOISTURE)
    micro.Init(cons,
               qc[lev],
               qv[lev],
               qi[lev],
               grids_to_evolve[lev],
               Geom(lev),
               dt_advance);

    micro.Cloud();
    micro.Diagnose();
    micro.IceFall();
    micro.Precip();
    micro.MicroPrecipFall();

    micro.Update(cons,
                 qv[lev],
                 qc[lev],
                 qi[lev],
                 qrain[lev],
                 qsnow[lev],
                 qgraup[lev]);
#endif
}
