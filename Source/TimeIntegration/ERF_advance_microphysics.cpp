#include <ERF.H>

using namespace amrex;

void ERF::advance_microphysics(const MultiFab& cons_in,
                               const MultiFab& qc_in,
                                     MultiFab& qv_in,
                               const MultiFab& qi_in,
                               const BoxArray& grids_to_evolve,
                               const Geometry& geom,
                               const Real& dt_advance)
{
#if defined(ERF_USE_MOISTURE)
    micro.Init(S_new,
               qc[lev],
               qv[lev],
               qi[lev],
               grids_to_evolve[lev],
               Geom(lev),
               dt_lev);
    micro.Cloud();
    micro.Diagnose();
    micro.IceFall();
    micro.Precip();
    micro.MicroPrecipFall();
    micro.Update(S_new,
                 qv[lev],
                 qc[lev],
                 qi[lev],
                 qrain[lev],
                 qsnow[lev],
                 qgraup[lev]);
#endif
}
