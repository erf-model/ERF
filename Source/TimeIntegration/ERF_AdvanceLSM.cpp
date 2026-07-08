#include <ERF.H>

using namespace amrex;

void ERF::advance_lsm (int lev,
                       MultiFab& cons_in,
                       MultiFab& xvel_in,
                       MultiFab& yvel_in,
                       const double& time,
                       const double& dt_advance)
{
    if (solverChoice.lsm_type != LandSurfaceType::None) {
        if (solverChoice.lsm_type == LandSurfaceType::NOAHMP) {
            // Cumulative surface precip accumulations (mm) from microphysics, used to
            // force Noah-MP (RAINBL). qmoist ordering is {rain_accum[, snow_accum,
            // graup_accum]} (MicVarMap); schemes without snow/graup (e.g. Kessler)
            // expose only rain_accum -> pass nullptr for the missing ones.
            const int nqm = (lev < int(qmoist.size())) ? int(qmoist[lev].size()) : 0;
            const MultiFab* rain_accum  = (nqm > 0) ? qmoist[lev][0] : nullptr;
            const MultiFab* snow_accum  = (nqm > 1) ? qmoist[lev][1] : nullptr;
            const MultiFab* graup_accum = (nqm > 2) ? qmoist[lev][2] : nullptr;
            lsm.Advance(lev, cons_in, xvel_in, yvel_in,
                        SFS_hfx3_lev[lev].get(), SFS_q1fx3_lev[lev].get(),
                        rain_accum, snow_accum, graup_accum,
                        time, dt_advance, istep[0], lsm.Get_LSM_Update_Status(0));
        } else {
            lsm.Advance(lev, dt_advance);
        }
    }
}
