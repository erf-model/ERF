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
            // Typed surface-precip accumulation sources from the active microphysics
            // scheme: borrowed const views of the scheme-native cumulative accumulators
            // plus their native->kg/m^2 (== water-equivalent mm) conversion factors.
            // Noah-MP forms the water-equivalent interval precip (RAINBL / MP_RAINNC /
            // MP_SNOW / MP_GRAUP / SR) from these, so schemes with differing accumulator
            // semantics or units (e.g. SAM's density-scaled snow/graupel) are handled
            // correctly. Empty when moisture/precip is off -> land model runs precip-free.
            const SurfacePrecipAccumulationSources precip_sources =
                micro ? micro->Get_Surface_Precip_Accumulation_Ptrs(lev)
                      : SurfacePrecipAccumulationSources{};
            lsm.Advance(lev, cons_in, xvel_in, yvel_in,
                        SFS_hfx3_lev[lev].get(), SFS_q1fx3_lev[lev].get(),
                        precip_sources,
                        time, dt_advance, istep[0], lsm.Get_LSM_Update_Status(0));
        } else {
            lsm.Advance(lev, dt_advance);
        }
    }
}
