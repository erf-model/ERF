#include <ERF.H>

using namespace amrex;

namespace {

void update_slm_precipitation (MultiFab& precip,
                               const SurfacePrecipAccumulationSources& sources,
                               const Real dt,
                               const bool store_accumulation)
{
    const bool use_total = surface_precip_has_total_source(sources);
    for (MFIter mfi(precip, TileNoZ()); mfi.isValid(); ++mfi) {
        const Box& box2d = mfi.tilebox();
        auto precip_arr = precip.array(mfi);

        Array4<const Real> total_arr{};
        Array4<const Real> rain_arr{};
        if (sources.total.accum)   { total_arr = sources.total.accum->const_array(mfi); }
        if (sources.rain.accum)    { rain_arr = sources.rain.accum->const_array(mfi); }

        const Real total_factor = sources.total.native_to_kg_m2;
        const Real rain_factor = sources.rain.native_to_kg_m2;

        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real accumulation = Real(0.0);
            if (use_total) {
                accumulation = total_arr(i,j,k) * total_factor;
            } else {
                // SLM receives liquid-water precipitation only.  When a
                // scheme does not provide the legacy total slot, use its
                // explicit liquid-rain accumulator; never add frozen species.
                if (rain_arr) { accumulation = rain_arr(i,j,k) * rain_factor; }
            }
            if (store_accumulation) {
                precip_arr(i,j,k) = accumulation;
            } else {
                // Convert cumulative total precipitation to the interval rate
                // expected by SLM.
                precip_arr(i,j,k) = amrex::max(Real(0.0),
                    (accumulation - precip_arr(i,j,k)) / dt);
            }
        });
    }
}

}

void ERF::advance_lsm (int lev,
                       MultiFab& cons_in,
                       MultiFab& xvel_in,
                       MultiFab& yvel_in,
                       const double& time,
                       const double& dt_advance)
{
    if (solverChoice.lsm_type != LandSurfaceType::None) {
        // Fill boundaries before getting cell-centered velocities
        xvel_in.FillBoundary(geom[lev].periodicity());
        yvel_in.FillBoundary(geom[lev].periodicity());

        lsm.Update_Lsm_Vars_Lev(lev, cons_in, xvel_in, yvel_in);
        const bool use_moist = solverChoice.moisture_type != MoistureType::None && solverChoice.moisture_type != MoistureType::Kessler_NoRain && solverChoice.moisture_type != MoistureType::SatAdj;
        const bool is_slm = solverChoice.lsm_type == LandSurfaceType::SLM;
        SurfacePrecipAccumulationSources slm_precip_sources{};
        if (use_moist && is_slm) {
            slm_precip_sources = micro ? micro->Get_Surface_Precip_Accumulation_Ptrs(lev)
                                       : SurfacePrecipAccumulationSources{};
            update_slm_precipitation(*precip[lev], slm_precip_sources, dt_advance, false);
            lsm.set_LSM_precip_input(lev, precip[lev].get());
        }

        lsm.set_LSM_terrain_inputs(lev, tsk_lev, lmask_lev);
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
            lsm.Advance(lev, dt_advance, t_new[lev], start_time);
        }
        lsm.Update_State_Vars_Lev(lev, cons_in);

        if (use_moist && is_slm) {
            // Retain the typed cumulative value for the next interval.
            update_slm_precipitation(*precip[lev], slm_precip_sources,
                                     dt_advance, true);
        }
    }
}
