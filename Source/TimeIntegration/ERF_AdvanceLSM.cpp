#include <ERF.H>

using namespace amrex;

void ERF::advance_lsm (int lev,
                       MultiFab& cons_in,
                       MultiFab& xvel_in,
                       MultiFab& yvel_in,
                       const Real& dt_advance)
{
    if (solverChoice.lsm_type != LandSurfaceType::None) {
        // Fill boundaries before getting cell-centered velocities
        xvel_in.FillBoundary(geom[lev].periodicity());
        yvel_in.FillBoundary(geom[lev].periodicity());

        lsm.Update_Lsm_Vars_Lev(lev, cons_in, xvel_in, yvel_in);
        if (solverChoice.rad_type != RadiationType::None) {
            lsm.set_LSM_flux_inputs(lev, sw_lw_fluxes[lev].get(), solar_zenith[lev].get());
        }
        const bool use_moist = solverChoice.moisture_type != MoistureType::None && solverChoice.moisture_type != MoistureType::Kessler_NoRain && solverChoice.moisture_type != MoistureType::SatAdj;
        int rain_comp = 0;
        if (use_moist) {
            if (solverChoice.moisture_type == MoistureType::Morrison || solverChoice.moisture_type == MoistureType::Morrison_NoIce)
            {
                rain_comp = 5;
            }

            // qmoist[0] is total rain accumulation in mm over entire simulation
            // Convert to a rate of mm/s for SLM
            for ( MFIter mfi(*precip[lev], TileNoZ()); mfi.isValid(); ++mfi) {
                const auto& box2d = mfi.tilebox();
                auto precip_array = precip[lev]->array(mfi);
                auto rain_accum_array   = qmoist[lev][rain_comp]->const_array(mfi);
                ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    precip_array(i, j, k) = std::max(0.0, (rain_accum_array(i, j, k) - precip_array(i, j, k)) / dt_advance);
                });
            }

            lsm.set_LSM_precip_input(lev, precip[lev].get());
        }

        lsm.set_LSM_terrain_inputs(lev, sst_lev, lmask_lev);
        if (solverChoice.lsm_type == LandSurfaceType::NOAHMP) {
            lsm.Advance(lev, cons_in, xvel_in, yvel_in, SFS_hfx3_lev[lev].get(), SFS_q1fx3_lev[lev].get(), dt_advance, istep[lev]);
        } else {
            lsm.Advance(lev, dt_advance);
        }
        lsm.Update_State_Vars_Lev(lev, cons_in);

        if (use_moist) {
            // copy current rain_accum into precip used to compute the rate for the next timestep
            for ( MFIter mfi(*precip[lev], TileNoZ()); mfi.isValid(); ++mfi) {
                const auto& box2d = mfi.tilebox();
                auto precip_array = precip[lev]->array(mfi);
                auto rain_accum_array   = qmoist[lev][rain_comp]->const_array(mfi);
                ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    precip_array(i, j, k) = rain_accum_array(i, j, k);
                });
            }
        }
    }
}
