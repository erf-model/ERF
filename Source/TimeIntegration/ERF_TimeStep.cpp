#include <ERF.H>
#include <ERF_Utils.H>
#include <ERF_ReadFromWRFBdy.H>

using namespace amrex;

/**
 * Function that coordinates the evolution across levels -- this calls Advance to do the
 * actual advance at this level,  then recursively calls itself at finer levels
 *
 * @param[in] lev level of refinement (coarsest level is 0)
 * @param[in] time start time for time advance
 * @param[in] iteration time step counter
 */

void
ERF::timeStep (int lev, Real time, int /*iteration*/)
{
    //
    // We need to FillPatch the coarse level before assessing whether to regrid
    // We have not done the swap yet so we fill the "new" which will become the "old"
    //
    MultiFab& S_new = vars_new[lev][Vars::cons];
    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

#ifdef ERF_USE_NETCDF
    //
    // Since we now only read in a subset of the time slices in wrfbdy and
    //     wrflowinp, we need to check whether it's time to read in more.
    //
    bool use_moist = (solverChoice.moisture_type != MoistureType::None);
    if (solverChoice.use_real_bcs && (lev==0)) {
        Real dT = bdy_time_interval;

        int n_time_old = static_cast<int>( (time        ) /  dT);
        int n_time_new = static_cast<int>( (time+dt[lev]) /  dT);

        int ntimes = bdy_data_xlo.size();
        for (int itime = 0; itime < ntimes; itime++)
        {
            //if (bdy_data_xlo[itime].size() > 0) {
            //    amrex::Print() << "HAVE  DATA AT TIME " << itime << std::endl;
            //} else {
            //    amrex::Print() << " NO   DATA AT TIME " << itime << std::endl;
            //}

            bool clear_itime = (itime < n_time_old);

            if (clear_itime && bdy_data_xlo[itime].size() > 0) {
                bdy_data_xlo[itime].clear();
                bdy_data_xhi[itime].clear();
                bdy_data_ylo[itime].clear();
                bdy_data_yhi[itime].clear();
                //amrex::Print() << "CLEAR BDY DATA AT TIME " << itime << std::endl;
            }

            bool need_itime = (itime >= n_time_old && itime <= n_time_new+1);
            //if (need_itime) amrex::Print()  << "NEED  BDY DATA AT TIME " << itime << std::endl;

            if (bdy_data_xlo[itime].size() == 0 && need_itime) {
                read_from_wrfbdy(itime,nc_bdy_file,geom[0].Domain(),
                                 bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                 real_width);

                convert_all_wrfbdy_data(itime, geom[0].Domain(), bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                    *mf_MUB, *mf_C1H, *mf_C2H,
                                    vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel], vars_new[lev][Vars::cons],
                                    geom[lev], use_moist);
           }
        } // itime
    } // use_real_bcs && lev == 0

    if (!nc_low_file.empty() && (lev==0)) {
        Real dT = low_time_interval;

        int n_time_old = static_cast<int>( (time        ) /  dT);
        int n_time_new = static_cast<int>( (time+dt[lev]) /  dT);

        int ntimes = bdy_data_xlo.size();
        for (int itime = 0; itime < ntimes; itime++)
        {
            bool clear_itime = (itime < n_time_old);

            if (clear_itime && low_data_zlo[itime].size() > 0) {
                low_data_zlo[itime].clear();
                //amrex::Print() << "CLEAR LOW DATA AT TIME " << itime << std::endl;
            }

            bool need_itime = (itime >= n_time_old && itime <= n_time_new+1);
            //if (need_itime) amrex::Print()  << "NEED  LOW DATA AT TIME " << itime << std::endl;

            if (low_data_zlo[itime].size() == 0 && need_itime) {
                read_from_wrflow(itime, nc_low_file, geom[lev].Domain(), low_data_zlo);

                update_sst_tsk(itime, geom[lev], ba2d[lev],
                               sst_lev[lev], tsk_lev[lev],
                               m_SurfaceLayer, low_data_zlo,
                               S_new, *mf_PSFC[lev],
                               solverChoice.rdOcp, lmask_lev[lev][0], use_moist);
            }
        } // itime
    } // have nc_low_file && lev == 0
#endif

    //
    // NOTE: the momenta here are not fillpatched (they are only used as scratch space)
    //
    if (lev == 0) {
        FillPatchCrseLevel(lev, time, {&S_new, &U_new, &V_new, &W_new});
    } else if (lev < finest_level) {
        FillPatchFineLevel(lev, time, {&S_new, &U_new, &V_new, &W_new},
                           {&S_new, &rU_new[lev], &rV_new[lev], &rW_new[lev]},
                           base_state[lev], base_state[lev]);
    }

    if (regrid_int > 0)  // We may need to regrid
    {
        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level)
        {
            if ( (istep[lev] % regrid_int == 0) && (istep[lev] > last_regrid_step[lev]) )
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level;

                regrid(lev, time);

#ifdef ERF_USE_PARTICLES
                if (finest_level != old_finest) {
                    particleData.Redistribute();
                }
#endif

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = dt[k-1] / static_cast<Real>(nsubsteps[k]);
                }
            } // if
        } // lev
    }

    // Update what we call "old" and "new" time
    t_old[lev] = t_new[lev];
    t_new[lev] += dt[lev];

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << std::setprecision(timeprecision)
                       << "ADVANCE from elapsed time = " << t_old[lev] << " to " << t_new[lev]
                       << " with dt = " << dt[lev] << std::endl;
    }

#ifdef ERF_USE_WW3_COUPLING
    amrex::Print() <<  " About to call send_to_ww3 from ERF_Timestep" << std::endl;
    send_to_ww3(lev);
    amrex::Print() <<  " About to call read_waves from ERF_Timestep"  << std::endl;
    read_waves(lev);
    //send_to_ww3(lev);
    //read_waves(lev);
    //send_to_ww3(lev);
#endif

    // Advance a single level for a single time step
    Advance(lev, time, dt[lev], istep[lev], nsubsteps[lev]);

    ++istep[lev];

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            Real strt_time_for_fine = time + (i-1)*dt[lev+1];
            timeStep(lev+1, strt_time_for_fine, i);
        }
    }

    if (verbose && lev == 0 && solverChoice.moisture_type != MoistureType::None) {
        amrex::Print() << "Cloud fraction " << time << "  " << cloud_fraction(time) << std::endl;
    }
}
