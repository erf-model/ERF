#include <ERF.H>

using namespace amrex;

/**
 * @brief Advance radiation diagnostics and heating rates for one time step.
 *
 * **Temporal semantics**
 *
 * This function is called exactly once per ERF::Advance() invocation, after
 * the SurfaceLayer and LSM updates and before the dycore slow and fast
 * substeps. It operates on the old state (t^n) at the beginning of the slow
 * step.
 *
 * - RRTMGP / Simple path (erf.radiation_model = RRTMGP or Simple): a full
 *   spectral model with its own time-centering and source-term semantics.
 *   Produces qheating_rates[lev].
 *
 * - Two-stream path (erf.radiation_model = TwoStream): a
 *   shortwave and longwave model that computes heating rates from the
 *   old-state atmosphere (t^n) with clear-sky and cloudy column algorithms.
 *   The heating rates go into qheating_rates[lev], a 2-component MultiFab
 *   holding shortwave and longwave.
 *
 * **Source-term application**
 *
 * The computed qheating_rates are injected into the RhoTheta source term in
 * ERF_MakeSources.cpp only while the slow RHS is being built (is_slow_step is
 * true), which ensures:
 * 1. Radiation tendencies are applied once per slow step, not per substep.
 * 2. The tendencies represent the old-state atmosphere throughout all fast
 *    substeps of the current slow step.
 * 3. There is no temporal aliasing from repeated calls to advance_radiation()
 *    within a slow step, since there is only one call per slow step.
 *
 * **Key contracts**
 *
 * - Radiation heating is an old-state forcing. The qheating_rates computed
 *   here are the radiative heating of the old-state atmosphere (t^n), applied
 *   as a source term while the slow RHS is built. That gives one radiative
 *   increment per slow step, consistent with the old state across every fast
 *   substep. Radiation does not adapt to the state within a slow step.
 *
 * - The two radiation paths are mutually exclusive. RRTMGP and two-stream
 *   never both run in one simulation; the if/else below selects one. Both
 *   produce qheating_rates in the same 2-component (SW, LW) format, and the
 *   source-term gate in ERF_MakeSources.cpp tests both, so exactly one
 *   matches in any given simulation.
 *
 * @param[in] lev Level of refinement (coarsest level is 0)
 * @param[in,out] cons Conservative quantities (Rho, RhoTheta, RhoQ*, RhoRE)
 * @param[in] dt_advance Time step for this slow-step stage [seconds]
 */
void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const double& dt_advance)
{
    BL_PROFILE("ERF::advance_radiation()");

    // Fill this level's radiation fields by interpolation from its parent.
    //
    // Two situations need this and need exactly the same work done:
    //
    //   * a nested patch -- a fine level that does not span the full column, for which
    //     RRTMGP cannot be run at all because it needs complete atmospheric columns; and
    //   * the first step of a level built by interp_atmos_from_coarse, whose atmospheric
    //     state came from FillCoarsePatch and is not yet thermodynamically consistent.
    //
    // NOTE: the coarse MultiFabs are handed to InterpFromCoarseLevel directly rather than
    //       being copied into a ghosted temporary first.  InterpFromCoarseLevel builds its
    //       own coarse patch and fills it with a ParallelCopy that takes the source's VALID
    //       region only (send_ghost = recv_ghost = 0), so a temporary carrying ghost cells
    //       contributes nothing to the result -- it only costs an allocation, a copy and,
    //       where FillBoundary was called on it, a round of communication per field per step.
    auto interp_rad_from_coarse = [&] ()
    {
        // Ensure the parent level's ghost cells are filled before interpolation.  This is
        // needed even when radiation did not run this step, especially with two-way coupling
        // where the grid structure may have changed.
        if (!rad[lev-1]->is_nested_patch()) {
            qheating_rates[lev-1]->FillBoundary(geom[lev-1].periodicity());
        }

        InterpFromCoarseLevel(*qheating_rates[lev], qheating_rates[lev]->nGrowVect(),
                              IntVect(0,0,0),
                              *qheating_rates[lev-1], 0, 0, 2,
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BCVars::cons_bc);

        // Radiation fluxes (needed for plotfiles and diagnostics)
        if (rad_fluxes[lev] && rad_fluxes[lev-1]) {
            const int nc = rad_fluxes[lev]->nComp();

            InterpFromCoarseLevel(*rad_fluxes[lev], rad_fluxes[lev]->nGrowVect(),
                                  IntVect(0,0,0),
                                  *rad_fluxes[lev-1], 0, 0, nc,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &cell_cons_interp,
                                  domain_bcs_type, BCVars::cons_bc);

            // The interpolation above cannot carry the top-of-atmosphere interface.
            // RRTMGP's layout puts the TOA fluxes in the z-ghost cell at k = khi+1 (see
            // ERF_MakeNewArrays.cpp), which is physical data rather than a halo, and it is
            // invisible to InterpFromCoarseLevel from both ends: the coarse TOA sits outside
            // the valid region its ParallelCopy reads, and on a level that spans the full
            // column the fine TOA sits outside the fine domain it writes.  A nested patch
            // does not reach the model top, so its top plane is an ordinary interior
            // interface that the call above has already filled.
            //
            // Move both planes into valid index space -- a single layer at the coarse domain
            // top -- so the ordinary machinery can interpolate them horizontally, then put
            // the result back in the fine level's ghost cell.
            if (!rad[lev]->is_nested_patch()) {
                const int khi_c = geom[lev-1].Domain().bigEnd(2);
                const int khi_f = geom[lev  ].Domain().bigEnd(2);

                auto top_slab = [&] (const BoxArray& ba_in) {
                    BoxList bl;
                    for (int i = 0, n = ba_in.size(); i < n; ++i) {
                        Box b = ba_in[i];
                        b.setSmall(2, khi_c); b.setBig(2, khi_c);
                        bl.push_back(b);
                    }
                    return BoxArray(std::move(bl));
                };

                MultiFab toa_crse(top_slab(rad_fluxes[lev-1]->boxArray()),
                                  rad_fluxes[lev-1]->DistributionMap(), nc, 0);
                MultiFab toa_fine(top_slab(rad_fluxes[lev  ]->boxArray()),
                                  rad_fluxes[lev  ]->DistributionMap(), nc, 0);

                for (MFIter mfi(toa_crse); mfi.isValid(); ++mfi) {
                    const Box& dbx = mfi.validbox();
                    Box sbx(dbx); sbx.shift(2, 1);   // the coarse TOA, one above the top layer
                    toa_crse[mfi].template copy<RunOn::Device>((*rad_fluxes[lev-1])[mfi],
                                                               sbx, 0, dbx, 0, nc);
                }

                // Both planes live at the same z index, so the ratio in z is one and the
                // interpolation is purely horizontal.
                IntVect rr2d(refRatio(lev-1)[0], refRatio(lev-1)[1], 1);
                InterpFromCoarseLevel(toa_fine, IntVect(0,0,0),
                                      IntVect(0,0,0),
                                      toa_crse, 0, 0, nc,
                                      geom[lev-1], geom[lev],
                                      rr2d, &cell_cons_interp,
                                      domain_bcs_type, BCVars::cons_bc);

                for (MFIter mfi(toa_fine); mfi.isValid(); ++mfi) {
                    const Box& sbx = mfi.validbox();
                    Box dbx(sbx); dbx.shift(2, khi_f + 1 - khi_c);
                    (*rad_fluxes[lev])[mfi].template copy<RunOn::Device>(toa_fine[mfi],
                                                                         sbx, 0, dbx, 0, nc);
                }
            }
        }

        // LSM radiation output fields (surface fluxes needed by NoahMP)
        if (solverChoice.lsm_type != LandSurfaceType::None) {
            Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();

            for (int i = 0; i < lsm_output_names.size(); ++i) {
                int varIdx_fine   = lsm.Get_DataIdx(lev  , lsm_output_names[i]);
                int varIdx_coarse = lsm.Get_DataIdx(lev-1, lsm_output_names[i]);
                if (varIdx_fine >= 0 && varIdx_coarse >= 0) {
                    MultiFab* lsm_fine   = lsm.Get_Data_Ptr(lev  , varIdx_fine);
                    MultiFab* lsm_coarse = lsm.Get_Data_Ptr(lev-1, varIdx_coarse);
                    if (lsm_fine && lsm_coarse && lsm_coarse->nComp() > 0) {
                        // LSM data are 2D surface fields: use 2D refinement ratio and pc_interp
                        // to avoid computing z-slopes from uninitialized ghost cells
                        IntVect rr2d(refRatio(lev-1)[0], refRatio(lev-1)[1], 1);
                        InterpFromCoarseLevel(*lsm_fine, IntVect(0,0,0),
                                              IntVect(0,0,0),
                                              *lsm_coarse, 0, 0, lsm_coarse->nComp(),
                                              geom[lev-1], geom[lev],
                                              rr2d, &pc_interp,
                                              domain_bcs_type, BCVars::cons_bc);
                    }
                }
            }
        }
    };

    if (solverChoice.rad_uses_interface()) {
        BL_PROFILE_VAR("ERF::advance_radiation():RRTMGP", rrtmgp_region);

        // On the first step of a level that has just been built by interp_atmos_from_coarse,
        // interpolate the heating rates and radiation fluxes from the parent instead of
        // computing them.  That level's atmospheric state came from FillCoarsePatch, which
        // interpolates rho, theta and qv independently and so leaves them thermodynamically
        // inconsistent; running RRTMGP on it produces NaNs (see MakeNewLevelFromCoarse).
        // Interpolating also gives the LSM the radiation fields it needs for that step.
        //
        // The flag is set by whichever routine built the level and is cleared here as soon as
        // it has been acted on, so exactly one step is skipped per level creation.  Levels
        // that read a full state of their own never have it set.
        if (lev > 0 &&
            lev < static_cast<int>(rad_interp_from_coarse_pending.size()) &&
            rad_interp_from_coarse_pending[lev]) {
            amrex::Print() << "Interpolating radiation heating rates and fluxes from level " << lev-1
                           << " to level " << lev << " on the first step after that level was built\n";
            interp_rad_from_coarse();
            rad_interp_from_coarse_pending[lev] = 0;
            return;
        }

#ifdef ERF_USE_NETCDF
        MultiFab *lat_ptr = lat_m[lev].get();
        MultiFab *lon_ptr = lon_m[lev].get();
#else
        MultiFab *lat_ptr = nullptr;
        MultiFab *lon_ptr = nullptr;
#endif
        // T surf from SurfaceLayer if we have it
        MultiFab* t_surf = (m_SurfaceLayer[Orientation(Direction::z, Orientation::low)]) ? m_SurfaceLayer[Orientation(Direction::z, Orientation::low)]->get_t_surf(lev) : nullptr;

        // RRTMGP inputs names and pointers
        Vector<std::string> lsm_input_names = rad[lev]->get_lsm_input_varnames();
        Vector<MultiFab*> lsm_input_ptrs(lsm_input_names.size(),nullptr);
        for (int i(0); i<lsm_input_ptrs.size(); ++i) {
            int varIdx = lsm.Get_DataIdx(lev,lsm_input_names[i]);
            if (varIdx >= 0) { lsm_input_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx); }
        }

        // RRTMGP output names and pointers
        Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();
        Vector<MultiFab*> lsm_output_ptrs(lsm_output_names.size(),nullptr);
        for (int i(0); i<lsm_output_ptrs.size(); ++i) {
            int varIdx = lsm.Get_DataIdx(lev,lsm_output_names[i]);
            if (varIdx >= 0) { lsm_output_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx); }
        }

        // Force radiation update to sync with lsm?
        bool lsm_updated = (lev==0 && max_level>0) ? lsm.Get_LSM_Update_Status(lev) : false;

        // Enter radiation class driver
        double time_for_rad = t_old[lev] + start_time;
        rad[lev]->Run(lev, istep[lev], time_for_rad, dt_advance,
                      cons.boxArray(), geom[lev], &(cons),
                      lmask_lev[lev][0].get(), t_surf,
                      lsm_input_ptrs, lsm_output_ptrs,
                      qheating_rates[lev].get(), rad_fluxes[lev].get(),
                      z_phys_nd[lev].get()     , lat_ptr, lon_ptr,
                      lsm_updated);

        // Fill ghost cells after radiation computes (needed for interpolation to finer levels)
        // This should be fast since it only fills this level's own ghost cells
        if (solverChoice.rad_type != RadiationType::None && !rad[lev]->is_nested_patch()) {
            qheating_rates[lev]->FillBoundary(geom[lev].periodicity());
        }

        // For nested patches (fine levels that don't reach model top), radiation
        // was skipped. Interpolate the radiation fields from the parent level.
        if (lev > 0 && rad[lev]->is_nested_patch()) {
            interp_rad_from_coarse();
        }
    }
    // Two-stream radiation driver, a separate path from the IRadiation
    // models above; erf.radiation_model selects exactly one of them.
    //
    // - The call happens exactly once per slow step (from ERF::Advance).
    // - The heating rates computed here are old-state based (t^n).
    // - They are injected into the RhoTheta source only on is_slow_step
    //   (see ERF_MakeSources.cpp), so there is no duplicate forcing.
    // - istep[lev] is the CSV row index, t_old[lev] the time logged with it,
    //   and dt_advance the step size (used by the surface-energy-balance
    //   update, which runs at the post-dycore call).
    // - The sun, the site and the surface temperature come from the same
    //   sources RRTMGP uses: start_time + t for the calendar, the lat_m/lon_m
    //   fields of a WRF or metgrid grid, and the surface layer's temperature.
    else if (solverChoice.rad_type == RadiationType::TwoStream) {
#ifdef ERF_USE_NETCDF
        const MultiFab* lat_ptr = lat_m[lev].get();
        const MultiFab* lon_ptr = lon_m[lev].get();
#else
        const MultiFab* lat_ptr = nullptr;
        const MultiFab* lon_ptr = nullptr;
#endif
        const MultiFab* t_surf = (m_SurfaceLayer[Orientation::zlo()])
                               ? m_SurfaceLayer[Orientation::zlo()]->get_t_surf(lev)
                               : nullptr;
        two_stream_rad.advance(lev, istep[lev], t_old[lev], dt_advance, "pre_dycore",
                               vars_old[lev][Vars::cons], z_phys_nd[lev].get(), geom[lev],
                               lsm, qheating_rates[lev].get(), rad_fluxes[lev].get(),
                               t_surf, lat_ptr, lon_ptr,
                               t_old[lev] + start_time, use_datetime);
    }
}
