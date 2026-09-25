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
    // NOTE: the coarse MultiFabs are handed to InterpFromCoarseLevel directly, with no
    //       intermediate ghosted temporary.  That is only correct because each of them
    //       carries -- and has filled -- the ghost cells this overload reads.  It builds a
    //       PhysBCFunctUseCoarseGhost, whose constructor derives the coarse halo the
    //       interpolation stencil needs and then asserts
    //
    //           cghost = min(cmf.nGrowVect(), src_ghost);
    //           AMREX_ALWAYS_ASSERT(cghost.allGE(src_ghost_outside_domain));
    //
    //       and the ParallelCopy that fills the coarse patch runs with send_ghost = cghost,
    //       i.e. it reads the coarse source's GHOST region, not only its valid region.  A
    //       coarse source with too few ghost cells therefore aborts on that assert, and one
    //       whose halo was never filled silently feeds garbage to the interpolation.  This
    //       is why qheating_rates and rad_fluxes are both defined with a (1,1,1) halo in
    //       ERF_MakeNewArrays.cpp, zeroed there, and FillBoundary'd below.
    // Can no column model run on this level, so that its radiation fields must be
    // interpolated from the parent instead?
    //
    // RRTMGP records its own predicate at Init and it stays authoritative for that model.
    // The two-stream model has no IRadiation object -- reading rad[] would dereference a
    // null pointer -- so derive it from the grids.
    //
    // The test is per box, not on the level's bounding box. A column sweep needs a whole
    // column inside ONE box, so the question is whether every box spans the domain in z,
    // not whether the boxes together do. A bounding-box test gets three layouts right and
    // one wrong: a level tagged at genuinely different heights in different horizontal
    // regions -- surface convection near klo here, cloud tops near khi there -- has a
    // bounding box that reaches both domain ends while no single box holds a column. That
    // level is no more sweepable than a shallow nest, and this sends it down the same path.
    auto level_needs_interpolation = [&] (int l) -> bool
    {
        if (l <= 0) { return false; }
        if (solverChoice.rad_uses_interface() && rad[l]) { return rad[l]->is_nested_patch(); }
        const Box& dom = geom[l].Domain();
        for (int ibox = 0; ibox < grids[l].size(); ++ibox) {
            const Box& b = grids[l][ibox];
            if (b.smallEnd(2) != dom.smallEnd(2) || b.bigEnd(2) != dom.bigEnd(2)) { return true; }
        }
        return false;
    };

    auto interp_rad_from_coarse = [&] ()
    {
        // Ensure the parent level's ghost cells are filled before interpolation.  This is
        // needed even when radiation did not run this step, especially with two-way coupling
        // where the grid structure may have changed.  A nested patch is exempt: it never
        // runs radiation itself, so its halo is already whatever its own pass through this
        // lambda interpolated into it.
        if (!level_needs_interpolation(lev-1)) {
            qheating_rates[lev-1]->FillBoundary(geom[lev-1].periodicity());
            if (rad_fluxes[lev-1]) {
                // The whole halo, as for qheating_rates -- the z ghosts matter too, since
                // grids split in z (amr.no_box_split_dir = -1) put interior z interfaces in
                // them -- but with z forced non-periodic.  The z-ghost at khi+1 holds the
                // top-of-atmosphere interface, physical data rather than a halo (see below);
                // it lies outside the domain, so no ordinary exchange reaches it, but a
                // z-periodic one would overwrite it with the bottom of the column.
                const IntVect per = geom[lev-1].periodicity().intVect();
                rad_fluxes[lev-1]->FillBoundary(Periodicity(IntVect(per[0],per[1],0)));
            }
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
            // ERF_MakeNewArrays.cpp), which is physical data rather than a halo.  On the
            // coarse side that cell does now fall inside the ghost region the ParallelCopy
            // reads, so it serves as the stencil's neighbor above the top layer -- the right
            // value for that role, being the next interface up in the same sequence -- but
            // it is only ever read there.  On the fine side, a level that spans the full
            // column has its TOA outside the fine domain, and the call above is asked for no
            // ghost cells outside the domain, so nothing is written there at all.  A nested
            // patch does not reach the model top, so its top plane is an ordinary interior
            // interface that the call above has already filled.
            //
            // Move both planes into valid index space -- a single layer at the coarse domain
            // top -- so the ordinary machinery can interpolate them horizontally, then put
            // the result back in the fine level's ghost cell.
            if (!level_needs_interpolation(lev)) {
                const int khi_c = geom[lev-1].Domain().bigEnd(2);
                const int khi_f = geom[lev  ].Domain().bigEnd(2);

                // Flatten a level's grids onto the single layer at k = khi_c, keeping only
                // the boxes that reach the top of their own domain.  A box that stops short
                // of it carries no TOA plane at all: flattening it would put a duplicate box
                // in the BoxArray, for which FillBoundary is ill-defined, and would ask the
                // copy below to read a z index outside the fab.  Every box spans the column
                // under the default amr.no_box_split_dir = 2, so this filter keeps all of
                // them there; it matters only for amr.no_box_split_dir = -1.  idx_out maps a
                // slab box back to the box it came from, since dropping boxes breaks the
                // one-to-one correspondence an MFIter would otherwise rely on.
                auto top_slab = [&] (const MultiFab& mf_in, int khi_in, BoxArray& ba_out,
                                     DistributionMapping& dm_out, Vector<int>& idx_out)
                {
                    const BoxArray&            ba_in = mf_in.boxArray();
                    const DistributionMapping& dm_in = mf_in.DistributionMap();
                    BoxList bl;
                    Vector<int> pmap;
                    for (int i = 0, n = int(ba_in.size()); i < n; ++i) {
                        if (ba_in[i].bigEnd(2) != khi_in) { continue; }
                        Box b = ba_in[i];
                        b.setSmall(2, khi_c); b.setBig(2, khi_c);
                        bl.push_back(b);
                        pmap.push_back(dm_in[i]);
                        idx_out.push_back(i);
                    }
                    ba_out = BoxArray(std::move(bl));
                    dm_out = DistributionMapping(std::move(pmap));
                };

                BoxArray            ba_toa_c, ba_toa_f;
                DistributionMapping dm_toa_c, dm_toa_f;
                Vector<int>         idx_toa_c, idx_toa_f;
                top_slab(*rad_fluxes[lev-1], khi_c, ba_toa_c, dm_toa_c, idx_toa_c);
                top_slab(*rad_fluxes[lev  ], khi_f, ba_toa_f, dm_toa_f, idx_toa_f);

                // toa_crse is the coarse source of an InterpFromCoarseLevel, so it must carry
                // the halo that stencil reads -- (1,1,0), the z ratio being one -- and that
                // halo must be defined.  Unlike rad_fluxes it is freshly allocated here and
                // the loop below writes valid boxes only, so zero it first: FillBoundary can
                // define only the ghosts backed by another box's valid data, and the ones
                // outside the domain are left to the setVal.
                MultiFab toa_crse(ba_toa_c, dm_toa_c, nc, IntVect(1,1,0));
                MultiFab toa_fine(ba_toa_f, dm_toa_f, nc, 0);
                toa_crse.setVal(Real(0.0));

                for (MFIter mfi(toa_crse); mfi.isValid(); ++mfi) {
                    const Box& dbx = mfi.validbox();
                    Box sbx(dbx); sbx.shift(2, 1);   // the coarse TOA, one above the top layer
                    toa_crse[mfi].template copy<RunOn::Device>(
                        (*rad_fluxes[lev-1])[idx_toa_c[mfi.index()]], sbx, 0, dbx, 0, nc);
                }
                // Unconditional, unlike the exchange at the top of this lambda: toa_crse is a
                // new MultiFab every step, so its halo is undefined here whatever the parent
                // level is -- being a nested patch says nothing about it.
                toa_crse.FillBoundary(geom[lev-1].periodicity());

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
                    (*rad_fluxes[lev])[idx_toa_f[mfi.index()]].template copy<RunOn::Device>(
                        toa_fine[mfi], sbx, 0, dbx, 0, nc);
                }
            }
        }

        // LSM radiation output fields (surface fluxes needed by NoahMP).  Gated on
        // rad_uses_interface() as well: the names come from rad[lev], which is null under
        // erf.radiation_model = TwoStream (see the Constructors), and that model reaches
        // this lambda too now that it runs on a hierarchy.  The two-stream model does not
        // feed the land surface these fields, so there is nothing to interpolate for it.
        if (solverChoice.rad_uses_interface() && solverChoice.lsm_type != LandSurfaceType::None) {
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

    // On the first step of a level that has just been built by interp_atmos_from_coarse,
    // interpolate the heating rates and radiation fluxes from the parent instead of
    // computing them.  That level's atmospheric state came from FillCoarsePatch, which
    // interpolates rho, theta and qv independently and so leaves them thermodynamically
    // inconsistent; running RRTMGP on it produces NaNs (see MakeNewLevelFromCoarse), and the
    // two-stream column sweep, which takes its optical properties from the same rho, theta
    // and qv, has no more claim on that state than RRTMGP does.  Interpolating also gives
    // the LSM the radiation fields it needs for that step.
    //
    // This sits ahead of the model dispatch because the flag is set for whichever model is
    // running (ERF_MakeNewLevel.cpp) and the reason for honouring it is a property of the
    // state, not of the model.  Skipping the pre-dycore sweep costs the two-stream model
    // nothing else: the flag is only ever set for lev > 0, and the surface energy balance
    // that the post-dycore call advances runs on level 0 alone.
    //
    // The flag is set by whichever routine built the level and is cleared here as soon as
    // it has been acted on, so exactly one step is skipped per level creation.  Levels
    // that read a full state of their own never have it set.
    if (lev > 0 &&
        (solverChoice.rad_uses_interface() || solverChoice.rad_type == RadiationType::TwoStream) &&
        lev < static_cast<int>(rad_interp_from_coarse_pending.size()) &&
        rad_interp_from_coarse_pending[lev]) {
        amrex::Print() << "Interpolating radiation heating rates and fluxes from level " << lev-1
                       << " to level " << lev << " on the first step after that level was built\n";
        interp_rad_from_coarse();
        rad_interp_from_coarse_pending[lev] = 0;
        return;
    }

    if (solverChoice.rad_uses_interface()) {
        BL_PROFILE_VAR("ERF::advance_radiation():RRTMGP", rrtmgp_region);

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
        if (solverChoice.rad_type != RadiationType::None && !level_needs_interpolation(lev)) {
            qheating_rates[lev]->FillBoundary(geom[lev].periodicity());
        }

        // For nested patches (fine levels that don't reach model top), radiation
        // was skipped. Interpolate the radiation fields from the parent level.
        if (lev > 0 && level_needs_interpolation(lev)) {
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
        // A nested patch has no complete column, so the sweep cannot run on it. Take the
        // same route RRTMGP takes: interpolate this level's heating rates and fluxes from
        // the parent rather than refusing the configuration. The two models now differ only
        // in how a level that *does* span the column is solved.
        if (lev > 0 && level_needs_interpolation(lev)) {
            interp_rad_from_coarse();
            return;
        }
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

        // Fill this level's halo so a finer level can interpolate from it. The
        // InterpFromCoarseLevel overload above reads the coarse source's ghost cells, not
        // only its valid region (see the note on interp_rad_from_coarse), so an unfilled
        // halo here would feed garbage to a nested child.
        qheating_rates[lev]->FillBoundary(geom[lev].periodicity());
        if (rad_fluxes[lev]) {
            // z left non-periodic: the ghost at khi+1 holds the top-of-atmosphere
            // interface, which is physical data rather than a halo.
            const IntVect per = geom[lev].periodicity().intVect();
            rad_fluxes[lev]->FillBoundary(Periodicity(IntVect(per[0],per[1],0)));
        }
    }
}
