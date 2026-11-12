#include "ERF_SrcHeaders.H"

#ifdef ERF_USE_FIRE
#include <ERF_FireCoupler.H>
#endif

/**
 *  Wrapper for calling the routine that creates the slow RHS
 */
    auto slow_rhs_fun_pre = [&,zero_d=zero,one_d=one](Vector<MultiFab>& S_rhs,
                                Vector<MultiFab>& S_old,
                                Vector<MultiFab>& S_data,
                                const double old_step_time,
                                const double old_stage_time,
                                const double new_stage_time,
                                const int nrk)
    {
        BL_PROFILE("slow_rhs_fun_pre");
        //
        // Define primitive variables for all later RK stages
        // (We have already done this for the first RK step)
        // Note that it is essential this happen before the call to make_mom_sources
        //    because some of the buoyancy routines use the primitive variables
        //
        if (nrk > 0) {
            int ng_cons = S_data[IntVars::cons].nGrow();
            cons_to_prim(S_data[IntVars::cons], S_prim, ng_cons);
            make_pi_stage(S_data[IntVars::cons]);
        }

        if (verbose) Print() << std::setprecision(timeprecision)
                             << "Making slow rhs at time " << old_stage_time
                             << " for fast variables advancing from " << old_step_time
                             << " to " << new_stage_time << std::endl;

        double slow_dt = new_stage_time - old_step_time;

        const GpuArray<Real, AMREX_SPACEDIM> dxInv = fine_geom.InvCellSizeArray();

        // *************************************************************************
        // Set up flux registers if using two_way coupling
        // *************************************************************************
        YAFluxRegister* fr_as_crse = nullptr;
        YAFluxRegister* fr_as_fine = nullptr;
        if (solverChoice.coupling_type == CouplingType::TwoWay && finest_level > 0) {
            if (level < finest_level) {
                fr_as_crse = getAdvFluxReg(level+1);
                // Once per step: the anelastic RK2 integrator registers fluxes at both stages
                // (RefluxStageWeight)
                if (nrk == 0) { fr_as_crse->reset(); }
            }
            if (level > 0) {
                fr_as_fine = getAdvFluxReg(level);
            }
        }

        // *************************************************************************
        // Get multifab pointers
        // *************************************************************************

        // Canopy data for mom sources
        MultiFab* forest_drag = (solverChoice.do_forest_drag) ?
            m_forest_drag[level]->get_drag_field() : nullptr;
        MultiFab* frontal_area = (solverChoice.do_forest_drag) ?
            m_forest_drag[level]->get_frontal_area() : nullptr;

        // Immersed Forcing
        MultiFab* terrain_blank = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                                   solverChoice.buildings_type == BuildingsType::ImmersedForcing) ?
            terrain_blanking[level].get() : nullptr;
        MultiFab* terrain_blank_xface = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                                         solverChoice.buildings_type == BuildingsType::ImmersedForcing) ?
            terrain_blanking_xface[level].get() : nullptr;
        MultiFab* terrain_blank_yface = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                                         solverChoice.buildings_type == BuildingsType::ImmersedForcing) ?
            terrain_blanking_yface[level].get() : nullptr;
        MultiFab* terrain_blank_zface = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                                         solverChoice.buildings_type == BuildingsType::ImmersedForcing) ?
            terrain_blanking_zface[level].get() : nullptr;

        // Update the total moisture variable *before* computing sources since this is used in
        //    the buoyancy calculation
        if (solverChoice.moisture_type != MoistureType::None) {
            int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
            make_qt(S_data[IntVars::cons], qt, n_qstate_into_total);
        }

        MultiFab p0_to_use, base_to_use;
        MultiFab *zpn_to_use, *zpc_to_use, *ax_to_use, *ay_to_use, *az_to_use, *dJ_to_use;

        // Moving terrain
        std::unique_ptr<MultiFab> z_t_pert;
        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
        {
            z_t_pert = std::make_unique<MultiFab>(S_data[IntVars::zmom].boxArray(), S_data[IntVars::zmom].DistributionMap(), 1, 1);
            update_terrain_stage(level, old_step_time, old_stage_time, new_stage_time, slow_dt);

              p0_to_use = MultiFab(base_state_new[level], make_alias, BaseState::p0_comp, 1);
            base_to_use = MultiFab(base_state_new[level], make_alias, 0, BaseState::num_comps);
             zpn_to_use =   z_phys_nd_src[level].get();
             zpc_to_use =   z_phys_cc_src[level].get();
              ax_to_use =      ax_src[level].get();
              ay_to_use =      ay_src[level].get();
              az_to_use =      az_src[level].get();
              dJ_to_use = detJ_cc_src[level].get();

        } else {
              p0_to_use = MultiFab(base_state[level], make_alias, BaseState::p0_comp, 1);
            base_to_use = MultiFab(base_state[level], make_alias, 0, BaseState::num_comps);
             zpn_to_use = z_phys_nd[level].get();
             zpc_to_use = z_phys_cc[level].get();
              ax_to_use =        ax[level].get();
              ay_to_use =        ay[level].get();
              az_to_use =        az[level].get();
              dJ_to_use =   detJ_cc[level].get();
        }

        // Get planar averages from persistent storage for immersed forcing
        bool l_use_IF = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                         solverChoice.buildings_type == BuildingsType::ImmersedForcing);
        Table1D<Real> r_avg_to_pass = l_use_IF ? r_plane_avg[level].table() : Table1D<Real>();
        Table1D<Real> t_avg_to_pass = l_use_IF ? t_plane_avg[level].table() : Table1D<Real>();

        // *****************************************************************************
        // Construct the source terms for the cell-centered (conserved) variables
        // *****************************************************************************
        make_sources(level, nrk, slow_dt, old_stage_time,
                     S_data, S_prim, cc_src, base_state[level], zpc_to_use,
                     xvel_new, yvel_new, zvel_new,
                     qheating_rates[level].get(),
                     terrain_blank, fine_geom, solverChoice,
                     mapfac[level],
                     rhotheta_src[level].get(), rhoqt_src[level].get(),
                     dptr_wbar_sub, d_rayleigh_ptrs_at_lev, d_sinesq_at_lev,
                     turbPert,
                     r_avg_to_pass, t_avg_to_pass,
                     forest_drag, &S_old[IntVars::cons],
                     true, ib_wall_hfx[level].get(), ib_wall_ktarget[level].get(), ib_wall_qfx[level].get());

#ifdef ERF_USE_DUST
        // *****************************************************************************
        // The dust model's emission, settling and deposition, added AFTER
        // make_sources, which resets cc_src at every RK stage
        // *****************************************************************************
        if (m_dust_model) {
            m_dust_model->add_sources(level, cc_src, S_data[IntVars::cons], *dJ_to_use, fine_geom, slow_dt);
        }
#endif

        // *****************************************************************************
        // Inject the fire-to-atmosphere heat flux AFTER make_sources, which resets
        // cc_src at every RK stage. The lagged fire flux (FireCoupler) is constant
        // across all RK stages and must be added here to survive the integrator.
        // *****************************************************************************
#ifdef ERF_USE_FIRE
        // on the level the fire grid refines; coarser levels get it by average-down.
        // The heat and vapour with the WRF-Fire profile, then the smoke tracer,
        // proportional to the fuel burnt (FireCoupler::add_sources).
        if (m_fire_coupler && level == m_fire_coupler->level()) {
            bool has_moisture = (solverChoice.moisture_type != MoistureType::None);
            m_fire_coupler->add_sources(cc_src,
                                        S_old[IntVars::cons],
                                        *zpc_to_use,
                                        zpn_to_use,
                                        fine_geom,
                                        has_moisture,
                                        istep[level]);
        }
#endif

        // Heat flux from the building faces (immersed-boundary surface energy
        // balance), added to the freshly rebuilt source; no-op unless enabled.
        if (ibseb_params.enable && ibseb_params.couple_heat &&
            level < static_cast<int>(m_ibseb.size()) && m_ibseb[level]) {
            m_ibseb[level]->add_heat_flux_to_source(cc_src, S_data[IntVars::cons], fine_geom,
                                                    solverChoice.c_p, solverChoice.rdOcp);
        }

        // Nudging of potential temperature towards observations at stations
        if (solverChoice.nudging_from_observations && obs_nudging) {
            obs_nudging->add_theta_source(level, old_stage_time, static_cast<Real>(slow_dt), fine_geom,
                                          S_data[IntVars::cons], cc_src, zpn_to_use, terrain_blank);
        }

        // *****************************************************************************
        // Add nudging terms on moist variables if we have a non-zero relaxation region
        // *****************************************************************************
#if defined(ERF_USE_NETCDF)
        if ( solverChoice.use_real_bcs && (level==0) &&
             (solverChoice.moisture_type != MoistureType::None) )
        {
            Real bdy_factor = solverChoice.bdy_nudge_factor;
            Real l_rdOcp    = solverChoice.rdOcp;
            Real l_c_p      = solverChoice.c_p;
            int moist_nudge_type = solverChoice.bdy_moist_nudge_type;
            int num_q = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
            AMREX_ALWAYS_ASSERT( (solverChoice.bdy_moist_nudge_type != 2) ||
                                 (solverChoice.moisture_type != MoistureType::Morrison_NoIce &&
                                  solverChoice.moisture_type != MoistureType::SAM_NoIce) );
            add_moist_nudging_terms(S_data[IntVars::cons], cc_src, num_q, slow_dt,
                                    old_stage_time, start_time,
                                    start_bdy_time, final_bdy_time, bdy_time_interval,
                                    bdy_factor, real_width, geom[level],
                                    bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                    m_r2d, l_c_p, l_rdOcp,
                                    solverChoice.moisture_indices,
                                    solverChoice.use_wrf_bdy_density, moist_nudge_type);
        }
#endif

        // Accumulate fixed-leaf-temperature canopy exchange in the staged
        // source storage.  The existing slow-pre and slow-post paths then
        // apply rho-theta and rho-qv with their normal RK/anelastic weighting.
        if (solverChoice.do_forest_drag && solverChoice.forest_biophysics &&
            solverChoice.forest_biophysics_heat && frontal_area) {
            AddCanopyBiophysicsHeatSources(cc_src,
                                           S_data[IntVars::cons],
                                           xvel_new, yvel_new, zvel_new,
                                           frontal_area, base_to_use,
                                           solverChoice);
        }

        // *****************************************************************************
        // Define the pressure gradient
        // *****************************************************************************
        make_gradp_pert(level, solverChoice, fine_geom, S_data,
                        p0_to_use, *zpn_to_use, *zpc_to_use,
                        mapfac[level],
                        get_eb(level), gradp[level]);

        // The stage pressure, in the stage whose update the step ends with, for the reflux of
        // the pressure impulse on the coarse/fine faces (ERF_MomentumReflux.cpp)
        if (IsRefluxStage(solverChoice, level, finest_level, nrk)) {
            if (MultiFab* P_reflux = momentum_reflux_pressure(level)) {
                momentum_reflux_add_slow(level, *P_reflux, S_data[IntVars::cons], p0_to_use, qt,
                                         solverChoice.moisture_type != MoistureType::None,
                                         RefluxStageWeight(solverChoice, level, finest_level, nrk) *
                                         static_cast<Real>(slow_dt));
            }
        }

        // *****************************************************************************
        // Define the buoyancy forcing term in the z-direction
        // *****************************************************************************
        int num_q = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
        make_buoyancy(level, S_data, S_prim, qt, buoyancy, fine_geom, solverChoice, base_to_use,
                      num_q, get_eb(level), solverChoice.compressibility_type[level] == CompressibilityType::Anelastic);

        // *****************************************************************************
        // Make remaining (not gradp or buoyancy) momentum sources
        // *****************************************************************************
        make_mom_sources(level, old_stage_time, slow_dt,
                         S_data, zpn_to_use, zpc_to_use, stretched_dz_h[level],
                         xvel_new, yvel_new, zvel_new,
                         xmom_src, ymom_src, zmom_src,
                         base_to_use, forest_drag, terrain_blank,
                         terrain_blank_xface, terrain_blank_yface, terrain_blank_zface,
                         cosPhi_m[level].get(), sinPhi_m[level].get(), fine_geom, solverChoice,
                         mapfac[level],
                         (solverChoice.have_geo_wind_profile) ? d_u_geos[level].data(): nullptr,
                         (solverChoice.have_geo_wind_profile) ? d_v_geos[level].data(): nullptr,
                         dptr_wbar_sub, d_rayleigh_ptrs_at_lev, d_sinesq_at_lev, d_sinesq_stag_at_lev,
                         d_sponge_ptrs_at_lev,
                         (solverChoice.hindcast_lateral_forcing? &forecast_state_interp[level] : nullptr),
                         sounding_at_level(level), lsf_at_level(level), lsf_data[level],
                         get_eb(level), true);

#ifdef ERF_USE_MOVING_BODIES
        // The structures' force on the air (erf.structures.model), fixed over the step
        if (m_structure_model && level == m_structure_model->level(finest_level)) {
            m_structure_model->add_momentum_sources(level, xmom_src, ymom_src, zmom_src);
        }
        // On a coarser level, the model's source restricted onto it (ERF::make_structure_coarse_sources)
        if (m_structure_model && level < m_structure_model->level(finest_level) &&
            level < static_cast<int>(m_structure_crse_src.size()) &&
            m_structure_crse_src[level][0].isDefined()) {
            MultiFab::Add(xmom_src, m_structure_crse_src[level][0], 0, 0, 1, 0);
            MultiFab::Add(ymom_src, m_structure_crse_src[level][1], 0, 0, 1, 0);
            MultiFab::Add(zmom_src, m_structure_crse_src[level][2], 0, 0, 1, 0);
        }
#endif

        // The impulse of the last momentum reflux onto this level, spread over its steps in
        // this coarse step (erf.momentum_reflux, ERF_MomentumReflux.cpp)
        if (momentum_reflux_enabled()) {
            momentum_reflux_add_sources(level, xmom_src, ymom_src);
        }

        // Nudging of the velocity towards observations at stations
        if (solverChoice.nudging_from_observations && obs_nudging) {
            obs_nudging->add_momentum_sources(level, old_stage_time, static_cast<Real>(slow_dt), fine_geom,
                                              S_data[IntVars::cons], xvel_new, yvel_new, zvel_new,
                                              xmom_src, ymom_src, zmom_src, zpn_to_use,
                                              terrain_blank_xface, terrain_blank_yface,
                                              terrain_blank_zface);
        }

        // *****************************************************************************
        // Add body sources if doing flow around a body
        // *****************************************************************************
        add_thin_body_sources(xmom_src, ymom_src, zmom_src,
                              xflux_imask[level], yflux_imask[level], zflux_imask[level],
                              thin_xforce[level], thin_yforce[level], thin_zforce[level]);

        // *****************************************************************************
        // Pseudo-incompressible: advect the cell-centered scalars of this stage with the upwind
        // direction of the projection that made S_data (kept for erf_slow_rhs_post, since the
        // projection at the end of this stage overwrites pi_upw_dir)
        // *****************************************************************************
        Array<const MultiFab*,AMREX_SPACEDIM> upw_dir{AMREX_D_DECL(nullptr,nullptr,nullptr)};
        const bool l_pi_src = (solverChoice.compressibility_type[level] == CompressibilityType::PseudoIncompressible &&
                               solverChoice.terrain_type != TerrainType::EB);
        if (l_pi_src) {
            for (auto* m : {&pi_adv_rt[level], &pi_src[level]}) {
                if (!*m || (*m)->boxArray() != grids[level] || (*m)->DistributionMap() != dmap[level]) {
                    *m = std::make_unique<MultiFab>(grids[level], dmap[level], 1, 0);
                    (*m)->setVal(zero_d);
                }
            }
        }
        if (l_pi_src) {
            define_pi_upw_dir(level);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                MultiFab::Copy(*pi_upw_dir_stage[level][idim], *pi_upw_dir[level][idim], 0, 0, 1, 0);
                upw_dir[idim] = pi_upw_dir_stage[level][idim].get();
            }
        }

        // *****************************************************************************
        // Define RHS for rho, rho_theta and momenta
        // *****************************************************************************
        erf_slow_rhs_pre(level, finest_level, nrk, slow_dt, S_rhs, S_old, S_data,
                         S_prim, qt, avg_xmom[level], avg_ymom[level], avg_zmom[level], upw_dir,
                         l_pi_src ? pi_adv_rt[level].get() : nullptr,
                         xvel_new, yvel_new, zvel_new,
                         z_t_rk[level], cc_src, xmom_src, ymom_src, zmom_src, buoyancy,
                         (level > 0) ? &zmom_crse_rhs[level] : nullptr,
                         Tau[level], Tau_corr[level], Tau_EB[level],
                         SmnSmn, eddyDiffs, Hfx1, Hfx2, Hfx3, Q1fx1, Q1fx2, Q1fx3, Q2fx3, Diss, mf_surface_fluxes,
                         Hfx3_EB, fine_geom, solverChoice, m_SurfaceLayer, domain_bcs_type_d, domain_bcs_type,
                         *zpn_to_use, *zpc_to_use, *ax_to_use, *ay_to_use, *az_to_use, *dJ_to_use,
                         stretched_dz_d[level], gradp[level],
                         mapfac[level], get_eb(level),
#ifdef ERF_USE_EAMXX_SHOC
                         eamxx_shoc_interface[level].get(),
#endif
                         native_shoc_driver[level].get(),
                         fr_as_crse, fr_as_fine, &base_to_use,
                         cloud_chamber_config.active ? &cloud_chamber_config : nullptr,
                         cloud_chamber_budget.get(),
                         ib_wall_face13[level].get(), ib_wall_face23[level].get(),
                         ib_wall_face33[level].get(), ib_wall_hfx[level].get());

        // Anelastic coarse/fine faces in time (ERF_CoarseFineFaces.cpp): a coarse level adds
        // this stage's mass fluxes to those of its step, a fine level registers the coarse
        // momenta its coarse/fine faces are filled from in this stage
        if (solverChoice.compressibility_type[level] != CompressibilityType::FullyCompressible) {
            if (level < finest_level) { cf_anelastic_record(level, nrk); }
            if (level > 0)            { cf_anelastic_stage_registration(level, nrk); }
        }

        if ((solverChoice.vert_implicit_fac[level][nrk] > zero_d) && solverChoice.implicit_before_substep) {
            const Real stage_dt = static_cast<Real>(slow_dt);
            /**
             * @brief Temporary MultiFab for conserved variable updates.
             */
            MultiFab scratch(S_data[IntVars::cons].boxArray(),S_data[IntVars::cons].DistributionMap(), 2,
                             S_data[IntVars::cons].nGrowVect());
            MultiFab::Copy(scratch, S_old[IntVars::cons], 0, 0, 2, S_data[IntVars::cons].nGrowVect()); // scratch := S_old               (for rho, rhotheta)
            MultiFab::Saxpy(scratch, stage_dt, S_rhs[IntVars::cons], 0, 0, 2, 0);                  // scratch := S_old + stage_dt*Src (for rho, rhotheta)
            scratch.FillBoundary(geom[level].periodicity());

            /**
             * @brief Temporary MultiFab for x-momentum updates.
             */
            MultiFab scratch_xmom(S_data[IntVars::xmom].boxArray(),
                                  S_data[IntVars::xmom].DistributionMap(), 1,
                                  S_data[IntVars::xmom].nGrowVect());
            /**
             * @brief Temporary MultiFab for y-momentum updates.
             */
            MultiFab scratch_ymom(S_data[IntVars::ymom].boxArray(),
                                  S_data[IntVars::ymom].DistributionMap(), 1,
                                  S_data[IntVars::ymom].nGrowVect());
            if (solverChoice.implicit_momentum_diffusion) {
                MultiFab::Copy(scratch_xmom, S_old[IntVars::xmom], 0, 0, 1,  S_data[IntVars::xmom].nGrowVect()); // scratch := S_old
                MultiFab::Saxpy(scratch_xmom, stage_dt, S_rhs[IntVars::xmom], 0, 0, 1, 0);                   // scratch := S_old + stage_dt*Src
                scratch_xmom.FillBoundary(geom[level].periodicity());

                MultiFab::Copy(scratch_ymom, S_old[IntVars::ymom], 0, 0, 1,  S_data[IntVars::ymom].nGrowVect()); // scratch := S_old
                MultiFab::Saxpy(scratch_ymom, stage_dt, S_rhs[IntVars::ymom], 0, 0, 1, 0);                   // scratch := S_old + stage_dt*Src
                scratch_ymom.FillBoundary(geom[level].periodicity());
            }

            // The ghost cells of scratch still hold S_old; at a coarse/fine boundary fill them with
            // the coarse data at the time we are advancing to, since a column that ends there is
            // closed with the value on the other side
            if (level > 0) {
                const bool l_mom = solverChoice.implicit_momentum_diffusion;
                Vector<MultiFab*> scratch_mom = {l_mom ? &scratch_xmom : nullptr,
                                                 l_mom ? &scratch_ymom : nullptr,
                                                 nullptr};
                FillCoarseFineGhostsForImplicit(level, new_stage_time, scratch, {RhoTheta_comp}, scratch_mom);
            }

#include "ERF_ImplicitPre.H"

            // NOTE: the implicit part of the vertical (rho theta) diffusion is not added to the flux
            //       registers.  (rho theta) is refluxed only on anelastic levels, and there no stage
            //       whose fluxes are registered solves implicitly: RK2 turns the solve off, and
            //       MidPoint solves in stage 0 but registers stage 1 (see SolverChoice and
            //       StageFluxWeight).

            MultiFab::Saxpy(scratch, -one_d, S_old[IntVars::cons], 1, 1, 1, 0); // scratch  := (S_new - S_old)           (for rhotheta only)
            scratch.mult(one_d / stage_dt);                                 // scratch  := (S_new - S_old) / stage_dt
            MultiFab::Copy(S_rhs[IntVars::cons], scratch, 1, 1, 1, 0);        // slow_rhs := (S_new - S_old) / stage_dt (for rhotheta only)

            if (solverChoice.implicit_momentum_diffusion) {
                MultiFab::Saxpy(scratch_xmom, -one_d, S_old[IntVars::xmom], 0, 0, 1, 0); // scratch  := (S_new - S_old)
                scratch_xmom.mult(one_d / stage_dt);                                 // scratch  := (S_new - S_old) / stage_dt
                MultiFab::Copy(S_rhs[IntVars::xmom], scratch_xmom, 0, 0, 1, 0);        // slow_rhs := (S_new - S_old) / stage_dt

                MultiFab::Saxpy(scratch_ymom, -one_d, S_old[IntVars::ymom], 0, 0, 1, 0); // scratch  := (S_new - S_old)
                scratch_ymom.mult(one_d / stage_dt);                                 // scratch  := (S_new - S_old) / stage_dt
                MultiFab::Copy(S_rhs[IntVars::ymom], scratch_ymom, 0, 0, 1, 0);        // slow_rhs := (S_new - S_old) / stage_dt
            }
        }

#ifdef ERF_USE_EAMXX_SHOC
        if (solverChoice.turbChoice[level].uses_eamxx_shoc() && eamxx_shoc_interface[level]) {
            eamxx_shoc_interface[level]->add_fast_tend(S_rhs);
        }
#endif
        if (solverChoice.turbChoice[level].uses_native_shoc() && native_shoc_driver[level]) {
            // Native SHOC applies its increment directly to state; ERF no
            // longer adds a separate fast RHS contribution.
        }

        // *****************************************************************************
        // Pseudo-incompressible: the non-advective tendency of (rho theta) in this stage, the
        // whole (rho theta) RHS (now including any implicit diffusion) less its advective part
        // *****************************************************************************
        if (l_pi_src) {
            MultiFab::Copy    (*pi_src[level], S_rhs[IntVars::cons], RhoTheta_comp, 0, 1, 0);
            MultiFab::Subtract(*pi_src[level], *pi_adv_rt[level], 0, 0, 1, 0);
        }

        // *****************************************************************************
        // Update for moving terrain
        // *****************************************************************************
        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
        {
            // The column update below works on valid cells, reading only the cells above and below in the
            // same column, so it needs no ghost-cell fill first.  The ghost planes below and above the
            // column keep the values the initialization gave them, which are finite (the base-state BCs
            // leave them alone on a moving mesh, see ERFPhysBCFunct_base): the bottom-face flux is
            // extrapolated from the three lowest cells, and the top-face flux is zero (z_t = 0 there).

            MultiFab  r_hse_new (base_state_new[level], make_alias, BaseState::r0_comp, 1);
            MultiFab  p_hse_new (base_state_new[level], make_alias, BaseState::p0_comp, 1);
            MultiFab pi_hse_new (base_state_new[level], make_alias, BaseState::pi0_comp, 1);
            MultiFab th_hse_new (base_state_new[level], make_alias, BaseState::th0_comp, 1);

            MultiFab* r0_new  = &r_hse_new;
            MultiFab* p0_new  = &p_hse_new;
            MultiFab* pi0_new = &pi_hse_new;
            MultiFab* th0_new = &th_hse_new;

            // We define and evolve (rho theta)_0 in order to re-create p_0 in a way that is consistent
            //    with our update of (rho theta) but does NOT maintain dp_0 / dz = -rho_0 g.  This is why
            //    we no longer discretize the vertical pressure gradient in perturbational form.
            MultiFab rt0(p0->boxArray(),p0->DistributionMap(),1,1);
            MultiFab rt0_new(p0->boxArray(),p0->DistributionMap(),1,1);
            MultiFab r0_temp(p0->boxArray(),p0->DistributionMap(),1,1);

            // Remember this does NOT maintain dp_0 / dz = -rho_0 g, so we can no longer
            //    discretize the vertical pressure gradient in perturbational form.
            AMREX_ALWAYS_ASSERT(solverChoice.advChoice.dycore_horiz_adv_type == AdvType::Centered_2nd);
            AMREX_ALWAYS_ASSERT(solverChoice.advChoice.dycore_vert_adv_type  == AdvType::Centered_2nd);

            double dt_base = new_stage_time - old_step_time;

            const Real l_rdOcp   = solverChoice.rdOcp;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
           // Tile only in x and y: the second kernel reads r0_tmp and rt0_tmp at k-1 and k+1,
           //    which the first kernel fills tile by tile, so a tile boundary in z would read
           //    values another thread may not have written yet
           for ( MFIter mfi(*p0,TileNoZ()); mfi.isValid(); ++mfi)
           {
               const Array4<Real      > rt0_arr     = rt0.array(mfi);
               const Array4<Real      > rt0_tmp_arr = rt0_new.array(mfi);

               const Array4<Real const>  r0_arr     = r0->const_array(mfi);
               const Array4<Real      >  r0_new_arr = r0_new->array(mfi);
               const Array4<Real      >  r0_tmp_arr = r0_temp.array(mfi);

               const Array4<Real const>  p0_arr     =  p0->const_array(mfi);
               const Array4<Real      >  p0_new_arr =  p0_new->array(mfi);
               const Array4<Real      > pi0_new_arr = pi0_new->array(mfi);
               const Array4<Real      > th0_new_arr = th0_new->array(mfi);

               const Array4<Real      >& z_t_arr    =  z_t_rk[level]->array(mfi);

               const Array4<Real const>& dJ_old_arr = detJ_cc[level]->const_array(mfi);
               const Array4<Real const>& dJ_new_arr = detJ_cc_new[level]->const_array(mfi);
               const Array4<Real const>& dJ_src_arr = detJ_cc_src[level]->const_array(mfi);

               // Valid cells and the planes above and below them, which the fluxes read
               Box gbx = mfi.tilebox(); gbx.grow(2,1);
               amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                   rt0_arr(i,j,k)     = getRhoThetagivenP(p0_arr(i,j,k));
                   rt0_tmp_arr(i,j,k) = getRhoThetagivenP(p0_new_arr(i,j,k));
                   r0_tmp_arr(i,j,k)  = r0_new_arr(i,j,k);
               });

               const int klo_base = fine_geom.Domain().smallEnd(2);
               // Valid cells only; the ghost cells are filled from them below.  Every grid spans the
               // full height of the domain (checked in ERF::ReadParameters), so the three lowest cells
               // the bottom-face extrapolation reads are in this box
               Box gbx2 = mfi.tilebox();
               amrex::ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                   Real zflux_r_lo = -z_t_arr(i,j,k  ) * myhalf * (r0_tmp_arr(i,j,k) + r0_tmp_arr(i,j,k-1));
                   Real zflux_r_hi = -z_t_arr(i,j,k+1) * myhalf * (r0_tmp_arr(i,j,k) + r0_tmp_arr(i,j,k+1));
                   if (k == klo_base) {
                       // The ghost cell below klo is not on the new mesh (the vertical ghost planes of the
                       // base state are not filled on a moving mesh), so the average above is not a
                       // bottom-face value; use the quadratic extrapolation to the bottom face from the
                       // three lowest cells, of equal spacing in the computational coordinate. At the
                       // top face z_t = 0, so the ghost cell above the column only has to be finite.
                       zflux_r_lo = -z_t_arr(i,j,k) *
                           ( Real(1.875)*r0_tmp_arr(i,j,k) - Real(1.25)*r0_tmp_arr(i,j,k+1)
                            + Real(0.375)*r0_tmp_arr(i,j,k+2) );
                   }

                   Real zflux_rt_lo = zflux_r_lo * myhalf * (rt0_tmp_arr(i,j,k)/r0_tmp_arr(i,j,k) + rt0_tmp_arr(i,j,k-1)/r0_tmp_arr(i,j,k-1));
                   Real zflux_rt_hi = zflux_r_hi * myhalf * (rt0_tmp_arr(i,j,k)/r0_tmp_arr(i,j,k) + rt0_tmp_arr(i,j,k+1)/r0_tmp_arr(i,j,k+1));
                   if (k == klo_base) {
                       // theta_0 at the bottom face by the same extrapolation
                       zflux_rt_lo = zflux_r_lo *
                           ( Real(1.875)*(rt0_tmp_arr(i,j,k  )/r0_tmp_arr(i,j,k  ))
                            - Real(1.25)*(rt0_tmp_arr(i,j,k+1)/r0_tmp_arr(i,j,k+1))
                           + Real(0.375)*(rt0_tmp_arr(i,j,k+2)/r0_tmp_arr(i,j,k+2)) );
                   }

                   Real invdetJ = one / dJ_src_arr(i,j,k);

                   Real src_r  = - invdetJ * ( zflux_r_hi  - zflux_r_lo  ) * dxInv[2];
                   Real src_rt = - invdetJ * ( zflux_rt_hi - zflux_rt_lo ) * dxInv[2];

                   Real dt_b = static_cast<Real>(dt_base);
                   Real rho0_new     = dJ_old_arr(i,j,k) *  r0_arr(i,j,k) + dt_b * dJ_src_arr(i,j,k) * src_r;
                   Real rt0_tmp_new  = dJ_old_arr(i,j,k) * rt0_arr(i,j,k) + dt_b * dJ_src_arr(i,j,k) * src_rt;

                   r0_new_arr(i,j,k) = rho0_new / dJ_new_arr(i,j,k);
                   rt0_tmp_new      /=            dJ_new_arr(i,j,k);

                    p0_new_arr(i,j,k) = getPgivenRTh(rt0_tmp_new);
                   pi0_new_arr(i,j,k) = getExnergivenRTh(rt0_tmp_new, l_rdOcp);
                   th0_new_arr(i,j,k) = rt0_tmp_new / r0_new_arr(i,j,k);
               });
            } // MFIter
            // Lateral and periodic ghost cells of the updated base state, on the new mesh (the functor
            // starts with the periodic and interior FillBoundary).  The vertical
            // ghost planes are again left alone: nothing that uses this base state on a moving mesh
            // changes the solution through them.  The pressure gradient, buoyancy and cell-centred
            // sources read only cells inside the domain; the momentum sponge reads the plane below at
            // the bottom face, where w is set by the boundary condition instead.
            (*physbcs_base_new[level])(base_state_new[level],0,BaseState::num_comps,base_state_new[level].nGrowVect());
        }

#ifdef ERF_USE_NETCDF
        // Populate RHS for relaxation zones if using real bcs
        if (solverChoice.use_real_bcs && (level == 0)) {
            const Real bdy_factor = solverChoice.bdy_nudge_factor;
            if (real_width>0) {
                //
                // Note that old_stage_time is elapsed time; start_bdy_time and final_bdy_time
                //           are total time
                //
                Real l_rdOcp    = solverChoice.rdOcp;
                Real l_c_p      = solverChoice.c_p;
                realbdy_compute_interior_ghost_rhs(old_stage_time, start_time, slow_dt,
                                                   start_bdy_time, final_bdy_time, bdy_time_interval,
                                                   bdy_factor, real_width, fine_geom,
                                                   S_rhs, S_data,
                                                   bdy_data_xlo, bdy_data_xhi,
                                                   bdy_data_ylo, bdy_data_yhi,
                                                   m_r2d, l_c_p, l_rdOcp,
                                                   solverChoice.use_wrf_bdy_density,
                                                   solverChoice.bdy_rho_nudge_factor);
            }
        }
#endif
    }; // end slow_rhs_fun_pre
