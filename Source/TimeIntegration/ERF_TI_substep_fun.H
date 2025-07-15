#include "ERF_SrcHeaders.H"

/**
 *  Wrapper for calling the routine that creates the fast RHS
 */
auto fast_rhs_fun = [&](int fast_step, int /*n_sub*/, int nrk,
                        Vector<MultiFab>& S_slow_rhs,
                        const Vector<MultiFab>& S_old,
                        Vector<MultiFab>& S_stage,
                        Vector<MultiFab>& S_data,
                        Vector<MultiFab>& S_scratch,
                        const Real dtau,
                        const Real inv_fac,
                        const Real old_substep_time,
                        const Real new_substep_time)
    {
        BL_PROFILE("fast_rhs_fun");
        if (verbose) amrex::Print() << "Fast time integration at level " << level
                                    << std::setprecision(timeprecision)
                                    << " from " << old_substep_time << " to " << new_substep_time
                                    << " with dt = " << dtau << std::endl;

        // Define beta_s here so that it is consistent between where we make the fast coefficients
        //    and where we use them
        // Per p2902 of Klemp-Skamarock-Dudhia-2007
        // beta_s = -1.0 : fully explicit
        // beta_s =  1.0 : fully implicit
        Real beta_s;
        if (solverChoice.substepping_type[level] == SubsteppingType::Implicit) {
            beta_s = 0.1;
        } else { // Fully explicit
            beta_s = -1.0;
        }

        // Canopy data for mom sources
        MultiFab* forest_drag = (solverChoice.do_forest_drag) ?
            m_forest_drag[level]->get_drag_field() : nullptr;

        // Immersed Forcing
        MultiFab* terrain_blank = (solverChoice.terrain_type == TerrainType::ImmersedForcing) ?
            terrain_blanking[level].get() : nullptr;

        Vector<MultiFab> Svec_to_use;
        if (fast_step == 0) {
            // If this is the first substep we pass in S_old as the previous step's solution
            Svec_to_use.push_back(MultiFab(S_old[IntVars::cons],make_alias,0,S_old[IntVars::cons].nComp()));
            Svec_to_use.push_back(MultiFab(S_old[IntVars::xmom],make_alias,0,1));
            Svec_to_use.push_back(MultiFab(S_old[IntVars::ymom],make_alias,0,1));
            Svec_to_use.push_back(MultiFab(S_old[IntVars::zmom],make_alias,0,1));
        } else {
            // If this is not the first substep we pass in S_data as the previous step's solution
            Svec_to_use.push_back(MultiFab(S_data[IntVars::cons],make_alias,0,S_data[IntVars::cons].nComp()));
            Svec_to_use.push_back(MultiFab(S_data[IntVars::xmom],make_alias,0,1));
            Svec_to_use.push_back(MultiFab(S_data[IntVars::ymom],make_alias,0,1));
            Svec_to_use.push_back(MultiFab(S_data[IntVars::zmom],make_alias,0,1));
        }

        make_sources(level, nrk, dtau, old_substep_time,
                     Svec_to_use, S_prim, cc_src, base_state[level], z_phys_cc[level].get(),
                     qheating_rates[level].get(),
                     terrain_blank, fine_geom, solverChoice,
                     mapfac[level],
                     dptr_rhotheta_src, dptr_rhoqt_src,
                     dptr_wbar_sub, d_rayleigh_ptrs_at_lev,
                     input_sounding_data, turbPert, false);

        // *************************************************************************
        // Set up flux registers if using two_way coupling
        // *************************************************************************
        const bool l_reflux = ( (solverChoice.coupling_type == CouplingType::TwoWay) && (nrk == 2) && (finest_level > 0) );

        YAFluxRegister* fr_as_crse = nullptr;
        YAFluxRegister* fr_as_fine = nullptr;
        if (l_reflux) {
            if (level < finest_level) {
                fr_as_crse = getAdvFluxReg(level+1);
            }
            if (level > 0) {
                fr_as_fine = getAdvFluxReg(level);
            }
        }

        MultiFab base_to_use;
        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh ) {
            base_to_use = MultiFab(base_state_new[level], make_alias, 0, BaseState::num_comps);
        } else {
            base_to_use = MultiFab(base_state[level], make_alias, 0, BaseState::num_comps);
        }

        make_mom_sources(old_substep_time, Svec_to_use,
                         *z_phys_nd[level].get(), *z_phys_cc[level].get(), stretched_dz_h[level],
                         xvel_new, yvel_new, zvel_new,
                         xmom_src, ymom_src, zmom_src,
                         base_to_use, forest_drag, terrain_blank,
                         cosPhi_m[level].get(), sinPhi_m[level].get(), fine_geom, solverChoice,
                         mapfac[level],
                         (solverChoice.have_geo_wind_profile) ? d_u_geos[level].data() : nullptr,
                         (solverChoice.have_geo_wind_profile) ? d_v_geos[level].data() : nullptr,
                         dptr_wbar_sub, d_rayleigh_ptrs_at_lev, d_sponge_ptrs_at_lev,
                         input_sounding_data, false);

        // Moving terrain
        std::unique_ptr<MultiFab> z_t_pert;
        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
        {
            z_t_pert = std::make_unique<MultiFab>(S_data[IntVars::zmom].boxArray(), S_data[IntVars::zmom].DistributionMap(), 1, 1);
            update_terrain_substep(level, old_substep_time, new_substep_time, dtau, S_data, z_t_pert);
        }

        if ( (fast_step == 0) || (solverChoice.terrain_type == TerrainType::MovingFittedMesh) )
        {
            make_fast_coeffs(level, fast_coeffs, S_stage, S_prim, pi_stage, fine_geom,
                             l_use_moisture, SolverChoice::mesh_type, solverChoice.gravity, solverChoice.c_p,
                             detJ_cc[level], r0, pi0, dtau, beta_s, phys_bc_type);
        }

        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
        {
            erf_fast_rhs_MT(fast_step, nrk, level, finest_level,
                            S_slow_rhs, Svec_to_use, S_stage, S_prim, qt, pi_stage, fast_coeffs,
                            S_data, S_scratch,
                            cc_src, xmom_src, ymom_src, zmom_src,
                            fine_geom,
                            solverChoice.gravity, solverChoice.use_lagged_delta_rt,
                            z_t_rk[level], z_t_pert.get(),
                            z_phys_nd[level], z_phys_nd_new[level], z_phys_nd_src[level],
                              detJ_cc[level],   detJ_cc_new[level],   detJ_cc_src[level],
                            dtau, beta_s, inv_fac, mapfac[level],
                            fr_as_crse, fr_as_fine, l_use_moisture, l_reflux, l_implicit_substepping);

        } else if (solverChoice.mesh_type != MeshType::ConstantDz) {

            erf_fast_rhs_T(fast_step, nrk, level, finest_level,
                           S_slow_rhs, Svec_to_use, S_stage, S_prim, qt, pi_stage, fast_coeffs,
                           S_data, S_scratch,
                           cc_src, xmom_src, ymom_src, zmom_src,
                           fine_geom, solverChoice.gravity,
                           z_phys_nd[level], detJ_cc[level], dtau, beta_s, inv_fac,
                           mapfac[level],
                           fr_as_crse, fr_as_fine, l_use_moisture, l_reflux, l_implicit_substepping);

        } else {

            erf_fast_rhs_N(fast_step, nrk, level, finest_level,
                           S_slow_rhs, Svec_to_use, S_stage, S_prim, qt, pi_stage, fast_coeffs,
                           S_data, S_scratch,
                           cc_src, xmom_src, ymom_src, zmom_src,
                           fine_geom, solverChoice.gravity,
                           dtau, beta_s, inv_fac, mapfac[level],
                           fr_as_crse, fr_as_fine, l_use_moisture, l_reflux, l_implicit_substepping);
        }

        // Even if we update all the conserved variables we don't need
        // to fillpatch the slow ones every acoustic substep

        // NOTE: Numerical diffusion is tested on in FillPatchIntermediate and dictates the size of the
        //       box over which VelocityToMomentum is computed. V2M requires one more ghost cell be
        //       filled for rho than velocity. This logical condition ensures we fill enough ghost cells
        //       when use_num_diff is true.
        int ng_cons = S_data[IntVars::cons].nGrowVect().max() - 1;
        int ng_vel  = S_data[IntVars::xmom].nGrowVect().max();
        if (!solverChoice.use_num_diff) {
            ng_cons = 1;
            ng_vel  = 1;
        }
        apply_bcs(S_data, new_substep_time, ng_cons, ng_vel, fast_only=true, vel_and_mom_synced=false);
    };
