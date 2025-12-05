#include "ERF_SrcHeaders.H"

/**
 *  Wrapper for calling the routine that creates the fast RHS
 */
auto acoustic_substepping_fun = [&](int fast_step, int n_sub, int nrk,
                                    Vector<MultiFab>& S_slow_rhs,
                                    const Vector<MultiFab>& S_old,
                                    Vector<MultiFab>& S_stage,
                                    Vector<MultiFab>& S_data,
                                    const Real dtau,
                                    const Real slow_dt,
                                    const Real inv_fac,
                                    const Real old_substep_time,
                                    const Real new_substep_time)
{
        BL_PROFILE("acoustic_substepping_fun");
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
            beta_s = solverChoice.beta_s;
        } else { // Fully explicit
            beta_s = -1.0;
        }

        const GpuArray<Real, AMREX_SPACEDIM> dxInv = fine_geom.InvCellSizeArray();

        // Canopy data for mom sources
        MultiFab* forest_drag = (solverChoice.do_forest_drag) ?
            m_forest_drag[level]->get_drag_field() : nullptr;

        // Immersed Forcing
        MultiFab* terrain_blank = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                                   solverChoice.buildings_type == BuildingsType::ImmersedForcing ) ?
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
                     xvel_new, yvel_new,
                     qheating_rates[level].get(),
                     terrain_blank, fine_geom, solverChoice,
                     mapfac[level],
                     dptr_rhotheta_src, dptr_rhoqt_src,
                     dptr_wbar_sub, d_rayleigh_ptrs_at_lev,
                     d_sinesq_at_lev,
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

        make_mom_sources(old_substep_time, dtau, Svec_to_use,
                         z_phys_nd[level].get(), z_phys_cc[level].get(), stretched_dz_h[level],
                         xvel_new, yvel_new, zvel_new,
                         xmom_src, ymom_src, zmom_src,
                         base_to_use, forest_drag, terrain_blank,
                         cosPhi_m[level].get(), sinPhi_m[level].get(), fine_geom, solverChoice,
                         mapfac[level],
                         (solverChoice.have_geo_wind_profile) ? d_u_geos[level].data() : nullptr,
                         (solverChoice.have_geo_wind_profile) ? d_v_geos[level].data() : nullptr,
                         dptr_wbar_sub, d_rayleigh_ptrs_at_lev,
                         d_sinesq_at_lev, d_sinesq_stag_at_lev, d_sponge_ptrs_at_lev,
                         (solverChoice.hindcast_lateral_forcing? &forecast_state_interp[level] : nullptr),
                         input_sounding_data, false);

        // Moving terrain
        std::unique_ptr<MultiFab> z_t_pert;
        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
        {
            z_t_pert = std::make_unique<MultiFab>(S_data[IntVars::zmom].boxArray(), S_data[IntVars::zmom].DistributionMap(), 1, 1);
            update_terrain_substep(level, old_substep_time, new_substep_time, dtau, S_data, z_t_pert);
        }

        bool l_use_moisture = ( solverChoice.moisture_type != MoistureType::None );

        if ( (fast_step == 0) || (solverChoice.terrain_type == TerrainType::MovingFittedMesh) )
        {
            make_fast_coeffs(level, fast_coeffs, S_stage, S_prim, pi_stage, fine_geom,
                             l_use_moisture, SolverChoice::mesh_type, solverChoice.gravity, solverChoice.c_p,
                             detJ_cc[level], r0, pi0, dtau, beta_s, phys_bc_type);
        }

        bool l_rayleigh_implicit = (solverChoice.dampingChoice.rayleigh_damping_type == RayleighDampingType::FastImplicit);
        Real l_damp_coef         = solverChoice.dampingChoice.rayleigh_dampcoef;

        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
        {
            erf_substep_MT(fast_step, nrk, level, finest_level,
                           S_slow_rhs, Svec_to_use, S_stage, S_prim, qt, pi_stage, fast_coeffs,
                           S_data, lagged_delta_rt[level], avg_xmom[level], avg_ymom[level], avg_zmom[level],
                           cc_src, xmom_src, ymom_src, zmom_src,
                           fine_geom,
                           solverChoice.gravity, solverChoice.use_lagged_delta_rt,
                           z_t_rk[level], z_t_pert.get(),
                           z_phys_nd[level], z_phys_nd_new[level], z_phys_nd_src[level],
                             detJ_cc[level],   detJ_cc_new[level],   detJ_cc_src[level],
                           dtau, beta_s, inv_fac, mapfac[level],
                           fr_as_crse, fr_as_fine, l_use_moisture, l_reflux);

        } else if (solverChoice.mesh_type == MeshType::VariableDz) {

            erf_substep_T(fast_step, nrk, level, finest_level,
                          S_slow_rhs, Svec_to_use, S_stage, S_prim, qt, pi_stage, fast_coeffs,
                          S_data, lagged_delta_rt[level], avg_xmom[level], avg_ymom[level], avg_zmom[level],
                          cc_src, xmom_src, ymom_src, zmom_src,
                          fine_geom, solverChoice.gravity,
                          z_phys_nd[level], detJ_cc[level], dtau, beta_s, inv_fac,
                          mapfac[level],
                          fr_as_crse, fr_as_fine, l_use_moisture, l_reflux,
                          (l_rayleigh_implicit) ? d_sinesq_ptrs[level].data() : nullptr, l_damp_coef);

        } else {

            erf_substep_NS(fast_step, nrk, level, finest_level,
                           S_slow_rhs, Svec_to_use, S_stage, S_prim, qt, pi_stage, fast_coeffs,
                           S_data, lagged_delta_rt[level], avg_xmom[level], avg_ymom[level], avg_zmom[level],
                           cc_src, xmom_src, ymom_src, zmom_src,
                           fine_geom, solverChoice.gravity, stretched_dz_d[level],
                           dtau, beta_s, inv_fac, mapfac[level],
                           fr_as_crse, fr_as_fine, l_use_moisture, l_reflux,
                           (l_rayleigh_implicit) ? d_sinesq_ptrs[level].data() : nullptr, l_damp_coef);
        }

        // Even if we update all the conserved variables we don't need
        // to fillpatch the slow ones every acoustic substep

        if ( (solverChoice.vert_implicit_fac[nrk] > 0.) && !solverChoice.implicit_before_substep &&
             (fast_step == n_sub-1) )
        {
            MultiFab scratch(S_data[IntVars::cons], make_alias, 0, 2);
            MultiFab scratch_xmom(S_data[IntVars::xmom], make_alias, 0, 1);
            MultiFab scratch_ymom(S_data[IntVars::ymom], make_alias, 0, 1);
#ifdef ERF_IMPLICIT_W
            MultiFab scratch_zmom(S_data[IntVars::zmom], make_alias, 0, 1);
#endif
#include "ERF_Implicit.H"
        }

        // Apply vertical velocity damping at the end of the timestep
        if ( (solverChoice.dampingChoice.w_damping) &&
             (nrk == 2) &&
             (fast_step == n_sub-1))
        {
            amrex::Real w_damping_const = solverChoice.dampingChoice.w_damping_const;
            amrex::Real w_damping_coeff = solverChoice.dampingChoice.w_damping_coeff;
            amrex::Real cflw_lim        = solverChoice.dampingChoice.w_damping_cfl;

            for ( MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi)
            {
                Box tbz = mfi.nodaltilebox(2);

                const Array4<const Real>& z_nd_arr  = z_phys_nd[level]->const_array(mfi);

                const Array4<const Real>& cell_data = S_data[IntVars::cons].const_array(mfi);
                const Array4<      Real>& rho_w     = S_data[IntVars::zmom].array(mfi);

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real dzInv;
                    if (z_nd_arr) {
                        // average z_nd to face then diff
                        if (k == domain.smallEnd(2)) {
                            dzInv = 4.0 / (
                                z_nd_arr(i  ,j  ,k+1) - z_nd_arr(i  ,j  ,k)
                              + z_nd_arr(i+1,j  ,k+1) - z_nd_arr(i+1,j  ,k)
                              + z_nd_arr(i  ,j+1,k+1) - z_nd_arr(i  ,j+1,k)
                              + z_nd_arr(i+1,j+1,k+1) - z_nd_arr(i+1,j+1,k) );
                        } else if (k == domain.bigEnd(2)+1) {
                            dzInv = 4.0 / (
                                z_nd_arr(i  ,j  ,k) - z_nd_arr(i  ,j  ,k-1)
                              + z_nd_arr(i+1,j  ,k) - z_nd_arr(i+1,j  ,k-1)
                              + z_nd_arr(i  ,j+1,k) - z_nd_arr(i  ,j+1,k-1)
                              + z_nd_arr(i+1,j+1,k) - z_nd_arr(i+1,j+1,k-1) );
                        } else {
                            // dz = 0.5 * (dz[i,j,k] + dz[i,j,k+1])
                            //    = 0.5 * (z_nd_face[i,j,k+1] - z_nd_face[i,j,k-1])
                            dzInv = 8.0 / (
                                z_nd_arr(i  ,j  ,k+1) - z_nd_arr(i  ,j  ,k-1)
                              + z_nd_arr(i+1,j  ,k+1) - z_nd_arr(i+1,j  ,k-1)
                              + z_nd_arr(i  ,j+1,k+1) - z_nd_arr(i  ,j+1,k-1)
                              + z_nd_arr(i+1,j+1,k+1) - z_nd_arr(i+1,j+1,k-1) );
                        }
                    } else {
                        dzInv = dxInv[2];
                    }

                    Real dampcoef = (w_damping_const > 0)
                                  ? w_damping_const
                                  : w_damping_coeff / (dzInv * dt_advance * dt_advance);

                    Real rho_on_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                    Real wmag = std::abs(rho_w(i,j,k) / rho_on_w_face);

                    // This is the dynamics timestep, not the slow or fast dt
                    Real cflw = wmag * dt_advance * dzInv;

                    if (cflw > cflw_lim) {
                        Real sgn_w = (rho_w(i,j,k) > 0) ? 1.0 : -1.0;
                        rho_w(i, j, k) -= rho_on_w_face * sgn_w * dampcoef * (cflw - cflw_lim) * dt_advance;
                        // Note: These prints are not immediately flushed to screen
#ifdef AMREX_USE_GPU
                        AMREX_DEVICE_PRINTF("t=%f w-damping applied at (%d,%d,%d) for w-CFL = %f > %f : %f --> %f\n",
                               old_time+dt_advance, i,j,k, cflw, cflw_lim, sgn_w*wmag, rho_w(i,j,k)/rho_on_w_face);
#else
                        printf("t=%f w-damping applied at (%d,%d,%d) for w-CFL = %f > %f : %f --> %f\n",
                               old_time+dt_advance, i,j,k, cflw, cflw_lim, sgn_w*wmag, rho_w(i,j,k)/rho_on_w_face);
#endif
                    }
                });
            }
        }

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
