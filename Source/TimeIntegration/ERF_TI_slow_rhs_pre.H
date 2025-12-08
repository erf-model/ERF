#include "ERF_SrcHeaders.H"

/**
 *  Wrapper for calling the routine that creates the slow RHS
 */
    auto slow_rhs_fun_pre = [&](Vector<MultiFab>& S_rhs,
                                Vector<MultiFab>& S_old,
                                Vector<MultiFab>& S_data,
                                const Real old_step_time,
                                const Real old_stage_time,
                                const Real new_stage_time,
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
            cons_to_prim(S_data[IntVars::cons], ng_cons);
        }

        if (verbose) Print() << std::setprecision(timeprecision)
                             << "Making slow rhs at time " << old_stage_time
                             << " for fast variables advancing from " << old_step_time
                             << " to " << new_stage_time << std::endl;

        Real slow_dt = new_stage_time - old_step_time;

        const GpuArray<Real, AMREX_SPACEDIM> dxInv = fine_geom.InvCellSizeArray();

        // *************************************************************************
        // Set up flux registers if using two_way coupling
        // *************************************************************************
        YAFluxRegister* fr_as_crse = nullptr;
        YAFluxRegister* fr_as_fine = nullptr;
        if (solverChoice.coupling_type == CouplingType::TwoWay && finest_level > 0) {
            if (level < finest_level) {
                fr_as_crse = getAdvFluxReg(level+1);
                fr_as_crse->reset();
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

        // Immersed Forcing
        MultiFab* terrain_blank = (solverChoice.terrain_type == TerrainType::ImmersedForcing ||
                                   solverChoice.buildings_type == BuildingsType::ImmersedForcing) ?
            terrain_blanking[level].get() : nullptr;

        // Update the total moisture variable *before* computing sources since this is used in
        //    the buoyancy calculation
        if (solverChoice.moisture_type != MoistureType::None) {
            make_qt(S_data[IntVars::cons], qt);
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

        // *****************************************************************************
        // Construct the source terms for the cell-centered (conserved) variables
        // *****************************************************************************
        make_sources(level, nrk, slow_dt, old_stage_time,
                     S_data, S_prim, cc_src, base_state[level], zpc_to_use,
                     xvel_new, yvel_new,
                     qheating_rates[level].get(),
                     terrain_blank, fine_geom, solverChoice,
                     mapfac[level],
                     dptr_rhotheta_src, dptr_rhoqt_src,
                     dptr_wbar_sub, d_rayleigh_ptrs_at_lev, d_sinesq_at_lev,
                     input_sounding_data, turbPert, true);

        // *****************************************************************************
        // Define the pressure gradient
        // *****************************************************************************
        make_gradp_pert(level, solverChoice, fine_geom, S_data,
                        p0_to_use, *zpn_to_use, *zpc_to_use,
                        mapfac[level],
                        get_eb(level), gradp[level]);

        // *****************************************************************************
        // Define the buoyancy forcing term in the z-direction
        // *****************************************************************************
        make_buoyancy(level, S_data, S_prim, qt, buoyancy, fine_geom, solverChoice, base_to_use,
                      micro->Get_Qstate_Moist_Size(), get_eb(level), solverChoice.anelastic[level]);

        // *****************************************************************************
        // Make remaining (not gradp or buoyancy) momentum sources
        // *****************************************************************************
        make_mom_sources(old_stage_time, slow_dt,
                         S_data, zpn_to_use, zpc_to_use, stretched_dz_h[level],
                         xvel_new, yvel_new, zvel_new,
                         xmom_src, ymom_src, zmom_src,
                         base_to_use, forest_drag, terrain_blank,
                         cosPhi_m[level].get(), sinPhi_m[level].get(), fine_geom, solverChoice,
                         mapfac[level],
                         (solverChoice.have_geo_wind_profile) ? d_u_geos[level].data(): nullptr,
                         (solverChoice.have_geo_wind_profile) ? d_v_geos[level].data(): nullptr,
                         dptr_wbar_sub, d_rayleigh_ptrs_at_lev, d_sinesq_at_lev, d_sinesq_stag_at_lev,
                         d_sponge_ptrs_at_lev,
                         (solverChoice.hindcast_lateral_forcing? &forecast_state_interp[level] : nullptr),
                         input_sounding_data, true);

        // *****************************************************************************
        // Add body sources if doing flow around a body
        // *****************************************************************************
        add_thin_body_sources(xmom_src, ymom_src, zmom_src,
                              xflux_imask[level], yflux_imask[level], zflux_imask[level],
                              thin_xforce[level], thin_yforce[level], thin_zforce[level]);

        // *****************************************************************************
        // Define RHS for rho, rho_theta and momenta
        // *****************************************************************************
        erf_slow_rhs_pre(level, finest_level, nrk, slow_dt, S_rhs, S_old, S_data,
                         S_prim, qt, avg_xmom[level], avg_ymom[level], avg_zmom[level],
                         xvel_new, yvel_new, zvel_new,
                         z_t_rk[level], cc_src, xmom_src, ymom_src, zmom_src, buoyancy,
                         (level > 0) ? &zmom_crse_rhs[level] : nullptr,
                         Tau[level], Tau_corr[level],
                         SmnSmn, eddyDiffs, Hfx1, Hfx2, Hfx3, Q1fx1, Q1fx2, Q1fx3, Q2fx3, Diss,
                         fine_geom, solverChoice, m_SurfaceLayer, domain_bcs_type_d, domain_bcs_type,
                         *zpn_to_use, *zpc_to_use, *ax_to_use, *ay_to_use, *az_to_use, *dJ_to_use,
                         stretched_dz_d[level], gradp[level],
                         mapfac[level], get_eb(level),
#ifdef ERF_USE_SHOC
                         shoc_interface[level],
#endif
                         fr_as_crse, fr_as_fine);

        if ((solverChoice.vert_implicit_fac[nrk] > 0.) && solverChoice.implicit_before_substep) {
            MultiFab scratch(S_data[IntVars::cons].boxArray(),S_data[IntVars::cons].DistributionMap(), 2,
                             S_data[IntVars::cons].nGrowVect());
            MultiFab::Copy(scratch, S_old[IntVars::cons], 0, 0, 2, S_data[IntVars::cons].nGrowVect()); // scratch := S_old               (for rho, rhotheta)
            MultiFab::Saxpy(scratch, slow_dt, S_rhs[IntVars::cons], 0, 0, 2, 0);                       // scratch := S_old + slow_dt*Src (for rho, rhotheta)
            scratch.FillBoundary(geom[level].periodicity());

            MultiFab scratch_xmom(S_data[IntVars::xmom].boxArray(),
                                  S_data[IntVars::xmom].DistributionMap(), 1,
                                  S_data[IntVars::xmom].nGrowVect());
            MultiFab scratch_ymom(S_data[IntVars::ymom].boxArray(),
                                  S_data[IntVars::ymom].DistributionMap(), 1,
                                  S_data[IntVars::ymom].nGrowVect());
#ifdef ERF_IMPLICIT_W
            MultiFab scratch_zmom(S_data[IntVars::zmom].boxArray(),
                                  S_data[IntVars::zmom].DistributionMap(), 1,
                                  S_data[IntVars::zmom].nGrowVect());
#endif
            if (solverChoice.implicit_momentum_diffusion) {
                MultiFab::Copy(scratch_xmom, S_old[IntVars::xmom], 0, 0, 1,  S_data[IntVars::xmom].nGrowVect()); // scratch := S_old
                MultiFab::Saxpy(scratch_xmom, slow_dt, S_rhs[IntVars::xmom], 0, 0, 1, 0);                        // scratch := S_old + slow_dt*Src
                scratch_xmom.FillBoundary(geom[level].periodicity());

                MultiFab::Copy(scratch_ymom, S_old[IntVars::ymom], 0, 0, 1,  S_data[IntVars::ymom].nGrowVect()); // scratch := S_old
                MultiFab::Saxpy(scratch_ymom, slow_dt, S_rhs[IntVars::ymom], 0, 0, 1, 0);                        // scratch := S_old + slow_dt*Src
                scratch_ymom.FillBoundary(geom[level].periodicity());
#ifdef ERF_IMPLICIT_W
                MultiFab::Copy(scratch_zmom, S_old[IntVars::zmom], 0, 0, 1,  S_data[IntVars::zmom].nGrowVect()); // scratch := S_old
                MultiFab::Saxpy(scratch_zmom, slow_dt, S_rhs[IntVars::zmom], 0, 0, 1, 0);                        // scratch := S_old + slow_dt*Src
                scratch_zmom.FillBoundary(geom[level].periodicity());
#endif
            }

#include "ERF_Implicit.H"

            MultiFab::Saxpy(scratch, -1.0, S_old[IntVars::cons], 1, 1, 1, 0); // scratch  := (S_new - S_old)           (for rhotheta only)
            scratch.mult(1.0 / slow_dt);                                      // scratch  := (S_new - S_old) / slow_dt
            MultiFab::Copy(S_rhs[IntVars::cons], scratch, 1, 1, 1, 0);        // slow_rhs := (S_new - S_old) / slow_dt (for rhotheta only)

            if (solverChoice.implicit_momentum_diffusion) {
                MultiFab::Saxpy(scratch_xmom, -1.0, S_old[IntVars::xmom], 0, 0, 1, 0); // scratch  := (S_new - S_old)
                scratch_xmom.mult(1.0 / slow_dt);                                      // scratch  := (S_new - S_old) / slow_dt
                MultiFab::Copy(S_rhs[IntVars::xmom], scratch_xmom, 0, 0, 1, 0);        // slow_rhs := (S_new - S_old) / slow_dt

                MultiFab::Saxpy(scratch_ymom, -1.0, S_old[IntVars::ymom], 0, 0, 1, 0); // scratch  := (S_new - S_old)
                scratch_ymom.mult(1.0 / slow_dt);                                      // scratch  := (S_new - S_old) / slow_dt
                MultiFab::Copy(S_rhs[IntVars::ymom], scratch_ymom, 0, 0, 1, 0);        // slow_rhs := (S_new - S_old) / slow_dt
#ifdef ERF_IMPLICIT_W
                MultiFab::Saxpy(scratch_zmom, -1.0, S_old[IntVars::zmom], 0, 0, 1, 0); // scratch  := (S_new - S_old)
                scratch_zmom.mult(1.0 / slow_dt);                                      // scratch  := (S_new - S_old) / slow_dt
                MultiFab::Copy(S_rhs[IntVars::zmom], scratch_zmom, 0, 0, 1, 0);        // slow_rhs := (S_new - S_old) / slow_dt
#endif
            }
        }

#ifdef ERF_USE_SHOC
        if (solverChoice.use_shoc) {
            shoc_interface[level]->add_fast_tend(S_rhs);
        }
#endif

        // *****************************************************************************
        // Update for moving terrain
        // *****************************************************************************
        if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
        {
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

            Real dt_base  = (new_stage_time - old_step_time);

            const Real l_rdOcp   = solverChoice.rdOcp;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
           for ( MFIter mfi(*p0,TilingIfNotGPU()); mfi.isValid(); ++mfi)
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

               Box gbx = mfi.growntilebox({1,1,1});
               amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                   rt0_arr(i,j,k)     = getRhoThetagivenP(p0_arr(i,j,k));
                   rt0_tmp_arr(i,j,k) = getRhoThetagivenP(p0_new_arr(i,j,k));
                   r0_tmp_arr(i,j,k)  = r0_new_arr(i,j,k);
               });

               Box gbx2 = mfi.growntilebox({1,1,0});
               amrex::ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                   Real zflux_r_lo = -z_t_arr(i,j,k  ) * 0.5 * (r0_tmp_arr(i,j,k) + r0_tmp_arr(i,j,k-1));
                   Real zflux_r_hi = -z_t_arr(i,j,k+1) * 0.5 * (r0_tmp_arr(i,j,k) + r0_tmp_arr(i,j,k+1));

                   Real zflux_rt_lo = zflux_r_lo * 0.5 * (rt0_tmp_arr(i,j,k)/r0_tmp_arr(i,j,k) + rt0_tmp_arr(i,j,k-1)/r0_tmp_arr(i,j,k-1));
                   Real zflux_rt_hi = zflux_r_hi * 0.5 * (rt0_tmp_arr(i,j,k)/r0_tmp_arr(i,j,k) + rt0_tmp_arr(i,j,k+1)/r0_tmp_arr(i,j,k+1));

                   Real invdetJ = 1.0 / dJ_src_arr(i,j,k);

                   Real src_r  = - invdetJ * ( zflux_r_hi  - zflux_r_lo  ) * dxInv[2];
                   Real src_rt = - invdetJ * ( zflux_rt_hi - zflux_rt_lo ) * dxInv[2];

                   Real rho0_new     = dJ_old_arr(i,j,k) *  r0_arr(i,j,k) + dt_base * dJ_src_arr(i,j,k) * src_r;
                   Real rt0_tmp_new  = dJ_old_arr(i,j,k) * rt0_arr(i,j,k) + dt_base * dJ_src_arr(i,j,k) * src_rt;

                   r0_new_arr(i,j,k) = rho0_new / dJ_new_arr(i,j,k);
                   rt0_tmp_new      /=            dJ_new_arr(i,j,k);

                    p0_new_arr(i,j,k) = getPgivenRTh(rt0_tmp_new);
                   pi0_new_arr(i,j,k) = getExnergivenRTh(rt0_tmp_new, l_rdOcp);
                   th0_new_arr(i,j,k) = rt0_tmp_new / r0_new_arr(i,j,k);
               });
            } // MFIter
            r0_new->FillBoundary(fine_geom.periodicity());
            p0_new->FillBoundary(fine_geom.periodicity());
            th0_new->FillBoundary(fine_geom.periodicity());
        }

#ifdef ERF_USE_NETCDF
        // Populate RHS for relaxation zones if using real bcs
        if (solverChoice.use_real_bcs && (level == 0)) {
            if (real_width>0) {
                realbdy_compute_interior_ghost_rhs(bdy_time_interval,
                                                   new_stage_time, slow_dt, stop_time-start_time,
                                                   real_width, real_set_width, fine_geom,
                                                   S_rhs, S_old, S_data,
                                                   bdy_data_xlo, bdy_data_xhi,
                                                   bdy_data_ylo, bdy_data_yhi);
            }
        }
#endif
    }; // end slow_rhs_fun_pre
