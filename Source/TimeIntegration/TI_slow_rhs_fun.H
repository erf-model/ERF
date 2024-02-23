#include "prob_common.H"

/**
 *  Wrapper for calling the routine that creates the slow RHS
 */
    auto slow_rhs_fun_pre = [&](Vector<MultiFab>& S_rhs,
                                Vector<MultiFab>& S_old,
                                Vector<MultiFab>& S_data,
                                Vector<MultiFab>& S_scratch,
                                const Real old_step_time,
                                const Real old_stage_time,
                                const Real new_stage_time,
                                const int nrk)
    {
        BL_PROFILE("slow_rhs_fun_pre");
        if (verbose) Print() << "Making slow rhs at time " << old_stage_time << " for fast variables advancing from " <<
                                old_step_time << " to " << new_stage_time << std::endl;

        Real slow_dt = new_stage_time - old_step_time;

        // *************************************************************************
        // Set up flux registers if using two_way coupling
        // *************************************************************************
        YAFluxRegister* fr_as_crse = nullptr;
        YAFluxRegister* fr_as_fine = nullptr;
        if (solverChoice.coupling_type == CouplingType::TwoWay) {
            if (level < finest_level) {
                fr_as_crse = getAdvFluxReg(level+1);
                fr_as_crse->reset();
            }
            if (level > 0) {
                fr_as_fine = getAdvFluxReg(level);
            }
        }

        // Moving terrain
        if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == TerrainType::Moving) )
        {
            // Note that the "old" and "new" metric terms correspond to
            // t^n and the RK stage (either t^*, t^** or t^{n+1} that this source
            // will be used to advance to

            // The "src" metric terms correspond to the time at which we are evaluating the source here,
            // aka old_stage_time

            if (verbose) Print() << "Re-making old geometry at old time   : " << old_step_time << std::endl;
            prob->init_custom_terrain(fine_geom,*z_phys_nd[level],old_step_time);
            init_terrain_grid  (fine_geom,*z_phys_nd[level], zlevels_stag);
            make_J             (fine_geom,*z_phys_nd[level], *detJ_cc[level]);

            if (verbose) Print() << "Making src geometry at old_stage_time:  " << old_stage_time << std::endl;
            prob->init_custom_terrain(fine_geom,*z_phys_nd_src[level],old_stage_time);
            init_terrain_grid  (fine_geom,*z_phys_nd_src[level], zlevels_stag);
            make_J             (fine_geom,*z_phys_nd_src[level], *detJ_cc_src[level]);

            if (verbose) Print() << "Making new geometry at new_stage_time: " << new_stage_time << std::endl;
            prob->init_custom_terrain(fine_geom,*z_phys_nd_new[level],new_stage_time);
            init_terrain_grid  (fine_geom,*z_phys_nd_new[level], zlevels_stag);
            make_J             (fine_geom,*z_phys_nd_new[level], *detJ_cc_new[level]);

            Real inv_dt  = 1./slow_dt;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*z_t_rk[level],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(IntVect(1,1,0));

                const Array4<Real      >& z_t_arr      =  z_t_rk[level]->array(mfi);
                const Array4<Real const>& z_nd_new_arr =  z_phys_nd_new[level]->const_array(mfi);
                const Array4<Real const>& z_nd_old_arr =  z_phys_nd[level]->const_array(mfi);

                // Loop over horizontal plane
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // Evaluate between RK stages assuming the geometry is linear between old and new time
                    z_t_arr(i,j,k) = 0.25 * inv_dt * (z_nd_new_arr(i+1,j+1,k) - z_nd_old_arr(i+1,j+1,k)
                                                     +z_nd_new_arr(i  ,j+1,k) - z_nd_old_arr(  i,j+1,k)
                                                     +z_nd_new_arr(i+1,j  ,k) - z_nd_old_arr(i+1,j  ,k)
                                                     +z_nd_new_arr(i  ,j  ,k) - z_nd_old_arr(i  ,j  ,k));
                });

            } // mfi

            MultiFab r_hse_new (base_state_new[level], make_alias, 0, 1); // r_0  is first  component
            MultiFab p_hse_new (base_state_new[level], make_alias, 1, 1); // p_0  is second component
            MultiFab pi_hse_new(base_state_new[level], make_alias, 2, 1); // pi_0 is second component

            MultiFab* r0_new = &r_hse_new;
            MultiFab* p0_new = &p_hse_new;

            make_buoyancy(S_data, S_prim, buoyancy, fine_geom,
                          solverChoice, r0_new, micro.Get_Qstate_Size());

            erf_slow_rhs_pre(level, finest_level, nrk, slow_dt, S_rhs, S_data, S_prim, S_scratch,
                             xvel_new, yvel_new, zvel_new,
                             z_t_rk[level], Omega, source, buoyancy,
#if defined(ERF_USE_RRTMGP)
                             qheating_rates[level],
#endif
                             Tau11, Tau22, Tau33, Tau12,
                             Tau13, Tau21,  Tau23, Tau31, Tau32, SmnSmn, eddyDiffs,
                             Hfx3, Diss,
                             fine_geom, solverChoice, m_most, domain_bcs_type_d, domain_bcs_type,
                             z_phys_nd_src[level], detJ_cc_src[level], p0_new,
                             mapfac_m[level], mapfac_u[level], mapfac_v[level],
                             fr_as_crse, fr_as_fine,
                             dptr_rhotheta_src, d_rayleigh_ptrs_at_lev);

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

            const GpuArray<Real, AMREX_SPACEDIM> dxInv = fine_geom.InvCellSizeArray();

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

               const Array4<Real const>  p0_arr     = p0->const_array(mfi);
               const Array4<Real      >  p0_new_arr = p0_new->array(mfi);

               const Array4<Real      >& z_t_arr      =  z_t_rk[level]->array(mfi);

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
               });
            } // MFIter
            r0_new->FillBoundary(fine_geom.periodicity());
            p0_new->FillBoundary(fine_geom.periodicity());

        } else {

            // If not moving_terrain
            make_buoyancy(S_data, S_prim, buoyancy,
                          fine_geom, solverChoice, r0, micro.Get_Qstate_Size());

            erf_slow_rhs_pre(level, finest_level, nrk, slow_dt, S_rhs, S_data, S_prim, S_scratch,
                             xvel_new, yvel_new, zvel_new,
                             z_t_rk[level], Omega, source, buoyancy,
#if defined(ERF_USE_RRTMGP)
                             qheating_rates[level],
#endif
                             Tau11, Tau22, Tau33, Tau12,
                             Tau13, Tau21,  Tau23, Tau31, Tau32, SmnSmn, eddyDiffs,
                             Hfx3, Diss,
                             fine_geom, solverChoice, m_most, domain_bcs_type_d, domain_bcs_type,
                             z_phys_nd[level], detJ_cc[level], p0,
                             mapfac_m[level], mapfac_u[level], mapfac_v[level],
                             fr_as_crse, fr_as_fine,
                             dptr_rhotheta_src, d_rayleigh_ptrs_at_lev);
        }

#ifdef ERF_USE_NETCDF
        // Populate RHS for relaxation zones if using real bcs
        if (use_real_bcs && (level == 0)) {
            if (real_width>0) {
                    realbdy_compute_interior_ghost_rhs(init_type,
                                                       bdy_time_interval, start_bdy_time, new_stage_time, slow_dt,
                                                       real_width, real_set_width, fine_geom,
                                                       S_rhs, S_old, S_data,
                                                       bdy_data_xlo, bdy_data_xhi,
                                                       bdy_data_ylo, bdy_data_yhi);
            }
        }
#endif

        // Compute RHS for fine interior ghost
        if (level > 0 && cf_width > 0) {
            fine_compute_interior_ghost_rhs(new_stage_time, slow_dt,
                                            cf_width, cf_set_width, fine_geom,
                                            &FPr_c[level-1], &FPr_u[level-1], &FPr_v[level-1], &FPr_w[level-1],
                                            domain_bcs_type, S_rhs, S_data);
        }

        // S_rhs[IntVar::cons].FillBoundary(fine_geom.periodicity());
        // S_rhs[IntVar::xmom].FillBoundary(fine_geom.periodicity());
        // S_rhs[IntVar::ymom].FillBoundary(fine_geom.periodicity());
        // S_rhs[IntVar::zmom].FillBoundary(fine_geom.periodicity());
    }; // end slow_rhs_fun_pre

    // *************************************************************
    // This called before RK stage
    // *************************************************************
    auto pre_update_fun = [&](Vector<MultiFab>& S_data, int ng_cons)
    {
        cons_to_prim(S_data[IntVar::cons], ng_cons);
    };

    // *************************************************************
    // This called after every RK stage -- from MRI or SRI
    // *************************************************************
    auto post_update_fun = [&](Vector<MultiFab>& S_data,
                               const Real time_for_fp, int ng_cons, int ng_vel)
    {
        bool fast_only = false;
        bool vel_and_mom_synced = false;
        apply_bcs(S_data, time_for_fp, ng_cons, ng_vel, fast_only, vel_and_mom_synced);
    };

    // *************************************************************
    // The "slow" integrator for MRI and the only integrator for SRI
    // *************************************************************
    auto slow_rhs_fun_post = [&](Vector<MultiFab>& S_rhs,
                                 Vector<MultiFab>& S_old,
                                 Vector<MultiFab>& S_new,
                                 Vector<MultiFab>& S_data,
                                 Vector<MultiFab>& S_scratch,
                                 const Real old_step_time,
                                 const Real old_stage_time,
                                 const Real new_stage_time,
                                 const int nrk)
    {
        amrex::ignore_unused(nrk);
        if (verbose) Print() << "Making slow rhs at time " << old_stage_time <<
                                " for slow variables advancing from " <<
                                old_step_time << " to " << new_stage_time << std::endl;

        // Note that the "old" and "new" metric terms correspond to
        // t^n and the RK stage (either t^*, t^** or t^{n+1} that this source
        // will be used to advance to
        Real slow_dt = new_stage_time - old_step_time;

#if defined(ERF_USE_NETCDF)
        bool moist_set_rhs = false;
        if ( use_real_bcs &&
             (level==0)   &&
             (real_set_width > 0) &&
             (solverChoice.moisture_type != MoistureType::None) )
        {
            moist_set_rhs = true;
        }
#endif

        // *************************************************************************
        // Set up flux registers if using two_way coupling
        // *************************************************************************
        YAFluxRegister* fr_as_crse = nullptr;
        YAFluxRegister* fr_as_fine = nullptr;
        if (solverChoice.coupling_type == CouplingType::TwoWay)
        {
            if (level < finest_level) {
                fr_as_crse = getAdvFluxReg(level+1);
            }
            if (level > 0) {
                fr_as_fine = getAdvFluxReg(level);
            }
        }

        // Moving terrain
        if ( solverChoice.use_terrain && (solverChoice.terrain_type == TerrainType::Moving) ) {
            erf_slow_rhs_post(level, finest_level, nrk, slow_dt,
                              S_rhs, S_old, S_new, S_data, S_prim, S_scratch,
                              xvel_new, yvel_new, zvel_new,
                              source, SmnSmn, eddyDiffs,
                              Hfx3, Diss,
                              fine_geom, solverChoice, m_most, domain_bcs_type_d,
                              z_phys_nd_src[level], detJ_cc[level], detJ_cc_new[level],
                              mapfac_m[level], mapfac_u[level], mapfac_v[level],
#if defined(ERF_USE_NETCDF)
                              moist_set_rhs, bdy_time_interval, start_bdy_time, new_stage_time,
                              real_width, real_set_width,
                              bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
#endif
                              fr_as_crse, fr_as_fine
                              );
        } else {
            erf_slow_rhs_post(level, finest_level, nrk, slow_dt,
                              S_rhs, S_old, S_new, S_data, S_prim, S_scratch,
                              xvel_new, yvel_new, zvel_new,
                              source, SmnSmn, eddyDiffs,
                              Hfx3, Diss,
                              fine_geom, solverChoice, m_most, domain_bcs_type_d,
                              z_phys_nd[level], detJ_cc[level], detJ_cc[level],
                              mapfac_m[level], mapfac_u[level], mapfac_v[level],
#if defined(ERF_USE_NETCDF)
                              moist_set_rhs, bdy_time_interval, start_bdy_time, new_stage_time,
                              real_width, real_set_width,
                              bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
#endif
                              fr_as_crse, fr_as_fine
                              );
        }
    }; // end slow_rhs_fun_post

#ifdef ERF_USE_POISSON_SOLVE
    auto slow_rhs_fun_inc = [&](Vector<MultiFab>& S_rhs,
                                Vector<MultiFab>& S_old,
                                Vector<MultiFab>& S_data,
                                Vector<MultiFab>& S_scratch,
                                const Real old_step_time,
                                const Real old_stage_time,
                                const Real new_stage_time,
                                const int nrk)
    {
        BL_PROFILE("slow_rhs_fun_inc");
        if (verbose) Print() << "Making slow rhs at time " << old_stage_time << " for fast variables advancing from " <<
                                old_step_time << " to " << new_stage_time << std::endl;

        Real slow_dt = new_stage_time - old_step_time;

        // If not moving_terrain
        make_buoyancy(S_data, S_prim, buoyancy,
                      fine_geom, solverChoice, r0, micro.Get_Qstate_Size());

        erf_slow_rhs_inc(level, nrk, slow_dt,
                         S_rhs, S_old, S_data, S_prim, S_scratch,
                         xvel_new, yvel_new, zvel_new,
                         z_t_rk[level], Omega, source, buoyancy, Tau11, Tau22, Tau33, Tau12,
                         Tau13, Tau21,  Tau23, Tau31, Tau32, SmnSmn, eddyDiffs,
                         Hfx3, Diss,
                         fine_geom, solverChoice, m_most, domain_bcs_type_d, domain_bcs_type,
                         z_phys_nd[level], detJ_cc[level], p0,
                         mapfac_m[level], mapfac_u[level], mapfac_v[level],
                         dptr_rhotheta_src, d_rayleigh_ptrs_at_lev);


#ifdef ERF_USE_NETCDF
        // Populate RHS for relaxation zones if using real bcs
        if (use_real_bcs && (level == 0)) {
            if (real_width>0) {
                    realbdy_compute_interior_ghost_rhs(init_type,
                                                       bdy_time_interval, start_bdy_time, new_stage_time, slow_dt,
                                                       real_width, real_set_width, fine_geom,
                                                       S_rhs, S_old, S_data,
                                                       bdy_data_xlo, bdy_data_xhi,
                                                       bdy_data_ylo, bdy_data_yhi);
            }
        }
#endif

        // Compute RHS for fine interior ghost
        if (level > 0 && cf_width>0) {
            fine_compute_interior_ghost_rhs(new_stage_time, slow_dt,
                                            cf_width, cf_set_width,
                                            &FPr_c[level-1], &FPr_u[level-1], &FPr_v[level-1], &FPr_w[level-1],
                                            boxes_at_level[level], domain_bcs_type,
                                            S_rhs, S_data);
        }
    }; // end slow_rhs_fun_inc
#endif
