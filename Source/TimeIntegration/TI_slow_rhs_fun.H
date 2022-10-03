#include "prob_common.H"

    // *************************************************************
    // The "slow" integrator for MRI and the only integrator for SRI
    // *************************************************************
    auto slow_rhs_fun_pre = [&](Vector<MultiFab>& S_rhs,
                                Vector<MultiFab>& S_data,
                                Vector<MultiFab>& S_scratch,
                                const Real old_step_time,
                                const Real old_stage_time,
                                const Real new_stage_time)
    {
        BL_PROFILE("slow_rhs_fun_pre");
        if (verbose) Print() << "Making slow rhs at time " << old_stage_time << " for fast variables advancing from " <<
                                old_step_time << " to " << new_stage_time << std::endl;

        // Moving terrain
        if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == 1) )
        {
            // Note that the "old" and "new" metric terms correspond to
            // t^n and the RK stage (either t^*, t^** or t^{n+1} that this source
            // will be used to advance to

            // The "src" metric terms correspond to the time at which we are evaluating the source here,
            // aka old_stage_time

            if (verbose) Print() << "Re-making old geometry at old time   : " << old_step_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd[level],old_step_time);
            init_terrain_grid  (fine_geom,*z_phys_nd[level]);
            make_metrics       (fine_geom,*z_phys_nd[level], *z_phys_cc[level], *detJ_cc[level]);

            if (verbose) Print() << "Making src geometry at old_stage_time:  " << old_stage_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd_src[level],old_stage_time);
            init_terrain_grid  (fine_geom,*z_phys_nd_src[level]);
            make_metrics       (fine_geom,*z_phys_nd_src[level], *z_phys_cc_src[level], *detJ_cc_src[level]);

            if (verbose) Print() << "Making new geometry at new_stage_time: " << new_stage_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd_new[level],new_stage_time);
            init_terrain_grid  (fine_geom,*z_phys_nd_new[level]);
            make_metrics       (fine_geom,*z_phys_nd_new[level], *z_phys_cc_new[level], *detJ_cc_new[level]);

            Real inv_dt  = 1./ (new_stage_time - old_step_time);

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

            erf_slow_rhs_pre(level, S_rhs, S_data, S_prim, S_scratch,
                             xvel_new, yvel_new, zvel_new, z_t_rk[level],
                             Omega, source, eddyDiffs,
                             fine_geom, ifr, solverChoice, m_most, domain_bcs_type_d, domain_bcs_type,
                             z_phys_nd_src[level], detJ_cc_src[level], r0, p0,
                             dptr_rayleigh_tau, dptr_rayleigh_ubar,
                             dptr_rayleigh_vbar, dptr_rayleigh_thetabar,
                             ingested_bcs);
        } else {
            erf_slow_rhs_pre(level, S_rhs, S_data, S_prim, S_scratch,
                             xvel_new, yvel_new, zvel_new, z_t_rk[level],
                             Omega, source, eddyDiffs,
                             fine_geom, ifr, solverChoice, m_most, domain_bcs_type_d, domain_bcs_type,
                             z_phys_nd[level], detJ_cc[level], r0, p0,
                             dptr_rayleigh_tau, dptr_rayleigh_ubar,
                             dptr_rayleigh_vbar, dptr_rayleigh_thetabar,
                             ingested_bcs);
        }

        // S_rhs[IntVar::cons].FillBoundary(fine_geom.periodicity());
        // S_rhs[IntVar::xmom].FillBoundary(fine_geom.periodicity());
        // S_rhs[IntVar::ymom].FillBoundary(fine_geom.periodicity());
        // S_rhs[IntVar::zmom].FillBoundary(fine_geom.periodicity());
    };

    // *************************************************************
    // This called after before RK stage
    // *************************************************************
    auto pre_update_fun = [&](Vector<MultiFab>& S_data, int ng_cons)
    {
        cons_to_prim(S_data[IntVar::cons], ng_cons);
    };

    // *************************************************************
    // This called after every RK stage -- from MRI or SRI
    // *************************************************************
    auto post_update_fun = [&](Vector<MultiFab>& S_data,
                               const Real time_for_fp,
                               int ng_cons, int ng_vel)
    {
        bool fast_only = false;
        apply_bcs(S_data, time_for_fp, ng_cons, ng_vel, fast_only);
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
                                 const Real new_stage_time)
    {
        if (verbose) Print() << "Making slow rhs at time " << old_stage_time <<
                                " for slow variables advancing from " <<
                                old_step_time << " to " << new_stage_time << std::endl;

        // Note that the "old" and "new" metric terms correspond to
        // t^n and the RK stage (either t^*, t^** or t^{n+1} that this source
        // will be used to advance to
        Real slow_dt = new_stage_time - old_step_time;

        // Moving terrain
        if ( solverChoice.use_terrain && (solverChoice.terrain_type == 1) )
        {
            erf_slow_rhs_post(level, slow_dt, S_rhs, S_old, S_new, S_data, S_prim, S_scratch,
                              xvel_new, yvel_new, zvel_new,
                              source, eddyDiffs,
                              fine_geom, solverChoice, m_most, domain_bcs_type_d,
                              detJ_cc[level], detJ_cc_new[level]);
        } else {
            erf_slow_rhs_post(level, slow_dt, S_rhs, S_old, S_new, S_data, S_prim, S_scratch,
                              xvel_new, yvel_new, zvel_new,
                              source, eddyDiffs,
                              fine_geom, solverChoice, m_most, domain_bcs_type_d,
                              detJ_cc[level], detJ_cc[level]);
        }
    };
