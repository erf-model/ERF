#include "prob_common.H"

    // *************************************************************
    // The "slow" integrator for MRI and the only integrator for SRI
    // *************************************************************
    auto slow_rhs_fun_pre = [&](Vector<MultiFab>& S_rhs,
                                Vector<MultiFab>& S_data,
                                Vector<MultiFab>& S_scratch,
                                const Real old_time,
                                const Real old_stage_time,
                                const Real new_stage_time)
    {
        BL_PROFILE("slow_rhs_fun_pre");
        if (verbose) Print() << "Making slow rhs at time " << old_stage_time <<
                                " for fast variables at level " << level <<
                                " for advancing from " << old_time << " to " << new_stage_time << std::endl;

        // Moving terrain
        if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == 1) )
        {

            const auto dz = fine_geom.CellSizeArray()[2];

            if (verbose) Print() << "Re-making src geometry at old time  " << old_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd[level],old_time);       // Defines h(i,j,k=0,t) at time at time t^n
            init_terrain_grid  (fine_geom,*z_phys_nd[level]);                // Defines z_phys       at time at time t^n
            make_metrics       (fine_geom,*z_phys_nd[level], *z_phys_cc[level], *detJ_cc[level]);  // Defines detJ, etc

            if (verbose) Print() << "Making src geometry at time  " << old_stage_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd_src[level],old_stage_time); // Defines h(i,j,k=0,t) at time at which src is evaluated
            init_terrain_grid  (fine_geom,*z_phys_nd_src[level]);                // Defines z_phys       at time at which src is evaluated
            make_metrics       (fine_geom,*z_phys_nd_src[level], *z_phys_cc_src[level], *detJ_cc_src[level]);  // Defines detJ, etc

            if (verbose) Print() << "Making new geometry at time " << new_stage_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd_new[level],new_stage_time); // Defines h(i,j,k=0,t) at new time at end of stage
            init_terrain_grid  (fine_geom,*z_phys_nd_new[level]);                // Defines z_phys       at new time at end of stage
            make_metrics       (fine_geom,*z_phys_nd_new[level], *z_phys_cc_new[level], *detJ_cc_new[level]);  // Defines detJ, etc

            const auto& dX = fine_geom.CellSizeArray();

            for (MFIter mfi(*z_t_rk[level],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(IntVect(1,1,0));

                // This is the top face
                int khi = gbx.bigEnd()[2];

                // Now make this a 2D box (after we extract khi)
                gbx.setRange(2,0);

                const Array4<const Real>& detJ_new  =   detJ_cc_new[level]->const_array(mfi);
                const Array4<const Real>& detJ_old  =       detJ_cc[level]->const_array(mfi);

                const Array4<      Real>& z_t_arr   =  z_t_rk[level]->array(mfi);

                // Note that delta_t here must correspond to the times used when making z_nd and detJ
                Real delta_t  = (new_stage_time - old_time);
                fill_dhdt(z_t_arr,gbx,dX,old_time,delta_t);

                // Loop over horizontal plane
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
                {
                    if (i == 3 and j == 3) amrex::Print() << "MAKING HT AT RK STAGE " << z_t_arr(i,j,0) << std::endl;
                    if (i == 3 and j == 3) amrex::Print() << "  WITH OLD TIME       " << old_time << " TO " << new_stage_time << std::endl;
                    // Evaluate between RK stages
                    for (int k = 1; k <= khi; k++)
                    {
                        // NOTE: J in this code is (1/J) in Sullivan et al (JAS 2014)
                        Real J_new = detJ_new(i,j,k-1);
                        Real J_old = detJ_old(i,j,k-1);
                        z_t_arr(i,j,k) = z_t_arr(i,j,k-1) + ( (J_new - J_old) / delta_t ) * dz;
                    }
                });
            } // mfi

            erf_slow_rhs_pre(level, S_rhs, S_data, S_prim, S_scratch,
                             xvel_new, yvel_new, zvel_new, z_t_rk[level],
                             Omega, source, eddyDiffs, diffflux,
                             fine_geom, ifr, solverChoice, m_most, domain_bcs_type_d, domain_bcs_type,
                             z_phys_nd_src[level], detJ_cc_src[level], r0, p0,
                             dptr_rayleigh_tau, dptr_rayleigh_ubar,
                             dptr_rayleigh_vbar, dptr_rayleigh_thetabar,
                             ingested_bcs);
        } else {
            erf_slow_rhs_pre(level, S_rhs, S_data, S_prim, S_scratch,
                             xvel_new, yvel_new, zvel_new, z_t_rk[level],
                             Omega, source, eddyDiffs, diffflux,
                             fine_geom, ifr, solverChoice, m_most, domain_bcs_type_d, domain_bcs_type,
                             z_phys_nd[level], detJ_cc[level], r0, p0,
                             dptr_rayleigh_tau, dptr_rayleigh_ubar,
                             dptr_rayleigh_vbar, dptr_rayleigh_thetabar,
                             ingested_bcs);
        }
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
                               const Real time_for_fp, const Real delta_t,
                               int ng_cons, int ng_vel)
    {
        // Moving terrain
        Real time_mt = time_for_fp - 0.5*delta_t;

        bool fast_only = false;
        apply_bcs(S_data, time_for_fp, time_mt, delta_t, ng_cons, ng_vel, fast_only);
    };

    // *************************************************************
    // The "slow" integrator for MRI and the only integrator for SRI
    // *************************************************************
    auto slow_rhs_fun_post = [&](Vector<MultiFab>& S_rhs,
                                 Vector<MultiFab>& S_old,
                                 Vector<MultiFab>& S_new,
                                 Vector<MultiFab>& S_data,
                                 Vector<MultiFab>& S_scratch,
                                 const Real old_time,
                                 const Real old_stage_time,
                                 const Real new_stage_time)
    {
        if (verbose) Print() << "Making slow rhs at time " << old_stage_time <<
                                " for slow variables at level " << level <<
                                " for advancing from " << old_time << " to " << new_stage_time << std::endl;

        Real slow_dt = new_stage_time - old_time;

        // Moving terrain
        if ( solverChoice.use_terrain && (solverChoice.terrain_type == 1) )
        {
            erf_slow_rhs_post(level, slow_dt, S_rhs, S_old, S_new, S_data, S_prim, S_scratch,
                              xvel_new, yvel_new, zvel_new,
                              source, eddyDiffs, diffflux,
                              fine_geom, solverChoice, m_most, domain_bcs_type_d,
                              detJ_cc[level], detJ_cc_src[level], detJ_cc_new[level]);
        } else {
            erf_slow_rhs_post(level, slow_dt, S_rhs, S_old, S_new, S_data, S_prim, S_scratch,
                              xvel_new, yvel_new, zvel_new,
                              source, eddyDiffs, diffflux,
                              fine_geom, solverChoice, m_most, domain_bcs_type_d,
                              detJ_cc[level], detJ_cc[level], detJ_cc[level]);
        }
    };
