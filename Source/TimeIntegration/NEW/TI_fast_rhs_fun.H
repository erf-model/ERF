    auto fast_rhs_fun = [&](int fast_step,
                            Vector<MultiFab>& S_slow_rhs,
                            const Vector<MultiFab>& S_old,
                            Vector<MultiFab>& S_stage,
                            Vector<MultiFab>& S_data,
                            Vector<MultiFab>& S_scratch,
                            const Real dtau,
                            const Real inv_fac,
                            const Real prev_time,
                            const Real stage_time,
                            const Real new_time)
    {
        BL_PROFILE("fast_rhs_fun");
        if (verbose) amrex::Print() << "Calling fast rhs at level " << level << " with dt = " << dtau << std::endl;

        // Moving terrain
        MultiFab* z_t_pert = nullptr;
        if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == 1) )
        {
            const auto dz = fine_geom.CellSizeArray()[2];

            // Make "old" fast geom -- store in z_phys_nd for convenience
            if (verbose) Print() << "Making geometry at start of substep time: " << prev_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd[level],prev_time);      // Defines h(i,j,k=0,t) at time at previous "fast" time
            init_terrain_grid  (fine_geom,*z_phys_nd[level]);                // Defines z_phys       at time at previous "fast" time
            make_metrics       (fine_geom,*z_phys_nd[level], *z_phys_cc[level], *detJ_cc[level]);  // Defines detJ, etc

            // Make geom at time that RK stage is defined
            if (verbose) Print() << "Making geometry at RK stage time        : " << stage_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd_src[level],stage_time);          // Defines h(i,j,k=0,t) at time at RK stage time
            init_terrain_grid  (fine_geom,*z_phys_nd_src[level]);                     // Defines z_phys       at time at RK stage time
            make_metrics       (fine_geom,*z_phys_nd_src[level], *z_phys_cc_src[level], *detJ_cc_src[level]);  // Defines detJ, etc

            // Make "new" fast geom
            if (verbose) Print() << "Making geometry for end of substep time :" << new_time << std::endl;
            init_custom_terrain(fine_geom,*z_phys_nd_new[level],new_time); // Defines h(i,j,k=0,t)
            init_terrain_grid  (fine_geom,*z_phys_nd_new[level]);            // Defines z_phys for all i,j,k
            make_metrics       (fine_geom,*z_phys_nd_new[level], *z_phys_cc_new[level], *detJ_cc_new[level]);  // Defines detJ

            z_t_pert = new MultiFab(S_data[IntVar::zmom].boxArray(), S_data[IntVar::zmom].DistributionMap(), 1, 1);

            for (MFIter mfi(*z_t_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(IntVect(1,1,0));

                // This is the top face
                int khi = gbx.bigEnd()[2];

                // Now make this a 2D box (after we extract khi)
                gbx.setRange(2,0);

                const Array4<const Real>& detJ_old_arr  =  detJ_cc[level]->const_array(mfi);
                const Array4<const Real>& detJ_new_arr  =  detJ_cc_new[level]->const_array(mfi);
                const Array4<      Real>& zpert_arr     =  z_t_pert->array(mfi);
                const Array4<      Real>& z_t_rk_arr    =  z_t_rk[level]->array(mfi);

                Real delta_t  = dtau;

                const auto& dX = fine_geom.CellSizeArray();

                amrex::Print() << " MAKING DHDT WITH TIME DT " << prev_time << " " << delta_t << std::endl;
                fill_dhdt( zpert_arr,gbx,dX,prev_time,delta_t);

                // Loop over horizontal plane
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
                {
                    // Evaluate from previous fast step to new fast step
                    for (int k = 1; k <= khi; k++)
                    {
                        // NOTE: J in this code is (1/J) in Sullivan et al (JAS 2014)
                        Real J_new = detJ_new_arr(i,j,k-1);
                        Real J_old = detJ_old_arr(i,j,k-1);
                        zpert_arr(i,j,k) = zpert_arr(i,j,k-1) + ( (J_new - J_old) / delta_t ) * dz;

                        // Convert to perturbation: z"_t(t) = z_t(t) - z_t^{RK}
                        zpert_arr(i,j,k) -= z_t_rk_arr(i,j,k);
                    }
                });
            } // mfi
        }


        if (solverChoice.use_terrain && solverChoice.terrain_type == 1) { // moving terrain

            // With moving terrain we remake this at ever time step since the "new" terrain changes
            make_fast_coeffs(level, fast_coeffs, S_stage, S_prim, pi_stage, fine_geom, ifr, solverChoice,
                             z_phys_nd_new[level], detJ_cc_new[level], r0, pi0, dtau, ingested_bcs);

            if (fast_step == 0) {

                // If this is the first substep we pass in S_old as the previous step's solution
                erf_fast_rhs_MT(fast_step, level, new_time, S_slow_rhs, S_old, S_stage, S_prim, pi_stage, fast_coeffs,
                                S_data, S_scratch, fine_geom, ifr, solverChoice, Omega, z_t_pert,
                                z_phys_nd[level], z_phys_nd_new[level],
                                detJ_cc[level], detJ_cc_new[level], dtau, inv_fac, ingested_bcs);
            } else {

                // If this is not the first substep we pass in S_data as the previous step's solution
                erf_fast_rhs_MT(fast_step, level, new_time, S_slow_rhs, S_data, S_stage, S_prim, pi_stage, fast_coeffs,
                                S_data, S_scratch, fine_geom, ifr, solverChoice, Omega, z_t_pert,
                                z_phys_nd[level], z_phys_nd_new[level],
                                detJ_cc[level], detJ_cc_new[level], dtau, inv_fac, ingested_bcs);
            }
        } else if (solverChoice.use_terrain && solverChoice.terrain_type == 0) { // fixed terrain
            if (fast_step == 0) {

                // If this is the first substep we make the coefficients since they are based only on stage data
                make_fast_coeffs(level, fast_coeffs, S_stage, S_prim, pi_stage, fine_geom, ifr, solverChoice,
                                 z_phys_nd[level], detJ_cc[level], r0, pi0, dtau, ingested_bcs);

                // If this is the first substep we pass in S_old as the previous step's solution
                erf_fast_rhs_T(fast_step, level, new_time, S_slow_rhs, S_old, S_stage, S_prim, pi_stage, fast_coeffs,
                               S_data, S_scratch, fine_geom, ifr, solverChoice, Omega,
                               z_phys_nd[level], detJ_cc[level], dtau, inv_fac, ingested_bcs);
            } else {
                // If this is not the first substep we pass in S_data as the previous step's solution
                erf_fast_rhs_T(fast_step, level, new_time, S_slow_rhs, S_data, S_stage, S_prim, pi_stage, fast_coeffs,
                               S_data, S_scratch, fine_geom, ifr, solverChoice, Omega,
                               z_phys_nd[level], detJ_cc[level], dtau, inv_fac, ingested_bcs);
            }
        } else { // no terrain
            if (fast_step == 0) {

                // If this is the first substep we make the coefficients since they are based only on stage data
                make_fast_coeffs(level, fast_coeffs, S_stage, S_prim, pi_stage, fine_geom, ifr, solverChoice,
                                 z_phys_nd[level], detJ_cc[level], r0, pi0, dtau, ingested_bcs);

                // If this is the first substep we pass in S_old as the previous step's solution
                erf_fast_rhs_N(fast_step, level, new_time, S_slow_rhs, S_old, S_stage, S_prim, pi_stage, fast_coeffs,
                               S_data, S_scratch, fine_geom, ifr, solverChoice, dtau, inv_fac, ingested_bcs);
            } else {
                // If this is not the first substep we pass in S_data as the previous step's solution
                erf_fast_rhs_N(fast_step, level, new_time, S_slow_rhs, S_data, S_stage, S_prim, pi_stage, fast_coeffs,
                               S_data, S_scratch, fine_geom, ifr, solverChoice, dtau, inv_fac, ingested_bcs);
            }
        }

        // Moving terrain
        Real time_mt = new_time - 0.5*dtau;

        // Even if we update all the conserved variables we don't need to fillpatch the slow ones every acoustic substep
        int ng_cons    = 1;
        int ng_vel     = 1;
        apply_bcs(S_data, new_time, time_mt, dtau, ng_cons, ng_vel, true);

        MultiFab::Copy(S_scratch[IntVar::cons], S_data[IntVar::cons], Cons::RhoTheta, Cons::RhoTheta, 1, ng_cons);
    };
