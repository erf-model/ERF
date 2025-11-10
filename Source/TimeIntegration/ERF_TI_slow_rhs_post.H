#include "ERF_SrcHeaders.H"

    auto slow_rhs_fun_post = [&](Vector<MultiFab>& S_rhs,
                                 Vector<MultiFab>& S_old,
                                 Vector<MultiFab>& S_new,
                                 Vector<MultiFab>& S_data,
                                 const Real old_step_time,
                                 const Real old_stage_time,
                                 const Real new_stage_time,
                                 const int nrk)
    {
        // Note that the "old" and "new" metric terms correspond to
        // t^n and the RK stage (either t^*, t^** or t^{n+1} that this source
        // will be used to advance to
        Real slow_dt = new_stage_time - old_step_time;

        if (verbose) amrex::Print() << "Time integration of scalars at level " << level
                                    << std::setprecision(timeprecision)
                                    << " from " << old_step_time << " to " << new_stage_time
                                    << " with dt = " << slow_dt
                                    << " using RHS created at " << old_stage_time << std::endl;

#if defined(ERF_USE_NETCDF)
        bool moist_set_rhs = false;
        if ( solverChoice.use_real_bcs && (level==0) && (real_set_width > 0) &&
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
        if (solverChoice.coupling_type == CouplingType::TwoWay && finest_level > 0)
        {
            if (level < finest_level) {
                fr_as_crse = getAdvFluxReg(level+1);
            }
            if (level > 0) {
                fr_as_fine = getAdvFluxReg(level);
            }
        }

        MultiFab* new_detJ =
           (solverChoice.terrain_type == TerrainType::MovingFittedMesh) ? (detJ_cc_new[level].get()) : (detJ_cc[level].get());

        erf_slow_rhs_post(level, finest_level, nrk, slow_dt, micro->Get_Qstate_Moist_Size(),
                          S_rhs, S_old, S_new, S_data, S_prim, avg_xmom[level], avg_ymom[level], avg_zmom[level],
                          xvel_new, yvel_new, zvel_new, cc_src, SmnSmn, eddyDiffs,
                          Hfx1, Hfx2, Hfx3, Q1fx1, Q1fx2, Q1fx3, Q2fx3, Diss,
                          fine_geom, solverChoice, m_SurfaceLayer, domain_bcs_type_d, domain_bcs_type,
                          z_phys_nd[level], z_phys_cc[level], ax[level], ay[level], az[level],
                          detJ_cc[level], new_detJ, stretched_dz_d[level], mapfac[level], EBFactory(level),
#if defined(ERF_USE_NETCDF)
                          moist_set_rhs, bdy_time_interval, new_stage_time, stop_time-start_time,
                          real_width, real_set_width,
                          bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
#endif
#ifdef ERF_USE_SHOC
                          shoc_interface[level],
#endif
                          fr_as_crse, fr_as_fine);

        // Apply state redistribution for cons states

        if (solverChoice.terrain_type == TerrainType::EB)
        {
            Vector<int> is_valid_slow_var; is_valid_slow_var.resize(RhoQ1_comp+1,0);
            if (solverChoice.turbChoice[level].use_tke) {is_valid_slow_var[    RhoKE_comp] = 1;}
                                                         is_valid_slow_var[RhoScalar_comp] = 1;
            if (solverChoice.moisture_type != MoistureType::None) {
                is_valid_slow_var[RhoQ1_comp] = 1;
            }
            const int num_comp_total = S_rhs[IntVars::cons].nComp();
            const int num_grow       = S_rhs[IntVars::cons].nGrow();
            const int nvars          = S_data[IntVars::cons].nComp();

            MultiFab dUdt_tmp(ba, dm, num_comp_total, num_grow, MFInfo(), EBFactory(level));
            dUdt_tmp.setVal(0.0, 0, num_comp_total, num_grow);

            int start_comp;
            int num_comp;

            for (int ivar(RhoKE_comp); ivar<= RhoQ1_comp; ++ivar)
            {
                if (is_valid_slow_var[ivar])
                {
                    start_comp = ivar;
                    num_comp   = 1;
                    if (ivar == RhoQ1_comp) {
                        num_comp = nvars - RhoQ1_comp;
                    } else if (ivar == RhoScalar_comp) {
                        num_comp = NSCALARS;
                    }
                    MultiFab::Copy(dUdt_tmp, S_rhs[IntVars::cons], start_comp, start_comp, num_comp, 0);
                }
            }
            dUdt_tmp.FillBoundary(fine_geom.periodicity());
            dUdt_tmp.setDomainBndry(1.234e10, 0, num_comp_total, fine_geom);

            const BCRec* bc_ptr_d = domain_bcs_type_d.data();

            // Update S_rhs by Redistribution.
            // To-do: Currently, redistributing all the scalar variables.
            //        This needs to be redistributed only for num_comp variables starting from ivar, for efficiency.
            redistribute_term ( num_comp_total, fine_geom, S_rhs[IntVars::cons], dUdt_tmp,
                                S_old[IntVars::cons], EBFactory(level), bc_ptr_d, slow_dt);

            // Update state using the updated S_rhs. (NOTE: redistribute_term returns RHS not state variables.)
            for ( MFIter mfi(S_new[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box tbx  = mfi.tilebox();
                const Array4<Real>& snew = S_new[IntVars::cons].array(mfi);
                const Array4<Real>& sold = S_old[IntVars::cons].array(mfi);
                const Array4<Real>& srhs = S_rhs[IntVars::cons].array(mfi);
                Array4<const Real> detJ_arr = EBFactory(level).getVolFrac().const_array(mfi);

                for (int ivar(RhoKE_comp); ivar<= RhoQ1_comp; ++ivar)
                {
                    if (is_valid_slow_var[ivar])
                    {
                        start_comp = ivar;
                        num_comp   = 1;
                        if (ivar == RhoQ1_comp) {
                            num_comp = nvars - RhoQ1_comp;
                        } else if (ivar == RhoScalar_comp) {
                            num_comp = NSCALARS;
                        }
                        ParallelFor(tbx, num_comp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn)
                        {
                            if (detJ_arr(i,j,k) > 0.0) {
                                const int n = start_comp + nn;
                                snew(i,j,k,n) = sold(i,j,k,n) + slow_dt * srhs(i,j,k,n);
                            }
                        });
                    }
                }
            }
        } // EB

        // Apply boundary conditions on all the state variables that have been updated
        // in both the fast and slow integrators
        apply_bcs(S_new, new_stage_time, S_new[IntVars::cons].nGrow(), S_new[IntVars::xmom].nGrow(),
                  fast_only=false, vel_and_mom_synced=false);

        if (solverChoice.moisture_tight_coupling) {
            // TODO: need iteration var for lagrangian microphysics
            // call signature in ERF::Advance() is
            //advance_microphysics(lev, S_new,    dt_lev,  iteration, time);
            advance_microphysics(level, S_new[0], slow_dt, 123456789, old_step_time);
        }
    }; // end slow_rhs_fun_post
