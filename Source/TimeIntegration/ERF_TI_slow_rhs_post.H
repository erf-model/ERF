#include "ERF_SrcHeaders.H"

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
                          S_rhs, S_old, S_new, S_data, S_prim, S_scratch,
                          xvel_new, yvel_new, zvel_new, cc_src, SmnSmn, eddyDiffs,
                          Hfx1, Hfx2, Hfx3, Q1fx1, Q1fx2, Q1fx3, Q2fx3, Diss,
                          fine_geom, solverChoice, m_SurfaceLayer, domain_bcs_type_d, domain_bcs_type,
                          z_phys_nd[level], z_phys_cc[level], ax[level], ay[level], az[level],
                          detJ_cc[level], new_detJ, stretched_dz_d[level], mapfac[level], EBFactory(level),
#if defined(ERF_USE_NETCDF)
                          moist_set_rhs, bdy_time_interval, start_bdy_time, new_stage_time,
                          real_width, real_set_width,
                          bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
#endif
                          fr_as_crse, fr_as_fine);

        // Apply boundary conditions on all the state variables that have been updated
        // in both the fast and slow integrators
        apply_bcs(S_new, new_stage_time, S_new[IntVars::cons].nGrow(), S_new[IntVars::xmom].nGrow(),
                  fast_only=false, vel_and_mom_synced=false);
    }; // end slow_rhs_fun_post
