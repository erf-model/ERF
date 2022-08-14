    auto fast_rhs_fun = [&](      Vector<MultiFab>& S_fast_rhs,
                                  Vector<MultiFab>& S_slow_rhs,
                                  Vector<MultiFab>& S_stage_data,
                            const Vector<MultiFab>& S_data,
                                  Vector<MultiFab>& S_scratch,
                            const Real fast_dt, const Real inv_fac)
    {
        if (verbose) amrex::Print() << "Calling fast rhs at level " << level << " with dt = " << fast_dt << std::endl;
        erf_fast_rhs(level, S_fast_rhs, S_slow_rhs, S_stage_data, S_prim,
                     S_data, S_scratch, advflux, fine_geom, ifr, solverChoice,
                     z_phys_nd[level], detJ_cc[level], r0, p0,
                     fast_dt, inv_fac);
    };


    auto post_substep_fun = [&](Vector<MultiFab>& S_data, const Real time_for_fp, int ng_cons, int ng_vel)
    {
        bool fast_only = true;
        apply_bcs(S_data, time_for_fp, ng_cons, ng_vel, fast_only);
    };
