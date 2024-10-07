/**
 *  Define the primitive variables by dividing the conserved variables by density
 */
    auto cons_to_prim = [&](const MultiFab& cons_state, int ng)
    {
        BL_PROFILE("cons_to_prim()");

        int ncomp_prim = S_prim.nComp();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(cons_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          const Box& gbx = mfi.growntilebox(ng);
          const Array4<const Real>& cons_arr     = cons_state.array(mfi);
          const Array4<      Real>& prim_arr     = S_prim.array(mfi);
          const Array4<      Real>& pi_stage_arr = pi_stage.array(mfi);
          const Real rdOcp = solverChoice.rdOcp;

          amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real rho       = cons_arr(i,j,k,Rho_comp);
            Real rho_theta = cons_arr(i,j,k,RhoTheta_comp);
            prim_arr(i,j,k,PrimTheta_comp) = rho_theta / rho;
            pi_stage_arr(i,j,k) = getExnergivenRTh(rho_theta, rdOcp);
            for (int n = 1; n < ncomp_prim; ++n) {
              prim_arr(i,j,k,PrimTheta_comp + n) = cons_arr(i,j,k,RhoTheta_comp + n) / rho;
            }
          });
      } // mfi
    };

/**
 *  This routine is called before the first step of the time integration, *and* in the case
 *  of a multi-stage method like RK3, this is called from "pre_update_fun" which is called
 *  before every subsequent stage.  Since we advance the variables in conservative form,
 *  we must convert momentum to velocity before imposing the bcs.
 */
    auto apply_bcs = [&](Vector<MultiFab>& S_data,
                         const Real time_for_fp, int ng_cons, int ng_vel,
                         bool fast_only, bool vel_and_mom_synced)
    {
        BL_PROFILE("apply_bcs()");

        int scomp_cons;
        int ncomp_cons;
        bool cons_only;

        int ng_cons_to_use;

        // **********************************************************************************
        // Because momentum is updated in the time-stepping routine, but boundary conditions
        //  are imposed on velocities, we must first update the velocity from the momentum
        //  before applying bcs.
        // **********************************************************************************
        if (!vel_and_mom_synced) {

            // **********************************************************************************
            // Call FillPatch routines for the density only because we need it to convert between
            //      momentum and velocity
            // This fills ghost cells/faces from
            //     1) coarser level if lev > 0
            //     2) physical boundaries
            //     3) other grids at the same level
            // **********************************************************************************

            // We must have at least one ghost cell of density to convert from momentum to velocity
            //    on the valid region
            AMREX_ALWAYS_ASSERT (ng_cons >= 1);

            // We must have at least one extra ghost cell of density to convert from velocity to momentum
            //    on the valid region
            ng_cons_to_use = std::max(ng_cons, ng_vel+1);

            scomp_cons = 0;
            ncomp_cons = 1;
            cons_only  = true;
            // **********************************************************************************
            // Fill ghost cells of density only
            // **********************************************************************************
            FillIntermediatePatch(level, time_for_fp,
                                  {&S_data[IntVars::cons], &xvel_new, &yvel_new, &zvel_new},
                                  {&S_data[IntVars::cons], &S_data[IntVars::xmom],
                                   &S_data[IntVars::ymom], &S_data[IntVars::zmom]},
                                  ng_cons_to_use, 0, cons_only, scomp_cons, ncomp_cons);
        }

        // ***************************************************************************************
        // Call FillPatch routines for all data except rho which was filled above
        // This fills ghost cells/faces from
        //     1) coarser level if lev > 0
        //     2) physical boundaries
        //     3) other grids at the same level
        // ***************************************************************************************
        if (vel_and_mom_synced) {
            if (fast_only) {
                scomp_cons = 0;
                ncomp_cons = 2; // rho and (rho theta) only
            } else {
                scomp_cons = 0;
                ncomp_cons = S_data[IntVars::cons].nComp();
            }
            // We must have at least one extra ghost cell of density to convert from velocity to momentum
            //    on the valid region
            ng_cons_to_use = std::max(ng_cons, ng_vel+1);

        } else {
            if (fast_only) {
                scomp_cons = 1;
                ncomp_cons = 1; // (rho theta) only since we filled rho above
            } else {
                scomp_cons = 1;
                ncomp_cons = S_data[IntVars::cons].nComp()-1; // since we filled rho above
            }
            ng_cons_to_use = ng_cons;
        }

        bool allow_most_bcs = true;
        if (fast_only) {
            allow_most_bcs = false;
        } else {
            allow_most_bcs = true;
        }

        // **********************************************************************************
        // NOTE: FillIntermediatePatch takes momenta at the new time, and returns
        //       BOTH updated velocities and momenta
        // **********************************************************************************
        cons_only = false;
        FillIntermediatePatch(level, time_for_fp,
                              {&S_data[IntVars::cons], &xvel_new, &yvel_new, &zvel_new},
                              {&S_data[IntVars::cons], &S_data[IntVars::xmom],
                               &S_data[IntVars::ymom], &S_data[IntVars::zmom]},
                              ng_cons_to_use, ng_vel, cons_only, scomp_cons, ncomp_cons,
                              allow_most_bcs);
    };
