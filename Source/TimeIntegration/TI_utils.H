/**
 *  Define the primitive variables by dividing the conserved variables by density
 */
    auto cons_to_prim = [&](const MultiFab& cons_state, int ng)
    {
        BL_PROFILE("cons_to_prim()");
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

          amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real rho       = cons_arr(i,j,k,Rho_comp);
            Real rho_theta = cons_arr(i,j,k,RhoTheta_comp);
            prim_arr(i,j,k,PrimTheta_comp) = rho_theta / rho;
            pi_stage_arr(i,j,k) = getExnergivenRTh(rho_theta, rdOcp);
            for (int n = 1; n < NUM_PRIM; ++n) {
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
        amrex::Array<const MultiFab*,3> cmf_const{&xmom_crse, &ymom_crse, &zmom_crse};
        amrex::Array<MultiFab*,3> fmf{&S_data[IntVar::xmom],
                                      &S_data[IntVar::ymom],
                                      &S_data[IntVar::zmom]};

        // ***************************************************************************************
        // Interpolate momentum from coarse faces to fine faces *only* on the coarse-fine boundary
        // ***************************************************************************************
        if (level > 0) {
            ifr->interp(fmf,cmf_const,0,1);
        }

        int scomp_cons;
        int ncomp_cons;
        bool cons_only;

        int ng_cons_to_use;

        if (!vel_and_mom_synced) {

            // **********************************************************************************
            // Call FillPatch routines for the density only because we need it to convert between
            //      momentum and velocity
            // This fills ghost cells/faces from
            //     1) coarser level if lev > 0
            //     2) physical boundaries
            //     3) other grids at the same level
            // **********************************************************************************
            scomp_cons = 0;
            ncomp_cons = 1;
            cons_only  = true;

            // We must have at least one ghost cell of density to convert from momentum to velocity
            //    on the valid region
            AMREX_ALWAYS_ASSERT (ng_cons >= 1);

            // We must have at least one extra ghost cell of density to convert from velocity to momentum
            //    on the valid region
            ng_cons_to_use = std::max(ng_cons, ng_vel+1);

            FillIntermediatePatch(level, time_for_fp,
                                  {&S_data[IntVar::cons], &xvel_new, &yvel_new, &zvel_new},
                                  ng_cons_to_use, 0, cons_only, scomp_cons, ncomp_cons);

            // Here we don't use include any of the ghost region because we have only updated
            //      momentum on valid faces
            MultiFab density(S_data[IntVar::cons], make_alias, Rho_comp, 1);
            MomentumToVelocity(xvel_new, yvel_new, zvel_new, density,
                               S_data[IntVar::xmom],
                               S_data[IntVar::ymom],
                               S_data[IntVar::zmom]);
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
                ncomp_cons = S_data[IntVar::cons].nComp();
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
                ncomp_cons = S_data[IntVar::cons].nComp()-1; // since we filled rho above
            }
            ng_cons_to_use = ng_cons;
        }

        bool allow_most_bcs = true;
        if (fast_only) {
            allow_most_bcs = false;
        } else {
            allow_most_bcs = true;
        }

        cons_only = false;
        FillIntermediatePatch(level, time_for_fp,
                              {&S_data[IntVar::cons], &xvel_new, &yvel_new, &zvel_new},
                              ng_cons_to_use, ng_vel, cons_only, scomp_cons, ncomp_cons,
                              allow_most_bcs);

        // Now we can convert back to momentum on valid+ghost since we have
        //     filled the ghost regions for both velocity and density
        MultiFab density(S_data[IntVar::cons], make_alias, Rho_comp, 1);
        VelocityToMomentum(xvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           yvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           zvel_new, IntVect(ng_vel,ng_vel,0),
                           density,
                           S_data[IntVar::xmom],
                           S_data[IntVar::ymom],
                           S_data[IntVar::zmom],
                           solverChoice.use_NumDiff);
    };
