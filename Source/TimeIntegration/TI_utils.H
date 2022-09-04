// #include <AMReX_BC_TYPES.H>
// #include <SpatialStencils.H>
// #include <prob_common.H>
// #include <AMReX_TimeIntegrator.H>
// #include <TerrainMetrics.H>
// #include <TimeIntegration.H>
// #include <ERF.H>

    // **************************************************************************************
    // Temporary array that we use to store primitive advected quantities for the RHS
    // **************************************************************************************
    auto cons_to_prim = [&](const MultiFab& cons_state, MultiFab& prim_state, int ng)
    {
        BL_PROFILE("cons_to_prim()");
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(cons_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          const Box& gbx = mfi.growntilebox(ng);
          const Array4<const Real>& cons_arr = cons_state.array(mfi);
          const Array4<Real>& prim_arr = prim_state.array(mfi);

          amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int n = 0; n < NUM_PRIM; ++n) {
              prim_arr(i,j,k,PrimTheta_comp + n) = cons_arr(i,j,k,RhoTheta_comp + n) / cons_arr(i,j,k,Rho_comp);
            }
          });
      } // mfi
    };

    // ***************************************************************************************
    // This routine is called before the first step of the time integration, *and* in the case
    //  of a multi-stage method like RK3, this is called from "post_update_fun" which is called
    //  before every subsequent stage.  Since we advance the variables in conservative form,
    //  we must convert momentum to velocity before imposing the bcs.
    // ***************************************************************************************
    auto apply_bcs = [&](Vector<MultiFab>& S_data, const Real time_for_fp,
                         const Real time_mt, const Real delta_t,
                         int ng_cons, int ng_vel, bool fast_only)
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

        // ***************************************************************************************
        // Call FillPatch routines for the density only because we need it to convert between
        //      momentum and velocity
        // This fills ghost cells/faces from
        //     1) coarser level if lev > 0
        //     2) physical boundaries
        //     3) other grids at the same level
        // ***************************************************************************************
        int scomp_cons = 0;
        int ncomp_cons = 1;
        bool cons_only = true;

        // We must have at least one ghost cell of density to convert from momentum to velocity
        //    on the valid region
        AMREX_ALWAYS_ASSERT (ng_cons >= 1);

        FillIntermediatePatch(level, time_for_fp, time_mt, delta_t,
                              {S_data[IntVar::cons], xvel_new, yvel_new, zvel_new},
                              ng_cons, 0, cons_only, scomp_cons, ncomp_cons);

        // Here we don't use include any of the ghost region because we have only updated
        //      momentum on valid faces
        MomentumToVelocity(xvel_new, IntVect::TheZeroVector(),
                           yvel_new, IntVect::TheZeroVector(),
                           zvel_new, IntVect::TheZeroVector(),
                           S_data[IntVar::cons],
                           S_data[IntVar::xmom],
                           S_data[IntVar::ymom],
                           S_data[IntVar::zmom]);

        // ***************************************************************************************
        // Call FillPatch routines for all data
        // This fills ghost cells/faces from
        //     1) coarser level if lev > 0
        //     2) physical boundaries
        //     3) other grids at the same level
        // ***************************************************************************************
        scomp_cons = 0;
        if (fast_only) {
            ncomp_cons = 2; // rho and (rho theta)
        } else {
            ncomp_cons = S_data[IntVar::cons].nComp();
        }

        bool allow_most_bcs = true;
        if (fast_only) {
            allow_most_bcs = false;
        } else {
            allow_most_bcs = true;
        }

        // We must have at least one extra ghost cell of density to convert from velocity to momentum
        //    on the valid region
        int ng_cons_to_use = std::max(ng_cons, ng_vel+1);

        cons_only = false;
        FillIntermediatePatch(level, time_for_fp, time_mt, delta_t,
                              {S_data[IntVar::cons], xvel_new, yvel_new, zvel_new},
                              ng_cons_to_use, ng_vel, cons_only, scomp_cons, ncomp_cons,
                              allow_most_bcs);

        // Now we can convert back to momentum on valid+ghost since we have
        //     filled the ghost regions for both velocity and density
        VelocityToMomentum(xvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           yvel_new, IntVect(ng_vel,ng_vel,ng_vel),
                           zvel_new, IntVect(ng_vel,ng_vel,0),
                           S_data[IntVar::cons],
                           S_data[IntVar::xmom],
                           S_data[IntVar::ymom],
                           S_data[IntVar::zmom]);
    };
