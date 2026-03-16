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

          //
          // We may need > one ghost cells of prim in order to compute higher order advective terms
          //
          amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              Real rho       = cons_arr(i,j,k,Rho_comp);
              Real rho_theta = cons_arr(i,j,k,RhoTheta_comp);
              prim_arr(i,j,k,PrimTheta_comp) = rho_theta / rho;
              for (int n = 1; n < ncomp_prim; ++n) {
                  prim_arr(i,j,k,PrimTheta_comp + n) = cons_arr(i,j,k,RhoTheta_comp + n) / rho;
              }
          });

          //
          // We only use one ghost cell of pi_stage so we only fill one here
          //
          const Box& gbx1 = mfi.growntilebox(1);

          const Array4<      Real>& pi_stage_arr = pi_stage.array(mfi);
          const Real rdOcp = solverChoice.rdOcp;
          const bool l_use_moisture  = (solverChoice.moisture_type != MoistureType::None);

          amrex::ParallelFor(gbx1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              auto qv_for_p = (l_use_moisture) ? cons_arr(i,j,k,RhoQ1_comp)/cons_arr(i,j,k,Rho_comp) : 0.0;
              pi_stage_arr(i,j,k) = getExnergivenRTh(cons_arr(i,j,k,RhoTheta_comp), rdOcp, qv_for_p);
          });
      } // mfi
    };

/**
 *  Define the total water content in each cell
 */
    auto make_qt = [&](const MultiFab& cons_state, MultiFab& qt)
    {
        BL_PROFILE("make_qt()");

        int n_qstate_moist = micro->Get_Qstate_Moist_Size();

        // All moisture models are guaranteed to have RhoQ1_comp
        MultiFab::Copy(qt, cons_state, RhoQ1_comp, 0, 1, qt.nGrowVect());

        // Sum up each component
        for (int n(1); n<n_qstate_moist; ++n) {
            MultiFab::Add(qt, cons_state, RhoQ1_comp+n, 0, 1, qt.nGrowVect());
        }

        // Divide out the dry density
        MultiFab::Divide(qt, cons_state, Rho_comp, 0, 1, qt.nGrowVect());
    };

/**
 *  This routine is called after the scalars are updated for each RK stage
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

        // **********************************************************************************
        // NOTE: FillIntermediatePatch takes momenta at the new time, and returns
        //       BOTH updated velocities and momenta
        // **********************************************************************************
        cons_only = false;
        FillIntermediatePatch(level, time_for_fp,
                              {&S_data[IntVars::cons], &xvel_new, &yvel_new, &zvel_new},
                              {&S_data[IntVars::cons], &S_data[IntVars::xmom],
                               &S_data[IntVars::ymom], &S_data[IntVars::zmom]},
                              ng_cons_to_use, ng_vel, cons_only, scomp_cons, ncomp_cons);
    };

    auto update_terrain_stage = [&](int lev, Real old_step_time, Real old_stage_time, Real new_stage_time, Real slow_dt)
    {
            // Note that the "old" and "new" metric terms correspond to
            // t^n and the RK stage (either t^*, t^** or t^{n+1} that this source
            // will be used to advance to

            // The "src" metric terms correspond to the time at which we are evaluating the source here,
            // aka old_stage_time

            if (verbose) Print() << "Re-making old geometry at old time   : " << old_step_time << std::endl;

            //
            // Make a temporary MF to fill the terrain
            //
            Box terrain_bx(surroundingNodes(Geom(0).Domain())); terrain_bx.grow(z_phys_nd[lev]->nGrow());
            FArrayBox terrain_fab_old(makeSlab(terrain_bx,2,0),1);
            FArrayBox terrain_fab_src(makeSlab(terrain_bx,2,0),1);
            FArrayBox terrain_fab_new(makeSlab(terrain_bx,2,0),1);

            prob->init_terrain_surface(fine_geom,terrain_fab_old,old_step_time);
            prob->init_terrain_surface(fine_geom,terrain_fab_src,old_stage_time);
            prob->init_terrain_surface(fine_geom,terrain_fab_new,new_stage_time);

            // Copy on intersection
            for (MFIter mfi(*z_phys_nd[lev],false); mfi.isValid(); ++mfi)
            {
                Box isect = terrain_fab_old.box() & (*z_phys_nd[lev])[mfi].box();
                (    *z_phys_nd[lev])[mfi].template copy<RunOn::Device>(terrain_fab_old,isect,0,isect,0,1);
                (*z_phys_nd_src[lev])[mfi].template copy<RunOn::Device>(terrain_fab_src,isect,0,isect,0,1);
                (*z_phys_nd_new[lev])[mfi].template copy<RunOn::Device>(terrain_fab_new,isect,0,isect,0,1);
            }

            make_terrain_fitted_coords(lev,fine_geom,*z_phys_nd[lev], zlevels_stag[lev], phys_bc_type);
            make_J             (fine_geom, *z_phys_nd[lev], *detJ_cc[lev]);
            make_areas         (fine_geom, *z_phys_nd[lev], *ax[lev], *ay[lev], *az[lev]);

            if (verbose) Print() << "Making src geometry at old_stage_time:  "
                                 << std::setprecision(timeprecision) << old_stage_time << std::endl;
            make_terrain_fitted_coords(lev,fine_geom,*z_phys_nd_src[lev], zlevels_stag[lev], phys_bc_type);
            make_J             (fine_geom, *z_phys_nd_src[lev], *detJ_cc_src[lev]);
            make_areas         (fine_geom, *z_phys_nd_src[lev], *ax_src[lev], *ay_src[lev], *az_src[lev]);
            make_zcc           (fine_geom, *z_phys_nd_src[lev], *z_phys_cc_src[lev]);

            if (verbose) Print() << "Making new geometry at new_stage_time: "
                                 << std::setprecision(timeprecision) << new_stage_time << std::endl;
            make_terrain_fitted_coords(lev,fine_geom,*z_phys_nd_new[lev], zlevels_stag[lev], phys_bc_type);
            make_J             (fine_geom, *z_phys_nd_new[lev], *detJ_cc_new[lev]);

            Real inv_dt  = 1./slow_dt;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*z_t_rk[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(IntVect(1,1,0));

                const Array4<Real      >& z_t_arr      =  z_t_rk[lev]->array(mfi);
                const Array4<Real const>& z_nd_new_arr =  z_phys_nd_new[lev]->const_array(mfi);
                const Array4<Real const>& z_nd_old_arr =  z_phys_nd[lev]->const_array(mfi);

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
    };

    auto update_terrain_substep = [&](int lev, Real old_substep_time, Real new_substep_time, Real dtau,
                                      Vector<MultiFab>& S_data, std::unique_ptr<MultiFab>& z_t_pert)
    {
        // Make "old" fast geom -- store in z_phys_nd for convenience
        Box terrain_bx(surroundingNodes(Geom(0).Domain())); terrain_bx.grow(z_phys_nd[lev]->nGrow());

        if (verbose) Print() << "Making geometry at start of substep time: "
                             << std::setprecision(timeprecision) << old_substep_time << std::endl;
        FArrayBox terrain_fab_old(makeSlab(terrain_bx,2,0),1);
        prob->init_terrain_surface(fine_geom,terrain_fab_old,old_substep_time);

        // Make "new" fast geom
        if (verbose) Print() << "Making geometry for end of substep time :"
                             << std::setprecision(timeprecision) << new_substep_time << std::endl;
        FArrayBox terrain_fab_new(makeSlab(terrain_bx,2,0),1);
        prob->init_terrain_surface(fine_geom,terrain_fab_new,new_substep_time);

        // Copy on intersection
        for (MFIter mfi(*z_phys_nd[lev],false); mfi.isValid(); ++mfi)
        {
            Box isect = terrain_fab_old.box() & (*z_phys_nd[lev])[mfi].box();
            (    *z_phys_nd[lev])[mfi].template copy<RunOn::Device>(terrain_fab_old,isect,0,isect,0,1);
            (*z_phys_nd_new[lev])[mfi].template copy<RunOn::Device>(terrain_fab_new,isect,0,isect,0,1);
        }

        make_terrain_fitted_coords(lev,fine_geom,*z_phys_nd[lev], zlevels_stag[lev], phys_bc_type);
        make_J             (fine_geom,*z_phys_nd[lev], *detJ_cc[lev]);

        make_terrain_fitted_coords(lev,fine_geom,*z_phys_nd_new[lev], zlevels_stag[lev], phys_bc_type);
        make_J             (fine_geom,*z_phys_nd_new[lev], *detJ_cc_new[lev]);

        Real inv_dt   = 1./dtau;

        z_t_pert = std::make_unique<MultiFab>(S_data[IntVars::zmom].boxArray(), S_data[IntVars::zmom].DistributionMap(), 1, 1);

        for (MFIter mfi(*z_t_rk[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.growntilebox(IntVect(1,1,0));

            const Array4<Real      >& z_t_arr      =  z_t_rk[lev]->array(mfi);
            const Array4<Real      >& zp_t_arr     =  z_t_pert->array(mfi);

            const Array4<Real const>& z_nd_new_arr =  z_phys_nd_new[lev]->const_array(mfi);
            const Array4<Real const>& z_nd_old_arr =  z_phys_nd[lev]->const_array(mfi);

            // Loop over horizontal plane
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Evaluate between RK stages assuming the geometry is linear between old and new time
                zp_t_arr(i,j,k) = 0.25 * inv_dt * (z_nd_new_arr(i+1,j+1,k) - z_nd_old_arr(i+1,j+1,k)
                                                  +z_nd_new_arr(i  ,j+1,k) - z_nd_old_arr(  i,j+1,k)
                                                  +z_nd_new_arr(i+1,j  ,k) - z_nd_old_arr(i+1,j  ,k)
                                                  +z_nd_new_arr(i  ,j  ,k) - z_nd_old_arr(i  ,j  ,k));
                // Convert to perturbation: z"_t(t) = z_t(t) - z_t^{RK}
                zp_t_arr(i,j,k) -= z_t_arr(i,j,k);
            });
        } // mfi
    };
