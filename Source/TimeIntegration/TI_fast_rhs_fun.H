    auto fast_rhs_fun = [&](      Vector<MultiFab>& S_fast_rhs,
                                  Vector<MultiFab>& S_slow_rhs,
                                  Vector<MultiFab>& S_stage_data,
                            const Vector<MultiFab>& S_data,
                                  Vector<MultiFab>& S_scratch,
                                  const Real dtau, const Real inv_fac, const Real time)
    {
        if (verbose) amrex::Print() << "Calling fast rhs at level " << level << " with dt = " << dtau << std::endl;

        // Moving terrain
        MultiFab* z_t_pert = nullptr;
        if ( solverChoice.use_terrain &&  (solverChoice.terrain_type == 1) )
        {
            if (verbose) Print() << "Making new geometry at substep " << time << std::endl;

            const auto dz = fine_geom.CellSizeArray()[2];

            init_custom_terrain(fine_geom,*z_phys_nd_new[level],time); // Defines h(i,j,k=0,t)
            init_terrain_grid  (fine_geom,*z_phys_nd_new[level]);            // Defines z_phys for all i,j,k
            make_metrics       (fine_geom,*z_phys_nd_new[level], *z_phys_cc_new[level], *detJ_cc_new[level]);  // Defines detJ

            z_t_pert = new MultiFab(S_data[IntVar::zmom].boxArray(), S_data[IntVar::zmom].DistributionMap(), 1, 1);

            for (MFIter mfi(*z_t_pert,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(IntVect(1,1,0));

                // This is the top face
                int khi = gbx.bigEnd()[2];

                // Now make this a 2D box (after we extract khi)
                gbx.setRange(2,0);

                const Array4<const Real>& detJ_new  =   detJ_cc_new[level]->const_array(mfi);
                const Array4<const Real>& detJ_old  =       detJ_cc[level]->const_array(mfi);

                const Array4<      Real>& zp_t_arr  =  z_t_pert->array(mfi);
                const Array4<      Real>& z_t_arr   =  z_t_rk[level]->array(mfi);

                Real delta_t  = time - 0.5*dtau;

                const auto& dX = fine_geom.CellSizeArray();

                // Loop over horizontal plane
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
                {
                    // Evaluate between RK stages
                    zp_t_arr(i,j,0) = dhdt(i,j,dX,time,delta_t);
                    for (int k = 1; k <= khi; k++)
                    {
                        // NOTE: J in this code is (1/J) in Sullivan et al (JAS 2014)
                        Real J_new = detJ_new(i,j,k-1);
                        Real J_old = detJ_old(i,j,k-1);
                        zp_t_arr(i,j,k) = zp_t_arr(i,j,k-1) + ( (J_new - J_old) / delta_t ) * dz;

                        // Convert to perturbation: z"_t(t) = z_t(t) - z_t^{RK}
                        zp_t_arr(i,j,k) -= z_t_arr(i,j,k);
                    }
                });
            } // mfi
        }

        erf_fast_rhs(level, time, S_fast_rhs, S_slow_rhs, S_stage_data, S_prim,
                     S_data, S_scratch, fine_geom, ifr, solverChoice, z_t_pert,
                     z_phys_nd[level], detJ_cc[level], r0, p0, dtau, inv_fac);
    };

    auto post_no_substep_fun = [&](Vector<MultiFab>& S_sum,
                                   Vector<MultiFab>& F_slow,
                                   const Real time_for_fp, const Real dt)
    {
        int n_data = IntVar::NumVars;

        const amrex::GpuArray<int, IntVar::NumVars> scomp_fast = {0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_fast = {2,1,1,1};

        // Update S_sum = S_stage only for the fast variables
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            for ( MFIter mfi(S_sum[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.tilebox();
                const Box tbx = mfi.nodaltilebox(0);
                const Box tby = mfi.nodaltilebox(1);
                const Box tbz = mfi.nodaltilebox(2);

                Vector<Array4<Real> >  ssum_h(n_data);
                Vector<Array4<Real> > fslow_h(n_data);

                for (int i = 0; i < n_data; ++i) {
                     ssum_h[i]  = S_sum[i].array(mfi);
                    fslow_h[i] = F_slow[i].array(mfi);
                }

                Gpu::AsyncVector<Array4<Real> >  ssum_d(n_data);
                Gpu::AsyncVector<Array4<Real> > fslow_d(n_data);

                Gpu::copy(Gpu::hostToDevice,  ssum_h.begin(),  ssum_h.end(),  ssum_d.begin());
                Gpu::copy(Gpu::hostToDevice, fslow_h.begin(), fslow_h.end(), fslow_d.begin());

                Array4<Real>*  ssum =  ssum_d.dataPtr();
                Array4<Real>* fslow = fslow_d.dataPtr();

                // Moving terrain
                if ( solverChoice.use_terrain && solverChoice.terrain_type == 1 )
                {
                    const Array4<const Real>& dJ_old =     detJ_cc[level]->const_array(mfi);
                    const Array4<const Real>& dJ_new = detJ_cc_new[level]->const_array(mfi);

                    ParallelFor(bx, ncomp_fast[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVar::cons] + nn;
                        ssum[IntVar::cons](i,j,k,n) += dt *
                           ( fslow[IntVar::cons](i,j,k,n) );
                        ssum[IntVar::cons](i,j,k,n) *= ( dJ_old(i,j,k) / dJ_new(i,j,k) );
                    });

                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::xmom](i,j,k) += dt *
                           ( fslow[IntVar::xmom](i,j,k) );
                        Real dJ_new_xface = 0.5 * (dJ_new(i,j,k) + dJ_new(i-1,j,k));
                        Real dJ_old_xface = 0.5 * (dJ_old(i,j,k) + dJ_old(i-1,j,k));
                        ssum[IntVar::xmom](i,j,k) *= ( dJ_old_xface / dJ_new_xface );
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::ymom](i,j,k) += dt *
                           ( fslow[IntVar::ymom](i,j,k) );
                        Real dJ_new_yface = 0.5 * (dJ_new(i,j,k) + dJ_new(i,j-1,k));
                        Real dJ_old_yface = 0.5 * (dJ_old(i,j,k) + dJ_old(i,j-1,k));
                        ssum[IntVar::ymom](i,j,k) *= ( dJ_old_yface / dJ_new_yface );
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::zmom](i,j,k) += dt *
                           ( fslow[IntVar::zmom](i,j,k) );
                        if (k > 0) {
                            Real dJ_new_zface = 0.5 * (dJ_new(i,j,k) + dJ_new(i,j,k-1));
                            Real dJ_old_zface = 0.5 * (dJ_old(i,j,k) + dJ_old(i,j,k-1));
                            ssum[IntVar::zmom](i,j,k) *= ( dJ_old_zface / dJ_new_zface );
                        }
                    });

                // No terrain, or terrain is static
                } else {

                    ParallelFor(bx, ncomp_fast[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVar::cons] + nn;
                        ssum[IntVar::cons](i,j,k,n) += dt *
                           ( fslow[IntVar::cons](i,j,k,n) );
                    });

                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::xmom](i,j,k) += dt *
                           ( fslow[IntVar::xmom](i,j,k) );
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::ymom](i,j,k) += dt *
                           ( fslow[IntVar::ymom](i,j,k) );
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::zmom](i,j,k) += dt *
                           ( fslow[IntVar::zmom](i,j,k) );
                    });
                }
            }
        }

        // Moving terrain
        Real time_mt = time_for_fp - 0.5*dt;

        // Even if we update all the conserved variables we don't need to fillpatch the slow ones every acoustic substep
        int ng_cons    = 1;
        int ng_vel     = 1;
        bool fast_only = true;
        apply_bcs(S_sum, time_for_fp, time_mt, dt, ng_cons, ng_vel, fast_only);
    };

    auto post_substep_fun = [&](Vector<MultiFab>& S_sum,
                                Vector<MultiFab>& F_slow, Vector<MultiFab>& F_pert,
                                const Real time_for_fp, const Real dtau)
    {
        int n_data = IntVar::NumVars;

        const amrex::GpuArray<int, IntVar::NumVars> scomp_fast = {0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_fast = {2,1,1,1};

        // Update S_sum = S_pert + S_stage only for the fast variables
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            for ( MFIter mfi(S_sum[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box bx = mfi.tilebox();
                const Box tbx = mfi.nodaltilebox(0);
                const Box tby = mfi.nodaltilebox(1);
                const Box tbz = mfi.nodaltilebox(2);

                Vector<Array4<Real> >  ssum_h(n_data);
                Vector<Array4<Real> > fslow_h(n_data);
                Vector<Array4<Real> > fpert_h(n_data);

                for (int i = 0; i < n_data; ++i) {
                     ssum_h[i]  = S_sum[i].array(mfi);
                    fslow_h[i] = F_slow[i].array(mfi);
                    fpert_h[i] = F_pert[i].array(mfi);
                }

                Gpu::AsyncVector<Array4<Real> >  ssum_d(n_data);
                Gpu::AsyncVector<Array4<Real> > fslow_d(n_data);
                Gpu::AsyncVector<Array4<Real> > fpert_d(n_data);

                Gpu::copy(Gpu::hostToDevice,  ssum_h.begin(),  ssum_h.end(),  ssum_d.begin());
                Gpu::copy(Gpu::hostToDevice, fslow_h.begin(), fslow_h.end(), fslow_d.begin());
                Gpu::copy(Gpu::hostToDevice, fpert_h.begin(), fpert_h.end(), fpert_d.begin());

                Array4<Real>*  ssum =  ssum_d.dataPtr();
                Array4<Real>* fslow = fslow_d.dataPtr();
                Array4<Real>* fpert = fpert_d.dataPtr();

                // Moving terrain
                if ( solverChoice.use_terrain && solverChoice.terrain_type == 1 )
                {
                    const Array4<const Real>& dJ_old =     detJ_cc[level]->const_array(mfi);
                    const Array4<const Real>& dJ_new = detJ_cc_new[level]->const_array(mfi);

                    ParallelFor(bx, ncomp_fast[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVar::cons] + nn;
                        ssum[IntVar::cons](i,j,k,n) += dtau *
                           ( fslow[IntVar::cons](i,j,k,n) + fpert[IntVar::cons](i,j,k,n) );
                        ssum[IntVar::cons](i,j,k,n) *= ( dJ_old(i,j,k) / dJ_new(i,j,k) );
                    });

                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::xmom](i,j,k) += dtau *
                           ( fslow[IntVar::xmom](i,j,k) + fpert[IntVar::xmom](i,j,k) );
                        Real dJ_new_xface = 0.5 * (dJ_new(i,j,k) + dJ_new(i-1,j,k));
                        Real dJ_old_xface = 0.5 * (dJ_old(i,j,k) + dJ_old(i-1,j,k));
                        ssum[IntVar::xmom](i,j,k) *= ( dJ_old_xface / dJ_new_xface );
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::ymom](i,j,k) += dtau *
                           ( fslow[IntVar::ymom](i,j,k) + fpert[IntVar::ymom](i,j,k));
                        Real dJ_new_yface = 0.5 * (dJ_new(i,j,k) + dJ_new(i,j-1,k));
                        Real dJ_old_yface = 0.5 * (dJ_old(i,j,k) + dJ_old(i,j-1,k));
                        ssum[IntVar::ymom](i,j,k) *= ( dJ_old_yface / dJ_new_yface );
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::zmom](i,j,k) += dtau *
                           ( fslow[IntVar::zmom](i,j,k) + fpert[IntVar::zmom](i,j,k) );
                        if (k > 0) {
                            Real dJ_new_zface = 0.5 * (dJ_new(i,j,k) + dJ_new(i,j,k-1));
                            Real dJ_old_zface = 0.5 * (dJ_old(i,j,k) + dJ_old(i,j,k-1));
                            ssum[IntVar::zmom](i,j,k) *= ( dJ_old_zface / dJ_new_zface );
                        }
                    });

                // No terrain, or terrain is static
                } else {

                    ParallelFor(bx, ncomp_fast[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVar::cons] + nn;
                        ssum[IntVar::cons](i,j,k,n) += dtau *
                           ( fslow[IntVar::cons](i,j,k,n) + fpert[IntVar::cons](i,j,k,n) );
                    });

                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::xmom](i,j,k) += dtau *
                           ( fslow[IntVar::xmom](i,j,k) + fpert[IntVar::xmom](i,j,k) );
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::ymom](i,j,k) += dtau *
                           ( fslow[IntVar::ymom](i,j,k) + fpert[IntVar::ymom](i,j,k));
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::zmom](i,j,k) += dtau *
                           ( fslow[IntVar::zmom](i,j,k) + fpert[IntVar::zmom](i,j,k) );
                    });
                }
            }
        }

        // Moving terrain
        Real time_mt = time_for_fp - 0.5*dtau;

        // Even if we update all the conserved variables we don't need to fillpatch the slow ones every acoustic substep
        int ng_cons    = 1;
        int ng_vel     = 1;
        bool fast_only = true;
        apply_bcs(S_sum, time_for_fp, time_mt, dtau, ng_cons, ng_vel, fast_only);
    };
