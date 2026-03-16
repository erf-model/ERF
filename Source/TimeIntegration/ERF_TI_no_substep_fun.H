/**
 *  Wrapper for advancing the solution with the slow RHS in the absence of acoustic substepping
 */
    auto no_substep_fun = [&](Vector<MultiFab>& S_sum,
                              Vector<MultiFab>& S_old,
                              Vector<MultiFab>& F_slow,
                              const Real time_for_fp, const Real slow_dt,
                              const int nrk)
    {
        BL_PROFILE("no_substep_fun");
        amrex::ignore_unused(nrk);
        int n_data = IntVars::NumTypes;

        const auto& dxInv = fine_geom.InvCellSizeArray();

        const amrex::GpuArray<int, IntVars::NumTypes> scomp_fast = {0,0,0,0};
        const amrex::GpuArray<int, IntVars::NumTypes> ncomp_fast = {2,1,1,1};

        if (verbose) amrex::Print() << " No-substepping time integration at level " << level
                                    << std::setprecision(timeprecision)
                                    << " to " << time_for_fp
                                    << " with dt = " << slow_dt << std::endl;

        // Update S_sum = S_stage only for the fast variables
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            for ( MFIter mfi(S_sum[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.tilebox();
                Box tbx = mfi.nodaltilebox(0);
                Box tby = mfi.nodaltilebox(1);
                Box tbz = mfi.nodaltilebox(2);

                Vector<Array4<Real> >  ssum_h(n_data);
                Vector<Array4<Real> >  sold_h(n_data);
                Vector<Array4<Real> > fslow_h(n_data);

                for (int i = 0; i < n_data; ++i) {
                     ssum_h[i]  = S_sum[i].array(mfi);
                     sold_h[i]  = S_old[i].array(mfi);
                    fslow_h[i] = F_slow[i].array(mfi);
                }

                Gpu::AsyncVector<Array4<Real> >  sold_d(n_data);
                Gpu::AsyncVector<Array4<Real> >  ssum_d(n_data);
                Gpu::AsyncVector<Array4<Real> > fslow_d(n_data);

                Gpu::copy(Gpu::hostToDevice,  sold_h.begin(),  sold_h.end(),  sold_d.begin());
                Gpu::copy(Gpu::hostToDevice,  ssum_h.begin(),  ssum_h.end(),  ssum_d.begin());
                Gpu::copy(Gpu::hostToDevice, fslow_h.begin(), fslow_h.end(), fslow_d.begin());

                Array4<Real>*  sold =  sold_d.dataPtr();
                Array4<Real>*  ssum =  ssum_d.dataPtr();
                Array4<Real>* fslow = fslow_d.dataPtr();

                // Moving terrain
                if ( solverChoice.terrain_type == TerrainType::MovingFittedMesh )
                {
                    const Array4<const Real>& dJ_old =     detJ_cc[level]->const_array(mfi);
                    const Array4<const Real>& dJ_new = detJ_cc_new[level]->const_array(mfi);
                    const Array4<const Real>& dJ_stg = detJ_cc_src[level]->const_array(mfi);

                    const Array4<const Real>& z_nd_old =     z_phys_nd[level]->const_array(mfi);
                    const Array4<const Real>& z_nd_new = z_phys_nd_new[level]->const_array(mfi);

                    const Array4<const Real>& mf_u = mapfac[level][MapFacType::u_x]->const_array(mfi);
                    const Array4<const Real>& mf_v = mapfac[level][MapFacType::v_y]->const_array(mfi);

                    const Array4<Real      >& z_t_arr  =  z_t_rk[level]->array(mfi);

                    // We have already scaled the slow source term to have the extra factor of dJ
                    ParallelFor(bx, ncomp_fast[IntVars::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVars::cons] + nn;
                        ssum[IntVars::cons](i,j,k,n) = dJ_old(i,j,k) *  sold[IntVars::cons](i,j,k,n)
                                           + slow_dt * dJ_stg(i,j,k) * fslow[IntVars::cons](i,j,k,n);
                        ssum[IntVars::cons](i,j,k,n) /= dJ_new(i,j,k);
                    });

                    // We have already scaled the slow source term to have the extra factor of dJ
                    ParallelFor(tbx, tby,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        Real h_zeta_old = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd_old);
                        Real h_zeta_new = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd_new);
                        ssum[IntVars::xmom](i,j,k) = ( h_zeta_old *  sold[IntVars::xmom](i,j,k)
                                                       + slow_dt * fslow[IntVars::xmom](i,j,k) ) / h_zeta_new;
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        Real h_zeta_old = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd_old);
                        Real h_zeta_new = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd_new);
                        ssum[IntVars::ymom](i,j,k) = ( h_zeta_old *  sold[IntVars::ymom](i,j,k)
                                                       + slow_dt * fslow[IntVars::ymom](i,j,k) ) / h_zeta_new;
                    });
                    ParallelFor(tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (k == 0) {
                            // Here we take advantage of the fact that moving terrain has a slip wall
                            // so we can just use the new value at (i,j,0).
                            Real rho_on_face =  ssum[IntVars::cons](i,j,k,Rho_comp);
                            ssum[IntVars::zmom](i,j,k) = WFromOmega(i,j,k,rho_on_face*z_t_arr(i,j,k),
                                                                    ssum[IntVars::xmom], ssum[IntVars::ymom],
                                                                    mf_u, mf_v, z_nd_new,dxInv);
                        } else {
                            Real dJ_old_kface = 0.5 * (dJ_old(i,j,k) + dJ_old(i,j,k-1));
                            Real dJ_new_kface = 0.5 * (dJ_new(i,j,k) + dJ_new(i,j,k-1));
                            ssum[IntVars::zmom](i,j,k) = ( dJ_old_kface *  sold[IntVars::zmom](i,j,k)
                                                             + slow_dt * fslow[IntVars::zmom](i,j,k) ) / dJ_new_kface;
                        }
                    });

                } else { // Fixed or no terrain
                    const Array4<const Real>& dJ_old = detJ_cc[level]->const_array(mfi);
                    ParallelFor(bx, ncomp_fast[IntVars::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVars::cons] + nn;
                        if (dJ_old(i,j,k) > 0.0) {
                            ssum[IntVars::cons](i,j,k,n) = sold[IntVars::cons](i,j,k,n) + slow_dt *
                               ( fslow[IntVars::cons](i,j,k,n) );
                        } else {
                            ssum[IntVars::cons](i,j,k,n) = sold[IntVars::cons](i,j,k,n);
                        }
                    });

                    // Commenting out the update is a HACK while developing the EB capability
                    if (solverChoice.terrain_type == TerrainType::EB)
                    {
                        const Array4<const Real>& vfrac_u = (get_eb(level).get_u_const_factory())->getVolFrac().const_array(mfi);
                        const Array4<const Real>& vfrac_v = (get_eb(level).get_v_const_factory())->getVolFrac().const_array(mfi);
                        const Array4<const Real>& vfrac_w = (get_eb(level).get_w_const_factory())->getVolFrac().const_array(mfi);

                        ParallelFor(tbx, tby, tbz,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            if (vfrac_u(i,j,k) > 0.0) {
                                // ssum[IntVars::xmom](i,j,k) = sold[IntVars::xmom](i,j,k);
                                ssum[IntVars::xmom](i,j,k) = sold[IntVars::xmom](i,j,k)
                                                           + slow_dt * fslow[IntVars::xmom](i,j,k);
                            } else {
                                ssum[IntVars::xmom](i,j,k) = sold[IntVars::xmom](i,j,k);
                            }
                        },
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            if (vfrac_v(i,j,k) > 0.0) {
                                // ssum[IntVars::ymom](i,j,k) = sold[IntVars::ymom](i,j,k);
                                ssum[IntVars::ymom](i,j,k) = sold[IntVars::ymom](i,j,k)
                                                           + slow_dt * fslow[IntVars::ymom](i,j,k);
                            } else {
                                ssum[IntVars::ymom](i,j,k) = sold[IntVars::ymom](i,j,k);
                            }
                        },
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            if (vfrac_w(i,j,k) > 0.0) {
                                // ssum[IntVars::zmom](i,j,k) = sold[IntVars::zmom](i,j,k);
                                ssum[IntVars::zmom](i,j,k) = sold[IntVars::zmom](i,j,k)
                                                           + slow_dt * fslow[IntVars::zmom](i,j,k);
                            } else {
                                ssum[IntVars::zmom](i,j,k) = sold[IntVars::zmom](i,j,k);
                            }
                        });

                    } else {
                        ParallelFor(tbx, tby, tbz,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            ssum[IntVars::xmom](i,j,k) = sold[IntVars::xmom](i,j,k)
                                                       + slow_dt * fslow[IntVars::xmom](i,j,k);
                        },
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            ssum[IntVars::ymom](i,j,k) = sold[IntVars::ymom](i,j,k)
                                                       + slow_dt * fslow[IntVars::ymom](i,j,k);
                        },
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            ssum[IntVars::zmom](i,j,k) = sold[IntVars::zmom](i,j,k)
                                                       + slow_dt * fslow[IntVars::zmom](i,j,k);
                        });
                    } // no EB
                } // not moving terrain
            } // mfi

            if (solverChoice.terrain_type == TerrainType::EB)
            {
                BL_PROFILE("no_substep_fun_redistribute");
                // Redistribute cons states (cell-centered)

                Vector<MultiFab> dUdt(IntVars::NumTypes);
                dUdt[IntVars::cons].define(ba, dm, F_slow[IntVars::cons].nComp(), F_slow[IntVars::cons].nGrow(), MFInfo(), EBFactory(level));
                for (int i = 1; i <= AMREX_SPACEDIM; ++i) {
                    dUdt[i].define(F_slow[i].boxArray(), F_slow[i].DistributionMap(), F_slow[i].nComp(), F_slow[i].nGrow(), MFInfo());
                }
                for (int i = 0; i <= AMREX_SPACEDIM; ++i) {
                    dUdt[i].setVal(0.0, 0, ncomp_fast[i], dUdt[i].nGrow());
                    MultiFab::Copy(dUdt[i], F_slow[i], 0, 0, F_slow[i].nComp(), 0);
                    dUdt[i].setDomainBndry(1.234e10, 0, ncomp_fast[i], fine_geom);
                    dUdt[i].FillBoundary(fine_geom.periodicity());
                }

                BCRec const* bc_ptr_d = domain_bcs_type_d.data();
                IntVect ngv_pert = S_old[IntVars::cons].nGrowVect();

                MultiFab S_pert_cons(ba, dm, S_old[IntVars::cons].nComp(), ngv_pert, MFInfo(), EBFactory(level));
                MultiFab rt0(ba, dm, 1, ngv_pert, MFInfo(), EBFactory(level));

                // Copy S_old into S_pert
                MultiFab::Copy(S_pert_cons, S_old[IntVars::cons],0,0,S_pert_cons.nComp(), ngv_pert);

                // Copy theta0 into rt0
                MultiFab::Copy(rt0,base_state[level],BaseState::th0_comp,0,1,ngv_pert);

                // Multiply theta0 by rho0 in rt0
                MultiFab::Multiply(rt0,base_state[level],BaseState::r0_comp,0,1,ngv_pert);

                MultiFab::Subtract(S_pert_cons, base_state[level], BaseState::r0_comp,0,1,ngv_pert);
                MultiFab::Subtract(S_pert_cons, rt0              ,                  0,1,1,ngv_pert);

                // Update F_slow for perturbational fast quantities (rho - rho0) and (rho_theta - rho_theta_0)
                redistribute_term(ncomp_fast[IntVars::cons], fine_geom, F_slow[IntVars::cons], dUdt[IntVars::cons],
                                  S_pert_cons, EBFactory(level), bc_ptr_d, slow_dt);

                // Update F_slow for momenta
                redistribute_term(ncomp_fast[IntVars::xmom], fine_geom, F_slow[IntVars::xmom], dUdt[IntVars::xmom],
                                  S_old[IntVars::xmom], *(get_eb(level).get_u_const_factory()), bc_ptr_d, slow_dt, IntVars::xmom);
                redistribute_term(ncomp_fast[IntVars::ymom], fine_geom, F_slow[IntVars::ymom], dUdt[IntVars::ymom],
                                  S_old[IntVars::ymom], *(get_eb(level).get_v_const_factory()), bc_ptr_d, slow_dt, IntVars::ymom);
                redistribute_term(ncomp_fast[IntVars::zmom], fine_geom, F_slow[IntVars::zmom], dUdt[IntVars::zmom],
                                  S_old[IntVars::zmom], *(get_eb(level).get_w_const_factory()), bc_ptr_d, slow_dt, IntVars::zmom);


                // Update state using the updated F_slow.
                for ( MFIter mfi(S_sum[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box bx = mfi.tilebox();
                    Box tbx = mfi.nodaltilebox(0);
                    Box tby = mfi.nodaltilebox(1);
                    Box tbz = mfi.nodaltilebox(2);

                    Vector<Array4<Real> >  ssum_h(n_data);
                    Vector<Array4<Real> >  sold_h(n_data);
                    Vector<Array4<Real> > fslow_h(n_data);

                    for (int i = 0; i < n_data; ++i) {
                         ssum_h[i] =  S_sum[i].array(mfi);
                         sold_h[i] =  S_old[i].array(mfi);
                        fslow_h[i] = F_slow[i].array(mfi);
                    }

                    Gpu::AsyncVector<Array4<Real> >  sold_d(n_data);
                    Gpu::AsyncVector<Array4<Real> >  ssum_d(n_data);
                    Gpu::AsyncVector<Array4<Real> > fslow_d(n_data);

                    Gpu::copy(Gpu::hostToDevice,  sold_h.begin(),  sold_h.end(),  sold_d.begin());
                    Gpu::copy(Gpu::hostToDevice,  ssum_h.begin(),  ssum_h.end(),  ssum_d.begin());
                    Gpu::copy(Gpu::hostToDevice, fslow_h.begin(), fslow_h.end(), fslow_d.begin());

                    Array4<Real>*  sold =  sold_d.dataPtr();
                    Array4<Real>*  ssum =  ssum_d.dataPtr();
                    Array4<Real>* fslow = fslow_d.dataPtr();

                    // const Array4<const Real>& vfrac_c = detJ_cc[level]->const_array(mfi);
                    const Array4<const Real>& vfrac_c = (get_eb(level).get_const_factory())->getVolFrac().const_array(mfi);

                    ParallelFor(bx, ncomp_fast[IntVars::cons], [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn)
                    {
                        const int n = scomp_fast[IntVars::cons] + nn;
                        if (vfrac_c(i,j,k) > 0.0) {
                            ssum[IntVars::cons](i,j,k,n) = sold[IntVars::cons](i,j,k,n) + slow_dt *
                                ( fslow[IntVars::cons](i,j,k,n) );
                        }
                    });

                    const Array4<const Real>& vfrac_u = (get_eb(level).get_u_const_factory())->getVolFrac().const_array(mfi);
                    const Array4<const Real>& vfrac_v = (get_eb(level).get_v_const_factory())->getVolFrac().const_array(mfi);
                    const Array4<const Real>& vfrac_w = (get_eb(level).get_w_const_factory())->getVolFrac().const_array(mfi);

                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (vfrac_u(i,j,k) > 0.0) {
                            ssum[IntVars::xmom](i,j,k) = sold[IntVars::xmom](i,j,k)
                            + slow_dt * fslow[IntVars::xmom](i,j,k);
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (vfrac_v(i,j,k) > 0.0) {
                            ssum[IntVars::ymom](i,j,k) = sold[IntVars::ymom](i,j,k)
                                                    + slow_dt * fslow[IntVars::ymom](i,j,k);
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (vfrac_w(i,j,k) > 0.0) {
                            ssum[IntVars::zmom](i,j,k) = sold[IntVars::zmom](i,j,k)
                                                    + slow_dt * fslow[IntVars::zmom](i,j,k);

                        }
                    });
                } // MFIter
            } // EB
        } // omp

        // Even if we update all the conserved variables we don't need
        // to fillpatch the slow ones every acoustic substep
        apply_bcs(S_sum, time_for_fp, S_sum[IntVars::cons].nGrow(), S_sum[IntVars::xmom].nGrow(),
                  fast_only=true, vel_and_mom_synced=false);

        if (solverChoice.anelastic[level]) {
            bool have_tb = (thin_xforce[0] || thin_yforce[0] || thin_zforce[0]);
            if (have_tb) {
                project_velocity_tb(level, slow_dt, S_sum);
            } else {
                project_momenta(level, time_for_fp, slow_dt, S_sum);
            }
        }
    };
