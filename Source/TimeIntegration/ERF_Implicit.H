        // *****************************************************************************
        // Do semi-implicit solve for diffusion: theta, u, and v
        // *****************************************************************************
        {
            const Real l_vert_implicit_fac = solverChoice.vert_implicit_fac[nrk];

            if (l_vert_implicit_fac > 0.) {

                const bool l_do_implicit_theta = solverChoice.implicit_thermal_diffusion;
                const bool l_do_implicit_mom = solverChoice.implicit_momentum_diffusion;

                // If we're doing an implicit solve for momenta (u and v only),
                // then we use the explicit solution at this point to derive
                // the surface velocity gradient -- to be consistent with the
                // hoextrap used with the surface layer BC. We then subtract
                // out the contribution from the gradient of tau_corr times
                // (1 - implicit_fac) and add in the implicit solution scaled
                // by implicit_fac.
                MultiFab* Tau13corr = (l_do_implicit_mom) ? Tau_corr[level][0].get() : nullptr;
                MultiFab* Tau23corr = (l_do_implicit_mom) ? Tau_corr[level][1].get() : nullptr;
#ifdef ERF_IMPLICIT_W
                MultiFab* Tau33corr = (l_do_implicit_mom) ? Tau_corr[level][2].get() : nullptr;
#endif

                // BCs
                const BCRec* bc_ptr_h = domain_bcs_type.data();
                GpuArray<Real, AMREX_SPACEDIM*2> l_bc_neumann_vals_d;
                for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
                    l_bc_neumann_vals_d[ori] = m_bc_neumann_vals[RhoTheta_comp][ori];
                }
                const bool l_use_SurfLayer = (m_SurfaceLayer != nullptr);

                const bool l_use_turb       = solverChoice.turbChoice[level].use_kturb;

                const bool l_use_stretched_dz = (solverChoice.mesh_type == MeshType::StretchedDz);

                for ( MFIter mfi(S_old[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
                {
                    Box bx  = mfi.tilebox();
                    Box tbx = mfi.nodaltilebox(0);
                    Box tby = mfi.nodaltilebox(1);

                    const Array4<      Real>& cell_data  = scratch.array(mfi);

                    const Array4<      Real>& rho_u      = (l_do_implicit_mom) ? scratch_xmom.array(mfi) : Array4<Real>{};
                    const Array4<      Real>& rho_v      = (l_do_implicit_mom) ? scratch_ymom.array(mfi) : Array4<Real>{};
                    const Array4<const Real>& tau13_corr = (l_do_implicit_mom) ? Tau13corr->array(mfi) : Array4<Real>{};
                    const Array4<const Real>& tau23_corr = (l_do_implicit_mom) ? Tau23corr->array(mfi) : Array4<Real>{};

#ifdef ERF_IMPLICIT_W
                    Box tbz = mfi.nodaltilebox(2);
                    const Array4<      Real>& rho_w      = (l_do_implicit_mom) ? scratch_zmom.array(mfi) : Array4<Real>{};
                    const Array4<const Real>& tau33_corr = (l_do_implicit_mom) ? Tau33corr->array(mfi) : Array4<Real>{};
#endif

                    const Array4<const Real>& z_nd_arr = z_phys_nd[level]->const_array(mfi);
                    const Array4<const Real>& detJ_arr = detJ_cc[level]->const_array(mfi);

                    const Array4<const Real>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

                    const Array4<const Real>& hfx_z   = Hfx3->const_array(mfi);

                    if (l_use_stretched_dz) {
                        if (l_do_implicit_theta) {
                            ImplicitDiffForState_S(bx, fine_geom.Domain(), level, slow_dt,
                                                   l_bc_neumann_vals_d,
                                                   cell_data,
                                                   stretched_dz_d[level], hfx_z,
                                                   mu_turb, solverChoice,
                                                   bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                        }
                        if (l_do_implicit_mom) {
                            ImplicitDiffForMom_S<0>(tbx, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_u, tau13_corr,
                                                    stretched_dz_d[level],
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);

                            ImplicitDiffForMom_S<1>(tby, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_v, tau23_corr,
                                                    stretched_dz_d[level],
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
#ifdef ERF_IMPLICIT_W
                            ImplicitDiffForMom_S<2>(tbz, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_w, tau33_corr,
                                                    stretched_dz_d[level],
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
#endif
                        }
                    } else if (l_use_terrain_fitted_coords) {
                        if (l_do_implicit_theta) {
                            ImplicitDiffForState_T(bx, fine_geom.Domain(), level, slow_dt,
                                                   l_bc_neumann_vals_d,
                                                   cell_data,
                                                   z_nd_arr, detJ_arr, dxInv, hfx_z,
                                                   mu_turb, solverChoice,
                                                   bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                        }
                        if (l_do_implicit_mom) {
                            ImplicitDiffForMom_T<0>(tbx, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_u, tau13_corr,
                                                    z_nd_arr, detJ_arr, dxInv,
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);

                            ImplicitDiffForMom_T<1>(tby, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_v, tau23_corr,
                                                    z_nd_arr, detJ_arr, dxInv,
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
#ifdef ERF_IMPLICIT_W
                            ImplicitDiffForMom_T<2>(tbz, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_w, tau33_corr,
                                                    z_nd_arr, detJ_arr, dxInv,
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
#endif
                        }
                    } else { // no terrain
                        if (l_do_implicit_theta) {
                            ImplicitDiffForState_N(bx, fine_geom.Domain(), level, slow_dt,
                                                   l_bc_neumann_vals_d,
                                                   cell_data,
                                                   dxInv, hfx_z,
                                                   mu_turb, solverChoice,
                                                   bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                        }
                        if (l_do_implicit_mom) {
                            ImplicitDiffForMom_N<0>(tbx, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_u, tau13_corr,
                                                    dxInv,
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);

                            ImplicitDiffForMom_N<1>(tby, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_v, tau23_corr,
                                                    dxInv,
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
#ifdef ERF_IMPLICIT_W
                            ImplicitDiffForMom_N<2>(tbz, fine_geom.Domain(), level, slow_dt,
                                                    cell_data, rho_w, tau33_corr,
                                                    dxInv,
                                                    mu_turb, solverChoice,
                                                    bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
#endif
                        }
                    }
                } // mfi

            } // if do implicit solve for diffusive contribution to (rho theta) update
        } // wrapper
