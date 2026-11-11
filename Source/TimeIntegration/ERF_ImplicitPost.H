        // *****************************************************************************
        // Do semi-implicit solve for vertical diffusion: KE and moisture
        // *****************************************************************************
        const bool l_do_implicit_ke    = solverChoice.implicit_ke_diffusion;
        const bool l_do_implicit_moist = solverChoice.implicit_moisture_diffusion;
        const Real l_vert_implicit_fac = solverChoice.vert_implicit_fac[level][nrk];
        if ( l_vert_implicit_fac > zero ) {

            const bool l_use_turb         = solverChoice.turbChoice[level].use_kturb;
            const bool l_use_stretched_dz = (solverChoice.mesh_type == MeshType::StretchedDz);

            if ( (l_do_implicit_ke) &&
                 (solverChoice.turbChoice[level].use_tke)) {

                // BCs for KE
                const BCRec* bc_ptr_h = domain_bcs_type.data();
                GpuArray<Real, AMREX_SPACEDIM*2> l_bc_neumann_vals_d;
                for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
                    l_bc_neumann_vals_d[ori] = m_bc_neumann_vals[RhoKE_comp][ori];
                }
                const bool l_use_SurfLayer = false;

                for ( MFIter mfi(S_new[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
                {
                    Box bx  = mfi.tilebox();

                    // NOTE: We have updated the state, use that directly
                    const Array4<      Real>& cell_data  = S_new[IntVars::cons].array(mfi);

                    const Array4<const Real>& z_nd_arr = z_phys_nd[level]->const_array(mfi);
                    const Array4<const Real>& detJ_arr = detJ_cc[level]->const_array(mfi);
                    const Array4<const Real> mf_mx = mapfac[level][MapFacType::m_x]->const_array(mfi);
                    const Array4<const Real> mf_my = mapfac[level][MapFacType::m_y]->const_array(mfi);
                    const bool l_implicit_metric = solverChoice.implicit_terrain_metric &&
                                                   l_use_terrain_fitted_coords;

                    const Array4<const Real>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) :
                                                                     Array4<const Real>{};

                    if (l_use_stretched_dz) {
                        ImplicitDiffForStateLU_S(bx, fine_geom.Domain(), level, RhoKE_comp,
                                                 slow_dt, l_bc_neumann_vals_d, cell_data,
                                                 stretched_dz_d[level], Array4<const Real>{},
                                                 mu_turb, solverChoice,
                                                 bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                    } else if (l_use_terrain_fitted_coords) {
                        ImplicitDiffForStateLU_T(bx, fine_geom.Domain(), level, RhoKE_comp,
                                                 slow_dt, l_bc_neumann_vals_d, cell_data,
                                                 z_nd_arr, detJ_arr, dxInv, Array4<const Real>{},
                                                 mu_turb, solverChoice,
                                                 bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac,
                                                 false, mf_mx, mf_my, l_implicit_metric);
                    } else { // no terrain
                        ImplicitDiffForStateLU_N(bx, fine_geom.Domain(), level, RhoKE_comp,
                                                 slow_dt, l_bc_neumann_vals_d, cell_data,
                                                 dxInv, Array4<const Real>{}, mu_turb, solverChoice,
                                                 bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac, false,
                                                 (ib_wall_face33[level]) ? ib_wall_face33[level]->const_array(mfi)
                                                                         : Array4<const Real>{},
                                                 (ib_wall_ktarget[level] && solverChoice.if_wall_tke_time_factor <= one)
                                                     ? ib_wall_ktarget[level]->const_array(mfi) : Array4<const Real>{});
                    }
                } // mfi
            } // do implicit tke and evolve tke

            if ( (l_do_implicit_moist) &&
                 (solverChoice.moisture_type != MoistureType::None) ) {

                // BCs for moisture (Neumann values are filled per species below)
                const BCRec* bc_ptr_h = domain_bcs_type.data();
                const bool l_use_SurfLayer = (m_SurfaceLayer[Orientation(Direction::z, Orientation::low)] != nullptr);

                const bool l_sbm = (solverChoice.moisture_type == MoistureType::SBM);

                for ( MFIter mfi(S_new[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
                {
                    Box bx  = mfi.tilebox();

                    // NOTE: We have updated the state, use that directly
                    const Array4<      Real>& cell_data  = S_new[IntVars::cons].array(mfi);

                    const Array4<const Real>& z_nd_arr = z_phys_nd[level]->const_array(mfi);
                    const Array4<const Real>& detJ_arr = detJ_cc[level]->const_array(mfi);
                    const Array4<const Real> mf_mx = mapfac[level][MapFacType::m_x]->const_array(mfi);
                    const Array4<const Real> mf_my = mapfac[level][MapFacType::m_y]->const_array(mfi);
                    const bool l_implicit_metric = solverChoice.implicit_terrain_metric &&
                                                   l_use_terrain_fitted_coords;

                    const Array4<const Real>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) :
                                                                     Array4<const Real>{};
                    for (int m_comp : implicit_post_comps) {
                        if (m_comp < RhoQ1_comp) { continue; }
                        AMREX_ALWAYS_ASSERT(m_comp < S_new[IntVars::cons].nComp());
                        erf_sbm::require_host_write_range_allowed(
                            l_sbm, m_comp, 1, erf_sbm::HostWritePath::Diffusion);
                        // Neumann values are stored by BC variable, not by state component
                        const int l_bc_comp =
                            ResolveNativeScalarDiffusionPolicy(m_comp, solverChoice.diffChoice).bc_comp;
                        GpuArray<Real, AMREX_SPACEDIM*2> l_bc_neumann_vals_d;
                        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
                            l_bc_neumann_vals_d[ori] = m_bc_neumann_vals[l_bc_comp][ori];
                        }

                        const bool is_q1 = (m_comp == RhoQ1_comp);
                        const Array4<const Real>& mfx_z = (is_q1 && Q1fx3) ? Q1fx3->const_array(mfi)
                                                                           : Array4<const Real> {};
                        const bool l_use_mrf_cg = solverChoice.turbChoice[level].enable_mrf_countergradient
                                                && (solverChoice.turbChoice[level].pbl_type == PBLType::MRF)
                                                && is_q1;
                        if (l_use_stretched_dz) {
                            ImplicitDiffForStateLU_S(bx, fine_geom.Domain(), level, m_comp,
                                                     slow_dt, l_bc_neumann_vals_d, cell_data,
                                                     stretched_dz_d[level], mfx_z,
                                                     mu_turb, solverChoice,
                                                     bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac,
                                                     l_use_mrf_cg);
                        } else if (l_use_terrain_fitted_coords) {
                            ImplicitDiffForStateLU_T(bx, fine_geom.Domain(), level, m_comp,
                                                     slow_dt, l_bc_neumann_vals_d, cell_data,
                                                     z_nd_arr, detJ_arr, dxInv, mfx_z,
                                                     mu_turb, solverChoice,
                                                     bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac,
                                                     l_use_mrf_cg,
                                                     mf_mx, mf_my, l_implicit_metric);
                        } else { // no terrain
                            // the fraction-stress wall face is closed to the moisture species as
                            // to theta and TKE (no wall value to pin: no wall_pin)
                            ImplicitDiffForStateLU_N(bx, fine_geom.Domain(), level, m_comp,
                                                     slow_dt, l_bc_neumann_vals_d, cell_data,
                                                     dxInv, mfx_z, mu_turb, solverChoice,
                                                     bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac,
                                                     l_use_mrf_cg,
                                                     (ib_wall_face33[level]) ? ib_wall_face33[level]->const_array(mfi)
                                                                             : Array4<const Real>{},
                                                     Array4<const Real>{});
                        }
                    } // m_comp
                } // mfi

                // The solves above ran on whole columns, the solid cells of an immersed terrain
                // included (the wall face decouples them from the air, the zlo surface layer
                // still feeds the bottom row): with the fraction-stress wall law hold the
                // moisture species of fully immersed cells at their step-start value again, as
                // erf_slow_rhs_post did before the solves
                if (solverChoice.if_fraction_stress && terrain_blank) {
                    const DiffChoice l_diff_choice = solverChoice.diffChoice;
                    for ( MFIter mfi(S_new[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
                    {
                        const Box bx = mfi.tilebox();
                        const Array4<      Real>& cell_data = S_new[IntVars::cons].array(mfi);
                        const Array4<const Real>& old_cons  = S_old[IntVars::cons].const_array(mfi);
                        const Array4<const Real>& t_blank   = terrain_blank->const_array(mfi);
                        for (int m_comp : implicit_post_comps) {
                            if (m_comp < RhoQ1_comp || !l_diff_choice.host_transports(m_comp)) { continue; }
                            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                            {
                                if (t_blank(i,j,k) == one) {
                                    cell_data(i,j,k,m_comp) = old_cons(i,j,k,m_comp);
                                }
                            });
                        }
                    }
                }
            } // if moist model and do implicit diff
        } // implicit vertical factor > 0
