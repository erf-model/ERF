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
                                                 bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                    }
                } // mfi
            } // do implicit tke and evolve tke

            if ( (l_do_implicit_moist) &&
                 (solverChoice.moisture_type != MoistureType::None) ) {

                // BCs for moisture (Neumann values are filled per species below)
                const BCRec* bc_ptr_h = domain_bcs_type.data();
                const bool l_use_SurfLayer = (m_SurfaceLayer[Orientation(Direction::z, Orientation::low)] != nullptr);

                // Same moisture range as NativeScalarDiffusionPolicy::is_moist (the water
                // species only), whose explicit vertical flux is scaled by
                // (1 - vert_implicit_fac).  Non-water species packed after them stay explicit.
                // With SBM only the host-owned qv is diffused; the projected liquid lanes are
                // owned by the spectral solver.
                const bool l_sbm = (solverChoice.moisture_type == MoistureType::SBM);
                const int l_n_moist = solverChoice.diffChoice.n_qstate_moist;
                AMREX_ALWAYS_ASSERT(RhoQ1_comp + l_n_moist <= S_new[IntVars::cons].nComp());
                const int l_moist_end = (l_sbm) ? RhoQ1_comp + amrex::min(l_n_moist, 1)
                                                : RhoQ1_comp + l_n_moist;
                erf_sbm::require_host_write_range_allowed(
                    l_sbm, RhoQ1_comp, l_moist_end - RhoQ1_comp, erf_sbm::HostWritePath::Diffusion);

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
                    for (int m_comp(RhoQ1_comp); m_comp<l_moist_end; ++m_comp) {
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
                            ImplicitDiffForStateLU_N(bx, fine_geom.Domain(), level, m_comp,
                                                     slow_dt, l_bc_neumann_vals_d, cell_data,
                                                     dxInv, mfx_z, mu_turb, solverChoice,
                                                     bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac,
                                                     l_use_mrf_cg);
                        }
                    } // m_comp
                } // mfi
            } // if moist model and do implicit diff
        } // implicit vertical factor > 0
