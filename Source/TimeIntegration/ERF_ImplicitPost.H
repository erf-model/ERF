        // *****************************************************************************
        // Do semi-implicit solve for diffusion: q1
        // *****************************************************************************
        const bool l_do_implicit_moist = solverChoice.implicit_moisture_diffusion;
        const Real l_vert_implicit_fac = solverChoice.vert_implicit_fac[nrk];
        if ( (l_do_implicit_moist) &&
             (l_vert_implicit_fac > zero) &&
             (solverChoice.moisture_type != MoistureType::None) ) {

            // BCs
            const BCRec* bc_ptr_h = domain_bcs_type.data();
            GpuArray<Real, AMREX_SPACEDIM*2> l_bc_neumann_vals_d;
            for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
                l_bc_neumann_vals_d[ori] = m_bc_neumann_vals[RhoQ1_comp][ori];
            }
            const bool l_use_SurfLayer = (m_SurfaceLayer != nullptr);

            const bool l_use_turb       = solverChoice.turbChoice[level].use_kturb;

            const bool l_use_stretched_dz = (solverChoice.mesh_type == MeshType::StretchedDz);

            for ( MFIter mfi(S_new[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                Box bx  = mfi.tilebox();

                // NOTE: We have updated the state, use that directly
                const Array4<      Real>& cell_data  = S_new[IntVars::cons].array(mfi);

                const Array4<const Real>& z_nd_arr = z_phys_nd[level]->const_array(mfi);
                const Array4<const Real>& detJ_arr = detJ_cc[level]->const_array(mfi);

                const Array4<const Real>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) :
                                                                     Array4<const Real>{};
                const Array4<const Real>& q1fx_z = Q1fx3->const_array(mfi);

                if (l_use_stretched_dz) {
                    ImplicitDiffForStateLU_S(bx, fine_geom.Domain(), level, RhoQ1_comp,
                                             slow_dt, l_bc_neumann_vals_d, cell_data,
                                             stretched_dz_d[level], q1fx_z,
                                             mu_turb, solverChoice,
                                             bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                } else if (l_use_terrain_fitted_coords) {
                    ImplicitDiffForStateLU_T(bx, fine_geom.Domain(), level, RhoQ1_comp,
                                             slow_dt, l_bc_neumann_vals_d, cell_data,
                                             z_nd_arr, detJ_arr, dxInv, q1fx_z,
                                             mu_turb, solverChoice,
                                             bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                } else { // no terrain
                    ImplicitDiffForStateLU_N(bx, fine_geom.Domain(), level, RhoQ1_comp,
                                             slow_dt, l_bc_neumann_vals_d, cell_data,
                                             dxInv, q1fx_z, mu_turb, solverChoice,
                                             bc_ptr_h, l_use_SurfLayer, l_vert_implicit_fac);
                }
            } // mfi
        } // if moist model and do implicit diff
