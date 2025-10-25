        // *****************************************************************************
        // Do semi-implicit solve for theta
        // *****************************************************************************

        const Real l_vert_implicit_fac = solverChoice.vert_implicit_fac[nrk];

        GpuArray<Real, AMREX_SPACEDIM*2> l_bc_neumann_vals_d;
        const int bc_comp = RhoTheta_comp;
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
            l_bc_neumann_vals_d[ori] = m_bc_neumann_vals[bc_comp][ori];
        }

        if (l_vert_implicit_fac > 0.) {
            const bool l_use_stretched_dz = (solverChoice.mesh_type == MeshType::StretchedDz);
            const bool l_use_SurfLayer = (m_SurfaceLayer != nullptr);
            const BCRec* bc_ptr_d = domain_bcs_type_d.data();

            const bool l_use_turb       = solverChoice.turbChoice[level].use_kturb;

            std::unique_ptr<MultiFab> dflux_z;
            dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, 1, 0);

            for ( MFIter mfi(S_old[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                Box bx  = mfi.tilebox();

              //const Array4<Real const> u_arr = xvel_new.const_array(mfi);
              //const Array4<Real const> v_arr = yvel_new.const_array(mfi);

                const Array4<      Real> & cell_data  = S_data[IntVars::cons].array(mfi);
              //const Array4<const Real> & cell_prim  = S_prim.array(mfi);

                const Array4<const Real>& z_nd_arr = z_phys_nd[level]->const_array(mfi);
                const Array4<const Real>& detJ_arr = detJ_cc[level]->const_array(mfi);

                const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

                const Array4<Real>& hfx_z   = Hfx3->array(mfi);

                if (l_use_stretched_dz) {
                    ImplicitDiffForState_S(bx, fine_geom.Domain(), level, dtau,
                                           l_bc_neumann_vals_d,
                                           cell_data,
                                           stretched_dz_d[level], hfx_z,
                                           mu_turb, solverChoice,
                                           bc_ptr_d, l_use_SurfLayer, l_vert_implicit_fac);
                } else if (l_use_terrain_fitted_coords) {
                    ImplicitDiffForState_T(bx, fine_geom.Domain(), level, dtau,
                                           l_bc_neumann_vals_d,
                                           cell_data,
                                           z_nd_arr, detJ_arr, dxInv, hfx_z,
                                           mu_turb, solverChoice,
                                           bc_ptr_d, l_use_SurfLayer, l_vert_implicit_fac);
                } else {
                    ImplicitDiffForState_N(bx, fine_geom.Domain(), level, dtau,
                                           l_bc_neumann_vals_d,
                                           cell_data,
                                           dxInv, hfx_z,
                                           mu_turb, solverChoice,
                                           bc_ptr_d, l_use_SurfLayer, l_vert_implicit_fac);
                }
            } // mfi
        } // if do implicit solve for diffusive contribution to (rho theta) update
