    // Using Deardorff (see Sullivan et al 1994)
    //    or k-eqn RANS (see Axell & Liungman 2001)
    //
    // Note: At this point, the thermal diffusivity ("Khv" field in ERF) and the
    //       subgrid dissipation ("diss" here) have been updated by the closure
    //       at the beginning of each timestep.
    //       The strain rate magnitude is updated at the beginning of the first
    //       RK stage only, therefore the shear production term also does not
    //       change between RK stages.
    //       hfx_z holds the face fluxes of the theta diffusion of this stage
    //       (the full flux, independent of the implicit/explicit split), so
    //       the buoyancy production averages the two faces of the cell for the
    //       cell-centred vertical heat flux. The surface heat flux is stored at
    //       hfx_z at the bottom face (overwritten by the surface layer at each
    //       RK stage), so with a surface layer the first cell uses that flux.
    if (l_use_keqn && (start_comp <= RhoKE_comp) && (end_comp >= RhoKE_comp)) {
        int qty_index = RhoKE_comp;
        const int  l_klo_dom       = domain.smallEnd(2);
        const bool l_surf_flux_klo = SurfLayer_zlo;
        ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Add Buoyancy Source
            // where the SGS buoyancy flux tau_{theta,i} = -KH * dtheta/dx_i,
            // such that for dtheta/dz < 0, there is a positive (upward) heat
            // flux; the TKE buoyancy production is then
            //   B = g/theta_0 * tau_{theta,w}
            // for a dry atmosphere.
            // TODO: To account for moisture, the Brunt-Vaisala frequency,
            //   N^2 = g[1/theta * dtheta/dz + ...]
            // **should** be a function of the water vapor and total water
            // mixing ratios, depending on whether conditions are saturated or
            // not (see the WRF model description, Skamarock et al 2019).
            // NOTE: reads hfx_z at k+1, one face above bx; Hfx3 is z-nodal and
            //       the theta diffusion writes every face of its valid box, top included.
            const bool l_surf_cell = l_surf_flux_klo && (k == l_klo_dom);
            const amrex::Real hfx_cc = (l_surf_cell)
                                     ? hfx_z(i,j,k)
                                     : amrex::Real(0.5) * ( hfx_z(i,j,k) + hfx_z(i,j,k+1) );
            amrex::Real buoyancy_src = l_abs_g * l_inv_theta0 * hfx_cc;
            if (!use_ref_theta) {
                // l_inv_theta0 == 1, divide by actual theta
                buoyancy_src /= cell_prim(i,j,k,PrimTheta_comp);
            }
            cell_rhs(i,j,k,qty_index) += buoyancy_src;

            // TKE shear production
            //   P = -tau_ij * S_ij = 2 * mu_turb * S_ij * S_ij
            // Note: This assumes that the horizontal and vertical diffusivities
            // of momentum are equal
            cell_rhs(i,j,k,qty_index) += two*mu_turb(i,j,k,EddyDiff::Mom_v) * SmnSmn_a(i,j,k);

            // TKE dissipation (explicit); with erf.implicit_tke_dissipation the
            // sink is folded into the update in erf_slow_rhs_post instead
            if (!l_implicit_diss) {
                cell_rhs(i,j,k,qty_index) -= diss(i,j,k);
            }
        });
    }
