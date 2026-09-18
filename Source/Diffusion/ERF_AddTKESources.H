    // Using Deardorff (see Sullivan et al 1994)
    //    or k-eqn RANS (see Axell & Liungman 2001)
    //
    // Note: At this point, the thermal diffusivity ("Khv" field in ERF) and the
    //       subgrid dissipation ("diss" here) have been updated by the closure
    //       at the beginning of each timestep.
    //       The strain rate magnitude is updated at the beginning of the first
    //       RK stage only, therefore the shear production term also does not
    //       change between RK stages.
    //       hfx_z (Hfx3) is z-nodal, (0,0,1): it lives at the w-faces, not at the
    //       cell centre, so a single hfx_z(i,j,k) is the flux at the LOWER face of
    //       the cell and using it as the cell-centred value is not formally correct.
    //       hfx_z holds the face fluxes of the theta diffusion of this stage
    //       (the full flux, independent of the implicit/explicit split); at the
    //       bottom face of the domain that is the surface-layer flux when a
    //       surface layer is used.  The buoyancy production therefore averages the
    //       two faces of every cell, the first cell included, for the cell-centred
    //       vertical heat flux.  The first cell is not special-cased to the surface
    //       flux alone: that would again place a w-face value at the cell centre
    //       (a MOST-profile gradient evaluated at the first cell centre would be
    //       the consistent alternative).
    //       On a terrain-fitted mesh hfx_z holds the flux along zeta only: the slope
    //       corrections and the Jacobian are applied where the fluxes enter the RHS, so
    //       over sloping terrain this source omits those terms (as does the hfx3 output).
    //       REQUIRED WRITER: the RhoTheta call of DiffusionSrcForState_{N,S,T}
    //       in erf_slow_rhs_pre, on the same RK stage, which writes every z-face
    //       of each valid box (zbx = surroundingNodes(bx,2), top face included).
    //       Running the RhoKE source without that same-stage theta pass would
    //       read stale faces.
    if (l_use_keqn && (start_comp <= RhoKE_comp) && (end_comp >= RhoKE_comp)) {
        int qty_index = RhoKE_comp;
        // The source reads hfx_z at k+1: the z-faces of bx must be in hfx_z
        AMREX_ASSERT(hfx_z.contains(bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)) &&
                     hfx_z.contains(bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)+1));
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
            const amrex::Real hfx_cc = amrex::Real(0.5) * ( hfx_z(i,j,k) + hfx_z(i,j,k+1) );
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
