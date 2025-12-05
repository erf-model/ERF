    // Using Deardorff (see Sullivan et al 1994)
    //    or k-eqn RANS (see Axell & Liungman 2001)
    //
    // Note: At this point, the thermal diffusivity ("Khv" field in ERF), the
    //       subgrid heat flux ("hfx_z" here), and the subgrid dissipation
    //       ("diss" here) have been updated by ComputeTurbulentViscosityLES --
    //       at the beginning of each timestep.
    //       The strain rate magnitude is updated at the beginning of the first
    //       RK stage only, therefore the shear production term also does not
    //       change between RK stages.
    //       The surface heat flux hfx_z(i,j,-1) is updated in MOSTStress at
    //       each RK stage, but that does not change the buoyancy production term here.
    if (l_use_keqn && (start_comp <= RhoKE_comp) && (end_comp >= RhoKE_comp)) {
        int qty_index = RhoKE_comp;
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
            cell_rhs(i,j,k,qty_index) += l_abs_g * l_inv_theta0 * hfx_z(i,j,k);

            // TKE shear production
            //   P = -tau_ij * S_ij = 2 * mu_turb * S_ij * S_ij
            // Note: This assumes that the horizontal and vertical diffusivities
            // of momentum are equal
            cell_rhs(i,j,k,qty_index) += 2.0*mu_turb(i,j,k,EddyDiff::Mom_v) * SmnSmn_a(i,j,k);

            // TKE dissipation
            cell_rhs(i,j,k,qty_index) -= diss(i,j,k);
        });
    }
