
    ParmParse pp_prob("prob");

    Real wbar_sub_max    = -0.65 ; pp_prob.query("wbar_sub_max", wbar_sub_max);
    Real wbar_cutoff_max = 1500.0; pp_prob.query("wbar_cutoff_max", wbar_cutoff_max);
    Real wbar_cutoff_min = 2100.0; pp_prob.query("wbar_cutoff_min", wbar_cutoff_min);
    AMREX_ASSERT_WITH_MESSAGE(wbar_cutoff_min > wbar_cutoff_max, "ERROR: wbar_cutoff_min < wbar_cutoff_max");

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.
    Vector<Real> zlevels;
    zlevels.resize(khi+1);
    reduce_to_max_per_height(zlevels, z_phys_nd);

    // Linearly increase wbar to the cutoff_max and then linearly decrease to cutoff_min
    Real z_0    = zlevels[0];

    Real slope1 =  wbar_sub_max / (wbar_cutoff_max - z_0);
    Real slope2 = -wbar_sub_max / (wbar_cutoff_min - wbar_cutoff_max);

    wbar[0]     = 0.0;

    for (int k = 1; k <= khi; k++) {
        const Real z_nd = zlevels[k];
        if (z_nd <= wbar_cutoff_max) {
            wbar[k] = slope1 * (z_nd - z_0);
        } else if (z_nd <= wbar_cutoff_min) {
            wbar[k] = slope2 * (z_nd - wbar_cutoff_max) + wbar_sub_max;
        } else {
            wbar[k] = 0.0;
        }
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
