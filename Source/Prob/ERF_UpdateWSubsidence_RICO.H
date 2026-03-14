
    ParmParse pp_prob("prob");

    Real wbar_cutoff_max = 1500.0; pp_prob.query("wbar_cutoff_max", wbar_cutoff_max);

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.
    Vector<Real> zlevels;
    zlevels.resize(khi+2);
    reduce_to_max_per_height(zlevels, z_phys_nd);

    wbar[0] = 0.0;
    for (int k = 1; k <= khi+1; k++) {
        const Real z = zlevels[k];
        if (z <= wbar_cutoff_max) {
            wbar[k] = -0.0000022*z;
        } else {
            wbar[k] = -0.0000022*2260.;
        }
    }

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
