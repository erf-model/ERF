
    ParmParse pp_prob("prob");

    Real wbar_cutoff_max = 1500.0; pp_prob.query("wbar_cutoff_max", wbar_cutoff_max);

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
