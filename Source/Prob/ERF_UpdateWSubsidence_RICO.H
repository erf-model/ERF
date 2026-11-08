
    ParmParse pp_prob("prob");

    // RICO: wbar decreases linearly up to 2260 m and is held constant above,
    // continuous at the cutoff
    Real wbar_cutoff_max = amrex::Real(2260.0); pp_prob.query("wbar_cutoff_max", wbar_cutoff_max);

    wbar[0] = zero;
    for (int k = 1; k <= khi+1; k++) {
        const Real z = zlevels[k];
        if (z <= wbar_cutoff_max) {
            wbar[k] = -amrex::Real(0.0000022)*z;
        } else {
            wbar[k] = -amrex::Real(0.0000022)*wbar_cutoff_max;
        }
    }

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
