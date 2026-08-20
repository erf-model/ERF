
    wbar[0] = zero;
    for (int k = 1; k <= khi+1; k++) {
        const Real z = zlevels[k];
        wbar[k] = -amrex::Real(0.00000375)*z;
    }

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
