
    // GATE: subsidence is disabled (all zeros)
    for (int k = 0; k <= khi+1; k++) {
        wbar[k] = 0.0;
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
