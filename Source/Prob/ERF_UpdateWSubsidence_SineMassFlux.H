    // Horizontally averaged density per level, summed over all boxes and ranks
    Gpu::DeviceVector<Real> rho_d(khi+1, zero);
    auto rho_arr = rho_d.data();
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto &box = mfi.tilebox();
        auto state_arr = state.const_array(mfi);
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Gpu::Atomic::AddNoRet(&rho_arr[k], state_arr(i, j, k, Rho_comp));
        });
    }

    Vector<Real> rho_h(khi+1, zero);
    Gpu::copy(Gpu::deviceToHost, rho_d.begin(), rho_d.end(), rho_h.begin());
    ParallelDescriptor::ReduceRealSum(rho_h.dataPtr(), khi+1);
    {
        const Box& dom = geom.Domain();
        const Real inv_ncell_xy = one / static_cast<Real>(dom.length(0)*dom.length(1));
        for (int k = 0; k <= khi; k++) { rho_h[k] *= inv_ncell_xy; }
    }

    // Staggered rho via linear interpolation
    Vector<Real> rho_stag(khi+2, zero);
    rho_stag[0]     = amrex::Real(1.5)*rho_h[0]   - myhalf*rho_h[1];
    rho_stag[khi+1] = amrex::Real(1.5)*rho_h[khi] - myhalf*rho_h[khi-1];
    for (int k = 1; k <= khi; k++) {
        rho_stag[k] = myhalf*(rho_h[k-1]+rho_h[k]);
    }

    // Sinusoidal mass flux: w = 2*sin(pi*t/amrex::Real(600))/rho for t < 600
    for (int k = 0; k <= khi+1; k++) {
        if (time < amrex::Real(600)) {
            Real real_time = static_cast<Real>(time) / amrex::Real(600);
            wbar[k] = two*std::sin(PI*real_time) / rho_stag[k];
        } else {
            wbar[k] = zero;
        }
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
