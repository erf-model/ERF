    const Real* prob_lo = geom.ProbLo();
    const auto dx = geom.CellSize();
    bool use_zlevels = (z_phys_cc != nullptr);

    for (MFIter mfi(*src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& src_arr = src->array(mfi);
        if (box.length(0) != 1)
        {
            src->setVal(0.0);
        } else {
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real z_cc = (use_zlevels) ? d_zlevels_arr[k] : prob_lo[2] + (k+0.5)*dx[2];
                if ((z_cc < 15000) && (z_cc > 380)) {
                    Real f_a[6] = {-1.22998448e-24,  5.86487977e-20, -1.11949167e-15,
                                    1.03417521e-11, -3.80858898e-08, -1.39371490e-05};
                    Real result = f_a[0];
                    for (int ind = 1; ind < 6; ++ind) {
                        result = result * z_cc + f_a[ind];
                    }
                    src_arr(i, j, k) = result;
                } else if (z_cc <= 380) {
                    Real f_a[2] = {1.73461284e-08, -3.00569096e-05};
                    Real result = f_a[0];
                    for (int ind = 1; ind < 2; ++ind) {
                        result = result * z_cc + f_a[ind];
                    }
                    src_arr(i, j, k) = result;
                } else {
                    src_arr(i, j, k) = 0.0;
                }
            });
        }
    }
