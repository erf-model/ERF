    const Real* prob_lo = geom.ProbLo();
    const auto dx = geom.CellSize();
    bool use_zlevels = (z_phys_cc != nullptr);

    for (MFIter mfi(*qsrc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& qsrc_arr = qsrc->array(mfi);
        if (box.length(0) != 1)
        {
            qsrc->setVal(0.0);
        } else {
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real z_cc = (use_zlevels) ? d_zlevels_arr[k] : prob_lo[2] + (k+0.5)*dx[2];
                Real f_a[9] = {-3.19860779e-39,  2.41718658e-34, -7.62715842e-30,  1.30171558e-25,
                               -1.30129938e-21,  7.71169160e-18, -2.59263914e-14,  4.05941658e-11,
                                2.56947279e-09};
                Real result = f_a[0];
                for (int ind = 1; ind < 9; ++ind) {
                    result = result * z_cc + f_a[ind];
                }
                qsrc_arr(i, j, k) = result;
            });
        }
    }
