    ParmParse pp_prob("prob");

    Real advection_moisture_rate      = 0.0; pp_prob.query("advection_moisture_rate", advection_moisture_rate);
    Real moisture_cutoff            = 500.0; pp_prob.query("moisture_source_cutoff", moisture_cutoff);

    // Only apply moisture source below nominal inversion height
    for (MFIter mfi(*qsrc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& qsrc_arr = qsrc->array(mfi);
        if (box.length(0) == 1)
        {
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real z_cc = d_zlevels_arr[k];
                if (z_cc < moisture_cutoff) {
                    qsrc_arr(i, j, k) =  5.226199353716133e-12*z_cc + advection_moisture_rate;
                } else {
                    qsrc_arr(i, j, k) = 5.226199353716133e-12*2980.0 + advection_moisture_rate;
                }
            });
        }
    }
