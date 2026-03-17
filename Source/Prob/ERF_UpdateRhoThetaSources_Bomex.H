    ParmParse pp_prob("prob");

    Real cutoff                 = 500.0; pp_prob.query("source_cutoff", cutoff);
    Real cutoff_transition      =  50.0; pp_prob.query("source_cutoff_transition", cutoff_transition);
    Real advection_heating_rate =   0.0; pp_prob.query("advection_heating_rate", advection_heating_rate);

    // Only apply temperature source below nominal inversion height
    for (MFIter mfi(*src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& src_arr = src->array(mfi);
        if (box.length(0) != 1)
        {
            src->setVal(0.0);
        } else {
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real z_cc = d_zlevels_arr[k];
                if (z_cc < cutoff) {
                    src_arr(i, j, k) = advection_heating_rate;
                } else if (z_cc < cutoff+cutoff_transition) {
                    Real slope = -advection_heating_rate / cutoff_transition;
                    src_arr(i, j, k) = (z_cc-cutoff) * slope + advection_heating_rate;
                } else {
                    src_arr(i, j, k) = 0.0;
                }
            });
        }
    }
