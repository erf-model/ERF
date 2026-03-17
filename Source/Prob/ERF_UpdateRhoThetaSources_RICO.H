    ParmParse pp_prob("prob");

    Real advection_heating_rate =   0.0; pp_prob.query("advection_heating_rate", advection_heating_rate);

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.

    // Only apply temperature source below nominal inversion height
    for ( amrex::MFIter mfi(*src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& src_arr = src->array(mfi);
        if (box.length(0) == 1) {
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                src_arr(i, j, k) = advection_heating_rate;
            });
        }
    }
