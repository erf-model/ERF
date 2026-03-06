    ParmParse pp_prob("prob");

    Real advection_moisture_rate      = 0.0; pp_prob.query("advection_moisture_rate", advection_moisture_rate);
    Real moisture_cutoff            = 500.0; pp_prob.query("moisture_source_cutoff", moisture_cutoff);
    Real moisture_cutoff_transition =  50.0; pp_prob.query("moisture_source_cutoff_transition", moisture_cutoff_transition);


    Vector<Real> zlevels;
    Gpu::DeviceVector<Real> d_zlevels;

    const int khi       = geom.Domain().bigEnd()[2];

    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.
    zlevels.resize(khi+1);
    d_zlevels.resize(khi+1);
    reduce_to_max_per_height(zlevels, z_phys_cc);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, zlevels.begin(), zlevels.end(), d_zlevels.begin());

    const Real* d_zlevels_arr = d_zlevels.dataPtr();

    // Only apply moisture source below nominal inversion height
    for ( amrex::MFIter mfi(*qsrc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& qsrc_arr = qsrc->array(mfi);
        if (box.length(0) != 1)
        {
            // spatially varying source not used in this problem
            qsrc->setVal(0.0);
        } else {
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real z_cc = d_zlevels_arr[k];
                if (z_cc < moisture_cutoff) {
                    qsrc_arr(i, j, k) = advection_moisture_rate;
                } else if (z_cc < moisture_cutoff+moisture_cutoff_transition) {
                    Real slope = -advection_moisture_rate / moisture_cutoff_transition;
                    qsrc_arr(i, j, k) = (z_cc-moisture_cutoff) * slope + advection_moisture_rate;
                } else {
                    qsrc_arr(i, j, k) = 0.0;
                }
            });
        }
    }
