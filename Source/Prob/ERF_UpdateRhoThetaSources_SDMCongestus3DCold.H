    ParmParse pp_prob("prob");

    Real advection_heating_rate = 5.4e-5; pp_prob.query("advection_heating_rate", advection_heating_rate);
    Real source_cutoff          = 500.0;  pp_prob.query("source_cutoff", source_cutoff);

    auto prob_lo = geom.ProbLoArray();
    auto prob_hi = geom.ProbHiArray();
    auto dx = geom.CellSizeArray();

    // Define a point (xc,yc) at the center of the domain
    Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
    Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);

    for (MFIter mfi(*src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& src_arr = src->array(mfi);
        if (box.length(0) != 1)
        {
            // Spatially varying source over x,y
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                auto x = prob_lo[0] + (i + 0.5) * dx[0];
                auto y = prob_lo[1] + (j + 0.5) * dx[1];
                auto z = prob_lo[2] + (k + 0.5) * dx[2];

                const auto c = 2890000.0;
                const auto r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
                const auto rsqr = r*r;
                src_arr(i, j, k) = (advection_heating_rate*exp(-rsqr / c))*exp(-z/source_cutoff);
            });
        } else {
            // src is a function over Z only
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                src_arr(i, j, k) = advection_heating_rate;
            });
        }
    }
