    ParmParse pp_prob("prob");

    Real advection_heating_rate      = amrex::Real(5.4e-5); pp_prob.query("advection_heating_rate"     , advection_heating_rate);
    Real advection_heating_rate_base = amrex::Real(5.4e-5); pp_prob.query("advection_heating_rate_base", advection_heating_rate_base);

    auto prob_lo = geom.ProbLoArray();
    auto prob_hi = geom.ProbHiArray();
    auto dx = geom.CellSizeArray();

    // Define a point (xc,yc,zc) at the center of the domain
    const Real xc = myhalf * (prob_lo[0] + prob_hi[0]);
    const Real yc = myhalf * (prob_lo[1] + prob_hi[1]);

    // Only apply temperature source below nominal inversion height
    for (MFIter mfi(*src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& src_arr = src->array(mfi);
        if (box.length(0) != 1)
        {
            // if z dimension size is 1, then src is a spatially varying function over x,y at k=0
            ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                auto x = prob_lo[0] + (i + myhalf) * dx[0];
                auto y = prob_lo[1] + (j + myhalf) * dx[1];
                auto z = prob_lo[2] + (k + myhalf) * dx[2];

                auto c = amrex::Real(2890000.0);
                auto r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
                auto rsqr = r*r;
                if (time < 3600) {
                   src_arr(i, j, k) = (advection_heating_rate_base)*exp(-z/100);
                } else {
                   src_arr(i, j, k) = (advection_heating_rate*exp(-rsqr / c))*exp(-z/100);
                }
            });
        } else {
            // src is a function over Z
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                src_arr(i, j, k) = advection_heating_rate;
            });
        }
    }
