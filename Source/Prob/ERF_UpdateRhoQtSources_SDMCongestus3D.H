    ParmParse pp_prob("prob");

    Real advection_moisture_rate      = amrex::Real(1.0e-8); pp_prob.query("advection_moisture_rate"     , advection_moisture_rate);
    Real advection_moisture_rate_base = amrex::Real(1.0e-8); pp_prob.query("advection_moisture_rate_base", advection_moisture_rate_base);

    auto prob_lo = geom.ProbLoArray();
    auto prob_hi = geom.ProbHiArray();
    auto dx = geom.CellSizeArray();

    // Define a point (xc,yc,zc) at the center of the domain
    auto xc = myhalf * (prob_lo[0] + prob_hi[0]);
    auto yc = myhalf * (prob_lo[1] + prob_hi[1]);

    // Only apply temperature source below nominal inversion height
    for (MFIter mfi(*qsrc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        const Array4<Real>& qsrc_arr = qsrc->array(mfi);
        if (box.length(0) != 1)
        {
            // if z dimension size is 1, then src is a spatially varying function over x,y at k=0
            ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                auto x = prob_lo[0] + (i + myhalf) * dx[0];
                auto y = prob_lo[1] + (j + myhalf) * dx[1];
                auto z = prob_lo[2] + (k + myhalf) * dx[2];

                //const Real c = two*dx[0];
                auto c = amrex::Real(1700.0)*amrex::Real(1700.0);
                auto r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
                auto rsqr = r*r;
                if (time < 3600) {
                   qsrc_arr(i, j, k) = (advection_moisture_rate_base)*exp(-z/100);
                } else {
                   qsrc_arr(i, j, k) = (advection_moisture_rate*exp(-rsqr / c))*exp(-z/100);
                }
            });
        } else {
            // src is a function over Z
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qsrc_arr(i, j, k) = zero;
            });
        }
    }
