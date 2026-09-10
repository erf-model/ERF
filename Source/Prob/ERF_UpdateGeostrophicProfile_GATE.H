    const Real* prob_lo = geom.ProbLo();
    const auto dx = geom.CellSize();

    for (int k = 0; k <= khi; k++) {
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)*dx[2];
        Real f_a[12] = {-1.27239999e-41,  1.11717520e-36, -4.33415278e-32,  9.86488742e-28,
                        -1.46662495e-23,  1.47887323e-19, -9.96258532e-16,  4.17832623e-12,
                        -9.41191724e-09,  8.63378170e-06, -4.29197266e-03, -8.18730365e-01};
        Real u_geo_wind = f_a[0];
        for (int ind = 1; ind < 12; ++ind) {
            u_geo_wind = u_geo_wind * z_cc + f_a[ind];
        }
        u_geos[k] = u_geo_wind;
        v_geos[k] = 0.0;
    }
