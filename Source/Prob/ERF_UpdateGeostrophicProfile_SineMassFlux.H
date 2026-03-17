    const Real* prob_lo = geom.ProbLo();
    const auto dx = geom.CellSize();

    for (int k = 0; k <= khi; k++) {
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)*dx[2];
        const Real u_geo_wind = -10.0 + z_cc * 0.0018;

        u_geos[k] = u_geo_wind;
        v_geos[k] = 0.0;
    }
