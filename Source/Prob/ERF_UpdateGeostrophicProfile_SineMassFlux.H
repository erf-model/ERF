    for (int k = 0; k <= khi; k++) {
        const Real z_cc = zlevels[k];
        const Real u_geo_wind = -10.0 + z_cc * 0.0018;

        u_geos[k] = u_geo_wind;
        v_geos[k] = 0.0;
    }
