
    // Only apply momentum source below nominal inversion height
    for (int k = 0; k <= khi; k++)
    {
        const Real z_cc = zlevels[k];
        const Real u_geo_wind = -amrex::Real(10.0) + z_cc * amrex::Real(0.0018);

        u_geos[k] =  u_geo_wind; // 0; // -coriolis_factor * v_geo_wind
        v_geos[k] =  0 ; // coriolis *  u_geo_wind;
    }
