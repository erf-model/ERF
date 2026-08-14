
    for (int k = 0; k <= khi; k++)
    {
        const Real z_cc = zlevels[k];
        const Real u_geo_wind = amrex::Real(3.0) + amrex::Real(0.0043)*z_cc;
        const Real v_geo_wind = -amrex::Real(9.0) + amrex::Real(0.0056)*z_cc;

        u_geos[k] = u_geo_wind;
        v_geos[k] = v_geo_wind;
    }
