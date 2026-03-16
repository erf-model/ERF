    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.
    Vector<Real> zlevels;
    zlevels.resize(khi+1);
    reduce_to_max_per_height(zlevels, z_phys_cc);

    // const Real coriolis = 2.0 * 2.0 * PI / 86400.0; // 0.376E-4;

    // Only apply momentum source below nominal inversion height
    for (int k = 0; k <= khi; k++)
    {
        const Real z_cc = zlevels[k];
        const Real u_geo_wind = -10.0 + z_cc * 0.0018;

        u_geos[k] =  u_geo_wind; // 0; // -coriolis_factor * v_geo_wind
        v_geos[k] =  0 ; // coriolis *  u_geo_wind;
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, u_geos.begin(), u_geos.end(), d_u_geos.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, v_geos.begin(), v_geos.end(), d_v_geos.begin());
