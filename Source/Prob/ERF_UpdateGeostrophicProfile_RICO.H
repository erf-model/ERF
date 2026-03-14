    // If the z coordinate varies in time and or space, then the the height
    // needs to be calculated at each time step. Here, we assume that only
    // grid stretching exists.
    Vector<Real> zlevels;
    zlevels.resize(khi+1);
    reduce_to_max_per_height(zlevels, z_phys_cc);

    // Only apply momentum source below nominal inversion height
    for (int k = 0; k <= khi; k++) {
        const Real z_cc = zlevels[k];
        const Real u_geo_wind = -9.9 + z_cc * 0.002;

        u_geos[k] = u_geo_wind;
        v_geos[k] = -3.8;
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, u_geos.begin(), u_geos.end(), d_u_geos.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, v_geos.begin(), v_geos.end(), d_v_geos.begin());
