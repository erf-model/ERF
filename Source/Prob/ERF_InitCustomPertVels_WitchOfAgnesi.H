
    amrex::Gpu::streamSynchronize();
    Real U_0 = 0.0;
    Real V_0 = 0.0;

    ParmParse pp_for_pert_vels("prob");
    pp_for_pert_vels.query("U_0", U_0);
    pp_for_pert_vels.query("V_0", V_0);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = U_0;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = V_0;
    });

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = 1. / dx[0]; dxInv[1] = 1. / dx[1]; dxInv[2] = 1. / dx[2];

    // Set the z-velocity from impenetrable condition
    if (sc.terrain_type == TerrainType::StaticFittedMesh) {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = WFromOmega(i, j, k, 0.0,
                                             x_vel_pert, y_vel_pert,
                                             mf_u, mf_v, z_nd, dxInv);
        });
    }
