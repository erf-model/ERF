
    amrex::Gpu::streamSynchronize();
    ParmParse pp_for_pert_vels("prob");
    Real U_0   = 0.0;  pp_for_pert_vels.query("U_0", U_0);
    Real V_0   = 0.0;  pp_for_pert_vels.query("V_0", V_0);
    Real W_0   = 0.0;  pp_for_pert_vels.query("W_0", W_0);

    Real z0    = 0.1;  pp_for_pert_vels.query("z0", z0);
    Real zRef  = 80.0; pp_for_pert_vels.query("zRef", zRef);
    Real uRef  = 0.0;  pp_for_pert_vels.query("uRef", uRef);

    int prob_type = -1; pp_for_pert_vels.query("prob_type", prob_type);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real*      dx = geomdata.CellSize();
        const Real        z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        if (uRef != 0.0) {
            x_vel_pert(i, j, k) = U_0 + uRef *
                                  std::log((z + z0)/z0) / std::log((zRef +z0)/z0);
        } else {
            x_vel_pert(i, j, k) = U_0;
        }
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = V_0;
    });

    // Set the z-velocity
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = W_0;
    });
