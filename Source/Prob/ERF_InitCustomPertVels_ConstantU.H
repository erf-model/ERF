
    ParmParse pp_prob("prob");
    Real U_0 = zero; pp_prob.query("U_0", U_0);
    Real V_0 = zero; pp_prob.query("V_0", V_0);
    Real W_0 = zero; pp_prob.query("W_0", W_0);

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

    // Set the z-velocity
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = W_0;
    });
