
    ParmParse pp("prob");
    Real V_0   =   1.0; pp.query("V_0", V_0);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.0) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = V_0 * sin(x) * cos(y) * cos(z);
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.0) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = - V_0 * cos(x) * sin(y) * cos(z);
    });
