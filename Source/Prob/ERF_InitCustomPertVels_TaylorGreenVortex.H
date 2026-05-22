
    ParmParse pp_for_pert_vels("prob");
    Real V_0   =   one; pp_for_pert_vels.query("V_0", V_0);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + zero) * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = prob_lo[2] + (k + myhalf) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = V_0 * std::sin(x) * std::cos(y) * std::cos(z);
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + (j + zero) * dx[1];
        const Real z = prob_lo[2] + (k + myhalf) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = - V_0 * std::cos(x) * std::sin(y) * std::cos(z);
    });
