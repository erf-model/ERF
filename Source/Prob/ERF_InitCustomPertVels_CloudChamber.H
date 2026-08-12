    ParmParse pp_prob("prob");
    Real U_0 = zero; pp_prob.query("U_0", U_0);

    const auto prob_lo = geomdata.ProbLo();
    const auto prob_hi = geomdata.ProbHi();
    const auto dx = geomdata.CellSize();
    const Real lx = prob_hi[0] - prob_lo[0];
    const Real ly = prob_hi[1] - prob_lo[1];
    const Real lz = prob_hi[2] - prob_lo[2];

    // Manufactured divergence-free field: u = d(psi)/dy, v = -d(psi)/dx,
    // with psi proportional to sin^2(pi*x/Lx) sin^2(pi*y/Ly)
    // sin^2(pi*z/Lz). It is zero at every no-slip face.
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real x = prob_lo[0] + Real(i) * dx[0];
        const Real y = prob_lo[1] + (Real(j) + myhalf) * dx[1];
        const Real z = prob_lo[2] + (Real(k) + myhalf) * dx[2];
        const Real sx = std::sin(PI * (x - prob_lo[0]) / lx);
        const Real sz = std::sin(PI * (z - prob_lo[2]) / lz);
        x_vel_pert(i, j, k) =
            U_0 * sx * sx *
            std::sin(two * PI * (y - prob_lo[1]) / ly) * sz * sz;
    });

    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real x = prob_lo[0] + (Real(i) + myhalf) * dx[0];
        const Real y = prob_lo[1] + Real(j) * dx[1];
        const Real z = prob_lo[2] + (Real(k) + myhalf) * dx[2];
        const Real sy = std::sin(PI * (y - prob_lo[1]) / ly);
        const Real sz = std::sin(PI * (z - prob_lo[2]) / lz);
        y_vel_pert(i, j, k) =
            -U_0 *
            std::sin(two * PI * (x - prob_lo[0]) / lx) * sy * sy * sz * sz;
    });

    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = zero;
    });
