    ParmParse pp_prob("prob");
    Real U_0 = zero; pp_prob.query("U_0", U_0);

    const GpuArray<Real, AMREX_SPACEDIM> prob_lo = {
        geomdata.ProbLo(0), geomdata.ProbLo(1), geomdata.ProbLo(2)};
    const auto prob_hi = geomdata.ProbHi();
    const auto dx = geomdata.CellSize();
    const Real lx = prob_hi[0] - prob_lo[0];
    const Real ly = prob_hi[1] - prob_lo[1];
    const Real lz = prob_hi[2] - prob_lo[2];
    const auto length = GpuArray<Real, AMREX_SPACEDIM>{lx, ly, lz};
    const int nx = geomdata.Domain().length(0);
    const int ny = geomdata.Domain().length(1);
    const Real v_scale = erf_cloud_chamber::staggered_velocity_v_scale(
        lx, ly, dx[0], dx[1], nx, ny);

    // Manufactured field with a discrete divergence-free staggered pattern.
    // The v prefactor accounts for the unequal x/y finite-volume derivative
    // factors; it is exactly one on square domains with equal resolution.
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real x = prob_lo[0] + Real(i) * dx[0];
        const Real y = prob_lo[1] + (Real(j) + myhalf) * dx[1];
        const Real z = prob_lo[2] + (Real(k) + myhalf) * dx[2];
        x_vel_pert(i, j, k) = erf_cloud_chamber::initial_u_velocity_perturbation(
            x, y, z, prob_lo, length, U_0, nx, ny);
    });

    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real x = prob_lo[0] + (Real(i) + myhalf) * dx[0];
        const Real y = prob_lo[1] + Real(j) * dx[1];
        const Real z = prob_lo[2] + (Real(k) + myhalf) * dx[2];
        y_vel_pert(i, j, k) = erf_cloud_chamber::initial_v_velocity_perturbation(
            x, y, z, prob_lo, length, U_0, v_scale, nx, ny);
    });

    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = erf_cloud_chamber::initial_w_velocity_perturbation();
    });
