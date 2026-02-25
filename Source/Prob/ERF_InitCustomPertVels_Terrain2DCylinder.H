
    amrex::Gpu::streamSynchronize();
    ParmParse pp("prob");
    Real U_0 = 10.0; pp.query("U_0", U_0);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = U_0;
    });
