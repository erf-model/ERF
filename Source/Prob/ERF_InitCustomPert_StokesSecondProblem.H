    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // This version perturbs rho but not p -- TODO: CHECK THIS
        state_pert(i, j, k, RhoTheta_comp) = std::pow(1.0,1.0/Gamma) * 101325.0 / 287.0 - p_hse(i,j,k);
    });
