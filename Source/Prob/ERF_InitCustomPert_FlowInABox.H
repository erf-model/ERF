    ParmParse pp_prob("prob");
    Real T_0_Pert_Mag = 0.0; pp_prob.query("T_0_Pert_Mag", T_0_Pert_Mag);

    // Add temperature perturbations
    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*T_0_Pert_Mag;
        state_pert(i, j, k, RhoTheta_comp) *= r_hse(i,j,k);
    });
