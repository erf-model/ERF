    ParmParse pp_prob("prob");
    Real T_0_Pert_Mag = zero; pp_prob.query("T_0_Pert_Mag", T_0_Pert_Mag);

    // Add temperature perturbations
    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        Real rand_double = amrex::Random(engine); // Between zero and one
        state_pert(i, j, k, RhoTheta_comp) = (rand_double*two - one)*T_0_Pert_Mag;
        state_pert(i, j, k, RhoTheta_comp) *= r_hse(i,j,k);
    });
