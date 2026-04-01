
    ParmParse pp("prob");

    Real KE_0  = amrex::Real(0.1); pp.query("KE_0", KE_0);
    Real KE_decay_height = -1; pp.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order  =  1; pp.query("KE_decay_order" , KE_decay_order );

    if (KE_decay_height > 0) {
        amrex::Print() << "Initial KE profile (order " << KE_decay_order
                       << ") will extend up to " << KE_decay_height
                       << std::endl;
    }

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Set an initial value for SGS KE
        const Real z = z_cc(i,j,k);
        if (state_pert.nComp() > RhoKE_comp) {
            state_pert(i, j, k, RhoKE_comp) = r_hse(i,j,k) * KE_0;
            // scale initial SGS kinetic energy with height
            if (KE_decay_height > 0) {
                state_pert(i, j, k, RhoKE_comp) *= max(
                    std::pow(1 - min(z/KE_decay_height,one), KE_decay_order),
                    1e-12);
            }
        }
    });
