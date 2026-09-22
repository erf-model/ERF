
    ParmParse pp("prob");

    Real KE_0  = zero; pp.query("KE_0", KE_0);
    Real KE_decay_height = -1; pp.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order  =  1; pp.query("KE_decay_order" , KE_decay_order );

    const bool use_eb = (sc.terrain_type == TerrainType::EB);
    const int klo = geomdata.Domain().smallEnd(2);
    if (KE_decay_height > 0 &&
        state_pert.nComp() > RhoKE_comp) {
        amrex::Print() << "Initial KE profile (order " << KE_decay_order
                       << ") will extend up to " << KE_decay_height
                       << std::endl;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(klo == bx.smallEnd(2),
        "Vertical decomposition with KE initialization is not allowed.");
    }

    if (state_pert.nComp() > RhoKE_comp) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(KE_0 >= zero,
            "KE_0 must be non-negative.");

        ParallelFor(bx, [=,one_d=one] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Set an initial value for SGS KE
            state_pert(i, j, k, RhoKE_comp) = r_hse(i,j,k) * KE_0;
            // scale initial SGS kinetic energy with height
            if (KE_decay_height > 0) {
                const Real z_sfc = (use_eb) ? zero :
                        fourth * ( z_nd(i  ,j  ,klo) + z_nd(i+1,j  ,klo)
                                 + z_nd(i  ,j+1,klo) + z_nd(i+1,j+1,klo) );
                const Real z_agl = z_cc(i,j,k) - z_sfc;
                state_pert(i, j, k, RhoKE_comp) *= max(
                    std::pow(1 - min(z_agl/KE_decay_height,one_d), KE_decay_order),
                    amrex::Real(1e-12));
            }
        });
    }
