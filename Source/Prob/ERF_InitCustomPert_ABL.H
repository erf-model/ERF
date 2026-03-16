
    ParmParse pp("prob");

    Real rho_0 = 0.0; pp.query("rho_0", rho_0);
    Real T_0   = 0.0; pp.query("T_0", T_0);
    Real A_0   = 1.0; pp.query("A_0", A_0);
    Real KE_0  = 0.1; pp.query("KE_0", KE_0);
    Real rhoKE_0 = -1; pp.query("rhoKE_0", rhoKE_0);

    Real KE_decay_height = -1; pp.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order  =  1; pp.query("KE_decay_order", KE_decay_order);

    // random initial perturbations (legacy code)
    Real T_0_Pert_Mag = 0.0; pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    bool pert_rhotheta = true; pp.query("pert_rhotheta", pert_rhotheta);

    // divergence-free initial perturbations
    Real pert_deltaU = 0.0; pp.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = 0.0; pp.query("pert_deltaV", pert_deltaV);
    Real pert_ref_height = 100.0; pp.query("pert_ref_height", pert_ref_height);

    // Capture by value for GPU
    const Real dx = geomdata.CellSize(0);
    const Real dy = geomdata.CellSize(1);
    const Real prob_lo_x = geomdata.ProbLo(0);
    const Real prob_lo_y = geomdata.ProbLo(1);
    const Real prob_lo_z = geomdata.ProbLo(2);
    const Real prob_hi_x = geomdata.ProbHi(0);
    const Real prob_hi_y = geomdata.ProbHi(1);
    const Real prob_hi_z = geomdata.ProbHi(2);

    // Define a point (xc,yc,zc) at the center of the domain
    const Real xc = 0.5 * (prob_lo_x + prob_hi_x);
    const Real yc = 0.5 * (prob_lo_y + prob_hi_y);
    const Real zc = 0.5 * (prob_lo_z + prob_hi_z);

    if (KE_decay_height > 0) {
        amrex::Print() << "Initial KE profile (order " << KE_decay_order
                       << ") will extend up to " << KE_decay_height
                       << std::endl;
    }

#ifndef AMREX_USE_GPU
    if (pert_ref_height > 0) {
        if (T_0_Pert_Mag != 0.0) {
            if (pert_rhotheta) {
                amrex::Print() << "Adding random rho*theta perturbations" << std::endl;
            } else {
                amrex::Print() << "Adding random theta perturbations" << std::endl;
            }
        }
    }
#endif

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real x = prob_lo_x + (i + 0.5) * dx;
        const Real y = prob_lo_y + (j + 0.5) * dy;
        const Real z = z_cc(i,j,k);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z <= pert_ref_height) && (T_0_Pert_Mag != 0.0)) {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*T_0_Pert_Mag;
            if (!pert_rhotheta) {
                // we're perturbing theta, not rho*theta
                state_pert(i, j, k, RhoTheta_comp) *= r_hse(i,j,k);
            }
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * exp(-10.*r*r);

        // Set an initial value for SGS KE
        if (state_pert.nComp() > RhoKE_comp) {
            // Deardorff
            if (rhoKE_0 > 0) {
                state_pert(i, j, k, RhoKE_comp) = rhoKE_0;
            } else {
                state_pert(i, j, k, RhoKE_comp) = r_hse(i,j,k) * KE_0;
            }
            if (KE_decay_height > 0) {
                // scale initial SGS kinetic energy with height
                state_pert(i, j, k, RhoKE_comp) *= amrex::max(
                    std::pow(1 - amrex::min(z/KE_decay_height,1.0), KE_decay_order),
                    1e-12);
            }
        }
    });

