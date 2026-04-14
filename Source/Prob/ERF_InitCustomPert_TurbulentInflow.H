
    ParmParse pp("prob");

    Real rho_0 = one; pp.query("rho_0", rho_0);
    Real T_0   = amrex::Real(300.0); pp.query("T_0", T_0);

    Real A_0   = one; pp.query("A_0", A_0);
    Real KE_0  = amrex::Real(0.1); pp.query("KE_0", KE_0);

    Real KE_decay_height = -1; pp.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order = 1; pp.query("KE_decay_order", KE_decay_order);

    Real U_0 = zero; pp.query("U_0", U_0);
    Real V_0 = zero; pp.query("V_0", V_0);
    Real W_0 = zero; pp.query("W_0", W_0);

    // random initial perturbations (legacy code)
    Real U_0_Pert_Mag = zero; pp.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = zero; pp.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = zero; pp.query("W_0_Pert_Mag", W_0_Pert_Mag);
    Real T_0_Pert_Mag = zero; pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    bool pert_rhotheta = true; pp.query("pert_rhotheta", pert_rhotheta);

    // divergence-free initial perturbations
    Real pert_deltaU = zero; pp.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = zero; pp.query("pert_deltaV", pert_deltaV);
    Real pert_periods_U = amrex::Real(5.0); pp.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V = amrex::Real(5.0); pp.query("pert_periods_V", pert_periods_V);
    Real pert_ref_height = amrex::Real(100.0); pp.query("pert_ref_height", pert_ref_height);

    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();

    if (KE_decay_height > 0) {
        amrex::Print() << "Initial KE profile (order " << KE_decay_order
                       << ") will extend up to " << KE_decay_height
                       << std::endl;
    }

    if (pert_ref_height > 0) {
        if ((pert_deltaU != zero) || (pert_deltaV != zero)) {
            amrex::Print() << "Adding divergence-free perturbations "
                           << pert_deltaU << " " << pert_deltaV
                           << std::endl;
        }
        if (U_0_Pert_Mag != zero) {
            amrex::Print() << "Adding random x-velocity perturbations" << std::endl;
        }
        if (V_0_Pert_Mag != zero) {
            amrex::Print() << "Adding random y-velocity perturbations" << std::endl;
        }
        if (T_0_Pert_Mag != zero) {
            if (pert_rhotheta) {
                amrex::Print() << "Adding random rho*theta perturbations" << std::endl;
            } else {
                amrex::Print() << "Adding random theta perturbations" << std::endl;
            }
        }
    }

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = z_cc(i,j,k);

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = myhalf * (prob_lo[0] + prob_hi[0]);
        const Real yc = myhalf * (prob_lo[1] + prob_hi[1]);
        const Real zc = myhalf * (prob_lo[2] + prob_hi[2]);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z <= pert_ref_height) && (T_0_Pert_Mag != zero)) {
            Real rand_double = amrex::Random(engine); // Between zero and one
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*two - one)*T_0_Pert_Mag;
            if (!pert_rhotheta) {
                // we're perturbing theta, not rho*theta
                state_pert(i, j, k, RhoTheta_comp) *= r_hse(i,j,k);
            }
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * exp(-amrex::Real(10.)*r*r);

        // Set an initial value for SGS KE
        if (state_pert.nComp() > RhoKE_comp) {
            // Deardorff
            state_pert(i, j, k, RhoKE_comp) = r_hse(i,j,k) * KE_0;
            if (KE_decay_height > 0) {
                // scale initial SGS kinetic energy with height
                state_pert(i, j, k, RhoKE_comp) *= max(
                    std::pow(1 - min(z/KE_decay_height,one), KE_decay_order), Real(1e-12));
            }
        }
    });
