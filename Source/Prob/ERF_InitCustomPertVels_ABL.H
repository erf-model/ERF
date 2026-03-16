
    ParmParse pp_prob("prob");

    Real U_0 = 0.0; pp_prob.query("U_0", U_0);
    Real V_0 = 0.0; pp_prob.query("V_0", V_0);
    Real W_0 = 0.0; pp_prob.query("W_0", W_0);

    // random initial perturbations (legacy code)
    Real U_0_Pert_Mag = 0.0; pp_prob.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = 0.0; pp_prob.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = 0.0; pp_prob.query("W_0_Pert_Mag", W_0_Pert_Mag);

    // divergence-free initial perturbations
    Real pert_deltaU = 0.0; pp_prob.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = 0.0; pp_prob.query("pert_deltaV", pert_deltaV);
    Real pert_periods_U = 5.0; pp_prob.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V = 5.0; pp_prob.query("pert_periods_V", pert_periods_V);
    Real pert_ref_height = 100.0; pp_prob.query("pert_ref_height", pert_ref_height);

    auto problo = geomdata.ProbLo();
    auto probhi = geomdata.ProbHi();

    Real aval = pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
    Real bval = pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
    Real ufac = pert_deltaU * std::exp(0.5) / pert_ref_height;
    Real vfac = pert_deltaV * std::exp(0.5) / pert_ref_height;

    if (pert_ref_height > 0) {
        if ((pert_deltaU != 0.0) || (pert_deltaV != 0.0)) {
            amrex::Print() << "Adding divergence-free perturbations "
                           << pert_deltaU << " " << pert_deltaV
                           << std::endl;
        }
        if (U_0_Pert_Mag != 0.0) {
            amrex::Print() << "Adding random x-velocity perturbations" << std::endl;
        }
        if (V_0_Pert_Mag != 0.0) {
            amrex::Print() << "Adding random y-velocity perturbations" << std::endl;
        }
    }

    // Set the x-velocity
    ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = (z_nd) ? 0.25*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                                     + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) )
                                   : prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = U_0;
        if ((z <= pert_ref_height) && (U_0_Pert_Mag != 0.0))
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real x_vel_prime = (rand_double*2.0 - 1.0)*U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
        if (pert_deltaU != 0.0)
        {
            const amrex::Real yl = y - prob_lo[1];
            const amrex::Real zl = z / pert_ref_height;
            const amrex::Real damp = std::exp(-0.5 * zl * zl);
            x_vel_pert(i, j, k) += ufac * damp * z * std::cos(aval * yl);
        }
    });

    // Set the y-velocity
    ParallelForRNG(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = (z_nd) ? 0.25*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                                     + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) )
                                   : prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = V_0;
        if ((z <= pert_ref_height) && (V_0_Pert_Mag != 0.0))
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real y_vel_prime = (rand_double*2.0 - 1.0)*V_0_Pert_Mag;
            y_vel_pert(i, j, k) += y_vel_prime;
        }
        if (pert_deltaV != 0.0)
        {
            const amrex::Real xl = x - prob_lo[0];
            const amrex::Real zl = z / pert_ref_height;
            const amrex::Real damp = std::exp(-0.5 * zl * zl);
            y_vel_pert(i, j, k) += vfac * damp * z * std::cos(bval * xl);
        }
    });

    // Set the z-velocity
    ParallelForRNG(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const int dom_lo_z = geomdata.Domain().smallEnd()[2];
        const int dom_hi_z = geomdata.Domain().bigEnd()[2];

        // Set the z-velocity
        if (k == dom_lo_z || k == dom_hi_z+1)
        {
            z_vel_pert(i, j, k) = 0.0;
        }
        else if (W_0_Pert_Mag != 0.0)
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real z_vel_prime = (rand_double*2.0 - 1.0)*W_0_Pert_Mag;
            z_vel_pert(i, j, k) = W_0 + z_vel_prime;
        }
    });
