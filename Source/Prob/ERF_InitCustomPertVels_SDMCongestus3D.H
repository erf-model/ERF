
    ParmParse pp_prob("prob");

    Real U_0 = zero; pp_prob.query("U_0", U_0);
    Real V_0 = zero; pp_prob.query("V_0", V_0);
    Real W_0 = zero; pp_prob.query("W_0", W_0);

    // random initial perturbations (legacy code)
    Real U_0_Pert_Mag = zero; pp_prob.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = zero; pp_prob.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = zero; pp_prob.query("W_0_Pert_Mag", W_0_Pert_Mag);

    Real pert_deltaU = zero; pp_prob.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = zero; pp_prob.query("pert_deltaV", pert_deltaV);

    Real pert_periods_U = amrex::Real(5.0); pp_prob.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V = amrex::Real(5.0); pp_prob.query("pert_periods_V", pert_periods_V);

    Real pert_ref_height = amrex::Real(100.0); pp_prob.query("pert_ref_height", pert_ref_height);

    Real prob_lo_x = geomdata.ProbLo(0);
    Real prob_lo_y = geomdata.ProbLo(1);
    Real prob_hi_x = geomdata.ProbHi(0);
    Real prob_hi_y = geomdata.ProbHi(1);

    Real aval = pert_periods_U  * two * PI / (prob_hi_y - prob_lo_y);
    Real bval = pert_periods_V  * two * PI / (prob_hi_x - prob_lo_x);

    Real ufac = pert_deltaU * std::exp(myhalf) / pert_ref_height;
    Real vfac = pert_deltaV * std::exp(myhalf) / pert_ref_height;

    // Set the x-velocity
    ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = prob_lo[2] + (k + myhalf) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = U_0;
        if ((z <= pert_ref_height) && (U_0_Pert_Mag != zero))
        {
            Real rand_double = amrex::Random(engine); // Between zero and one
            Real x_vel_prime = (rand_double*two - one)*U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
        if (pert_deltaU != zero)
        {
            const amrex::Real yl = y - prob_lo[1];
            const amrex::Real zl = z / pert_ref_height;
            const amrex::Real damp = std::exp(-myhalf * zl * zl);
            x_vel_pert(i, j, k) += ufac * damp * z * std::cos(aval * yl);
        }
    });

    // Set the y-velocity
    ParallelForRNG(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real z = prob_lo[2] + (k + myhalf) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = V_0;
        if ((z <= pert_ref_height) && (V_0_Pert_Mag != zero))
        {
            Real rand_double = amrex::Random(engine); // Between zero and one
            Real y_vel_prime = (rand_double*two - one)*V_0_Pert_Mag;
            y_vel_pert(i, j, k) += y_vel_prime;
        }
        if (pert_deltaV != zero)
        {
            const amrex::Real xl = x - prob_lo[0];
            const amrex::Real zl = z / pert_ref_height;
            const amrex::Real damp = std::exp(-myhalf * zl * zl);
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
            z_vel_pert(i, j, k) = zero;
        }
        else if (W_0_Pert_Mag != zero)
        {
            Real rand_double = amrex::Random(engine); // Between zero and one
            Real z_vel_prime = (rand_double*two - one)*W_0_Pert_Mag;
            z_vel_pert(i, j, k) = W_0 + z_vel_prime;
        }
    });
