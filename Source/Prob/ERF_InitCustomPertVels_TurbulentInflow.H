    ParmParse pp_for_pert_vels("prob");

    Real rho_0 = one; pp_for_pert_vels.query("rho_0", rho_0);
    Real T_0   = amrex::Real(300.0); pp_for_pert_vels.query("T_0", T_0);

    Real A_0   = one; pp_for_pert_vels.query("A_0", A_0);
    Real KE_0  = amrex::Real(0.1); pp_for_pert_vels.query("KE_0", KE_0);

    Real KE_decay_height = -1; pp_for_pert_vels.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order = 1; pp_for_pert_vels.query("KE_decay_order", KE_decay_order);

    Real U_0 = zero; pp_for_pert_vels.query("U_0", U_0);
    Real V_0 = zero; pp_for_pert_vels.query("V_0", V_0);
    Real W_0 = zero; pp_for_pert_vels.query("W_0", W_0);

    // random initial perturbations (legacy code)
    Real U_0_Pert_Mag = zero; pp_for_pert_vels.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = zero; pp_for_pert_vels.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = zero; pp_for_pert_vels.query("W_0_Pert_Mag", W_0_Pert_Mag);
    Real T_0_Pert_Mag = zero; pp_for_pert_vels.query("T_0_Pert_Mag", T_0_Pert_Mag);
    bool pert_rhotheta = true; pp_for_pert_vels.query("pert_rhotheta", pert_rhotheta);

    // divergence-free initial perturbations
    Real pert_deltaU = zero; pp_for_pert_vels.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = zero; pp_for_pert_vels.query("pert_deltaV", pert_deltaV);
    Real pert_periods_U = amrex::Real(5.0); pp_for_pert_vels.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V = amrex::Real(5.0); pp_for_pert_vels.query("pert_periods_V", pert_periods_V);
    Real pert_ref_height = amrex::Real(100.0); pp_for_pert_vels.query("pert_ref_height", pert_ref_height);

    // helper vars
    const Real prob_lo_x = geomdata.ProbLo(0);
    const Real prob_lo_y = geomdata.ProbLo(1);
    const Real prob_hi_x = geomdata.ProbHi(0);
    const Real prob_hi_y = geomdata.ProbHi(1);
    Real aval = pert_periods_U * two * PI / (prob_hi_y - prob_lo_y);
    Real bval = pert_periods_V * two * PI / (prob_hi_x - prob_lo_x);
    Real ufac = pert_deltaU * std::exp(myhalf) / pert_ref_height;
    Real vfac = pert_deltaV * std::exp(myhalf) / pert_ref_height;

    const bool use_terrain  = (SolverChoice::terrain_type != TerrainType::None);

    // Set the x-velocity
    ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = use_terrain ? fourth*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                                          + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) )
                                   : prob_lo[2] + (k + myhalf) * dx[2];

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
            const Real yl = y - prob_lo[1];
            const Real zl = z / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);
            x_vel_pert(i, j, k) += ufac * damp * z * std::cos(aval * yl);
        }
    });

  // Set the y-velocity
  ParallelForRNG(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + myhalf) * dx[0];
    const Real z = use_terrain ? fourth*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                                      + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) )
                               : prob_lo[2] + (k + myhalf) * dx[2];

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
        const Real xl = x - prob_lo[0];
        const Real zl = z / pert_ref_height;
        const Real damp = std::exp(-myhalf * zl * zl);
        y_vel_pert(i, j, k) += vfac * damp * z * std::cos(bval * xl);
    }
  });

  // Set the z-velocity
  ParallelForRNG(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
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
