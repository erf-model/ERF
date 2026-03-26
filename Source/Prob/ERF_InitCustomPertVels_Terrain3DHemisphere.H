
    amrex::Gpu::streamSynchronize();
    ParmParse pp_for_pert_vels("prob");
    Real U_0 = zero; pp_for_pert_vels.query("U_0", U_0);
    Real V_0 = zero; pp_for_pert_vels.query("V_0", V_0);
    Real U_0_Pert_Mag = zero; pp_for_pert_vels.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = zero; pp_for_pert_vels.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = zero; pp_for_pert_vels.query("W_0_Pert_Mag", W_0_Pert_Mag);
    Real pert_deltaU = zero; pp_for_pert_vels.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = zero; pp_for_pert_vels.query("pert_deltaV", pert_deltaV);
    Real pert_periods_U  = amrex::Real(5.0); pp_for_pert_vels.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V  = amrex::Real(5.0); pp_for_pert_vels.query("pert_periods_V", pert_periods_V);;
    Real pert_ref_height = one; pp_for_pert_vels.query("pert_ref_height", pert_ref_height);

    pp_for_pert_vels.query("pert_ref_height", pert_ref_height);

    // Real aval = pert_periods_U * two * PI / (probhi[1] - problo[1]);
    // Real ufac = pert_deltaU * std::exp(myhalf) / pert_ref_height;

    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real bval = pert_periods_V * two * PI / (probhi[0] - problo[0]);
    Real vfac = pert_deltaV * std::exp(myhalf) / pert_ref_height;

    const bool use_terrain  = (SolverChoice::terrain_type != TerrainType::None);

    // Set the x-velocity
    ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real z = use_terrain ? fourth*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                                          + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) )
            : prob_lo[2] + (k + myhalf) * dx[2];

        x_vel_pert(i, j, k) = zero;
        if ((z <= pert_ref_height) && (U_0_Pert_Mag != zero))
        {
            Real rand_double = amrex::Random(engine); // Between zero and one
            Real x_vel_prime = (rand_double*two - one)*U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
    });

    // Set the y-velocity
    ParallelForRNG(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real z = use_terrain ? fourth*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                                          + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) )
            : prob_lo[2] + (k + myhalf) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = zero;
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

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = one / dx[0];
    dxInv[1] = one / dx[1];
    dxInv[2] = one / dx[2];

    // Set the z-velocity from impenetrable condition
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
            z_vel_pert(i, j, k) = z_vel_prime;
        }
    });
