    ParmParse pp_prob("prob");
    Real Ampl = zero;        pp_prob.query("Ampl", Ampl);
    Real wavelength = one; pp_prob.get("wavelength", wavelength);

    Real H           = geomdata.ProbHi()[2];
    Real kp          = two * PI / wavelength;
    Real g           = CONST_GRAV;
    Real omega       = std::sqrt(g * kp);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();

        const Real x = prob_lo[0] + i * dx[0];
        Real z = fourth * (z_nd(i,j,k) + z_nd(i,j+1,k) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1));

        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;

        Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);

        x_vel_pert(i, j, k) = -Ampl * omega * fac * std::sin(kp * x);
    });

    // Set the z-velocity from impenetrable condition
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto dx         = geomdata.CellSize();

        Real x   = (i + myhalf) * dx[0];
        Real z   = fourth * ( z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k) );
        Real z_base = Ampl * std::sin(kp * x);

        z -= z_base;
        Real fac = std::sinh( kp * (z - H) ) / std::sinh(kp * H);

        z_vel_pert(i, j, k) = Ampl * omega * fac * std::cos(kp * x);
    });
