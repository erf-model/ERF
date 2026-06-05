
    // Parse params
    ParmParse pp("prob");
    Real Ampl = zero       ; pp.query("Ampl", Ampl);
    Real wavelength = one; pp.get("wavelength", wavelength);

    Real H           = geomdata.ProbHi()[2];
    Real kp          = two * PI / wavelength;
    Real g           = CONST_GRAV;
    Real omega       = std::sqrt(g * kp);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];

        Real z = z_cc(i,j,k);
        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;

        Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
        Real p_prime = -(Ampl * omega * omega / kp) * fac * std::sin(kp * x) * r_hse(i,j,k);
        Real p_total = p_hse(i,j,k) + p_prime;

        // Define (rho theta) given pprime
        state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p_total) - getRhoThetagivenP(p_hse(i,j,k));

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = state_pert(i,j,k,Rho_comp);
    });
