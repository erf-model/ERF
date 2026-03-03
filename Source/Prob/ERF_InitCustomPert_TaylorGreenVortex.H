
    ParmParse pp("prob");
    Real rho_0 =   1.0; pp.query("rho_0", rho_0);
    Real   T_0 = 300.0; pp.query("T_0", T_0);
    Real M_0   =   1.0; pp.query("M_0", M_0);
    Real V_0   =   1.0; pp.query("V_0", V_0);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Initial potential temperature (actually rho*theta) perturbation
        const Real p = rho_0 * V_0*V_0*
                              ( 1.0 / (Gamma * M_0 * M_0)
                              + (1.0 / 16.0) * (cos(2 * x) + cos(2 * y)) * (cos(2 * z) + 2));
        state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p) - getRhoThetagivenP(p_hse(i,j,k));

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 1.0 * rho_0;
    });
