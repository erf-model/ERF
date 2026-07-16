
    ParmParse pp("prob");
    Real rho_0 =   one; pp.query("rho_0", rho_0);
    Real   T_0 = amrex::Real(300.0); pp.query("T_0", T_0);
    Real M_0   =   one; pp.query("M_0", M_0);
    Real V_0   =   one; pp.query("V_0", V_0);

    ParallelFor(bx, [=,myhalf_d=myhalf,one_d=one,Gamma_d=Gamma] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf_d) * dx[0];
        const Real y = prob_lo[1] + (j + myhalf_d) * dx[1];
        const Real z = prob_lo[2] + (k + myhalf_d) * dx[2];

        // Initial potential temperature (actually rho*theta) perturbation
        const Real p = rho_0 * V_0*V_0*
                              ( one_d / (Gamma_d * M_0 * M_0)
                              + (one_d / amrex::Real(16.0)) * (std::cos(2 * x) + std::cos(2 * y)) * (std::cos(2 * z) + 2));
        state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p) - getRhoThetagivenP(p_hse(i,j,k));

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = one_d * rho_0;
    });
