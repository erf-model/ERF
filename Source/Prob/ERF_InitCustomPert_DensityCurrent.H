    // Parse params
    ParmParse pp("prob");
    Real x_c =      0.; pp.query("x_c", x_c);
    Real z_c =   3000.; pp.query("z_c", z_c);
    Real x_r =   4000.; pp.query("x_r", x_r);
    Real z_r =   2000.; pp.query("z_r", z_r);
    Real T_pert = -15.; pp.query("T_pert", T_pert);

    const bool const_rho    = (sc.fixed_density[lev] == 1);

    const Real rdOcp = sc.rdOcp;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo = geomdata.ProbLo();
        const auto dx      = geomdata.CellSize();

        const Real x = ( prob_lo[0] + (i + 0.5) * dx[0] ) / mf_m(i,j,0);
        const Real z = z_cc(i,j,k);

        Real L = std::sqrt(
            std::pow((x - x_c)/x_r, 2) +
            std::pow((z - z_c)/z_r, 2));

        if (L <= 1.0)
        {
            Real dT = T_pert * (std::cos(PI*L) + 1.0)/2.0;
            Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

            // Note: dT is a perturbation in temperature, theta_perturbed is base state + perturbation
            Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), rdOcp);
            Real theta_0         = (Tbar_hse   )*std::pow(p_0/p_hse(i,j,k), rdOcp);

            if (const_rho) {
                state_pert(i, j, k, RhoTheta_comp) = r_hse(i,j,k) * (theta_perturbed  - theta_0);
            } else {
                state_pert(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
            }
        }
    });
