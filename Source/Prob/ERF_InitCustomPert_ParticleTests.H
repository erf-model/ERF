  ParmParse pp("prob");
  Real x_c =    zero; pp.query("x_c", x_c);
  Real z_c = amrex::Real(3000.0); pp.query("z_c", z_c);
  Real x_r = amrex::Real(4000.0); pp.query("x_r", x_r);
  Real z_r = amrex::Real(2000.0); pp.query("z_r", z_r);
  Real T_pert = -amrex::Real(15.0); pp.query("T_pert", T_pert);

  // Overridden physical constants
  Real C_p = amrex::Real(1004.0);

  ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + myhalf) * dx[0];
    const Real z = z_cc(i,j,k);

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

    Real L = std::sqrt(
        std::pow((x - x_c)/x_r, 2) +
        std::pow((z - z_c)/z_r, 2)
    );
    if (L <= one) {
        Real dT = T_pert * (std::cos(PI*L) + one)/two;

        // Note: dT is a temperature perturbation, theta_perturbed is base state + perturbation in theta
        Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), R_d/C_p);

        // This version perturbs rho but not p
        state_pert(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
    }
  });
