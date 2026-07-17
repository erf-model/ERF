
    ParmParse pp("prob");
    Real rho_0 = one;  pp.query("rho_0", rho_0);
    Real A_0   = one;  pp.query("A_0", A_0);
    Real B_0   = zero;  pp.query("B_0", B_0);
    Real rad_0 = zero;  pp.query("rad_0", rad_0);

    // Location of "center" of scalar (multiplies domain length)
    Real xc_frac = myhalf; pp.query("xc_frac", xc_frac);
    Real yc_frac = myhalf; pp.query("yc_frac", yc_frac);
    Real zc_frac = myhalf; pp.query("zc_frac", zc_frac);

    int prob_type = -1; pp.query("prob_type", prob_type);
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = prob_lo[2] + (k + myhalf) * dx[2];

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = xc_frac * (prob_lo[0] + prob_hi[0]);
        const Real yc = yc_frac * (prob_lo[1] + prob_hi[1]);
        const Real zc = zc_frac * (prob_lo[2] + prob_hi[2]);

        // Define ellipse parameters
        const Real r0   = rad_0 * (prob_hi[0] - prob_lo[0]);

        const Real r3d    = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
        const Real r2d_xy = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
        const Real r2d_xz = std::sqrt((x-xc)*(x-xc) + (z-zc)*(z-zc));

        if (prob_type == 10)
        {
            // Set scalar = A_0*std::exp(-10r^2), where r is distance from center of domain,
            //            + B_0*std::sin(x)
            state_pert(i, j, k, RhoScalar_comp) = A_0 * std::exp(-amrex::Real(10.)*r3d*r3d) + B_0*std::sin(x);

        } else if (prob_type == 11) {
            state_pert(i, j, k, RhoScalar_comp) = A_0 * fourth * (one + std::cos(PI * std::min(r2d_xy, r0) / r0));
        } else if (prob_type == 12) {
            state_pert(i, j, k, RhoScalar_comp) = A_0 * fourth * (one + std::cos(PI * std::min(r2d_xz, r0) / r0));
        } else if (prob_type == 13) {
            const Real r0_z = rad_0 * (prob_hi[2] - prob_lo[2]);
            //ellipse for mapfac shear validation
            const Real r2d_xz_ell = std::sqrt((x-xc)*(x-xc)/(r0*r0) + (z-zc)*(z-zc)/(r0_z*r0_z));
            state_pert(i, j, k, RhoScalar_comp) = A_0 * fourth * (one + std::cos(PI * std::min(r2d_xz_ell, r0_z) / r0_z));
        } else if (prob_type == 14) {
            state_pert(i, j, k, RhoScalar_comp) = std::cos(PI*x);
        } else {
            // Set scalar = A_0 in a ball of radius r0 and 0 elsewhere
            if (r3d < r0) {
               state_pert(i, j, k, RhoScalar_comp) = A_0;
            } else {
               state_pert(i, j, k, RhoScalar_comp) = zero;
            }
        }

        state_pert(i, j, k, RhoScalar_comp) *= rho_0;
  });
