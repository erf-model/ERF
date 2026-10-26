
    // On a mesh other than ConstantDz, heights come from z_cc and the vertical levels,
    //    not from geomdata's prob_lo[2], prob_hi[2] and CellSize(2)
    const bool use_z_phys = use_physical_z();
    const Real phys_z_lo  = phys_zlo(geomdata, lev);
    const Real phys_z_hi  = phys_zhi(geomdata, lev);

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

    // prob_type 11, 12 and 13 divide by the radius
    if (prob_type == 11 || prob_type == 12 || prob_type == 13) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rad_0 > zero,
            "prob.rad_0 must be > 0 for prob.prob_type = 11, 12 and 13");
    }
    ParallelFor(bx, [=,one_d=one] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = use_z_phys ? z_cc(i,j,k) : prob_lo[2] + (k + myhalf) * dx[2];

        // Define a point (xc,yc,zc) at the given fractions of the domain extent
        const Real xc = prob_lo[0] + xc_frac * (prob_hi[0] - prob_lo[0]);
        const Real yc = prob_lo[1] + yc_frac * (prob_hi[1] - prob_lo[1]);
        const Real zc = phys_z_lo  + zc_frac * (phys_z_hi  - phys_z_lo);

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
            state_pert(i, j, k, RhoScalar_comp) = A_0 * fourth * (one_d + std::cos(PI * std::min(r2d_xy, r0) / r0));
        } else if (prob_type == 12) {
            state_pert(i, j, k, RhoScalar_comp) = A_0 * fourth * (one_d + std::cos(PI * std::min(r2d_xz, r0) / r0));
        } else if (prob_type == 13) {
            const Real r0_z = rad_0 * (phys_z_hi - phys_z_lo);
            //ellipse for mapfac shear validation
            //
            // NOTE: unlike r2d_xy and r2d_xz above, this radius is already normalized by the
            //       semi-axes, so it is dimensionless and equals one *on* the ellipse.  We
            //       therefore clip it at one rather than at a length -- dividing it by r0_z
            //       would scale the blob as 1/r0_z instead of giving the intended ellipse.
            //       With r0_z = r0 this reduces exactly to prob_type 12.
            //
            const Real r2d_xz_ell = std::sqrt((x-xc)*(x-xc)/(r0*r0) + (z-zc)*(z-zc)/(r0_z*r0_z));
            state_pert(i, j, k, RhoScalar_comp) = A_0 * fourth * (one_d + std::cos(PI * std::min(r2d_xz_ell, one_d)));
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
