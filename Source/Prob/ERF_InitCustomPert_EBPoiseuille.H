    ParmParse pp("prob");
    Real A_0 = one  ; pp.query("A_0", A_0);
    Real B_0 = zero  ; pp.query("B_0", B_0);
    Real rad_0 = zero; pp.query("rad_0", rad_0);
    Real rho_0 = one; pp.query("rho_0", rho_0);

    // Location of "center" of scalar (multiplies domain length)
    Real xc_frac = myhalf; pp.query("xc_frac", xc_frac);
    Real yc_frac = myhalf; pp.query("yc_frac", yc_frac);
    Real zc_frac = myhalf; pp.query("zc_frac", zc_frac);

    int prob_type = -1;
    pp.query("prob_type", prob_type);

    Real xradius = amrex::Real(10.0); pp.query("xradius", xradius); // x-radius of scalar bubble
    Real zradius = amrex::Real(10.0); pp.query("zradius", zradius); // z-radius of scalar bubble

    AMREX_ALWAYS_ASSERT(sc.terrain_type == TerrainType::EB);

    // Set the state_pert
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
        const Real r2d_xz = std::sqrt((x-xc)*(x-xc) + (z-zc)*(z-zc));
        const Real r2d_xz_nd   = std::sqrt((x-xc)*(x-xc)/xradius/xradius
                                          +(z-zc)*(z-zc)/zradius/zradius);

        if (prob_type == 10)
        {
            // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain,
            //            + B_0*sin(x)
            // state_pert(i, j, k, RhoScalar_comp) = A_0 * std::exp(-amrex::Real(10.)*r3d*r3d) + B_0*std::sin(x);
            state_pert(i, j, k, RhoScalar_comp) = A_0 * std::exp(-amrex::Real(0.1)*r2d_xz*r2d_xz) + B_0*std::sin(x);

        } else if (prob_type == 11) {
            if (r2d_xz_nd < one)
            {
                state_pert(i, j, k, RhoScalar_comp) = myhalf * A_0 * (one + std::cos(PI*r2d_xz_nd));
            }
        } else {
            // Set scalar = A_0 in a ball of radius r0 and 0 elsewhere
            if (r3d < r0) {
                state_pert(i, j, k, RhoScalar_comp) = A_0;
            }
        }

        state_pert(i, j, k, RhoScalar_comp) *= rho_0;
    });
