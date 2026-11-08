
    // On a mesh other than ConstantDz, heights come from z_nd and the vertical levels,
    //    not from geomdata's prob_lo[2], prob_hi[2] and CellSize(2)
    const bool use_z_phys = use_physical_z();

    amrex::Gpu::streamSynchronize();
    ParmParse pp_for_pert_vels("prob");
    Real U_0   = zero;  pp_for_pert_vels.query("U_0", U_0);
    Real V_0   = zero;  pp_for_pert_vels.query("V_0", V_0);
    Real W_0   = zero;  pp_for_pert_vels.query("W_0", W_0);

    Real z0    = amrex::Real(0.1);  pp_for_pert_vels.query("z0", z0);
    Real zRef  = amrex::Real(80.0); pp_for_pert_vels.query("zRef", zRef);
    Real uRef  = zero;  pp_for_pert_vels.query("uRef", uRef);

    int prob_type = -1; pp_for_pert_vels.query("prob_type", prob_type);

    // An optional smooth test flow added to the velocity above, divergence-free in the
    //    continuum and varying in all three directions: the x-z roll of the streamfunction
    //    psi = A (Lx / 2 pi) sin(2 pi x / Lx) sin(pi z / H), whose w vanishes at the bottom and
    //    top, and the horizontal shears A sin(2 pi y / Ly) in u and A sin(2 pi x / Lx) in v
    //    (periodic in x and y).  The initial projection makes it discretely divergence-free.
    Real test_amp = zero; pp_for_pert_vels.query("test_flow_amp", test_amp);
    const Real tf_xlo = geomdata.ProbLo()[0], tf_ylo = geomdata.ProbLo()[1];
    const Real tf_Lx  = geomdata.ProbHi()[0] - tf_xlo;
    const Real tf_Ly  = geomdata.ProbHi()[1] - tf_ylo;
    const Real tf_zlo = geomdata.ProbLo()[2];
    const Real tf_H   = geomdata.ProbHi()[2] - tf_zlo;
    const Real tf_kx  = Real(2.0) * PI / tf_Lx;
    const Real tf_ky  = Real(2.0) * PI / tf_Ly;
    const Real tf_kz  = PI / tf_H;

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real*      dx = geomdata.CellSize();
        const Real        z = use_z_phys ? erf_prob_heights::z_xface(i,j,k,z_nd) : prob_lo[2] + (k + myhalf) * dx[2];

        // Set the x-velocity
        if (uRef != zero) {
            x_vel_pert(i, j, k) = U_0 + uRef *
                                  std::log((z + z0)/z0) / std::log((zRef +z0)/z0);
        } else {
            x_vel_pert(i, j, k) = U_0;
        }
        if (test_amp != zero) {
            const Real x = prob_lo[0] + i * dx[0];
            const Real y = prob_lo[1] + (j + myhalf) * dx[1];
            x_vel_pert(i, j, k) += test_amp * (tf_kz / tf_kx) * std::sin(tf_kx * (x - tf_xlo))
                                                              * std::cos(tf_kz * (z - tf_zlo))
                                 + test_amp * std::sin(tf_ky * (y - tf_ylo));
        }
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = V_0;
        if (test_amp != zero) {
            const Real* prob_lo = geomdata.ProbLo();
            const Real*      dx = geomdata.CellSize();
            const Real x = prob_lo[0] + (i + myhalf) * dx[0];
            y_vel_pert(i, j, k) += test_amp * std::sin(tf_kx * (x - tf_xlo));
        }
    });

    // Set the z-velocity
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = W_0;
        if (test_amp != zero) {
            const Real* prob_lo = geomdata.ProbLo();
            const Real*      dx = geomdata.CellSize();
            const Real x = prob_lo[0] + (i + myhalf) * dx[0];
            const Real z = use_z_phys ? erf_prob_heights::z_zface(i,j,k,z_nd) : prob_lo[2] + k * dx[2];
            z_vel_pert(i, j, k) -= test_amp * std::cos(tf_kx * (x - tf_xlo)) * std::sin(tf_kz * (z - tf_zlo));
        }
    });
