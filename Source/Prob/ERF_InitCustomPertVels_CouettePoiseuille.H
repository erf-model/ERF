    Real U_0   = 0.0;
    Real V_0   = 0.0;
    Real W_0   = 0.0;

    ParmParse pp("prob");

    pp.query("U_0", U_0);
    pp.query("V_0", V_0);
    pp.query("W_0", W_0);

    int prob_type;
    pp.get("prob_type", prob_type);

    AMREX_ALWAYS_ASSERT (prob_type == 1 ||
                         prob_type == 10 || prob_type == 11 ||
                         prob_type == 20 || prob_type == 21);


    Real pert_periods_u = 5.0; pp.query("pert_periods_u", pert_periods_u);
    Real pert_periods_v = 5.0; pp.query("pert_periods_v", pert_periods_v);

    Real pert_delta_u = 0.0; pp.query("pert_delta_u", pert_delta_u);
    Real pert_delta_v = 0.0; pp.query("pert_delta_u", pert_delta_v);

    Real pert_lo = -1e34; pp.query("pert_lo", pert_lo);
    Real pert_hi =  1e34; pp.query("pert_hi", pert_hi);

    auto problo = geomdata.ProbLo();
    auto probhi = geomdata.ProbHi();

    Real aval = pert_periods_u * 2.0 * PI / (probhi[1] - problo[1]);
    Real bval = pert_periods_v * 2.0 * PI / (probhi[0] - problo[0]);

    // Couette flow
    if (prob_type == 1) {

        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const auto *const prob_hi  = geomdata.ProbHi();
            const auto *const dx       = geomdata.CellSize();
            const Real z = (k + 0.5) * dx[2];
            x_vel_pert(i, j, k) = U_0 * z / prob_hi[2];
        });

        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const auto *const prob_hi  = geomdata.ProbHi();
            const auto *const dx       = geomdata.CellSize();
            const Real z = (k + 0.5) * dx[2];
            y_vel_pert(i, j, k) = V_0 * z / prob_hi[2];
        });

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            z_vel_pert(i, j, k) = W_0;
        });

    // Poiseuille flow
    } else if (prob_type == 10 || prob_type == 11) {

        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real* prob_lo = geomdata.ProbLo();
            const Real* dx      = geomdata.CellSize();
            const Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

            // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
            if (prob_type == 10) {
                x_vel_pert(i, j, k) = 1.0 - z_h * z_h;
            } else {
                x_vel_pert(i, j, k) = 0.0;
            }
        });

        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real* prob_lo = geomdata.ProbLo();
            const Real* dx      = geomdata.CellSize();
            const Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

            // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
            if (prob_type == 11) {
               y_vel_pert(i, j, k) = 1.0 - z_h * z_h;
            } else {
               y_vel_pert(i, j, k) = 0.0;
            }
        });

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = 0.0;
        });

    // plane channel flow initialization
    } else if (prob_type == 20 || prob_type == 21) {

        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real* prob_lo = geomdata.ProbLo();
            const Real* prob_hi = geomdata.ProbHi();
            const Real* dx      = geomdata.CellSize();


            const Real z = 0.25*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1) + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) );

            // Normalized wall-normal dist between -1 and 1
            Real y_h = (prob_type == 20) ? 2.0 * (j + 0.5) * dx[1] / (prob_hi[1] - prob_lo[1]) - 1.0
                                                 : 2.0 * (z - prob_lo[2])  / (prob_hi[2] - prob_lo[2]) - 1.0;

            x_vel_pert(i, j, k) = U_0 * (1.0 - y_h * y_h);

            if (pert_delta_u != 0.0) {
                const Real yl = (j + 0.5) * dx[1];
                const Real scaling = std::cos(PI/2.0 * y_h);
                x_vel_pert(i, j, k) += pert_delta_u * scaling * std::cos(aval * yl);
            }
        });

        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real* prob_lo = geomdata.ProbLo();
            const Real* prob_hi = geomdata.ProbHi();
            const Real* dx      = geomdata.CellSize();

            // Normalized wall-normal dist between -1 and 1
            y_vel_pert(i, j, k) = 0.0;

            const Real z = 0.25*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1) + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) );

            if (pert_delta_u != 0.0) {
                const Real xl = (i + 0.5) * dx[0];
                Real y_h = (prob_type == 20) ? 2.0 * (j + 0.5) * dx[1] / (prob_hi[1] - prob_lo[1]) - 1.0
                                                     : 2.0 * (z - prob_lo[2])  / (prob_hi[2] - prob_lo[2]) - 1.0;
                const Real scaling = std::cos(PI/2.0 * y_h);
                y_vel_pert(i, j, k) += pert_delta_v * scaling * std::cos(bval * xl);
            }
        });

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = 0.0;
        });
    } // prob_type
