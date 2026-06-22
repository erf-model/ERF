
    // Parse params
    ParmParse pp_for_pert_vels("prob");

    amrex::Real M_inf   = amrex::Real(0.2);   pp_for_pert_vels.query("M_inf", M_inf);//freestream Mach number [-]
    amrex::Real T_inf   = amrex::Real(300.0); pp_for_pert_vels.query("T_inf", T_inf);//freestream temperature [K]
    amrex::Real alpha   = zero;               pp_for_pert_vels.query("alpha", alpha);//inflow angle, 0 --> x-aligned [rad]
    amrex::Real gamma   = Gamma;              pp_for_pert_vels.query("gamma", gamma);//specific heat ratio [-]
    amrex::Real beta    = amrex::Real(0.01);  pp_for_pert_vels.query("beta", beta);//non-dimensional max perturbation strength [-]
    amrex::Real sigma   = one;                pp_for_pert_vels.query("sigma", sigma);//Gaussian standard deviation, i.e., spreading parameter [-]
    amrex::Real R       = two;                pp_for_pert_vels.query("R", R);//characteristic length scale for grid [m]
    amrex::Real xc_frac = myhalf;             pp_for_pert_vels.query("xc", xc_frac);//normalized x-location of vortex center [-]
    amrex::Real yc_frac = myhalf;                         pp_for_pert_vels.query("yc", yc_frac);//normalized y-location of vortex center [-]

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();

    amrex::Real xc = xc_frac * (prob_hi[0] - prob_lo[0]);
    amrex::Real yc = yc_frac * (prob_hi[1] - prob_lo[1]);

    amrex::Real a_inf = std::sqrt(gamma * R_d * T_inf);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* dx = geomdata.CellSize();

        // Note that since we only use (x-xc) and (y-yc), we neglect the prob_lo offset in these
        const Real x =  i           * dx[0]; // face center
        const Real y = (j + myhalf) * dx[1]; // cell center
        const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);

        x_vel_pert(i, j, k) = (M_inf * std::cos(alpha)
                            - (y - yc)/R * Omg) * a_inf;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* dx = geomdata.CellSize();

        // Note that since we only use (x-xc) and (y-yc), we neglect the prob_lo offset in these
        const Real x = (i + myhalf) * dx[0]; // cell center
        const Real y =  j           * dx[1]; // face center
        const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);

        y_vel_pert(i, j, k) = (M_inf * std::sin(alpha)
                            + (x - xc)/R * Omg) * a_inf;
    });
