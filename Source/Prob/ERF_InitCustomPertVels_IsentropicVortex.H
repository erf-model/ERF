
    // Parse params
    ParmParse pp("prob");

    amrex::Real M_inf = 0.2;   pp.query("M_inf", M_inf);//freestream Mach number [-]
    amrex::Real T_inf = 300.0; pp.query("T_inf", T_inf);//freestream temperature [K]
    amrex::Real alpha = 0.0;   pp.query("alpha", alpha);//inflow angle, 0 --> x-aligned [rad]
    amrex::Real gamma = Gamma; pp.query("gamma", gamma);//specific heat ratio [-]
    amrex::Real beta = 0.01;   pp.query("beta", beta);//non-dimensional max perturbation strength [-]
    amrex::Real sigma = 1.0;   pp.query("sigma", sigma);//Gaussian standard deviation, i.e., spreading parameter [-]
    amrex::Real R = 2.0;       pp.query("R", R);//characteristic length scale for grid [m]
    amrex::Real xc_frac = 0.5; pp.query("xc", xc_frac);//normalized x-location of vortex center [-]
    amrex::Real yc_frac = 0.5; pp.query("yc", yc_frac);//normalized y-location of vortex center [-]

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
        const Real x =  i        * dx[0]; // face center
        const Real y = (j + 0.5) * dx[1]; // cell center
        const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);

        x_vel_pert(i, j, k) = (M_inf * std::cos(alpha)
                            - (y - yc)/R * Omg) * a_inf;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* dx = geomdata.CellSize();

        // Note that since we only use (x-xc) and (y-yc), we neglect the prob_lo offset in these
        const Real x = (i + 0.5) * dx[0]; // cell center
        const Real y =  j        * dx[1]; // face center
        const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);

        y_vel_pert(i, j, k) = (M_inf * std::sin(alpha)
                            + (x - xc)/R * Omg) * a_inf;
    });
