
    // Parse params
    ParmParse pp("prob");

    amrex::Real p_inf = p_0;   pp.query("p_inf", p_inf); //freestream pressure [Pa]
    amrex::Real T_inf = 300.0; pp.query("T_inf", T_inf);//freestream temperature [K]
    amrex::Real M_inf = 0.2;   pp.query("M_inf", M_inf);//freestream Mach number [-]
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

    // amrex::Print() << "  vortex initialized at ("
    //                << parms.xc << ", "
    //                << parms.yc << ")"
    //                << std::endl;

    amrex::Real inv_gm1 = 1.0 / (gamma - 1.0);

    // amrex::Print() << "  reference pressure = " << parms.p_inf << " Pa" << std::endl;
    // amrex::Print() << "  reference temperature = " << parms.T_inf << " K" << std::endl;
    // amrex::Print() << "  reference potential temperature (not used) = " << parms.T_0 << " K" << std::endl;

    amrex::Real rho_0 = p_inf / (R_d * T_inf);
    // amrex::Print() << "  calculated freestream air density = "
    //                << parms.rho_0 << " kg/m^3"
    //                << std::endl;

    amrex::Real a_inf = std::sqrt(gamma * R_d * T_inf);
    // amrex::Print() << "  calculated speed of sound, a = "
    //                << parms.a_inf << " m/s"
    //                << std::endl;

    // amrex::Print() << "  freestream u/a = "
    //                << parms.M_inf * std::cos(parms.alpha)
    //                << std::endl;
    // amrex::Print() << "  freestream v/a = "
    //                << parms.M_inf * std::sin(parms.alpha)
    //                << std::endl;

    const amrex::Real rdOcp = sc.rdOcp;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real* dx = geomdata.CellSize();

        // Note that since we only use (x-xc) and (y-yc), we neglect the prob_lo offset in these
        const Real x = (i + 0.5) * dx[0]; // cell center
        const Real y = (j + 0.5) * dx[1]; // cell center

        // Calculate perturbation temperature
        const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);
        const Real deltaT = -(gamma - 1.0)/(2.0*sigma*sigma) * Omg*Omg;

        // Set the perturbation density
        const Real rho_norm = std::pow(1.0 + deltaT, inv_gm1);
        state_pert(i, j, k, Rho_comp) = rho_norm * rho_0 - r_hse(i,j,k);

        // Initial _potential_ temperature
        const Real T = (1.0 + deltaT) * T_inf;
        const Real p = std::pow(rho_norm, Gamma) / Gamma  // isentropic relation
                              * rho_0*a_inf*a_inf;
        const Real rho_theta = rho_0 * rho_norm * (T * std::pow(p_0 / p, rdOcp)); // T --> theta
        state_pert(i, j, k, RhoTheta_comp) = rho_theta - getRhoThetagivenP(p_hse(i,j,k)); // Set the perturbation rho*theta

        const Real r2d_xy = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
        state_pert(i, j, k, RhoScalar_comp) = 0.25 * (1.0 + std::cos(PI * std::min(r2d_xy, R) / R));
    });
