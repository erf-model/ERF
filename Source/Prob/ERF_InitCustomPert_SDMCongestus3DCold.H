
    ParmParse pp("prob");

    amrex::Real rho_0 = 0.0; pp.query("rho_0", rho_0);
    amrex::Real T_0   = 0.0; pp.query("T_0"  , T_0);
    amrex::Real A_0   = 1.0; pp.query("A_0"  , A_0);
    amrex::Real KE_0  = 0.1; pp.query("KE_0" , KE_0);

    Real  T_0_Pert_Mag   =   0.0; pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    Real  Q_0_Pert_Mag   =   0.0; pp.query("Q_0_Pert_Mag", Q_0_Pert_Mag);
    Real pert_ref_height = 100.0; pp.query("pert_ref_height", pert_ref_height);

    Real pert_deltaT     =   0.0; pp.query("pert_deltaT", pert_deltaT);
    Real pert_deltaQV    =   0.0; pp.query("pert_deltaQV", pert_deltaQV);
    Real pert_periods_T  =   0.0; pp.query("pert_periods_T", pert_periods_T);
    Real pert_periods_QV =   0.0; pp.query("pert_periods_QV", pert_periods_QV);

    auto problo_arr = geomdata.ProbLo();
    auto probhi_arr = geomdata.ProbHi();

    Real cval  = pert_periods_T  * 2.0 * PI / (probhi_arr[0] - problo_arr[0]);
    Real dval  = pert_periods_QV * 2.0 * PI / (probhi_arr[0] - problo_arr[0]);
    Real tfac  = pert_deltaT  * std::exp(0.5) / pert_ref_height;
    Real qvfac = pert_deltaQV * std::exp(0.5) / pert_ref_height;

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const Real rdOcp   = sc.rdOcp;

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
        const Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
        const Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z <= pert_ref_height) && (T_0_Pert_Mag != 0.0)) {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*T_0_Pert_Mag;
        }
        if ((z <= pert_ref_height) && (Q_0_Pert_Mag != 0.0)) {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            state_pert(i, j, k, RhoQ1_comp) = (rand_double*2.0 - 1.0)*Q_0_Pert_Mag;
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * exp(-10.*r*r);

        // Set an initial value for KE
        state_pert(i, j, k, RhoKE_comp) = KE_0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }

        // sinusoidal variation of theta and qv for regression tests
        if (pert_deltaT != 0.0)
        {
            const Real zl = z / pert_ref_height;
            const Real damp = std::exp(-0.5 * zl * zl);

            Real rhotheta  = state(i,j,k,RhoTheta_comp);
            Real rho       = state(i,j,k,Rho_comp);
            Real qv        = state(i,j,k,RhoQ1_comp) / rho;
            Real Told      = getTgivenRandRTh(rho,rhotheta,qv);
            Real P         = getPgivenRTh(rhotheta,qv);

            Real Tpert    = tfac * damp * z * std::cos(cval * (x - xc));
            Real Tnew     = Told + Tpert;

            Real theta_new = getThgivenTandP(Tnew,P,rdOcp);
            Real rhonew    = getRhogivenThetaPress(theta_new,P,rdOcp,qv);

            state_pert(i, j, k, Rho_comp) = rhonew - rho;
            state_pert(i, j, k, RhoTheta_comp) = 0.0;
        }

        if (use_moisture && pert_deltaQV != 0.0)
        {
            const Real zl = z / pert_ref_height;
            const Real damp = std::exp(-0.5 * zl * zl);

            Real rhoold = state(i,j,k,Rho_comp);
            Real rhonew = rhoold + state_pert(i,j,k,Rho_comp);

            Real qvold = state(i,j,k,RhoQ1_comp) / rhoold;
            Real qvnew = qvold + qvfac * damp * z * std::cos(dval * (x - xc));

            state_pert(i, j, k, RhoQ1_comp) = rhonew * qvnew - rhoold * qvold;
        }
    });
