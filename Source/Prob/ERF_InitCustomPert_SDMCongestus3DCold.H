
    // On a mesh other than ConstantDz, heights come from z_cc and the vertical levels,
    //    not from geomdata's prob_lo[2], prob_hi[2] and CellSize(2)
    const bool use_z_phys = use_physical_z();
    const Real phys_z_lo  = phys_zlo(geomdata, lev);
    const Real phys_z_hi  = phys_zhi(geomdata, lev);

    ParmParse pp("prob");

    amrex::Real rho_0 = zero; pp.query("rho_0", rho_0);
    amrex::Real T_0   = zero; pp.query("T_0"  , T_0);
    amrex::Real A_0   = one; pp.query("A_0"  , A_0);
    amrex::Real KE_0  = amrex::Real(0.1); pp.query("KE_0" , KE_0);

    Real  T_0_Pert_Mag   =   zero; pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    Real  Q_0_Pert_Mag   =   zero; pp.query("Q_0_Pert_Mag", Q_0_Pert_Mag);
    Real pert_ref_height = amrex::Real(100.0); pp.query("pert_ref_height", pert_ref_height);

    Real pert_deltaT     =   zero; pp.query("pert_deltaT", pert_deltaT);
    Real pert_deltaQV    =   zero; pp.query("pert_deltaQV", pert_deltaQV);
    Real pert_periods_T  =   zero; pp.query("pert_periods_T", pert_periods_T);
    Real pert_periods_QV =   zero; pp.query("pert_periods_QV", pert_periods_QV);

    auto problo_arr = geomdata.ProbLo();
    auto probhi_arr = geomdata.ProbHi();

    Real cval  = pert_periods_T  * two * PI / (probhi_arr[0] - problo_arr[0]);
    Real dval  = pert_periods_QV * two * PI / (probhi_arr[0] - problo_arr[0]);
    Real tfac  = pert_deltaT  * std::exp(myhalf) / pert_ref_height;
    Real qvfac = pert_deltaQV * std::exp(myhalf) / pert_ref_height;

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const Real rdOcp   = sc.rdOcp;

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = use_z_phys ? z_cc(i,j,k) : prob_lo[2] + (k + myhalf) * dx[2];
        const Real z_agl = z - erf_prob_heights::z_ground_cc(i,j,z_sfc_nd); // height above the ground

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = myhalf * (prob_lo[0] + prob_hi[0]);
        const Real yc = myhalf * (prob_lo[1] + prob_hi[1]);
        const Real zc = myhalf * (phys_z_lo + phys_z_hi);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z_agl <= pert_ref_height) && (T_0_Pert_Mag != zero)) {
            Real rand_double = amrex::Random(engine); // Between zero and one
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*two - one)*T_0_Pert_Mag;
        }
        if (use_moisture && (z_agl <= pert_ref_height) && (Q_0_Pert_Mag != zero)) {
            Real rand_double = amrex::Random(engine); // Between zero and one
            state_pert(i, j, k, RhoQ1_comp) = (rand_double*two - one)*Q_0_Pert_Mag;
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * exp(-amrex::Real(10.)*r*r);

        // Set an initial value for KE
        state_pert(i, j, k, RhoKE_comp) = KE_0;

        // sinusoidal variation of theta and qv for regression tests
        if (pert_deltaT != zero)
        {
            const Real zl = z_agl / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);

            Real rhotheta  = state(i,j,k,RhoTheta_comp);
            Real rho       = state(i,j,k,Rho_comp);
            Real qv        = (use_moisture) ? state(i,j,k,RhoQ1_comp) / rho : zero;
            Real Told      = getTgivenRandRTh(rho,rhotheta,qv);
            Real P         = getPgivenRTh(rhotheta,qv);

            Real Tpert    = tfac * damp * z_agl * std::cos(cval * (x - xc));
            Real Tnew     = Told + Tpert;

            Real theta_new = getThgivenTandP(Tnew,P,rdOcp);
            Real rhonew    = getRhogivenThetaPress(theta_new,P,rdOcp,qv);

            state_pert(i, j, k, Rho_comp) = rhonew - rho;
            state_pert(i, j, k, RhoTheta_comp) = zero;
        }

        if (use_moisture && pert_deltaQV != zero)
        {
            const Real zl = z_agl / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);

            Real rhoold = state(i,j,k,Rho_comp);
            Real rhonew = rhoold + state_pert(i,j,k,Rho_comp);

            Real qvold = state(i,j,k,RhoQ1_comp) / rhoold;
            Real qvnew = qvold + qvfac * damp * z_agl * std::cos(dval * (x - xc));

            state_pert(i, j, k, RhoQ1_comp) = rhonew * qvnew - rhoold * qvold;
        }
    });
