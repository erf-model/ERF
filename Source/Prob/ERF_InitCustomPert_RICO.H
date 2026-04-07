    ParmParse pp_prob("prob");

    Real A_0  = one; pp_prob.query("A_0", A_0);
    Real KE_0 = amrex::Real(0.1); pp_prob.query("KE_0", KE_0);

    // sinusoidal perturbations for regression tests
    Real pert_deltaT = zero; pp_prob.query("pert_deltaT", pert_deltaT);
    Real pert_deltaQV = zero; pp_prob.query("pert_deltaQV", pert_deltaQV);
    Real pert_periods_T = zero; pp_prob.query("pert_periods_T", pert_periods_T);
    Real pert_periods_QV = zero; pp_prob.query("pert_periods_QV", pert_periods_QV);

    Real T_0_Pert_Mag = zero; pp_prob.query("T_0_Pert_Mag", T_0_Pert_Mag);
    Real qv_0_Pert_Mag = zero; pp_prob.query("qv_0_Pert_Mag", qv_0_Pert_Mag);

    Real wbar_sub_max    = -amrex::Real(0.65); pp_prob.query("wbar_sub_max", wbar_sub_max);

    Real pert_ref_height = amrex::Real(100.0); pp_prob.query("pert_ref_height", pert_ref_height);

    bool custom_TKE = false; pp_prob.query("custom_TKE", custom_TKE);

    const Real prob_lo_x = geomdata.ProbLo()[0];
    const Real prob_hi_x = geomdata.ProbHi()[0];

    Real cval = pert_periods_T  * two * PI / (prob_hi_x - prob_lo_x);
    Real dval = pert_periods_QV * two * PI / (prob_hi_x - prob_lo_x);

    Real tfac  = pert_deltaT  * std::exp(myhalf) / pert_ref_height;
    Real qvfac = pert_deltaQV * std::exp(myhalf) / pert_ref_height;

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const Real rdOcp   = sc.rdOcp;

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        const Real z = prob_lo[2] + (k + myhalf) * dx[2];

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = myhalf * (prob_lo[0] + prob_hi[0]);
        const Real yc = myhalf * (prob_lo[1] + prob_hi[1]);
        const Real zc = myhalf * (prob_lo[2] + prob_hi[2]);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
        //
        // Add temperature perturbations -- we want to keep pressure constant
        //     so these effectively end up as density perturbations
        //
        if ((z <= pert_ref_height) && (T_0_Pert_Mag != zero)) {

            Real rhotheta  = state(i,j,k,RhoTheta_comp);
            Real rho       = state(i,j,k,Rho_comp);
            Real qv        = state(i,j,k,RhoQ1_comp) / rho;
            Real Told      = getTgivenRandRTh(rho,rhotheta,qv);
            Real P         = getPgivenRTh(rhotheta,qv);

            Real rand_double = amrex::Random(engine); // Between zero and one
            Real Tpert    = (rand_double*two - one)*T_0_Pert_Mag;
            Real Tnew     = Told + Tpert;

            Real theta_new = getThgivenTandP(Tnew,P,rdOcp);
            Real rhonew    = getRhogivenThetaPress(theta_new,P,rdOcp,qv);
            state_pert(i, j, k, Rho_comp) = rhonew - rho;

            // Note we do not perturb this
            state_pert(i, j, k, RhoTheta_comp) = zero;

            //  Instead of perturbing (rho theta) we perturb T and hold (rho theta) fixed,
            //  which ends up being stored as a perturbation in rho
            // state_pert(i, j, k, RhoTheta_comp) = (rand_double*two - one)*T_0_Pert_Mag;
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * exp(-amrex::Real(10.)*r*r);

        // Set an initial value for KE
        if (custom_TKE) {
            state_pert(i, j, k, RhoKE_comp) = (one - z/prob_hi[2]) * r_hse(i,j,k);
        } else {
            state_pert(i, j, k, RhoKE_comp) = KE_0;
        }

        if (use_moisture) {
            if ((z <= pert_ref_height) && (qv_0_Pert_Mag != zero))
            {
                Real rhoold = state(i,j,k,Rho_comp);
                Real rhonew = rhoold + state_pert(i,j,k,Rho_comp);

                Real qvold = state(i,j,k,RhoQ1_comp) / rhoold;

                Real rand_double = amrex::Random(engine); // Between zero and one
                Real qvnew = qvold + (rand_double*two - one)*qv_0_Pert_Mag;

                state_pert(i, j, k, RhoQ1_comp) = rhonew * qvnew - rhoold * qvold;
            }
        }

        // sinusoidal variation of theta and qv for regression tests
        if (pert_deltaT != zero)
        {
            const Real zl = z / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);

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
            state_pert(i, j, k, RhoTheta_comp) = zero;
        }

        if (use_moisture && pert_deltaQV != zero)
        {
            const Real zl = z / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);

            Real rhoold = state(i,j,k,Rho_comp);
            Real rhonew = rhoold + state_pert(i,j,k,Rho_comp);

            Real qvold = state(i,j,k,RhoQ1_comp) / rhoold;
            Real qvnew = qvold + qvfac * damp * z * std::cos(dval * (x - xc));

            state_pert(i, j, k, RhoQ1_comp) = rhonew * qvnew - rhoold * qvold;
        }
    });
