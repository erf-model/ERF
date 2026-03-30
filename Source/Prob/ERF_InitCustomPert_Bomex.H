
    ParmParse pp("prob");

    amrex::Real rho_0 = zero; pp.query("rho_0", rho_0);
    amrex::Real T_0   = zero; pp.query("T_0"  , T_0);
    amrex::Real A_0   = one; pp.query("A_0"  , A_0);
    amrex::Real KE_0  = amrex::Real(0.1); pp.query("KE_0" , KE_0);

    Real  T_0_Pert_Mag   =   zero; pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    Real qv_0_Pert_Mag   =   zero; pp.query("qv_0_Pert_Mag", qv_0_Pert_Mag);
    Real pert_ref_height = amrex::Real(100.0); pp.query("pert_ref_height", pert_ref_height);
    bool custom_TKE      = false; pp.query("custom_TKE", custom_TKE);

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
            state_pert(i, j, k, RhoQ1_comp) = zero;
            state_pert(i, j, k, RhoQ2_comp) = zero;
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
    });
