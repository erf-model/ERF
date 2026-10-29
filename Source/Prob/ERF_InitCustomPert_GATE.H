
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
    Real qv_0_Pert_Mag   =   zero; pp.query("qv_0_Pert_Mag", qv_0_Pert_Mag);
    Real pert_ref_height = amrex::Real(100.0); pp.query("pert_ref_height", pert_ref_height);
    bool custom_TKE      = false; pp.query("custom_TKE", custom_TKE);

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

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
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

        // Temperature and qv perturbations: the random ones (T_0_Pert_Mag, qv_0_Pert_Mag)
        //     and the sinusoidal ones for regression tests (pert_deltaT, pert_deltaQV) are
        //     summed, then applied once, so that turning on one does not discard the other
        const bool rand_T  = (z_agl <= pert_ref_height) && (T_0_Pert_Mag != zero);
        const bool rand_qv = use_moisture && (z_agl <= pert_ref_height) && (qv_0_Pert_Mag != zero);
        const bool sin_T   = (pert_deltaT != zero);
        const bool sin_qv  = use_moisture && (pert_deltaQV != zero);

        // The random numbers are drawn in the same order as when each was applied on its own
        Real dT  = zero;
        Real dqv = zero;
        if (rand_T) {
            Real rand_double = amrex::Random(engine); // Between zero and one
            dT += (rand_double*two - one)*T_0_Pert_Mag;
        }
        if (rand_qv) {
            Real rand_double = amrex::Random(engine); // Between zero and one
            dqv += (rand_double*two - one)*qv_0_Pert_Mag;
        }
        if (sin_T || sin_qv) {
            const Real zl = z_agl / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);
            if (sin_T)  { dT  += tfac  * damp * z_agl * std::cos(cval * (x - xc)); }
            if (sin_qv) { dqv += qvfac * damp * z_agl * std::cos(dval * (x - xc)); }
        }

        // Apply both at once (ERF_MoistPerturbation.H): the temperature perturbation at
        //     fixed pressure ends up as a density perturbation (rho*theta is held), and qv
        //     becomes qv + dqv on the new density, so a temperature perturbation alone
        //     leaves qv, and the pressure, as they were
        if (rand_T || sin_T || rand_qv || sin_qv) {
            Real drho, drhoqv;
            erf_moist_perturbation::perturb_temperature_and_vapour(
                state(i,j,k,Rho_comp), state(i,j,k,RhoTheta_comp),
                use_moisture ? state(i,j,k,RhoQ1_comp) : zero, dT, dqv,
                rand_T || sin_T, use_moisture, rdOcp, drho, drhoqv);
            if (rand_T || sin_T) {
                state_pert(i, j, k, Rho_comp) = drho;
                // Note we do not perturb this
                state_pert(i, j, k, RhoTheta_comp) = zero;
            }
            if (use_moisture) { state_pert(i, j, k, RhoQ1_comp) = drhoqv; }
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * std::exp(-amrex::Real(10.)*r*r);

        // Set an initial value for KE
        if (custom_TKE) {
            state_pert(i, j, k, RhoKE_comp) = (one - z/phys_z_hi) * r_hse(i,j,k);
        } else {
            state_pert(i, j, k, RhoKE_comp) = KE_0;
        }

        if (use_moisture) {
            // The density change above would change qc with it; rescale rho*qc so the
            //     cloud-water mixing ratio is held, as in Aaron Lattanzi's #4772 for BOMEX
            state_pert(i, j, k, RhoQ2_comp) = state_pert(i,j,k,Rho_comp)
                                            * state(i,j,k,RhoQ2_comp) / state(i,j,k,Rho_comp);
        }
    });
