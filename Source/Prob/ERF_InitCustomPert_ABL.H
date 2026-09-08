
    ParmParse pp_prob("prob");

    int deterministic_ic_pert = 0;
    int fix_random_seed = 0;
    long random_seed = -1;
    {
        ParmParse pp_erf("erf");
        pp_erf.query("deterministic_ic_pert", deterministic_ic_pert);
        pp_erf.query("fix_random_seed", fix_random_seed);
        pp_erf.query("random_seed", random_seed);
    }
    const std::uint64_t ic_pert_seed = fix_random_seed
        ? std::uint64_t(1024)
        : ((random_seed >= 0) ? static_cast<std::uint64_t>(random_seed) : std::uint64_t(0));
    constexpr int temperature_hash_comp = 0;

    Real rho_0 = amrex::Real(0); pp_prob.query("rho_0", rho_0);
    Real T_0   = amrex::Real(0); pp_prob.query("T_0", T_0);
    Real A_0   = amrex::Real(1); pp_prob.query("A_0", A_0);
    Real KE_0  = amrex::Real(0.1); pp_prob.query("KE_0", KE_0);
    Real rhoKE_0 = -1; pp_prob.query("rhoKE_0", rhoKE_0);

    Real KE_decay_height = -1; pp_prob.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order  =  1; pp_prob.query("KE_decay_order", KE_decay_order);

    // random initial perturbations (legacy code)
    Real T_0_Pert_Mag = amrex::Real(0); pp_prob.query("T_0_Pert_Mag", T_0_Pert_Mag);
    bool pert_rhotheta = true; pp_prob.query("pert_rhotheta", pert_rhotheta);

    // divergence-free initial perturbations
    Real pert_deltaU = amrex::Real(0); pp_prob.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = amrex::Real(0); pp_prob.query("pert_deltaV", pert_deltaV);
    Real pert_ref_height = amrex::Real(100.0); pp_prob.query("pert_ref_height", pert_ref_height);

    // sinusoidal perturbations for regression tests.  pert_periods_T defaults to
    // 5 to match pert_periods_U/V in ERF_InitCustomPertVels_ABL.H: a zero period
    // would make cos(cval*(x-xc)) identically one, i.e. a uniform warm bias
    // rather than a zero-mean perturbation.
    Real pert_deltaT    = amrex::Real(0); pp_prob.query("pert_deltaT", pert_deltaT);
    Real pert_periods_T = amrex::Real(5.0); pp_prob.query("pert_periods_T", pert_periods_T);

    // Capture by value for GPU
    const Real dx = geomdata.CellSize(0);
    const Real dy = geomdata.CellSize(1);
    const Real prob_lo_x = geomdata.ProbLo(0);
    const Real prob_lo_y = geomdata.ProbLo(1);
    const Real prob_lo_z = geomdata.ProbLo(2);
    const Real prob_hi_x = geomdata.ProbHi(0);
    const Real prob_hi_y = geomdata.ProbHi(1);
    const Real prob_hi_z = geomdata.ProbHi(2);

    // Define a point (xc,yc,zc) at the center of the domain
    const Real xc = myhalf * (prob_lo_x + prob_hi_x);
    const Real yc = myhalf * (prob_lo_y + prob_hi_y);
    const Real zc = myhalf * (prob_lo_z + prob_hi_z);

    // Sinusoidal perturbation wavenumber and amplitude factor.  These mirror
    // the divergence-free velocity perturbations in ERF_InitCustomPertVels_ABL.H
    // so that tfac*damp*z peaks at pert_deltaT when z = pert_ref_height.
    // pert_ref_height <= 0 is used elsewhere in this file to mean "no
    // perturbations", so guard the division: an infinite tfac would meet a zero
    // damp inside the kernel and yield NaN under amrex.fpe_trap_invalid.
    const bool add_sinusoidal_T = (pert_deltaT != amrex::Real(0)) &&
                                  (pert_ref_height > amrex::Real(0));
    const Real cval = pert_periods_T * two * PI / (prob_hi_x - prob_lo_x);
    const Real tfac = add_sinusoidal_T
                    ? pert_deltaT * std::exp(myhalf) / pert_ref_height
                    : amrex::Real(0);

    const bool use_moisture = (sc.moisture_type != MoistureType::None);
    const Real rdOcp        = sc.rdOcp;

    if (KE_decay_height > 0) {
        amrex::Print() << "Initial KE profile (order " << KE_decay_order
                       << ") will extend up to " << KE_decay_height
                       << std::endl;
    }

#ifndef AMREX_USE_GPU
    if (pert_ref_height > 0) {
        if (T_0_Pert_Mag != amrex::Real(0)) {
            if (pert_rhotheta) {
                amrex::Print() << "Adding random rho*theta perturbations" << std::endl;
            } else {
                amrex::Print() << "Adding random theta perturbations" << std::endl;
            }
        }
        if (add_sinusoidal_T) {
            amrex::Print() << "Adding sinusoidal temperature perturbations "
                           << pert_deltaT << " over " << pert_periods_T
                           << " periods" << std::endl;
        }
    }
    // The sinusoidal temperature perturbation is carried entirely by the density
    // perturbation, which ERF_InitCustomPertState.cpp drops when the level uses a
    // fixed-in-time density (implied by erf.anelastic).  Say so rather than fail
    // silently, since the companion velocity perturbations still apply.
    if (add_sinusoidal_T && sc.fixed_density[lev]) {
        amrex::Print() << "WARNING: prob.pert_deltaT has no effect at level " << lev
                       << " because fixed_density is set there" << std::endl;
    }
#endif

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real x = prob_lo_x + (i + myhalf) * dx;
        const Real y = prob_lo_y + (j + myhalf) * dy;
        const Real z = z_cc(i,j,k);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z <= pert_ref_height) && (T_0_Pert_Mag != amrex::Real(0))) {
            Real rand_double;
            if (deterministic_ic_pert) {
                rand_double = erf_hash_rng::hash_uniform(
                    i, j, k, temperature_hash_comp, lev, ic_pert_seed);
            } else {
                rand_double = amrex::Random(engine); // Between zero and one
            }
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*amrex::Real(2) - amrex::Real(1))*T_0_Pert_Mag;
            if (!pert_rhotheta) {
                // we're perturbing theta, not rho*theta
                state_pert(i, j, k, RhoTheta_comp) *= r_hse(i,j,k);
            }
        }

        // Sinusoidal temperature perturbation for regression tests.  Unlike the
        // random branch above, this is a pure function of position, so the initial
        // state does not depend on the box decomposition, the MPI rank count, or
        // the OpenMP tile size -- which is what a bitwise gold-file comparison
        // requires.
        //
        // The perturbation is isobaric: T is perturbed at fixed pressure, so it is
        // carried by the density perturbation while rho*theta stays unchanged.
        // Because P = P(rho*theta, qv), holding P fixed also requires qv to be
        // unchanged, so RhoQ1 has to be rescaled to track the new density --
        // otherwise qv = RhoQ1/rho drifts and drags the pressure along with it.
        //
        // This deliberately differs from ERF_InitCustomPert_SDMCongestus3D.H, which
        // assigns state_pert(RhoTheta_comp) = 0 at this point: that would discard
        // any random perturbation written by the T_0_Pert_Mag branch above.  Leaving
        // RhoTheta_comp alone is equivalent whenever T_0_Pert_Mag == 0 (state_pert is
        // zero-initialized by the caller) and non-destructive when it is not -- but
        // setting both knobs does perturb rho*theta, and the result is then no
        // longer isobaric.
        if (add_sinusoidal_T) {
            const Real zl   = z / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);

            const Real rho      = state(i,j,k,Rho_comp);
            const Real rhotheta = state(i,j,k,RhoTheta_comp);
            const Real qv       = (use_moisture) ? state(i,j,k,RhoQ1_comp) / rho : amrex::Real(0);

            const Real Told = getTgivenRandRTh(rho,rhotheta,qv);
            const Real P    = getPgivenRTh(rhotheta,qv);

            const Real Tnew      = Told + tfac * damp * z * std::cos(cval * (x - xc));
            const Real theta_new = getThgivenTandP(Tnew,P,rdOcp);
            const Real rho_new   = getRhogivenThetaPress(theta_new,P,rdOcp,qv);

            state_pert(i, j, k, Rho_comp) += rho_new - rho;
            if (use_moisture) {
                // Keep qv, and therefore P, fixed as the density changes.
                state_pert(i, j, k, RhoQ1_comp) += (rho_new - rho) * qv;
            }
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * std::exp(-amrex::Real(10.)*r*r);

        // Set an initial value for SGS KE
        if (state_pert.nComp() > RhoKE_comp) {
            // Deardorff
            if (rhoKE_0 > 0) {
                state_pert(i, j, k, RhoKE_comp) = rhoKE_0;
            } else {
                state_pert(i, j, k, RhoKE_comp) = r_hse(i,j,k) * KE_0;
            }
            if (KE_decay_height > 0) {
                // scale initial SGS kinetic energy with height
                state_pert(i, j, k, RhoKE_comp) *= amrex::max(
                    std::pow(1 - amrex::min(z/KE_decay_height,amrex::Real(1)), KE_decay_order), Real(1e-12));
            }
        }
    });
