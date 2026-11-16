
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

    // Optional warm bubble, made as WRF's idealized fire initialisation makes it
    // (module_initialize_fire.F: delt_perturbation, x/y/zrad_perturbation, hght_perturbation):
    // theta += bubble_dT cos^2(pi r / 2) for r = |((x-x_c)/x_r, (y-y_c)/y_r, (z-z_c)/z_r)| <= 1,
    // with z the cell-centre height z_cc (the same height the other perturbations here use).
    // bubble_dT = 0 (the default) adds no bubble.
    Real bubble_dT  = amrex::Real(0);
    const bool bubble_dT_given = pp_prob.query("bubble_dT", bubble_dT);
    Real bubble_x_c = amrex::Real(0); pp_prob.query("bubble_x_c", bubble_x_c);
    Real bubble_y_c = amrex::Real(0); pp_prob.query("bubble_y_c", bubble_y_c);
    Real bubble_z_c = amrex::Real(0); pp_prob.query("bubble_z_c", bubble_z_c);
    Real bubble_x_r = amrex::Real(0); pp_prob.query("bubble_x_r", bubble_x_r);
    Real bubble_y_r = amrex::Real(0); pp_prob.query("bubble_y_r", bubble_y_r);
    Real bubble_z_r = amrex::Real(0); pp_prob.query("bubble_z_r", bubble_z_r);
    const bool add_bubble = (bubble_dT != amrex::Real(0));
    // Report once per level: the tile holding the level's lowest domain corner prints it.
    // A refined level whose grids do not cover that corner gets the bubble without the line
    // (this block sees only its own tile, not the level's BoxArray).
    const bool bubble_reporter = bx.contains(geomdata.Domain().smallEnd());
    if (add_bubble) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(bubble_dT) &&
            std::isfinite(bubble_x_c) && std::isfinite(bubble_y_c) && std::isfinite(bubble_z_c),
            "prob.bubble_dT and prob.bubble_x_c, bubble_y_c, bubble_z_c must be finite numbers");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bubble_x_r > 0 && bubble_y_r > 0 && bubble_z_r > 0 &&
            std::isfinite(bubble_x_r) && std::isfinite(bubble_y_r) && std::isfinite(bubble_z_r),
            "prob.bubble_dT needs positive, finite prob.bubble_x_r, bubble_y_r and bubble_z_r");
        if (bubble_reporter) {
            amrex::AllPrint() << "Level " << lev << ": adding a warm bubble prob.bubble_dT = " << bubble_dT
                              << " K at (" << bubble_x_c << ", " << bubble_y_c << ", " << bubble_z_c
                              << ") m with radii (" << bubble_x_r << ", " << bubble_y_r << ", "
                              << bubble_z_r << ") m, "
                              << (sc.fixed_density[lev] ? "at the fixed density" : "at constant pressure")
                              << std::endl;
        }
    } else {
        // The bubble keys mean nothing without bubble_dT: stop if it was never set, warn if
        // it was set to zero (the bubble switched off on purpose).
        const bool bubble_keys_given =
            pp_prob.contains("bubble_x_c") || pp_prob.contains("bubble_y_c") ||
            pp_prob.contains("bubble_z_c") || pp_prob.contains("bubble_x_r") ||
            pp_prob.contains("bubble_y_r") || pp_prob.contains("bubble_z_r");
        if (bubble_keys_given) {
            if (!bubble_dT_given) {
                amrex::Abort("prob.bubble_x_c/y_c/z_c/x_r/y_r/z_r are set but prob.bubble_dT is not: "
                             "set prob.bubble_dT (the bubble's theta excess in K) or remove them");
            } else if (bubble_reporter) {
                amrex::AllPrint() << "WARNING: prob.bubble_dT = 0, so the prob.bubble_* centre and "
                                  << "radii are not used" << std::endl;
            }
        }
    }
    // Under a fixed-in-time density (anelastic) only rho*theta is perturbed, so the bubble
    // is a theta excess at that density. Otherwise it is made at constant pressure, as WRF
    // makes it (theta raised, density recomputed from the unchanged pressure).
    const bool bubble_fixed_density = (sc.fixed_density[lev] != 0);
    Real pert_periods_T = amrex::Real(5.0); pp_prob.query("pert_periods_T", pert_periods_T);

    // Capture by value for GPU
    const Real dx = geomdata.CellSize(0);
    const Real dy = geomdata.CellSize(1);
    const Real prob_lo_x = geomdata.ProbLo(0);
    const Real prob_lo_y = geomdata.ProbLo(1);
    const Real prob_lo_z = phys_zlo(geomdata, lev);   // physical, not the computational, bottom
    const Real prob_hi_x = geomdata.ProbHi(0);
    const Real prob_hi_y = geomdata.ProbHi(1);
    const Real prob_hi_z = phys_zhi(geomdata, lev);   // physical, not the computational, top

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

    if (KE_decay_height > 0 &&
        state_pert.nComp() > RhoKE_comp) {
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
    // fixed-in-time density (implied by a non-FullyCompressible erf.compressibility_type).  Say so rather than fail
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
        const Real z_agl = z - erf_prob_heights::z_ground_cc(i,j,z_sfc_nd); // height above the ground

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z_agl <= pert_ref_height) && (T_0_Pert_Mag != amrex::Real(0))) {
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

        // Warm bubble (added to any random perturbation above)
        if (add_bubble) {
            const Real rx = (x - bubble_x_c) / bubble_x_r;
            const Real ry = (y - bubble_y_c) / bubble_y_r;
            const Real rz = (z - bubble_z_c) / bubble_z_r;
            const Real rr = std::sqrt(rx*rx + ry*ry + rz*rz);
            if (rr <= amrex::Real(1)) {
                const Real cb     = std::cos(myhalf * PI * rr);
                const Real dtheta = bubble_dT * cb * cb;
                const Real rho      = state(i,j,k,Rho_comp);
                const Real rhotheta = state(i,j,k,RhoTheta_comp);
                if (bubble_fixed_density) {
                    state_pert(i, j, k, RhoTheta_comp) += rho * dtheta;
                } else {
                    // The pressure is a function of rho*theta and qv only, so keeping both
                    // keeps the pressure: the density drops to rho*theta / (theta + dtheta).
                    const Real rho_new = rhotheta / (rhotheta / rho + dtheta);
                    state_pert(i, j, k, Rho_comp) += rho_new - rho;
                    if (use_moisture) {
                        // keep qv = RhoQ1/rho, and so the pressure, as the density changes
                        state_pert(i, j, k, RhoQ1_comp) += (rho_new - rho) * state(i,j,k,RhoQ1_comp) / rho;
                    }
                }
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
            const Real zl   = z_agl / pert_ref_height;
            const Real damp = std::exp(-myhalf * zl * zl);

            const Real rho      = state(i,j,k,Rho_comp);
            const Real rhotheta = state(i,j,k,RhoTheta_comp);
            const Real qv       = (use_moisture) ? state(i,j,k,RhoQ1_comp) / rho : amrex::Real(0);

            const Real Told = getTgivenRandRTh(rho,rhotheta,qv);
            const Real P    = getPgivenRTh(rhotheta,qv);

            const Real Tnew      = Told + tfac * damp * z_agl * std::cos(cval * (x - xc));
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
        //
        // KE_decay_height is a depth of turbulence above the ground, so the profile
        // is a function of height above the terrain surface, not of absolute z.  Over
        // terrain higher than KE_decay_height the latter would clamp the decay factor
        // to its floor at the ground -- exactly backwards.  On a flat mesh this is the
        // height above geometry.prob_lo[2].  Under EB, z_cc is already measured
        // from the EB surface, so the ground is at zero there (see z_sfc_nd).  Matches
        // the surface height used by SurfaceLayer::init_tke_from_ustar.
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
                    std::pow(1 - amrex::min(z_agl/KE_decay_height,amrex::Real(1)), KE_decay_order),
                    amrex::Real(1e-12));
            }
        }
    });
