
    ParmParse pp("prob");

    Real A_0 = amrex::Real(1.0); pp.query("A_0", A_0);
    Real KE_0 = amrex::Real(0.1); pp.query("KE_0", KE_0);

    // Random temperature perturbation amplitude and taper/cutoff height.
    Real T_0_Pert_Mag = amrex::Real(0.0); pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    Real qv_0_Pert_Mag = amrex::Real(0.0); pp.query("qv_0_Pert_Mag", qv_0_Pert_Mag);
    Real pert_ref_height = amrex::Real(200.0); pp.query("pert_ref_height", pert_ref_height);
    bool custom_TKE = false; pp.query("custom_TKE", custom_TKE);

    const bool use_moisture = (sc.moisture_type != MoistureType::None);
    const bool use_terrain = (sc.terrain_type != TerrainType::None);
    const bool const_rho = (sc.fixed_density[lev] == 1);
    const Real rdOcp = sc.rdOcp;

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k,
                                             const RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = use_terrain ? z_cc(i, j, k) : prob_lo[2] + (k + 0.5) * dx[2];

        const Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
        const Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
        const Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);
        const Real r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Apply a random temperature perturbation that tapers linearly from
        // T_0_Pert_Mag at the surface to zero at pert_ref_height.
        if (pert_ref_height > 0.0 && z <= pert_ref_height &&
            T_0_Pert_Mag != 0.0) {
            const Real rhotheta = state(i,j,k,RhoTheta_comp);
            const Real rho = state(i,j,k,Rho_comp);
            const Real qv = state(i,j,k,RhoQ1_comp) / rho;
            const Real Told = getTgivenRandRTh(rho, rhotheta, qv);
            const Real P = getPgivenRTh(rhotheta, qv);

            const Real rand_double = amrex::Random(engine);
            const Real Tnew = Told + T_0_Pert_Mag * (1.0 - 2.0 * rand_double)
                              * (1.0 - z / pert_ref_height);
            const Real theta_new = getThgivenTandP(Tnew, P, rdOcp);
            const Real rhonew = getRhogivenThetaPress(theta_new, P, rdOcp, qv);

            if (const_rho) {
                state_pert(i, j, k, Rho_comp) = 0.0;
                state_pert(i, j, k, RhoTheta_comp) =
                    r_hse(i, j, k) * (theta_new - (rhotheta / rho));
            } else {
                state_pert(i, j, k, Rho_comp) = rhonew - rho;
                state_pert(i, j, k, RhoTheta_comp) = 0.0;
            }
        }

        state_pert(i, j, k, RhoScalar_comp) = A_0 * exp(-10.0*r*r);

        if (custom_TKE) {
            state_pert(i, j, k, RhoKE_comp) =
                (1.0 - z / prob_hi[2]) * r_hse(i,j,k);
        } else {
            state_pert(i, j, k, RhoKE_comp) = KE_0;
        }

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
            if (z <= pert_ref_height && qv_0_Pert_Mag != 0.0) {
                const Real rhoold = state(i,j,k,Rho_comp);
                const Real rhonew = rhoold + state_pert(i,j,k,Rho_comp);
                const Real qvold = state(i,j,k,RhoQ1_comp) / rhoold;
                const Real rand_double = amrex::Random(engine);
                const Real qvnew = qvold +
                    (2.0 * rand_double - 1.0) * qv_0_Pert_Mag;
                state_pert(i, j, k, RhoQ1_comp) =
                    rhonew * qvnew - rhoold * qvold;
            }
        }
    });
