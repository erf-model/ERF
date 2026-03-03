    ParmParse pp_prob("prob");

    Real T_0 = 300.0; pp_prob.query("T_0", T_0);
    Real x_c = 0.0; pp_prob.query("x_c", x_c);
    Real y_c = 0.0; pp_prob.query("y_c", y_c);
    Real z_c = 0.0; pp_prob.query("z_c", z_c);
    Real x_r = 0.0; pp_prob.query("x_r", x_r);
    Real y_r = 0.0; pp_prob.query("y_r", y_r);
    Real z_r = 0.0; pp_prob.query("z_r", z_r);

    Real T_pert = -15.0; pp_prob.query("T_pert", T_pert);

    bool T_pert_is_airtemp = true; pp_prob.query("T_pert_is_airtemp", T_pert_is_airtemp);
    bool perturb_rho = true; pp_prob.query("perturb_rho", perturb_rho);

    bool do_moist_bubble = false; pp_prob.query("do_moist_bubble", do_moist_bubble);
    Real theta_pert  = 2.0; pp_prob.query("theta_pert", theta_pert);

    const int khi = geomdata.Domain().bigEnd()[2];

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    const Real dz        = geomdata.CellSize()[2];
    const Real rdOcp     = sc.rdOcp;

    amrex::Print() << "Bubble delta T = " << T_pert << " K" << std::endl;
    amrex::Print() << "  centered at ("
                   << x_c << " " << y_c << " " << z_c << ")" << std::endl;
    amrex::Print() << "  with extent ("
                   << x_r << " " << y_r << " " << z_r << ")" << std::endl;

    if (T_0 <= 0)
    {
        amrex::Print() << "Ignoring T_0 = " << T_0
                       << ", background fields should have been initialized with erf.init_type"
                       << std::endl;
    }

    Real qt_init = 0.02; pp_prob.query("qt_init", qt_init);

    Real eq_pot_temp = 320.0; pp_prob.query("eq_pot_temp",eq_pot_temp);

    bool use_empirical = false;
    pp_prob.query("use_empircal_psat",use_empirical);

    if (do_moist_bubble) {
        Vector<Real> h_r(khi+2);
        Vector<Real> h_p(khi+2);
        Vector<Real> h_t(khi+2);
        Vector<Real> h_q_v(khi+2);

        Gpu::DeviceVector<Real> d_r(khi+2);
        Gpu::DeviceVector<Real> d_p(khi+2);
        Gpu::DeviceVector<Real> d_t(khi+2);
        Gpu::DeviceVector<Real> d_q_v(khi+2);

        HSEutils::init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(),h_q_v.data(),dz,khi,
                                                 qt_init,eq_pot_temp,use_empirical,false);

        Gpu::copyAsync(Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
        Gpu::copyAsync(Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
        Gpu::copyAsync(Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
        Gpu::copyAsync(Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

        Real* theta_back   = d_t.data();
        Real* p_back   = d_p.data();
        Real* q_v_back = d_q_v.data();

        int moisture_type = 1;

        if(sc.moisture_type == MoistureType::SAM) {
            moisture_type = 1;
        } else if(sc.moisture_type == MoistureType::SAM_NoIce ||
                  sc.moisture_type == MoistureType::SAM_NoPrecip_NoIce) {
            moisture_type = 2;
        }

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();
            const amrex::Real x        = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real y        = prob_lo[1] + (j + 0.5) * dx[1];
            const amrex::Real z        = prob_lo[2] + (k + 0.5) * dx[2];

            Real rad, delta_theta, theta_total, rho, RH;

            // Introduce the warm bubble. Assume that the bubble is pressure matched with the background
            rad = 0.0;
            if (x_r > 0) rad += std::pow((x - x_c)/x_r, 2);
            if (y_r > 0) rad += std::pow((y - y_c)/y_r, 2);
            if (z_r > 0) rad += std::pow((z - z_c)/z_r, 2);
            rad = std::sqrt(rad);

            if (rad <= 1.0) {
                delta_theta = theta_pert*std::pow(cos(PI*rad/2.0),2);
            } else {
                 delta_theta = 0.0;
            }

            theta_total  = theta_back[k]*(delta_theta/300.0 + 1);
            Real T = getTgivenPandTh(p_back[k], theta_total, (R_d/Cp_d));
            rho    = p_back[k]/(R_d*T*(1.0 + (R_v/R_d)*q_v_back[k]));
            RH     = compute_relative_humidity(use_empirical);
            Real q_v_hot = vapor_mixing_ratio(p_back[k], T, RH);

            // Compute background quantities
            Real T_back   = getTgivenPandTh(p_back[k], theta_back[k], (R_d/Cp_d));
            Real rho_back = p_back[k]/(R_d*T_back*(1.0 + (R_v/R_d)*q_v_back[k]));

            // This version perturbs rho but not p
            state_pert(i, j, k, RhoTheta_comp) = rho*theta_total - rho_back*theta_back[k]*(1.0 + (R_v/R_d)*q_v_back[k]);
            state_pert(i, j, k, Rho_comp)      = rho - rho_back*(1.0 + qt_init);

            // Set scalar = 0 everywhere
            state_pert(i, j, k, RhoScalar_comp) = 0.0;

            // mean states
            state_pert(i, j, k, RhoQ1_comp) = rho*q_v_hot;
            state_pert(i, j, k, RhoQ2_comp) = rho*(qt_init - q_v_hot);

            // Cold microphysics are present
            int nstate = state_pert.nComp();
            if (nstate == NVAR_max) {
                Real omn;
                if(moisture_type == 1) {
                    omn = std::max(0.0,std::min(1.0,(T-tbgmin)*a_bg));
                } else if(moisture_type == 2) {
                    omn = 1.0;
                } else {
                    omn = 0.0;
                    Abort("Invalid moisture type specified");
                }
                Real qn  = state_pert(i, j, k, RhoQ2_comp);
                state_pert(i, j, k, RhoQ2_comp) = qn * omn;
                state_pert(i, j, k, RhoQ3_comp) = qn * (1.0-omn);
            }
        });
    } else {
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const Real z = prob_lo[2] + (k + 0.5) * dx[2];

            //perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k), rdOcp,
            //                  state_pert(i, j, k, Rho_comp),
            //                  state_pert(i, j, k, RhoTheta_comp));

            // Perturbation air temperature
            // - The bubble is either cylindrical (for 2-D problems, if two
            //   radial extents are specified) or an ellipsoid (if all three
            //   radial extents are specified).
            amrex::Real L = 0.0;
            if (x_r > 0) L += std::pow((x - x_c)/x_r, 2);
            if (y_r > 0) L += std::pow((y - y_c)/y_r, 2);
            if (z_r > 0) L += std::pow((z - z_c)/z_r, 2);
            L = std::sqrt(L);
            amrex::Real dT;
            if (L > 1.0) {
                dT = 0.0;
            }
            else {
                dT = T_pert * std::pow(cos(PI*L/2.0),2);
            }

            // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
            const amrex::Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

            // Note: theta_perturbed is theta PLUS perturbation in theta
            amrex::Real theta_perturbed;
            if (T_pert_is_airtemp) {
                // dT is air temperature
                theta_perturbed = (Tbar_hse + dT)*std::pow(p_0/p_hse(i,j,k), rdOcp);
            } else {
                // dT is potential temperature
                theta_perturbed = Tbar_hse*std::pow(p_0/p_hse(i,j,k), rdOcp) + dT;
            }

            if (perturb_rho)
            {
                // this version perturbs rho but not p (i.e., rho*theta)
                // - hydrostatic rebalance is needed (TODO: is this true?)
                // - this is the approach taken in the density current problem
                state_pert(i,j,k,Rho_comp)      = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
                state_pert(i,j,k,RhoTheta_comp) = 0.0; // i.e., hydrostatically balanced pressure stays const
            }
            else
            {
                // this version perturbs rho*theta (i.e., p) but not rho
                state_pert(i,j,k,Rho_comp)      = 0.0; // i.e., hydrostatically balanced density stays const
                state_pert(i,j,k,RhoTheta_comp) = r_hse(i,j,k) * theta_perturbed - getRhoThetagivenP(p_hse(i,j,k));
            }
        });
    } // do_moist_bubble
