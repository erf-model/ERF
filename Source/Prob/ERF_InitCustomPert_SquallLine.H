
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
    const Real& dz        = geomdata.CellSize()[2];

    // Call the routine to calculate the 1d background condition

    Vector<Real> h_r(khi+2);
    Vector<Real> h_p(khi+2);
    Vector<Real> h_t(khi+2);
    Vector<Real> h_q_v(khi+2);

    Gpu::DeviceVector<Real> d_r(khi+2);
    Gpu::DeviceVector<Real> d_p(khi+2);
    Gpu::DeviceVector<Real> d_t(khi+2);
    Gpu::DeviceVector<Real> d_q_v(khi+2);

    ParmParse pp_prob("prob");
    Real q_t = amrex::Real(0.02);
    pp_prob.query("qt_init",q_t);

    Real eq_pot_temp = amrex::Real(320.0);
    pp_prob.query("eq_pot_temp",eq_pot_temp);

    bool use_empirical = false;
    pp_prob.query("use_empirical_psat",use_empirical);

    Real height = amrex::Real(1200.0);
    pp_prob.query("height",height);

    Real z_tr = amrex::Real(12000.0);
    pp_prob.query("z_tr",z_tr);

    Real x_c  = zero  ; pp_prob.query("x_c", x_c);
    Real z_c  = amrex::Real(1.5e3); pp_prob.query("z_c", z_c);
    Real x_r  = amrex::Real(4.0e3); pp_prob.query("x_r", x_r);
    Real z_r  = amrex::Real(1.5e3); pp_prob.query("z_r", z_r);
    Real T_tr = amrex::Real(213.0); pp_prob.query("T_tr", T_tr);
    Real theta_c = three; pp_prob.query("theta_c", theta_c);

    Real theta_0  = amrex::Real(300.0); pp_prob.query("theta_0", theta_0);
    Real theta_tr = amrex::Real(343.0); pp_prob.query("theta_tr", theta_tr);

    Real r_c = one;

    HSEutils::init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(),h_q_v.data(),dz,khi,
                                              q_t,eq_pot_temp,use_empirical,true,height,z_tr,
                                              theta_0, theta_tr, T_tr);

     amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
     amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
     amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
     amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

    Real* t   = d_t.data();
    Real* p   = d_p.data();

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real z        = prob_lo[2] + (k + myhalf) * dx[2];
        const Real x        = prob_lo[0] + (i + myhalf) * dx[0];
        Real rad, delta_theta, theta_total, temperature, rho, RH, scalar;

        Real scaled_height = z / z_tr;
        int which_zone = -1;
        if (z <= height) {
            which_zone = 1;
        } else if (z <= z_tr) {
            which_zone = 2;
        } else {
            which_zone = 3;
        }

        // Introduce the warm bubble. Assume that the bubble is pressure matche with the background

        rad = std::pow(std::pow((x - x_c)/x_r,2) + std::pow((z - z_c)/z_r,2), myhalf);

    if(rad <= r_c){
        delta_theta = theta_c*std::pow(cos(PI*rad/two),2);
        scalar = std::pow(cos(PI*rad/two),2);
    }
    else{
        delta_theta = zero;
        scalar = zero;
    }

    theta_total     = t[k] + delta_theta;

    temperature     = getTgivenPandTh(p[k], theta_total, (R_d/Cp_d));
    Real T_b        = getTgivenPandTh(p[k], t[k]       , (R_d/Cp_d));

    RH = HSEutils::compute_relative_humidity(p[k], T_b, use_empirical, which_zone, scaled_height);

    Real q_v_hot;
    if (z < height) {
        q_v_hot = amrex::Real(0.014);
    } else {
        q_v_hot = HSEutils::vapor_mixing_ratio(p[k], T_b, RH, use_empirical, which_zone);
    }

    rho = p[k]/(R_d*temperature*(one + (R_v/R_d)*q_v_hot));

    // Compute background quantities
    Real temperature_back = getTgivenPandTh(p[k], t[k], (R_d/Cp_d));
    Real T_back           = getTgivenPandTh(p[k], t[k], (R_d/Cp_d));

    Real RH_back          = HSEutils::compute_relative_humidity(p[k], T_back, use_empirical, which_zone, scaled_height);

    Real q_v_back;
    if (z < height) {
        q_v_back = amrex::Real(0.014);
    } else {
        q_v_back = HSEutils::vapor_mixing_ratio(p[k], T_back, RH_back, use_empirical, which_zone);
    }
    Real rho_back         = p[k]/(R_d*temperature_back*(one + (R_v/R_d)*q_v_back));

    // This version perturbs rho but not p
    state_pert(i, j, k, RhoTheta_comp) = rho*theta_total - rho_back*t[k]*(one + (R_v/R_d)*q_v_back);// rho*d_t[k]*(one + R_v_by_R_d*q_v_hot);
    state_pert(i, j, k, Rho_comp)      = rho - rho_back*(one + q_v_back);

    // Set scalar = 0 everywhere
    state_pert(i, j, k, RhoScalar_comp) = rho*scalar;

    // mean states
    if (use_moisture) {
        state_pert(i, j, k, RhoQ1_comp) = rho*q_v_hot;
        state_pert(i, j, k, RhoQ2_comp) = zero;
    }

  });

