#include "prob.H"
#include "Microphysics_Utils.H"
#include "ERF_Constants.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  // Parse params
  ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("x_c", parms.x_c);
  pp.query("y_c", parms.y_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("y_r", parms.y_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
  pp.query("T_pert_is_airtemp", parms.T_pert_is_airtemp);
  pp.query("perturb_rho", parms.perturb_rho);
  pp.query("dampcoef", parms.dampcoef);
  pp.query("zdamp", parms.zdamp);

  pp.query("do_moist_bubble", parms.do_moist_bubble);
  pp.query("theta_pert", parms.theta_pert);
  pp.query("qt_init", parms.qt_init);
  pp.query("eq_pot_temp", parms.eq_pot_temp);

  init_base_parms(parms.rho_0, parms.T_0);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_saturation_pressure (const Real T_b)
{
    return erf_esatw(T_b)*100.0;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_relative_humidity ()
{
    return 1.0;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_vapor_pressure (const Real p_s, const Real RH)
{
    return p_s*RH;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real vapor_mixing_ratio (const Real p_b, const Real T_b, const Real RH)
{
    Real p_s = compute_saturation_pressure(T_b);
    Real p_v = compute_vapor_pressure(p_s, RH);
    Real q_v = Rd_on_Rv*p_v/(p_b - p_v);
    return q_v;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real Problem::compute_F_for_temp(const Real T_b, const Real p_b)
{
    Real q_t = parms.qt_init;
    Real fac = Cp_d + Cp_l*q_t;
    Real RH  = compute_relative_humidity();
    Real q_v = vapor_mixing_ratio(p_b, T_b, RH);
    Real p_s = compute_saturation_pressure(T_b);
    Real p_v = compute_vapor_pressure(p_s, RH);
    return  parms.eq_pot_temp - T_b*std::pow((p_b - p_v)/p_0, -R_d/fac)*std::exp(L_v*q_v/(fac*T_b));
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real Problem::compute_temperature (const Real p_b)
{
    Real T_b = 200.0, delta_T; // Initial guess
    for(int iter=0; iter<20; iter++)
    {
        Real F         = compute_F_for_temp(T_b, p_b);
        Real F_plus_dF = compute_F_for_temp(T_b+1e-10, p_b);
        Real F_prime   = (F_plus_dF - F)/1e-10;
        delta_T = -F/F_prime;
        T_b = T_b + delta_T;
    }

    if(std::fabs(delta_T) > 1e-8){
         amrex::Abort("Newton Raphson for temperature could not converge");
    }

    return T_b;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_dewpoint_temperature (const Real T_b, const Real RH)
{
    Real T_dp, gamma, T;
    T = T_b - 273.15;

    Real b = 18.678, c = 257.14, d = 234.5;
    gamma = log(RH*exp((b - T/d)*T/(c + T)));

    T_dp = c*gamma/(b - gamma);

    return T_dp;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real Problem::compute_theta(const Real T_b, const Real p_b)
{
    return T_b*std::pow(p_0/p_b, R_d/Cp_d);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_temperature_from_theta(const Real theta, const Real p)
{
    return theta*std::pow(p/p_0, R_d/Cp_d);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
void Problem::compute_rho (const Real& pressure, Real& theta, Real& rho, Real& q_v, Real& T_dp, Real& T_b)
{
    T_b     = compute_temperature(pressure);
    theta   = compute_theta(T_b, pressure);
    Real RH = compute_relative_humidity();
    q_v     = vapor_mixing_ratio(pressure, T_b, RH);
    rho     = pressure/(R_d*T_b*(1.0 + (R_v/R_d)*q_v));
    rho     = rho*(1.0 + parms.qt_init); // q_t = 0.02 a given constant for this test case
    T_dp    = compute_dewpoint_temperature(T_b, RH);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real Problem::compute_F (const Real& p_k, const Real& p_k_minus_1, Real &theta_k, Real& rho_k, Real& q_v_k,
                         Real& T_dp, Real& T_b, const Real& dz, const Real& rho_k_minus_1)
{
    Real F;

    if(rho_k_minus_1 == 0.0) // This loop is for the first point above the ground
    {
        compute_rho(p_k, theta_k, rho_k, q_v_k, T_dp, T_b);
        F = p_k - p_k_minus_1 + rho_k*CONST_GRAV*dz/2.0;
    }
    else
    {
        compute_rho(p_k, theta_k, rho_k, q_v_k, T_dp, T_b);
        F = p_k - p_k_minus_1 + rho_k_minus_1*CONST_GRAV*dz/2.0 + rho_k*CONST_GRAV*dz/2.0;
    }

    return F;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real Problem::compute_p_k (Real &p_k, const Real p_k_minus_1, Real &theta_k, Real& rho_k, Real& q_v_k,
                           Real& T_dp, Real& T_b, const Real dz, const Real rho_k_minus_1)
{
    Real delta_p_k;

    for(int iter=0; iter<20; iter++)
    {
        Real F         = compute_F(p_k, p_k_minus_1, theta_k, rho_k, q_v_k, T_dp, T_b, dz, rho_k_minus_1);
        Real F_plus_dF = compute_F(p_k+1e-10, p_k_minus_1, theta_k, rho_k, q_v_k, T_dp, T_b, dz, rho_k_minus_1);
        Real F_prime   = (F_plus_dF - F)/1e-10;
        delta_p_k = -F/F_prime;
        p_k            = p_k + delta_p_k;
    }

    if(std::fabs(delta_p_k) > 1e-8){
        amrex::Abort("Newton Raphson for pressure could not converge");
    }

    return p_k;
}

void
Problem::init_isentropic_hse_no_terrain(Real *theta, Real* r, Real* p, Real *q_v,
                                        const Real& dz,
                                        const int& khi)
{
    Real T_b, T_dp;

    // Compute the quantities at z = 0.5*dz (first cell center)
    p[0] = p_0;
    compute_p_k(p[0], p_0, theta[0], r[0], q_v[0], T_dp, T_b, dz, 0.0);

    for (int k=1;k<=khi;k++){
        p[k] = p[k-1];
        compute_p_k(p[k], p[k-1], theta[k], r[k], q_v[k], T_dp, T_b, dz, r[k-1]);
    }

    r[khi+1] = r[khi];
}

void
Problem::erf_init_dens_hse_moist (MultiFab& rho_hse,
                                  std::unique_ptr<MultiFab>& z_phys_nd,
                                  Geometry const& geom)
{
    const Real dz        = geom.CellSize()[2];
    const int khi        = geom.Domain().bigEnd()[2];

    // use_terrain = 1
    if (z_phys_nd) {
        amrex::Abort("Terrain not implemented for moist Bubble case!");
    } else { // use_terrain = 0

        // These are at cell centers (unstaggered)
        Vector<Real> h_r(khi+2);
        Vector<Real> h_p(khi+2);
        Vector<Real> h_t(khi+2);
        Vector<Real> h_q_v(khi+2);

        amrex::Gpu::DeviceVector<Real> d_r(khi+2);
        amrex::Gpu::DeviceVector<Real> d_p(khi+2);
        amrex::Gpu::DeviceVector<Real> d_t(khi+2);
        amrex::Gpu::DeviceVector<Real> d_q_v(khi+2);

        init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(), h_q_v.data(), dz, khi);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

        Real* r = d_r.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
          for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              const Box& bx = mfi.growntilebox(1);
              const Array4<Real> rho_hse_arr   = rho_hse[mfi].array();
              ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
              {
                  int kk = std::max(k,0);
                  rho_hse_arr(i,j,k) = r[kk];
              });
          } // mfi
    } // no terrain
}


void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const int khi = geomdata.Domain().bigEnd()[2];

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    const Real dz        = geomdata.CellSize()[2];
    const Real rdOcp     = sc.rdOcp;

    amrex::Print() << "Bubble delta T = " << parms.T_pert << " K" << std::endl;
    amrex::Print() << "  centered at ("
                   << parms.x_c << " " << parms.y_c << " " << parms.z_c << ")" << std::endl;
    amrex::Print() << "  with extent ("
                   << parms.x_r << " " << parms.y_r << " " << parms.z_r << ")" << std::endl;

    if (parms.T_0 <= 0)
    {
        amrex::Print() << "Ignoring parms.T_0 = " << parms.T_0
                       << ", background fields should have been initialized with erf.init_type"
                       << std::endl;
    }

    if (z_cc) { // nonflat terrain
        if (parms.do_moist_bubble) {
            amrex::Abort("Terrain not implemented for moist Bubble case!");
        } else {
            amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Geometry (note we must include these here to get the data on device)
                const auto prob_lo         = geomdata.ProbLo();
                const auto dx              = geomdata.CellSize();

                const Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const Real y = prob_lo[1] + (j + 0.5) * dx[1];
                const Real z = z_cc(i,j,k);

                perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k),
                                  parms, rdOcp,
                                  state_pert(i, j, k, Rho_comp),
                                  state_pert(i, j, k, RhoTheta_comp));

                state_pert(i, j, k, RhoScalar_comp) = 0.0;

                if (use_moisture) {
                    state_pert(i, j, k, RhoQ1_comp) = 0.0;
                    state_pert(i, j, k, RhoQ2_comp) = 0.0;
                    state_pert(i, j, k, RhoQ3_comp) = 0.0;
                }

            });
        } // do_moist_bubble
    } else {
        if (parms.do_moist_bubble) {
            Vector<Real> h_r(khi+2);
            Vector<Real> h_p(khi+2);
            Vector<Real> h_t(khi+2);
            Vector<Real> h_q_v(khi+2);

            amrex::Gpu::DeviceVector<Real> d_r(khi+2);
            amrex::Gpu::DeviceVector<Real> d_p(khi+2);
            amrex::Gpu::DeviceVector<Real> d_t(khi+2);
            amrex::Gpu::DeviceVector<Real> d_q_v(khi+2);

            init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(),h_q_v.data(),dz, khi);

            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

            Real* theta_back   = d_t.data();
            Real* p_back   = d_p.data();
            Real* q_v_back = d_q_v.data();


            const Real x_c = parms.x_c, y_c = parms.y_c, z_c = parms.z_c;
            const Real x_r = parms.x_r, y_r = parms.y_r, z_r = parms.z_r;
            const Real theta_pert = parms.theta_pert;

             int moisture_type = 1;

            if(sc.moisture_type == MoistureType::SAM) {
                moisture_type = 1;
            } else if(sc.moisture_type == MoistureType::SAM_NoIce ||
                      sc.moisture_type == MoistureType::SAM_NoPrecip_NoIce) {
                moisture_type = 2;
            }

            amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k)
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

                if(rad <= 1.0){
                    delta_theta = theta_pert*std::pow(cos(PI*rad/2.0),2);
                }
                else{
                    delta_theta = 0.0;
                }

                theta_total  = theta_back[k]*(delta_theta/300.0 + 1);
                Real T = compute_temperature_from_theta(theta_total, p_back[k]);
                rho    = p_back[k]/(R_d*T*(1.0 + (R_v/R_d)*q_v_back[k]));
                RH     = compute_relative_humidity();
                Real q_v_hot = vapor_mixing_ratio(p_back[k], T, RH);

                // Compute background quantities
                Real T_back   = compute_temperature_from_theta(theta_back[k], p_back[k]);
                Real rho_back = p_back[k]/(R_d*T_back*(1.0 + (R_v/R_d)*q_v_back[k]));

                // This version perturbs rho but not p
                state_pert(i, j, k, RhoTheta_comp) = rho*theta_total - rho_back*theta_back[k]*(1.0 + (R_v/R_d)*q_v_back[k]);
                state_pert(i, j, k, Rho_comp)      = rho - rho_back*(1.0 + parms.qt_init);

                // Set scalar = 0 everywhere
                state_pert(i, j, k, RhoScalar_comp) = 0.0;

                // mean states
                state_pert(i, j, k, RhoQ1_comp) = rho*q_v_hot;
                state_pert(i, j, k, RhoQ2_comp) = rho*(parms.qt_init - q_v_hot);
                state_pert(i, j, k, RhoQ3_comp) = 0.0;

                // Cold microphysics are present
                int nstate = state_pert.ncomp;
                if (nstate == NVAR_max) {
                    Real omn;
                    if(moisture_type == 1) {
                        omn = std::max(0.0,std::min(1.0,(T-tbgmin)*a_bg));
                    } else if(moisture_type == 2) {
                        omn = 1.0;
                    }
                    Real qn  = state_pert(i, j, k, RhoQ2_comp);
                    state_pert(i, j, k, RhoQ2_comp) = qn * omn;
                    state_pert(i, j, k, RhoQ3_comp) = qn * (1.0-omn);
                }
            });
        } else {
            amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Geometry (note we must include these here to get the data on device)
                const auto prob_lo         = geomdata.ProbLo();
                const auto dx              = geomdata.CellSize();

                const Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const Real y = prob_lo[1] + (j + 0.5) * dx[1];
                const Real z = prob_lo[2] + (k + 0.5) * dx[2];

                perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k),
                                  parms, rdOcp,
                                  state_pert(i, j, k, Rho_comp),
                                  state_pert(i, j, k, RhoTheta_comp));

                state_pert(i, j, k, RhoScalar_comp) = 0.0;

                if (use_moisture) {
                    state_pert(i, j, k, RhoQ1_comp) = 0.0;
                    state_pert(i, j, k, RhoQ2_comp) = 0.0;
                    state_pert(i, j, k, RhoQ3_comp) = 0.0;
                }
            });
        } // do_moist_bubble
    } // use_terrain

    const Real u0 = parms.U_0;
    const Real v0 = parms.V_0;
    const Real w0 = parms.W_0;

    // Set the x-velocity
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = u0;
    });

    // Set the y-velocity
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = v0;
    });

    // Set the z-velocity
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = w0;
    });

    amrex::Gpu::streamSynchronize();
}
