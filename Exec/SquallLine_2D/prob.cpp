#include "prob.H"
#include <ERF_Constants.H>

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem ()
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("z_tr", parms.z_tr);
  pp.query("height", parms.height);
  pp.query("theta_0", parms.theta_0);
  pp.query("theta_tr", parms.theta_tr);
  pp.query("T_tr", parms.T_tr);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("theta_c", parms.theta_c);
}

Real Problem::compute_theta (const Real z)
{
    Real theta_0 = parms.theta_0, theta_tr = parms.theta_tr, T_tr = parms.T_tr, z_tr = parms.z_tr;
    if(z <= z_tr) {
        return theta_0 + (theta_tr - theta_0)*std::pow(z/z_tr,1.25);
    } else {
        return theta_tr*exp(CONST_GRAV/(Cp_d*T_tr)*(z - z_tr));
    }

}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_saturation_pressure (const Real T_b)
{

    Real p_s = exp(34.494 - 4924.99/(T_b - 273.15 + 237.1))/std::pow(T_b - 273.15 + 105.0,1.57);

    //Real T = T_b - 273.15;

    //Real p_s = 0.61121e3*exp((18.678 - T/234.5)*(T/(257.14 + T)));

    return p_s;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_relative_humidity (const Real z, const Real height, const Real z_tr, const Real p_b, const Real T_b)
{
    Real p_s = compute_saturation_pressure(T_b);

    Real q_s = Rd_on_Rv*p_s/(p_b - p_s);

    if(z <= height){
        return 0.014/q_s;
    }else if(z <= z_tr){
        return 1.0 - 0.75*std::pow(z/z_tr,1.25);
    }else{
        return 0.25;
    }
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_vapor_pressure (const Real p_s, const Real RH)
{
    return p_s*RH;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real vapor_mixing_ratio (const Real z, const Real height, const Real p_b, const Real T_b, const Real RH)
{
    Real p_s = compute_saturation_pressure(T_b);
    Real p_v = compute_vapor_pressure(p_s, RH);

    Real q_v = Rd_on_Rv*p_v/(p_b - p_v);

    if(z < height){
        return 0.014;
    }else{
        return q_v;
    }
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_temperature (const Real p_b, const Real theta_b)
{
    return theta_b*std::pow(p_b/p_0,R_d/Cp_d);
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

void Problem::compute_rho (const Real& z, const Real& pressure, Real& theta, Real& rho, Real& q_v, Real& T_dp, Real& T_b)
{

    theta   = compute_theta(z);
    T_b     = compute_temperature(pressure, theta);
    Real RH = compute_relative_humidity(z, parms.height, parms.z_tr, pressure, T_b);
    q_v     = vapor_mixing_ratio(z, parms.height, pressure, T_b, RH);
    rho     = pressure/(R_d*T_b*(1.0 + (R_v/R_d)*q_v));
    rho     = rho*(1.0 + q_v);
    T_dp    = compute_dewpoint_temperature(T_b, RH);
}

Real Problem::compute_F (const Real& p_k, const Real& p_k_minus_1, Real &theta_k, Real& rho_k, Real& q_v_k,
                  Real& T_dp, Real& T_b, const Real& dz, const Real& z, const Real& rho_k_minus_1)
{

    Real F;

    if(rho_k_minus_1 == 0.0) // This loop is for the first point above the ground
    {
        compute_rho(z, p_k, theta_k, rho_k, q_v_k, T_dp, T_b);
        F = p_k - p_k_minus_1 + rho_k*CONST_GRAV*dz/2.0;
    }
    else
    {
        compute_rho(z, p_k, theta_k, rho_k, q_v_k, T_dp, T_b);
        F = p_k - p_k_minus_1 + rho_k_minus_1*CONST_GRAV*dz/2.0 + rho_k*CONST_GRAV*dz/2.0;
    }

    return F;
}

Real Problem::compute_p_k (Real &p_k, const Real p_k_minus_1, Real &theta_k, Real& rho_k, Real& q_v_k,
                  Real& T_dp, Real& T_b, const Real dz, const Real z, const Real rho_k_minus_1)
{

    for(int iter=0;iter<20;iter++)
    {
        Real F         = compute_F(p_k, p_k_minus_1, theta_k, rho_k, q_v_k, T_dp, T_b, dz, z, rho_k_minus_1);
        Real F_plus_dF = compute_F(p_k+1e-10, p_k_minus_1, theta_k, rho_k, q_v_k, T_dp, T_b, dz, z, rho_k_minus_1);
        Real F_prime   = (F_plus_dF - F)/1e-10;
        Real delta_p_k = -F/F_prime;
        p_k            = p_k + delta_p_k;
    }

    return p_k;
}

void
Problem::init_isentropic_hse_no_terrain(Real *theta, Real* r, Real* p, Real *q_v,
                               const Real& dz, const Real&  prob_lo_z,
                               const int& khi)
{

    //FILE *file_IC;
    //file_IC = fopen("input_sounding_probcpp.txt","w");
    Real z, T_b, T_dp;

    // Compute the quantities at z = 0.5*dz (first cell center)
    z = prob_lo_z + 0.5*dz;
    p[0] = p_0;
    compute_p_k(p[0], p_0, theta[0], r[0], q_v[0], T_dp, T_b, dz, z, 0.0);
    //fprintf(file_IC, "%0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g\n", z, T_b-273.15, T_dp, p[0], r[0], theta[0], q_v[0]);


    for (int k=1;k<=khi;k++){
        z = prob_lo_z + (k+0.5)*dz;
        p[k] = p[k-1];
        compute_p_k(p[k], p[k-1], theta[k], r[k], q_v[k], T_dp, T_b, dz, z, r[k-1]);
        //fprintf(file_IC, "%0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g\n", z, T_b-273.15, T_dp, p[k], r[k], theta[k], q_v[k]);
    }
    //fclose(file_IC);


    r[khi+1] = r[khi];
}

void
Problem::erf_init_dens_hse_moist (MultiFab& rho_hse,
                                  std::unique_ptr<MultiFab>& z_phys_nd,
                                  Geometry const& geom)
{


    const Real prob_lo_z = geom.ProbLo()[2];
    const Real dz        = geom.CellSize()[2];
    const int khi        = geom.Domain().bigEnd()[2];

    // use_terrain = 1
    if (z_phys_nd) {

        if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

        for ( MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > rho_arr  = rho_hse.array(mfi);
            //Array4<Real const> z_cc_arr = z_phys_cc->const_array(mfi);

            // Create a flat box with same horizontal extent but only one cell in vertical
            const Box& tbz = mfi.nodaltilebox(2);
            Box b2d = tbz; // Copy constructor
            b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
            b2d.setRange(2,0);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
              Array1D<Real,0,255> r;

              //init_isentropic_hse_terrain(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

              for (int k = 0; k <= khi; k++) {
                 rho_arr(i,j,k) = r(k);
              }
              rho_arr(i,j,   -1) = rho_arr(i,j,0);
              rho_arr(i,j,khi+1) = rho_arr(i,j,khi);
            });
        } // mfi
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

        init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(), h_q_v.data(), dz,prob_lo_z,khi);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

        Real* r     = d_r.data();

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
Problem::init_custom_pert (
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& state,
    Array4<Real      > const& x_vel,
    Array4<Real      > const& y_vel,
    Array4<Real      > const& z_vel,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];

  // Call the routine to calculate the 1d background condition

   Vector<Real> h_r(khi+2);
   Vector<Real> h_p(khi+2);
   Vector<Real> h_t(khi+2);
   Vector<Real> h_q_v(khi+2);

   amrex::Gpu::DeviceVector<Real> d_r(khi+2);
   amrex::Gpu::DeviceVector<Real> d_p(khi+2);
   amrex::Gpu::DeviceVector<Real> d_t(khi+2);
   amrex::Gpu::DeviceVector<Real> d_q_v(khi+2);

   init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(),h_q_v.data(),dz,prob_lo_z,khi);

   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

    Real* t   = d_t.data();
    Real* p   = d_p.data();


    const Real x_c = parms.x_c, z_c = parms.z_c, x_r = parms.x_r, z_r = parms.z_r, theta_c = parms.theta_c, r_c = 1.0;
    //const Real x_c = 0.0, z_c = 2.0e3, x_r = 10.0e3, z_r = 1.5e3, r_c = 1.0, theta_c = 3.0;

    Real Rd_by_Cp = sc.rdOcp;
    Rd_by_Cp = Rd_by_Cp;

    Real height = parms.height;
    Real z_tr = parms.z_tr;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo         = geomdata.ProbLo();
    const auto dx              = geomdata.CellSize();
    const amrex::Real z        = prob_lo[2] + (k + 0.5) * dx[2];
    const amrex::Real x        = prob_lo[0] + (i + 0.5) * dx[0];
    Real rad, delta_theta, theta_total, temperature, rho, RH, scalar;

    // Introduce the warm bubble. Assume that the bubble is pressure matche with the background

    rad = std::pow(std::pow((x - x_c)/x_r,2) + std::pow((z - z_c)/z_r,2), 0.5);

    if(rad <= r_c){
        delta_theta = theta_c*std::pow(cos(PI*rad/2.0),2);
        scalar = std::pow(cos(PI*rad/2.0),2);
    }
    else{
        delta_theta = 0.0;
        scalar = 0.0;
    }

    theta_total     = t[k] + delta_theta;
    temperature     = compute_temperature(p[k], theta_total);
    Real T_b        = compute_temperature(p[k], t[k]);
    RH              = compute_relative_humidity(z, height, z_tr, p[k], T_b);
    Real q_v_hot    = vapor_mixing_ratio(z, height, p[k], T_b, RH);
    rho             = p[k]/(R_d*temperature*(1.0 + (R_v/R_d)*q_v_hot));

    // Compute background quantities
    Real temperature_back = compute_temperature(p[k], t[k]);
    Real T_back           = compute_temperature(p[k], t[k]);
    Real RH_back          = compute_relative_humidity(z, height, z_tr, p[k], T_back);
    Real q_v_back         = vapor_mixing_ratio(z, height, p[k], T_back, RH_back);
    Real rho_back         = p[k]/(R_d*temperature_back*(1.0 + (R_v/R_d)*q_v_back));

    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = rho*theta_total - rho_back*t[k]*(1.0 + (R_v/R_d)*q_v_back);// rho*d_t[k]*(1.0 + R_v_by_R_d*q_v_hot);
    state(i, j, k, Rho_comp)      = rho - rho_back*(1.0 + q_v_back);

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = rho*scalar;

    // mean states
    if (use_moisture) {
        state(i, j, k, RhoQ1_comp) = rho*q_v_hot;
        state(i, j, k, RhoQ2_comp) = 0.0;
    }

  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real z = prob_lo_z + (k+0.5) * dz;
    x_vel(i,j,k) = -12.0*std::max(0.0, (2.5e3 - z)/2.5e3);
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
}

void
Problem::erf_init_rayleigh (amrex::Vector<amrex::Real>& tau,
                            amrex::Vector<amrex::Real>& ubar,
                            amrex::Vector<amrex::Real>& vbar,
                            amrex::Vector<amrex::Real>& wbar,
                            amrex::Vector<amrex::Real>& thetabar,
                            amrex::Geometry      const& geom,
                            std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      tau[k]  = 1.0;
      ubar[k] = 2.0;
      vbar[k] = 1.0;
      wbar[k] = 0.0;
      thetabar[k] = parms.Theta_0;
  }
}


