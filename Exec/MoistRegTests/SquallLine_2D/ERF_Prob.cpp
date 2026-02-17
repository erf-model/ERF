#include "ERF_Prob.H"
#include <ERF_Constants.H>
#include <ERF_HSEUtils.H>

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
  pp.query("theta_0", parms.theta_0);
  pp.query("theta_tr", parms.theta_tr);
  pp.query("T_tr", parms.T_tr);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("theta_c", parms.theta_c);
}

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& sc,
    const int /*lev*/)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real& dz        = geomdata.CellSize()[2];

  // Call the routine to calculate the 1d background condition

   Vector<Real> h_r(khi+2);
   Vector<Real> h_p(khi+2);
   Vector<Real> h_t(khi+2);
   Vector<Real> h_q_v(khi+2);

   amrex::Gpu::DeviceVector<Real> d_r(khi+2);
   amrex::Gpu::DeviceVector<Real> d_p(khi+2);
   amrex::Gpu::DeviceVector<Real> d_t(khi+2);
   amrex::Gpu::DeviceVector<Real> d_q_v(khi+2);

   ParmParse pp("prob");
   Real q_t = 0.02;
   pp.query("qt_init",q_t);

   Real eq_pot_temp = 320.0;
   pp.query("eq_pot_temp",eq_pot_temp);

   bool use_empirical = false;
   pp.query("use_empirical_psat",use_empirical);

   Real height = -1.;
   pp.query("height",height);

   Real z_tr = -1.;
   pp.query("z_tr",z_tr);

   Real theta_tr = 0.;
   pp.query("theta_tr",theta_tr);

   Real theta_0 = 0.;
   pp.query("theta_0",theta_0);

   Real T_tr = 0.;
   pp.query("T_tr",T_tr);

   HSEutils::init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(),h_q_v.data(),dz,khi,
                                            q_t,eq_pot_temp,use_empirical,true,height,z_tr,
                                            theta_0, theta_tr, T_tr);

   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

    Real* t   = d_t.data();
    Real* p   = d_p.data();

    const Real x_c = parms.x_c, z_c = parms.z_c, x_r = parms.x_r, z_r = parms.z_r, theta_c = parms.theta_c, r_c = 1.0;

  ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real z        = prob_lo[2] + (k + 0.5) * dx[2];
    const Real x        = prob_lo[0] + (i + 0.5) * dx[0];
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

    temperature     = getTgivenPandTh(p[k], theta_total, (R_d/Cp_d));
    Real T_b        = getTgivenPandTh(p[k], t[k]       , (R_d/Cp_d));

    RH = HSEutils::compute_relative_humidity(p[k], T_b, use_empirical, which_zone, scaled_height);

    Real q_v_hot;
    if (z < height) {
        q_v_hot = 0.014;
    } else {
        q_v_hot = HSEutils::vapor_mixing_ratio(p[k], T_b, RH, use_empirical, which_zone);
    }

    rho = p[k]/(R_d*temperature*(1.0 + (R_v/R_d)*q_v_hot));

    // Compute background quantities
    Real temperature_back = getTgivenPandTh(p[k], t[k], (R_d/Cp_d));
    Real T_back           = getTgivenPandTh(p[k], t[k], (R_d/Cp_d));

    Real RH_back          = HSEutils::compute_relative_humidity(p[k], T_back, use_empirical, which_zone, scaled_height);

    Real q_v_back;
    if (z < height) {
        q_v_back = 0.014;
    } else {
        q_v_back = HSEutils::vapor_mixing_ratio(p[k], T_back, RH_back, use_empirical, which_zone);
    }
    Real rho_back         = p[k]/(R_d*temperature_back*(1.0 + (R_v/R_d)*q_v_back));

    // This version perturbs rho but not p
    state_pert(i, j, k, RhoTheta_comp) = rho*theta_total - rho_back*t[k]*(1.0 + (R_v/R_d)*q_v_back);// rho*d_t[k]*(1.0 + R_v_by_R_d*q_v_hot);
    state_pert(i, j, k, Rho_comp)      = rho - rho_back*(1.0 + q_v_back);

    // Set scalar = 0 everywhere
    state_pert(i, j, k, RhoScalar_comp) = rho*scalar;

    // mean states
    if (use_moisture) {
        state_pert(i, j, k, RhoQ1_comp) = rho*q_v_hot;
        state_pert(i, j, k, RhoQ2_comp) = 0.0;
    }

  });

  Gpu::streamSynchronize();
}

void
Problem::init_custom_pert_vels (
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real const> const& /*z_nd*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];

  // Set the x-velocity
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real z = prob_lo_z + (k+0.5) * dz;
      x_vel_pert(i,j,k) = -12.0*std::max(0.0, (2.5e3 - z)/2.5e3);
  });

  // Set the y-velocity
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      y_vel_pert(i, j, k) = 0.0;
  });

  // Set the z-velocity
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      z_vel_pert(i, j, k) = 0.0;
  });

  Gpu::streamSynchronize();
}

void
Problem::erf_init_rayleigh(
    amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
    amrex::Geometry const& geom,
    std::unique_ptr<MultiFab>& /*z_phys_nd*/,
    amrex::Real /*zdamp*/)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      rayleigh_ptrs[Rayleigh::ubar][k]     = 2.0;
      rayleigh_ptrs[Rayleigh::vbar][k]     = 1.0;
      rayleigh_ptrs[Rayleigh::wbar][k]     = 0.0;
      rayleigh_ptrs[Rayleigh::thetabar][k] = parms.T_0;
  }
}
