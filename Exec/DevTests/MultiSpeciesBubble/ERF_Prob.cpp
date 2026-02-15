#include "ERF_Prob.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_HSEUtils.H"

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

  pp.query("do_moist_bubble", parms.do_moist_bubble);
  pp.query("theta_pert", parms.theta_pert);

  // check for SDM species
  m_species.clear();
  ParmParse pp_sdm("super_droplets_moisture");
  std::string species_input = "species";
  if (pp_sdm.contains(species_input.c_str())) {
      int num_species = pp_sdm.countval(species_input.c_str());
      std::string sp_name;
      for (int i = 0; i < num_species; i++) {
          pp_sdm.get(species_input.c_str(), sp_name, i);
          m_species.push_back(sp_name);
      }
  }
  m_num_species = m_species.size();

  m_qv_init_species = std::vector<amrex::Real>(m_num_species,0.0);
  for (int i = 0; i < m_num_species; i++) {
      std::string key_str = "qv_init_" + m_species[i];
      pp.query(key_str.c_str(), m_qv_init_species[i]);
  }

  init_base_parms(parms.rho_0, parms.T_0);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real compute_relative_humidity ()
{
    return 1.0;
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

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& state,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& sc,
    const int /*lev*/ )
{
    AMREX_ALWAYS_ASSERT(sc.moisture_type == MoistureType::SuperDroplets);

    const auto khi     = geomdata.Domain().bigEnd()[2];
    const Real dz      = geomdata.CellSize()[2];
    const Real rdOcp   = sc.rdOcp;

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

    ParmParse pp("prob");
    Real q_t = 0.02;
    pp.query("qt_init", q_t);

    Real eq_pot_temp = 320.0;
    pp.query("eq_pot_temp",eq_pot_temp);

    bool use_empirical = false;
    pp.query("use_empircal_psat",use_empirical);

    if (parms.do_moist_bubble) {
        Vector<Real> h_r(khi+2);
        Vector<Real> h_p(khi+2);
        Vector<Real> h_t(khi+2);
        Vector<Real> h_q_v(khi+2);

        Gpu::DeviceVector<Real> d_r(khi+2);
        Gpu::DeviceVector<Real> d_p(khi+2);
        Gpu::DeviceVector<Real> d_t(khi+2);
        Gpu::DeviceVector<Real> d_q_v(khi+2);

        HSEutils::init_isentropic_hse_no_terrain(h_t.data(), h_r.data(), h_p.data(),
                                                 h_q_v.data(), dz, khi, q_t, eq_pot_temp, use_empirical);

        Gpu::copy(Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
        Gpu::copy(Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
        Gpu::copy(Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
        Gpu::copy(Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

        Real* theta_back   = d_t.data();
        Real* p_back   = d_p.data();
        Real* q_v_back = d_q_v.data();

        const auto x_c = parms.x_c, y_c = parms.y_c, z_c = parms.z_c;
        const auto x_r = parms.x_r, y_r = parms.y_r, z_r = parms.z_r;
        const auto theta_pert = parms.theta_pert;

        int which_zone = -1;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const auto x = prob_lo[0] + (i + 0.5) * dx[0];
            const auto y = prob_lo[1] + (j + 0.5) * dx[1];
            const auto z = prob_lo[2] + (k + 0.5) * dx[2];

            Real rad, delta_theta, theta_total, rho, RH;

            // Introduce the warm bubble.
            // Assume that the bubble is pressure matched with the background
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
            Real T = getTgivenPandTh(p_back[k], theta_total, (R_d/Cp_d));
            rho    = p_back[k]/(R_d*T*(1.0 + (R_v/R_d)*q_v_back[k]));
            RH     = compute_relative_humidity();
            Real q_v_hot = HSEutils::vapor_mixing_ratio(p_back[k], T, RH, use_empirical, which_zone);

            // Compute background quantities
            Real T_back   = getTgivenPandTh(p_back[k], theta_back[k], (R_d/Cp_d));
            Real rho_back = p_back[k]/(R_d*T_back*(1.0 + (R_v/R_d)*q_v_back[k]));

            // This version perturbs rho but not p
            state_pert(i, j, k, RhoTheta_comp) = rho*theta_total - rho_back*theta_back[k]*(1.0 + (R_v/R_d)*q_v_back[k]);
            state_pert(i, j, k, Rho_comp)      = rho - rho_back*(1.0 + q_t);

            // mean states
            state_pert(i, j, k, RhoQ1_comp) = rho*q_v_hot;
        });
    } else {
        ParallelFor(bx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const Real z = prob_lo[2] + (k + 0.5) * dx[2];

            perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k),
                              parms_d, rdOcp,
                              state_pert(i, j, k, Rho_comp),
                              state_pert(i, j, k, RhoTheta_comp));
        });
    } // do_moist_bubble

    // if other species are present
    auto n_sp = m_num_species;
    auto nstart = m_nstart_sp;
    int ncomp_species = state_pert.nComp() - nstart;
    AMREX_ALWAYS_ASSERT(n_sp*m_ncomp_per_species == ncomp_species);
    Gpu::DeviceVector<Real> qv_init_d(n_sp);
    Gpu::copy( Gpu::hostToDevice,
               m_qv_init_species.begin(),
               m_qv_init_species.end(),
               qv_init_d.begin() );
    auto qv_arr = qv_init_d.data();

    ParallelFor(bx, [=, parms_d=parms]
                AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        amrex::Real L = 0.0;
        if (parms_d.x_r > 0) L += std::pow((x - parms_d.x_c)/parms_d.x_r, 2);
        if (parms_d.y_r > 0) L += std::pow((y - parms_d.y_c)/parms_d.y_r, 2);
        if (parms_d.z_r > 0) L += std::pow((z - parms_d.z_c)/parms_d.z_r, 2);
        L = std::sqrt(L);

        if (L < 1.0) {
            auto rho = state(i,j,k,Rho_comp) + state_pert(i,j,k,Rho_comp);
            for (int ns = 0; ns < n_sp; ns++) {
                state_pert(i, j, k, nstart+2*ns)   = rho*qv_arr[ns];
                state_pert(i, j, k, nstart+2*ns+1) = 0.0;
            }
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
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    const Real u0 = parms.U_0;
    const Real v0 = parms.V_0;
    const Real w0 = parms.W_0;

    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = u0;
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = v0;
    });
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = w0;
    });

    Gpu::streamSynchronize();
}
