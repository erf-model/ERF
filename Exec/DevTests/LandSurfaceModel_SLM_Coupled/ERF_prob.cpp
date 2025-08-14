#include "ERF_prob.H"
#include "AMReX_Random.H"
#include <Utils/ERF_ParFunctions.H>
#include <Utils/ERF_Interpolation_1D.H>

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem (const Real* problo, const Real* probhi)
{
    // Parse params
    ParmParse pp("prob");
    pp.query("rho_0", parms.rho_0);
    pp.query("T_0", parms.T_0);
    pp.query("A_0", parms.A_0);
    pp.query("QKE_0", parms.KE_0);

    pp.query("U_0", parms.U_0);
    pp.query("V_0", parms.V_0);
    pp.query("W_0", parms.W_0);
    pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
    pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
    pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
    pp.query("T_0_Pert_Mag", parms.T_0_Pert_Mag);
    pp.query("qv_0_Pert_Mag", parms.qv_0_Pert_Mag);

    pp.query("pert_deltaU", parms.pert_deltaU);
    pp.query("pert_deltaV", parms.pert_deltaV);
    pp.query("pert_periods_U", parms.pert_periods_U);
    pp.query("pert_periods_V", parms.pert_periods_V);
    pp.query("pert_ref_height", parms.pert_ref_height);
    parms.aval = parms.pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
    parms.bval = parms.pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
    parms.ufac = parms.pert_deltaU * std::exp(0.5) / parms.pert_ref_height;
    parms.vfac = parms.pert_deltaV * std::exp(0.5) / parms.pert_ref_height;

    //===========================================================================
    // READ USER-DEFINED INPUTS
    pp.query("advection_heating_rate", parms.advection_heating_rate);
    pp.query("restart_time",parms.restart_time);
    pp.query("source_cutoff", parms.cutoff);
    pp.query("source_cutoff_transition", parms.cutoff_transition);
    pp.query("advection_moisture_rate", parms.advection_moisture_rate);
    pp.query("moisture_source_cutoff", parms.moisture_cutoff);
    pp.query("moisture_source_cutoff_transition", parms.moisture_cutoff_transition);
    pp.query("wbar_sub_max", parms.wbar_sub_max);
    pp.query("wbar_cutoff_max", parms.wbar_cutoff_max);
    pp.query("wbar_cutoff_min", parms.wbar_cutoff_min);
    AMREX_ASSERT_WITH_MESSAGE(parms.wbar_cutoff_min > parms.wbar_cutoff_max, "ERROR: wbar_cutoff_min < wbar_cutoff_max");
    pp.query("custom_TKE", parms.custom_TKE);


    // Forcing files for SLM test
    /*
    pp.query("large_scale_forcing_file", parms.lsf_file);
    pp.query("sounding_file", parms.snd_file);
    pp.query("forcing_timescale", parms.lsf_scale);

    if (parms.lsf_file != "")
    {
        // read the lsf file
        read_forcing_file(parms.lsf_file, nt_lsf, nz_lsf, times_lsf, z_lsf, p_lsf, t_lsf, q_lsf, u_lsf, v_lsf, w_lsf);
    }
    */
    //===========================================================================

    init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert (
    const Box&  bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& state,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& z_nd,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);
    const bool use_terrain = (sc.terrain_type != TerrainType::None);

    const Real rdOcp   = sc.rdOcp;

    // interpolate the LSF data now using the geometry
    // z_cc is on device!
    amrex::Vector<amrex::Real> zlevels_stag;
    /*
    if (use_terrain) {
        int khi = geomdata.Domain().bigEnd()[2] + 2;
        amrex::Gpu::DeviceVector<Real> d_zlevels_stag_vec(khi);
        amrex::Real *d_zlevels_stag = d_zlevels_stag_vec.data();
        ParallelFor(khi, [=] AMREX_GPU_DEVICE (int k) noexcept {
            int i = bx.smallEnd(0);
            int j = bx.smallEnd(1);
            d_zlevels_stag[k] = z_cc(i, j, k);
        });
        zlevels_stag.resize(khi);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_zlevels_stag_vec.begin(), d_zlevels_stag_vec.end(), zlevels_stag.begin());
    }
    interp_forcing(geomdata, zlevels_stag, times_lsf, z_lsf, t_lsf, q_lsf, u_lsf, v_lsf, w_lsf);
    */
   
    ParallelForRNG(bx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = (use_terrain) ? z_cc(i, j, k) : prob_lo[2] + (k + 0.5) * dx[2];

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
        const Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
        const Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
        //
        // Add temperature perturbations -- we want to keep pressure constant
        //     so these effectively end up as density perturbations
        //
        //if ((z <= parms_d.pert_ref_height) && (parms_d.T_0_Pert_Mag != 0.0)) {
        if ((z <= parms_d.pert_ref_height)) {
            Real rhotheta  = state(i,j,k,RhoTheta_comp);
            Real rho       = state(i,j,k,Rho_comp);
            Real qv        = state(i,j,k,RhoQ1_comp) / rho;
            Real Told      = getTgivenRandRTh(rho,rhotheta,qv);
            Real P         = getPgivenRTh(rhotheta,qv);

            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            //Real Tpert    = (rand_double*2.0 - 1.0)*parms_d.T_0_Pert_Mag;
            Real Tpert    = (1.0 - 2.0*rand_double);
            Real Tnew     = Told + 0.1*Tpert*(1.0 - z / parms_d.pert_ref_height);

            Real theta_new = getThgivenTandP(Tnew,P,rdOcp);
            Real rhonew    = getRhogivenThetaPress(theta_new,P,rdOcp,qv);
            //state_pert(i, j, k, Rho_comp) = rhonew - rho;
            state_pert(i, j, k, Rho_comp) = 0.0;

            // Note we do not perturb this
            //state_pert(i, j, k, RhoTheta_comp) = 0.0;
            state_pert(i, j, k, RhoTheta_comp) = rhotheta - (theta_new * rhonew);

            //  Instead of perturbing (rho theta) we perturb T and hold (rho theta) fixed,
            //  which ends up being stored as a perturbation in rho
            // state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*parms_d.T_0_Pert_Mag;
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = parms_d.A_0 * exp(-10.*r*r);

        // Set an initial value for QKE
        if (parms_d.custom_TKE) {
            state_pert(i, j, k, RhoKE_comp) = (1.0 - z/prob_hi[2]) * r_hse(i,j,k);
        } else {
            state_pert(i, j, k, RhoKE_comp) = parms_d.KE_0;
        }

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
            /*
            if ((z <= parms_d.pert_ref_height) && (parms_d.qv_0_Pert_Mag != 0.0))
            {
                Real rhoold = state(i,j,k,Rho_comp);
                Real rhonew = rhoold + state_pert(i,j,k,Rho_comp);

                Real qvold = state(i,j,k,RhoQ1_comp) / rhoold;

                Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
                Real qvnew = qvold + (rand_double*2.0 - 1.0)*parms_d.qv_0_Pert_Mag;

                state_pert(i, j, k, RhoQ1_comp) = rhonew * qvnew - rhoold * qvold;
            }
            */
        }
    });

    // Set the x-velocity
    ParallelForRNG(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = parms_d.U_0;
        if ((z <= parms_d.pert_ref_height) && (parms_d.U_0_Pert_Mag != 0.0))
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real x_vel_prime = (rand_double*2.0 - 1.0)*parms_d.U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
        if (parms_d.pert_deltaU != 0.0)
        {
            const amrex::Real yl = y - prob_lo[1];
            const amrex::Real zl = z / parms_d.pert_ref_height;
            const amrex::Real damp = std::exp(-0.5 * zl * zl);
            x_vel_pert(i, j, k) += parms_d.ufac * damp * z * std::cos(parms_d.aval * yl);
        }
    });

  // Set the y-velocity
  ParallelForRNG(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the y-velocity
      y_vel_pert(i, j, k) = parms_d.V_0;
      if ((z <= parms_d.pert_ref_height) && (parms_d.V_0_Pert_Mag != 0.0))
      {
          Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
          Real y_vel_prime = (rand_double*2.0 - 1.0)*parms_d.V_0_Pert_Mag;
          y_vel_pert(i, j, k) += y_vel_prime;
      }
      if (parms_d.pert_deltaV != 0.0)
      {
          const amrex::Real xl = x - prob_lo[0];
          const amrex::Real zl = z / parms_d.pert_ref_height;
          const amrex::Real damp = std::exp(-0.5 * zl * zl);
          y_vel_pert(i, j, k) += parms_d.vfac * damp * z * std::cos(parms_d.bval * xl);
      }
  });

  // Set the z-velocity
  ParallelForRNG(zbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
  {
      const int dom_lo_z = geomdata.Domain().smallEnd()[2];
      const int dom_hi_z = geomdata.Domain().bigEnd()[2];

      // Set the z-velocity
      if (k == dom_lo_z || k == dom_hi_z+1)
      {
          z_vel_pert(i, j, k) = 0.0;
      }
      else if (parms_d.W_0_Pert_Mag != 0.0)
      {
          Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
          Real z_vel_prime = (rand_double*2.0 - 1.0)*parms_d.W_0_Pert_Mag;
          z_vel_pert(i, j, k) = parms_d.W_0 + z_vel_prime;
      }
  });
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_rhotheta_sources (const Real& time,
                                  Vector<Real>& src,
                                  Gpu::DeviceVector<Real>& d_src,
                                  const Geometry& geom,
                                  std::unique_ptr<MultiFab>& z_phys_cc)
{
    if (src.empty()) return;

    const int khi       = geom.Domain().bigEnd()[2];
    const Real* prob_lo = geom.ProbLo();
    const auto dx       = geom.CellSize();

    // Note: If z_phys_cc, then use_terrain=1 was set. If the z coordinate
    // varies in time and or space, then the the height needs to be
    // calculated at each time step. Here, we assume that only grid
    // stretching exists.
    if (z_phys_cc && zlevels.empty()) {
        Print() << "Initializing z levels on stretched grid" << std::endl;
        zlevels.resize(khi+1);
        reduce_to_max_per_height(zlevels, z_phys_cc);
    }

    int itime_curr = 0;
    int itime_next = 0;
    amrex::Real coeff_curr = 1.0;
    amrex::Real coeff_next = 0.0;
    const bool lsf_forcing = times_lsf.size() > 0;
    amrex::Real inv_scale = 1.0 / parms.lsf_scale;
    if (lsf_forcing)
    {
       get_forcing_time_coeffs(time, itime_curr, itime_next, coeff_curr, coeff_next);
    }

    // Only apply temperature source below nominal inversion height
    for (int k = 0; k <= khi; k++) {
        /*
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)* dx[2];

        if (lsf_forcing)
        {
            // apply tls forcing
            src[k] = t_int_lsf[itime_curr][k] * coeff_curr + t_int_lsf[itime_next][k] * coeff_next;
            amrex::Print() << " rhotheta forcing: k = " << k << ": tls_curr = " << t_int_lsf[itime_curr][k] << " tls_next = " << t_int_lsf[itime_next][k] << " src = " << src[k] << " src_scaled = " << src[k] * inv_scale << std::endl;
            src[k] *= inv_scale;
        } else {
            src[k] = 0.0;
        }
        */
        src[k] = 0.0;
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, src.begin(), src.end(), d_src.begin());
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_rhoqt_sources (const Real& time,
                               Vector<Real>& qsrc,
                               Gpu::DeviceVector<Real>& d_qsrc,
                               const Geometry& geom,
                               std::unique_ptr<MultiFab>& z_phys_cc)
{
    if (qsrc.empty()) return;

    const int khi       = geom.Domain().bigEnd()[2];
    const Real* prob_lo = geom.ProbLo();
    const auto dx       = geom.CellSize();

    // Note: If z_phys_cc, then use_terrain=1 was set. If the z coordinate
    // varies in time and or space, then the the height needs to be
    // calculated at each time step. Here, we assume that only grid
    // stretching exists.
    if (z_phys_cc && zlevels.empty()) {
        Print() << "Initializing z levels on stretched grid" << std::endl;
        zlevels.resize(khi+1);
        reduce_to_max_per_height(zlevels, z_phys_cc);
    }

    int itime_curr = 0;
    int itime_next = 0;
    amrex::Real coeff_curr = 1.0;
    amrex::Real coeff_next = 0.0;
    const bool lsf_forcing = times_lsf.size() > 0;
    amrex::Real inv_scale = 1.0 / parms.lsf_scale;
    if (lsf_forcing)
    {
       get_forcing_time_coeffs(time, itime_curr, itime_next, coeff_curr, coeff_next);
    }

    // Only apply temperature source below nominal inversion height
    for (int k = 0; k <= khi; k++) {
        /*
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)* dx[2];

        if (lsf_forcing)
        {
            // apply qls forcing
            qsrc[k] = q_int_lsf[itime_curr][k] * coeff_curr + q_int_lsf[itime_next][k] * coeff_next;
            //amrex::Print() << " rhoqt forcing: k = " << k << ": qls_curr = " << q_int_lsf[itime_curr][k] << " qls_next = " << q_int_lsf[itime_next][k] << " src = " << qsrc[k] << " src_scaled = " << qsrc[k] * inv_scale << std::endl;
            qsrc[k] *= inv_scale;
        } else {
            qsrc[k] = 0.0;
        }
        */
        qsrc[k] = 0.0;
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, qsrc.begin(), qsrc.end(), d_qsrc.begin());
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_w_subsidence (const Real& time,
                              Vector<Real>& wbar,
                              Gpu::DeviceVector<Real>& d_wbar,
                              const Geometry& geom,
                              std::unique_ptr<MultiFab>& z_phys_cc)
{
    if (wbar.empty()) return;

    const int khi       = geom.Domain().bigEnd()[2] + 1; // lives on z-faces
    const Real* prob_lo = geom.ProbLo();
    const auto dx       = geom.CellSize();

    // Note: If z_phys_cc, then use_terrain=1 was set. If the z coordinate
    // varies in time and or space, then the the height needs to be
    // calculated at each time step. Here, we assume that only grid
    // stretching exists.

    if (z_phys_cc && (zlevels.empty() || zlevels.size() != khi+1)) {
        Print() << "Initializing z levels on stretched grid: w sub" << std::endl;
        zlevels.resize(khi+1);
        reduce_to_max_per_height(zlevels, z_phys_cc);
    }

    /*
    for (int k = 0; k < zlevels.size(); k++)
    {
        Print() << " k = " << k << " zlevel = " << zlevels[k] << std::endl;
    }
    */

    int itime_curr = 0;
    int itime_next = 0;
    amrex::Real coeff_curr = 1.0;
    amrex::Real coeff_next = 0.0;
    const bool lsf_forcing = times_lsf.size() > 0;
    amrex::Real inv_scale = 1.0 / parms.lsf_scale;
    if (lsf_forcing)
    {
       get_forcing_time_coeffs(time, itime_curr, itime_next, coeff_curr, coeff_next);
    }

    for (int k = 0; k <= khi; k++) {
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + k*dx[2];

        /*
        if (lsf_forcing)
        {
            // apply wls forcing
            wbar[k] = w_int_lsf[itime_curr][k] * coeff_curr + w_int_lsf[itime_next][k] * coeff_next;
            //amrex::Print() << " wsub forcing: k = " << k << ": wls_curr = " << w_int_lsf[itime_curr][k] << " wls_next = " << w_int_lsf[itime_next][k] << " src = " << wbar[k] << " src_scaled = " << wbar[k] * inv_scale << std::endl;
            wbar[k] *= inv_scale;
        } else {
            wbar[k] = 0.0;
        }
        */
        wbar[k] = 0.0;
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_geostrophic_profile (const Real& time,
                                     Vector<Real>& u_geos,
                                     Gpu::DeviceVector<Real>& d_u_geos,
                                     Vector<Real>& v_geos,
                                     Gpu::DeviceVector<Real>& d_v_geos,
                                     const Geometry& geom,
                                     std::unique_ptr<MultiFab>& z_phys_cc)
{
    if (u_geos.empty()) return;

    const int khi       = geom.Domain().bigEnd()[2];
    const Real* prob_lo = geom.ProbLo();
    const auto dx       = geom.CellSize();

    // Note: If z_phys_cc, then use_terrain=1 was set. If the z coordinate
    // varies in time and or space, then the the height needs to be
    // calculated at each time step. Here, we assume that only grid
    // stretching exists.
    if (z_phys_cc && zlevels.empty()) {
        Print() << "Initializing z levels on stretched grid" << std::endl;
        zlevels.resize(khi+1);
        reduce_to_max_per_height(zlevels, z_phys_cc);
    }

    // const Real coriolis = 2.0 * 2.0 * PI / 86400.0; // 0.376E-4;
    int itime_curr = 0;
    int itime_next = 0;
    amrex::Real coeff_curr = 1.0;
    amrex::Real coeff_next = 0.0;
    const bool lsf_forcing = times_lsf.size() > 0;
    amrex::Real inv_scale = 1.0 / parms.lsf_scale;
    if (lsf_forcing)
    {
       get_forcing_time_coeffs(time, itime_curr, itime_next, coeff_curr, coeff_next);
    }

    // Only apply temperature source below nominal inversion height
    for (int k = 0; k <= khi; k++) {
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)* dx[2];
        //const Real u_geo_wind = -10.0 + z_cc * 0.0018;

        //u_geos[k] =  u_geo_wind; // 0; // -coriolis_factor * v_geo_wind
        //v_geos[k] =  0 ; // coriolis *  u_geo_wind;

        /*
        if (lsf_forcing)
        {
            // apply uls and vls forcing
            u_geos[k] = u_int_lsf[itime_curr][k] * coeff_curr + u_int_lsf[itime_next][k] * coeff_next;
            //amrex::Print() << " u forcing: k = " << k << ": uls_curr = " << u_int_lsf[itime_curr][k] << " uls_next = " << u_int_lsf[itime_next][k] << " src = " << u_geos[k] << " src_scaled = " << u_geos[k] * inv_scale << std::endl;
            //u_geos[k] *= inv_scale;

            v_geos[k] = v_int_lsf[itime_curr][k] * coeff_curr + v_int_lsf[itime_next][k] * coeff_next;
            //amrex::Print() << " v forcing: k = " << k << ": vls_curr = " << v_int_lsf[itime_curr][k] << " vls_next = " << v_int_lsf[itime_next][k] << " src = " << v_geos[k] << " src_scaled = " << v_geos[k] * inv_scale << std::endl;
            //v_geos[k] *= inv_scale;
        } else {
            u_geos[k] = 0.0;
            v_geos[k] = 0.0;
        }
        */
       u_geos[k] = 0.0;
       v_geos[k] = 0.0;
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, u_geos.begin(), u_geos.end(), d_u_geos.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, v_geos.begin(), v_geos.end(), d_v_geos.begin());
}

int get_columns_in_str(std::string &str)
{
    std::istringstream iss(str);
    std::string tmp;
    // Get the number of columns in the file
    int ncols = 0;
    while (iss >> tmp) {
        ncols+= 1;
    }
    return ncols;
}

void
Problem::read_forcing_file (const std::string lsf_file,
                            int& nt,
                            int& nz,
                            amrex::Vector<amrex::Real>& times,
                            amrex::Vector<amrex::Vector<amrex::Real>>& z_lsf,
                            amrex::Vector<amrex::Vector<amrex::Real>>& p_lsf,
                            amrex::Vector<amrex::Vector<amrex::Real>>& t_lsf,
                            amrex::Vector<amrex::Vector<amrex::Real>>& q_lsf,
                            amrex::Vector<amrex::Vector<amrex::Real>>& u_lsf,
                            amrex::Vector<amrex::Vector<amrex::Real>>& v_lsf,
                            amrex::Vector<amrex::Vector<amrex::Real>>& w_lsf)
{
    amrex::Print() << "  Opening LSF file '" << lsf_file << "'" << std::endl;
    std::ifstream ifs(lsf_file);
    if (!ifs.is_open())
    {
        amrex::Error("Error opening input forcing file " + lsf_file);
    }

    std::string line;

    // temporary vectors for storing the values at each time
    amrex::Vector<amrex::Real> z_in, p_in, t_in, q_in, u_in, v_in, w_in;

    //   skip first header line
    std::getline(ifs, line);
    // header should have 7 columns: z, p, tls, qls, uls, vls, wls
    AMREX_ALWAYS_ASSERT(get_columns_in_str(line) == 7);

    nt = 0;
    nz = 0;

    amrex::Real curr_time = 0.0;
    amrex::Real start_time = 0.0; // first time in file, in days, other inputs are relative to this time

    amrex::Real pres0;
    while(std::getline(ifs, line)) {
        std::istringstream iss(line);

        // The start of each day has 6 columns: day, levels, pres0, "day, levels, pres0"
        int ncol = get_columns_in_str(line);

        if (ncol == 6)
        {
            // this is a new set of inputs at a given time

            // have to read these as strings first to split on commas
            std::string s_time, s_lev, s_pres0;
            iss >> s_time >> s_lev >> s_pres0;
            amrex::Real time = std::stod(s_time);
            int nlev = std::stoi(s_lev);
            pres0 = std::stod(s_pres0);

            // convert time to seconds (relative to first time), and pres from mb to Pa
            if (nt == 0)
            {
                start_time = time;
            }
            time = (time - start_time) * 86400.0;

            if (times.size() > 0 && time != times.back())
            {
                // if we started a new time, store data for previous time and reset tmp arrays
                AMREX_ALWAYS_ASSERT(z_in.size() == nz + 1);
                z_lsf.push_back(z_in);
                p_lsf.push_back(p_in);
                t_lsf.push_back(t_in);
                q_lsf.push_back(q_in);
                u_lsf.push_back(u_in);
                v_lsf.push_back(v_in);
                w_lsf.push_back(w_in);
                
                z_in.clear();
                p_in.clear();
                t_in.clear();
                q_in.clear();
                u_in.clear();
                v_in.clear();
                w_in.clear();
            }

            times.push_back(time);
            nt++;

            pres0 *= 100.0;
            if (nz == 0)
            {
                nz = nlev;
            } else {
                // make sure this has the same number of inputs.
                // TODO: this probably doesn't have to be true since we are interpolating onto the grid?
                AMREX_ALWAYS_ASSERT(nlev == nz);
            }

            z_in.push_back(0.0); // assuming zbot = 0.0 for now, this will get interpolated later once we have grid info
            p_in.push_back(pres0);
            t_in.push_back(0.0); // assuming temperature change at surf is 0
            q_in.push_back(0.0);
            u_in.push_back(0.0);
            v_in.push_back(0.0);
            w_in.push_back(0.0);
        } 
        else if (ncol == 7)
        {
            // this is defining a level at the current time
            amrex::Real z, p, tls, qls, uls, vls, wls;
            iss >> z >> p >> tls >> qls >> uls >> vls >> wls;

            if (z_in.size() == 1)
	        {
                u_in[0] = uls;
                v_in[0] = vls;
                w_in[0] = wls;
                //t_in[0] = getThgivenPandT(tls, p*100.0, R_d/Cp_d);
                t_in[0] = tls*std::pow(pres0/(p*100.0), R_d/Cp_d);
                q_in[0] = qls;
            }

            z_in.push_back(z);
            p_in.push_back(p * 100.0);
            //t_in.push_back(getThgivenPandT(tls, p*100.0, R_d/Cp_d));
            t_in.push_back(tls*std::pow(pres0/(p*100.0), R_d/Cp_d));
            q_in.push_back(qls);
            u_in.push_back(uls);
            v_in.push_back(vls);
            w_in.push_back(wls);
        }
        else {
            amrex::Error("unexpected line in input file: " + line);
        }
    }

    ifs.close();

    // store the last set from the EOF
    AMREX_ALWAYS_ASSERT(z_in.size() == nz + 1);
    z_lsf.push_back(z_in);
    p_lsf.push_back(p_in);
    t_lsf.push_back(t_in);
    q_lsf.push_back(q_in);
    u_lsf.push_back(u_in);
    v_lsf.push_back(v_in);
    w_lsf.push_back(w_in);

    nt = times.size();
    amrex::Print() << "  Read " << nz << " levels at " << nt << " times" << std::endl;
    
    // debugging print
    for (int i = 0; i < times.size(); i++)
    {
        amrex::Print() << " LSF at time = " << times[i] << std::endl;
        for (int lev = 0; lev < nz + 1; lev++)
        {
            amrex::Print() << "   level " << lev << ": z = " << z_lsf[i][lev] << " p = " << p_lsf[i][lev] << " tls = " << t_lsf[i][lev] << " qls = " << q_lsf[i][lev] << " uls = " << u_lsf[i][lev] << " vls = " << v_lsf[i][lev] << " wls = " << w_lsf[i][lev] << std::endl;
        }
        break;
    }
}


void
Problem::interp_forcing(const amrex::GeometryData &geom,
                        const amrex::Vector<amrex::Real>& zlevels_stag,
                        const amrex::Vector<amrex::Real>& times,
                        amrex::Vector<amrex::Vector<amrex::Real>>& z_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& t_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& q_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& u_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& v_lsf,
                        const amrex::Vector<amrex::Vector<amrex::Real>>& w_lsf)
{
    const int klo = 0;
    const int khi = geom.Domain().bigEnd()[AMREX_SPACEDIM-1];
    const int Nz = geom.Domain().size()[AMREX_SPACEDIM-1];
    const amrex::Real dz = geom.CellSize()[AMREX_SPACEDIM-1];

    const bool use_terrain = (zlevels_stag.size() > 0);
    const amrex::Real zbot = (use_terrain) ? zlevels_stag[klo]   : geom.ProbLo(AMREX_SPACEDIM-1);
    const amrex::Real ztop = (use_terrain) ? zlevels_stag[khi+1] : geom.ProbHi(AMREX_SPACEDIM-1);

    int num_times = times.size();
    AMREX_ALWAYS_ASSERT(num_times == nt_lsf);

    z_int_lsf.resize(num_times);
    t_int_lsf.resize(num_times);
    q_int_lsf.resize(num_times);
    u_int_lsf.resize(num_times);
    v_int_lsf.resize(num_times);
    w_int_lsf.resize(num_times);
    for (int itime = 0; itime < times.size(); itime++)
    {
        z_int_lsf[itime].resize(Nz+2);
        t_int_lsf[itime].resize(Nz+2);
        q_int_lsf[itime].resize(Nz+2);
        u_int_lsf[itime].resize(Nz+2);
        v_int_lsf[itime].resize(Nz+2);
        w_int_lsf[itime].resize(Nz+2);

        z_lsf[itime][0] = zbot; // update with zbot now that we have grid info

        // surface level
        z_int_lsf[itime][0] = z_lsf[itime][0];
        t_int_lsf[itime][0] = t_lsf[itime][0];
        q_int_lsf[itime][0] = q_lsf[itime][0];
        u_int_lsf[itime][0] = u_lsf[itime][0];
        v_int_lsf[itime][0] = v_lsf[itime][0];
        w_int_lsf[itime][0] = w_lsf[itime][0];

        for (int k = 0; k < Nz; k++)
        {
            z_int_lsf[itime][k+1] = (use_terrain) ? 0.5 * (zlevels_stag[k] + zlevels_stag[k+1])
                                                 : zbot + (k + 0.5) * dz;
            t_int_lsf[itime][k+1] = interpolate_1d(z_lsf[itime].dataPtr(), t_lsf[itime].dataPtr(), z_int_lsf[itime][k+1], nz_lsf);
            q_int_lsf[itime][k+1] = interpolate_1d(z_lsf[itime].dataPtr(), q_lsf[itime].dataPtr(), z_int_lsf[itime][k+1], nz_lsf);
            u_int_lsf[itime][k+1] = interpolate_1d(z_lsf[itime].dataPtr(), u_lsf[itime].dataPtr(), z_int_lsf[itime][k+1], nz_lsf);
            v_int_lsf[itime][k+1] = interpolate_1d(z_lsf[itime].dataPtr(), v_lsf[itime].dataPtr(), z_int_lsf[itime][k+1], nz_lsf);
            w_int_lsf[itime][k+1] = interpolate_1d(z_lsf[itime].dataPtr(), w_lsf[itime].dataPtr(), z_int_lsf[itime][k+1], nz_lsf);
        }

        z_int_lsf[itime][Nz+1] = ztop;
        t_int_lsf[itime][Nz+1] = interpolate_1d(z_lsf[itime].dataPtr(), t_lsf[itime].dataPtr(), ztop, nz_lsf);
        q_int_lsf[itime][Nz+1] = interpolate_1d(z_lsf[itime].dataPtr(), q_lsf[itime].dataPtr(), ztop, nz_lsf);
        u_int_lsf[itime][Nz+1] = interpolate_1d(z_lsf[itime].dataPtr(), u_lsf[itime].dataPtr(), ztop, nz_lsf);
        v_int_lsf[itime][Nz+1] = interpolate_1d(z_lsf[itime].dataPtr(), v_lsf[itime].dataPtr(), ztop, nz_lsf);
        w_int_lsf[itime][Nz+1] = interpolate_1d(z_lsf[itime].dataPtr(), w_lsf[itime].dataPtr(), ztop, nz_lsf);
    }

    // debugging print
    for (int i = 0; i < times.size(); i++)
    {
        amrex::Print() << " INTERPOLATED LSF at time = " << times[i] << std::endl;
        for (int k = 0; k <= Nz + 1; k++)
        {
            amrex::Print() << "   level " << k << ": z = " << z_int_lsf[i][k] << " tls = " << t_int_lsf[i][k] << " qls = " << q_int_lsf[i][k] << " uls = " << u_int_lsf[i][k] << " vls = " << v_int_lsf[i][k] << " wls = " << w_int_lsf[i][k] << std::endl;
        }
        break;
    }
}

void
Problem::get_forcing_time_coeffs(const amrex::Real& time,
                                 int& curr,
                                 int& next,
                                 amrex::Real& coeff_curr,
                                 amrex::Real& coeff_next)
{
    // helper function to get the current time indices of forcing data to use and
    // coefficients based on the current simulation time
    curr = 0;
    next = 0;
    coeff_curr = 1.0;
    coeff_next = 0.0;

    // find the time index of forcing data to use
    for (int nt = 1; nt < nt_lsf; nt++)
    {
        if (time > times_lsf[nt])
        {
            curr = nt;
        }
    }

    if (curr == nt_lsf - 1)
    {
        // use last set if time > last forcing time
        next = nt_lsf;
    } else {
        next = curr + 1;
        coeff_next = (time - times_lsf[curr]) / (times_lsf[next] - times_lsf[curr]);
        coeff_curr = (1.0 - coeff_next);
    }
}
