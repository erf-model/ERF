#include "prob.H"
#include "AMReX_Random.H"
#include <Utils/ParFunctions.H>

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem (const amrex::Real* problo, const amrex::Real* probhi)
{
    // Parse params
    ParmParse pp("prob");
    pp.query("rho_0", parms.rho_0);
    pp.query("T_0", parms.T_0);
    pp.query("A_0", parms.A_0);
    pp.query("QKE_0", parms.QKE_0);

    pp.query("U_0", parms.U_0);
    pp.query("V_0", parms.V_0);
    pp.query("W_0", parms.W_0);
    pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
    pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
    pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
    pp.query("T_0_Pert_Mag", parms.T_0_Pert_Mag);

    pp.query("pert_deltaU", parms.pert_deltaU);
    pp.query("pert_deltaV", parms.pert_deltaV);
    pp.query("pert_periods_U", parms.pert_periods_U);
    pp.query("pert_periods_V", parms.pert_periods_V);
    pp.query("pert_ref_height", parms.pert_ref_height);
    parms.aval = parms.pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
    parms.bval = parms.pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
    parms.ufac = parms.pert_deltaU * std::exp(0.5) / parms.pert_ref_height;
    parms.vfac = parms.pert_deltaV * std::exp(0.5) / parms.pert_ref_height;

    pp.query("dampcoef", parms.dampcoef);
    pp.query("zdamp", parms.zdamp);

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
    //===========================================================================

    init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert (
    const amrex::Box&  bx,
    const amrex::Box& xbx,
    const amrex::Box& ybx,
    const amrex::Box& zbx,
    amrex::Array4<amrex::Real const> const& /*state*/,
    amrex::Array4<amrex::Real      > const& state_pert,
    amrex::Array4<amrex::Real      > const& x_vel_pert,
    amrex::Array4<amrex::Real      > const& y_vel_pert,
    amrex::Array4<amrex::Real      > const& z_vel_pert,
    amrex::Array4<amrex::Real      > const& /*r_hse*/,
    amrex::Array4<amrex::Real      > const& /*p_hse*/,
    amrex::Array4<amrex::Real const> const& /*z_nd*/,
    amrex::Array4<amrex::Real const> const& /*z_cc*/,
    amrex::GeometryData const& geomdata,
    amrex::Array4<amrex::Real const> const& /*mf_m*/,
    amrex::Array4<amrex::Real const> const& /*mf_u*/,
    amrex::Array4<amrex::Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    ParallelForRNG(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
        const Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
        const Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z <= parms.pert_ref_height) && (parms.T_0_Pert_Mag != 0.0)) {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*parms.T_0_Pert_Mag;
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = parms.A_0 * exp(-10.*r*r);

        // Set an initial value for QKE
        state_pert(i, j, k, RhoQKE_comp) = parms.QKE_0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    // Set the x-velocity
    ParallelForRNG(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = parms.U_0;
        if ((z <= parms.pert_ref_height) && (parms.U_0_Pert_Mag != 0.0))
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real x_vel_prime = (rand_double*2.0 - 1.0)*parms.U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
        if (parms.pert_deltaU != 0.0)
        {
            const amrex::Real yl = y - prob_lo[1];
            const amrex::Real zl = z / parms.pert_ref_height;
            const amrex::Real damp = std::exp(-0.5 * zl * zl);
            x_vel_pert(i, j, k) += parms.ufac * damp * z * std::cos(parms.aval * yl);
        }
    });

  // Set the y-velocity
  ParallelForRNG(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the y-velocity
      y_vel_pert(i, j, k) = parms.V_0;
      if ((z <= parms.pert_ref_height) && (parms.V_0_Pert_Mag != 0.0))
      {
          Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
          Real y_vel_prime = (rand_double*2.0 - 1.0)*parms.V_0_Pert_Mag;
          y_vel_pert(i, j, k) += y_vel_prime;
      }
      if (parms.pert_deltaV != 0.0)
      {
          const amrex::Real xl = x - prob_lo[0];
          const amrex::Real zl = z / parms.pert_ref_height;
          const amrex::Real damp = std::exp(-0.5 * zl * zl);
          y_vel_pert(i, j, k) += parms.vfac * damp * z * std::cos(parms.bval * xl);
      }
  });

  // Set the z-velocity
  ParallelForRNG(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
  {
      const int dom_lo_z = geomdata.Domain().smallEnd()[2];
      const int dom_hi_z = geomdata.Domain().bigEnd()[2];

      // Set the z-velocity
      if (k == dom_lo_z || k == dom_hi_z+1)
      {
          z_vel_pert(i, j, k) = 0.0;
      }
      else if (parms.W_0_Pert_Mag != 0.0)
      {
          Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
          Real z_vel_prime = (rand_double*2.0 - 1.0)*parms.W_0_Pert_Mag;
          z_vel_pert(i, j, k) = parms.W_0 + z_vel_prime;
      }
  });
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_rhotheta_sources (const amrex::Real& /*time*/,
                                  amrex::Vector<amrex::Real>& src,
                                  amrex::Gpu::DeviceVector<amrex::Real>& d_src,
                                  const amrex::Geometry& geom,
                                  std::unique_ptr<amrex::MultiFab>& z_phys_cc)
{
    if (src.empty()) return;

    const int khi              = geom.Domain().bigEnd()[2];
    const amrex::Real* prob_lo = geom.ProbLo();
    const auto dx              = geom.CellSize();

    // Note: If z_phys_cc, then use_terrain=1 was set. If the z coordinate
    // varies in time and or space, then the the height needs to be
    // calculated at each time step. Here, we assume that only grid
    // stretching exists.
    if (z_phys_cc && zlevels.empty()) {
        amrex::Print() << "Initializing z levels on stretched grid" << std::endl;
        zlevels.resize(khi+1);
        reduce_to_max_per_level(zlevels, z_phys_cc);
    }

    // Only apply temperature source below nominal inversion height
    for (int k = 0; k <= khi; k++) {
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)* dx[2];
        if (z_cc < parms.cutoff) {
            src[k] = parms.advection_heating_rate;
        } else if (z_cc < parms.cutoff+parms.cutoff_transition) {
            src[k] = parms.advection_heating_rate * (z_cc-parms.cutoff)/parms.cutoff_transition;
        } else {
            src[k] = 0.0;
        }
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, src.begin(), src.end(), d_src.begin());
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_rhoqt_sources (const amrex::Real& /*time*/,
                               amrex::Vector<amrex::Real>& qsrc,
                               amrex::Gpu::DeviceVector<amrex::Real>& d_qsrc,
                               const amrex::Geometry& geom,
                               std::unique_ptr<amrex::MultiFab>& z_phys_cc)
{
    if (qsrc.empty()) return;

    const int khi              = geom.Domain().bigEnd()[2];
    const amrex::Real* prob_lo = geom.ProbLo();
    const auto dx              = geom.CellSize();

    // Note: If z_phys_cc, then use_terrain=1 was set. If the z coordinate
    // varies in time and or space, then the the height needs to be
    // calculated at each time step. Here, we assume that only grid
    // stretching exists.
    if (z_phys_cc && zlevels.empty()) {
        amrex::Print() << "Initializing z levels on stretched grid" << std::endl;
        zlevels.resize(khi+1);
        reduce_to_max_per_level(zlevels, z_phys_cc);
    }

    // Only apply temperature source below nominal inversion height
    for (int k = 0; k <= khi; k++) {
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)* dx[2];
        if (z_cc < parms.moisture_cutoff) {
            qsrc[k] = parms.advection_moisture_rate;
        } else if (z_cc < parms.moisture_cutoff+parms.moisture_cutoff_transition) {
            qsrc[k] = parms.advection_moisture_rate * (z_cc-parms.moisture_cutoff)/parms.moisture_cutoff_transition;
        } else {
            qsrc[k] = 0.0;
        }
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, qsrc.begin(), qsrc.end(), d_qsrc.begin());
}

//=============================================================================
// USER-DEFINED FUNCTION
//=============================================================================
void
Problem::update_w_subsidence (const amrex::Real& /*time*/,
                               amrex::Vector<amrex::Real>& wbar,
                               amrex::Gpu::DeviceVector<amrex::Real>& d_wbar,
                               const amrex::Geometry& geom,
                               std::unique_ptr<amrex::MultiFab>& z_phys_cc)
{
    if (wbar.empty()) return;

    const int khi              = geom.Domain().bigEnd()[2];
    const amrex::Real* prob_lo = geom.ProbLo();
    const auto dx              = geom.CellSize();

    // Note: If z_phys_cc, then use_terrain=1 was set. If the z coordinate
    // varies in time and or space, then the the height needs to be
    // calculated at each time step. Here, we assume that only grid
    // stretching exists.
    if (z_phys_cc && zlevels.empty()) {
        amrex::Print() << "Initializing z levels on stretched grid" << std::endl;
        zlevels.resize(khi+1);
        reduce_to_max_per_level(zlevels, z_phys_cc);
    }

    // Linearly increase wbar to the cutoff_max and then linearly decrease to cutoff_min
    Real z_0    = (z_phys_cc) ? zlevels[0] : prob_lo[2] + 0.5 * dx[2];
    Real slope1 =  parms.wbar_sub_max / (parms.wbar_cutoff_max - z_0);
    Real slope2 = -parms.wbar_sub_max / (parms.wbar_cutoff_min - parms.wbar_cutoff_max);
    wbar[0]     = 0.0;
    for (int k = 1; k <= khi; k++) {
        const Real z_cc = (z_phys_cc) ? zlevels[k] : prob_lo[2] + (k+0.5)* dx[2];
        if (z_cc <= parms.wbar_cutoff_max) {
            wbar[k] = slope1 * (z_cc - z_0);
        } else if (z_cc <= parms.wbar_cutoff_min) {
            wbar[k] = slope2 * (z_cc - parms.wbar_cutoff_max) + parms.wbar_sub_max;
        } else {
            wbar[k] = 0.0;
        }
    }

    // Copy from host version to device version
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, wbar.begin(), wbar.end(), d_wbar.begin());
}
