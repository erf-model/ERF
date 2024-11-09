#include "ERF_prob.H"
#include "AMReX_Random.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);
  pp.query("A_0", parms.A_0);
  pp.query("KE_0", parms.KE_0);
  pp.query("KE_decay_height", parms.KE_decay_height);
  pp.query("KE_decay_order", parms.KE_decay_order);

  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
  pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
  pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
  pp.query("T_0_Pert_Mag", parms.T_0_Pert_Mag);
  pp.query("pert_rhotheta", parms.pert_rhotheta);

  pp.query("pert_deltaU", parms.pert_deltaU);
  pp.query("pert_deltaV", parms.pert_deltaV);
  pp.query("pert_periods_U", parms.pert_periods_U);
  pp.query("pert_periods_V", parms.pert_periods_V);
  pp.query("pert_ref_height", parms.pert_ref_height);
  parms.aval = parms.pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
  parms.bval = parms.pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
  parms.ufac = parms.pert_deltaU * std::exp(0.5) / parms.pert_ref_height;
  parms.vfac = parms.pert_deltaV * std::exp(0.5) / parms.pert_ref_height;

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const amrex::Box&  bx,
    const amrex::Box& xbx,
    const amrex::Box& ybx,
    const amrex::Box& zbx,
    amrex::Array4<amrex::Real const> const& /*state*/,
    amrex::Array4<amrex::Real      > const& state_pert,
    amrex::Array4<amrex::Real      > const& x_vel_pert,
    amrex::Array4<amrex::Real      > const& y_vel_pert,
    amrex::Array4<amrex::Real      > const& z_vel_pert,
    amrex::Array4<amrex::Real      > const& r_hse,
    amrex::Array4<amrex::Real      > const& /*p_hse*/,
    amrex::Array4<amrex::Real const> const& z_nd,
    amrex::Array4<amrex::Real const> const& z_cc,
    amrex::GeometryData const& geomdata,
    amrex::Array4<amrex::Real const> const& /*mf_m*/,
    amrex::Array4<amrex::Real const> const& /*mf_u*/,
    amrex::Array4<amrex::Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const bool use_terrain = sc.use_terrain;

    if (parms.KE_decay_height > 0) {
        amrex::Print() << "Initial KE profile (order " << parms.KE_decay_order
                       << ") will extend up to " << parms.KE_decay_height
                       << std::endl;
    }

    if (parms.pert_ref_height > 0) {
        if ((parms.pert_deltaU != 0.0) || (parms.pert_deltaV != 0.0)) {
            amrex::Print() << "Adding divergence-free perturbations "
                           << parms.pert_deltaU << " " << parms.pert_deltaV
                           << std::endl;
        }
        if (parms.U_0_Pert_Mag != 0.0) {
            amrex::Print() << "Adding random x-velocity perturbations" << std::endl;
        }
        if (parms.V_0_Pert_Mag != 0.0) {
            amrex::Print() << "Adding random y-velocity perturbations" << std::endl;
        }
        if (parms.T_0_Pert_Mag != 0.0) {
            if (parms.pert_rhotheta) {
                amrex::Print() << "Adding random rho*theta perturbations" << std::endl;
            } else {
                amrex::Print() << "Adding random theta perturbations" << std::endl;
            }
        }
    }

  ParallelForRNG(bx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = use_terrain ? z_cc(i,j,k) : prob_lo[2] + (k + 0.5) * dx[2];

    // Define a point (xc,yc,zc) at the center of the domain
    const Real xc = 0.5 * (prob_lo[0] + prob_hi[0]);
    const Real yc = 0.5 * (prob_lo[1] + prob_hi[1]);
    const Real zc = 0.5 * (prob_lo[2] + prob_hi[2]);

    const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

    // Add temperature perturbations
    if ((z <= parms_d.pert_ref_height) && (parms_d.T_0_Pert_Mag != 0.0)) {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*parms_d.T_0_Pert_Mag;
        if (!parms_d.pert_rhotheta) {
            // we're perturbing theta, not rho*theta
            state_pert(i, j, k, RhoTheta_comp) *= r_hse(i,j,k);
        }
    }

    // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
    state_pert(i, j, k, RhoScalar_comp) = parms_d.A_0 * exp(-10.*r*r);

    // Set an initial value for SGS KE
    if (state_pert.nComp() > RhoKE_comp) {
        // Deardorff
        state_pert(i, j, k, RhoKE_comp) = r_hse(i,j,k) * parms_d.KE_0;
        if (parms_d.KE_decay_height > 0) {
            // scale initial SGS kinetic energy with height
            state_pert(i, j, k, RhoKE_comp) *= max(
                std::pow(1 - min(z/parms_d.KE_decay_height,1.0), parms_d.KE_decay_order),
                1e-12);
        }
    }

    if (use_moisture) {
        state_pert(i, j, k, RhoQ1_comp) = 0.0;
        state_pert(i, j, k, RhoQ2_comp) = 0.0;
    }
  });

  // Set the x-velocity
  ParallelForRNG(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = use_terrain ? 0.25*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                                      + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) )
                               : prob_lo[2] + (k + 0.5) * dx[2];

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
  ParallelForRNG(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = use_terrain ? 0.25*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                                      + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) )
                               : prob_lo[2] + (k + 0.5) * dx[2];

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
  ParallelForRNG(zbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
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
