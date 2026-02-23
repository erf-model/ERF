#include "ERF_Prob.H"
#include "AMReX_Random.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const Real* problo, const Real* probhi)
{
  ParmParse pp("prob");

  Real rho_0 = 1.0; pp.query("rho_0", rho_0);
  Real T_0   = 300.0; pp.query("T_0", T_0);

}

void
Problem::init_custom_pert (
    const amrex::Box&  bx,
    amrex::Array4<Real const> const& /*state*/,
    amrex::Array4<Real      > const& state_pert,
    amrex::Array4<Real      > const& r_hse,
    amrex::Array4<Real      > const& /*p_hse*/,
    amrex::Array4<Real const> const& /*z_nd*/,
    amrex::Array4<Real const> const& z_cc,
    amrex::GeometryData const& geomdata,
    amrex::Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp("prob");

    Real rho_0 = 1.0; pp.query("rho_0", rho_0);
    Real T_0   = 300.0; pp.query("T_0", T_0);

    Real A_0   = 1.0; pp.query("A_0", A_0);
    Real KE_0  = 0.1; pp.query("KE_0", KE_0);

    Real KE_decay_height = -1; pp.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order = 1; pp.query("KE_decay_order", KE_decay_order);

    Real U_0 = 0.0; pp.query("U_0", U_0);
    Real V_0 = 0.0; pp.query("V_0", V_0);
    Real W_0 = 0.0; pp.query("W_0", W_0);

    // random initial perturbations (legacy code)
    Real U_0_Pert_Mag = 0.0; pp.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = 0.0; pp.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = 0.0; pp.query("W_0_Pert_Mag", W_0_Pert_Mag);
    Real T_0_Pert_Mag = 0.0; pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    bool pert_rhotheta = true; pp.query("pert_rhotheta", pert_rhotheta);

    // divergence-free initial perturbations
    Real pert_deltaU = 0.0; pp.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = 0.0; pp.query("pert_deltaV", pert_deltaV);
    Real pert_periods_U = 5.0; pp.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V = 5.0; pp.query("pert_periods_V", pert_periods_V);
    Real pert_ref_height = 100.0; pp.query("pert_ref_height", pert_ref_height);

    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();

    // helper vars
    Real aval = pert_periods_U * 2.0 * PI / (prob_hi[1] - prob_lo[1]);
    Real bval = pert_periods_V * 2.0 * PI / (prob_hi[0] - prob_lo[0]);
    Real ufac = pert_deltaU * std::exp(0.5) / pert_ref_height;
    Real vfac = pert_deltaV * std::exp(0.5) / pert_ref_height;

    const bool use_terrain  = (SolverChoice::terrain_type != TerrainType::None);

    if (KE_decay_height > 0) {
        amrex::Print() << "Initial KE profile (order " << KE_decay_order
                       << ") will extend up to " << KE_decay_height
                       << std::endl;
    }

    if (pert_ref_height > 0) {
        if ((pert_deltaU != 0.0) || (pert_deltaV != 0.0)) {
            amrex::Print() << "Adding divergence-free perturbations "
                           << pert_deltaU << " " << pert_deltaV
                           << std::endl;
        }
        if (U_0_Pert_Mag != 0.0) {
            amrex::Print() << "Adding random x-velocity perturbations" << std::endl;
        }
        if (V_0_Pert_Mag != 0.0) {
            amrex::Print() << "Adding random y-velocity perturbations" << std::endl;
        }
        if (T_0_Pert_Mag != 0.0) {
            if (pert_rhotheta) {
                amrex::Print() << "Adding random rho*theta perturbations" << std::endl;
            } else {
                amrex::Print() << "Adding random theta perturbations" << std::endl;
            }
        }
    }

    ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
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
        if ((z <= pert_ref_height) && (T_0_Pert_Mag != 0.0)) {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*2.0 - 1.0)*T_0_Pert_Mag;
            if (!pert_rhotheta) {
                // we're perturbing theta, not rho*theta
                state_pert(i, j, k, RhoTheta_comp) *= r_hse(i,j,k);
            }
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = A_0 * exp(-10.*r*r);

        // Set an initial value for SGS KE
        if (state_pert.nComp() > RhoKE_comp) {
            // Deardorff
            state_pert(i, j, k, RhoKE_comp) = r_hse(i,j,k) * KE_0;
            if (KE_decay_height > 0) {
                // scale initial SGS kinetic energy with height
                state_pert(i, j, k, RhoKE_comp) *= max(
                    std::pow(1 - min(z/KE_decay_height,1.0), KE_decay_order),
                    1e-12);
            }
        }
    });
}

void
Problem::init_custom_pert_vels (
    const amrex::Box& xbx,
    const amrex::Box& ybx,
    const amrex::Box& zbx,
    amrex::Array4<Real      > const& x_vel_pert,
    amrex::Array4<Real      > const& y_vel_pert,
    amrex::Array4<Real      > const& z_vel_pert,
    amrex::Array4<Real const> const& z_nd,
    amrex::GeometryData const& geomdata,
    amrex::Array4<Real const> const& /*mf_u*/,
    amrex::Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp("prob");

    Real rho_0 = 1.0; pp.query("rho_0", rho_0);
    Real T_0   = 300.0; pp.query("T_0", T_0);

    Real A_0   = 1.0; pp.query("A_0", A_0);
    Real KE_0  = 0.1; pp.query("KE_0", KE_0);

    Real KE_decay_height = -1; pp.query("KE_decay_height", KE_decay_height);
    Real KE_decay_order = 1; pp.query("KE_decay_order", KE_decay_order);

    Real U_0 = 0.0; pp.query("U_0", U_0);
    Real V_0 = 0.0; pp.query("V_0", V_0);
    Real W_0 = 0.0; pp.query("W_0", W_0);

    // random initial perturbations (legacy code)
    Real U_0_Pert_Mag = 0.0; pp.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = 0.0; pp.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = 0.0; pp.query("W_0_Pert_Mag", W_0_Pert_Mag);
    Real T_0_Pert_Mag = 0.0; pp.query("T_0_Pert_Mag", T_0_Pert_Mag);
    bool pert_rhotheta = true; pp.query("pert_rhotheta", pert_rhotheta);

    // divergence-free initial perturbations
    Real pert_deltaU = 0.0; pp.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = 0.0; pp.query("pert_deltaV", pert_deltaV);
    Real pert_periods_U = 5.0; pp.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V = 5.0; pp.query("pert_periods_V", pert_periods_V);
    Real pert_ref_height = 100.0; pp.query("pert_ref_height", pert_ref_height);

    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();

    // helper vars
    Real aval = pert_periods_U * 2.0 * PI / (prob_hi[1] - prob_lo[1]);
    Real bval = pert_periods_V * 2.0 * PI / (prob_hi[0] - prob_lo[0]);
    Real ufac = pert_deltaU * std::exp(0.5) / pert_ref_height;
    Real vfac = pert_deltaV * std::exp(0.5) / pert_ref_height;

    const bool use_terrain  = (SolverChoice::terrain_type != TerrainType::None);

    // Set the x-velocity
    ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = use_terrain ? 0.25*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                                          + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) )
                                   : prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = U_0;
        if ((z <= pert_ref_height) && (U_0_Pert_Mag != 0.0))
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real x_vel_prime = (rand_double*2.0 - 1.0)*U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
        if (pert_deltaU != 0.0)
        {
            const Real yl = y - prob_lo[1];
            const Real zl = z / pert_ref_height;
            const Real damp = std::exp(-0.5 * zl * zl);
            x_vel_pert(i, j, k) += ufac * damp * z * std::cos(aval * yl);
        }
    });

  // Set the y-velocity
  ParallelForRNG(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = use_terrain ? 0.25*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                                      + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) )
                               : prob_lo[2] + (k + 0.5) * dx[2];

    // Set the y-velocity
    y_vel_pert(i, j, k) = V_0;
    if ((z <= pert_ref_height) && (V_0_Pert_Mag != 0.0))
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real y_vel_prime = (rand_double*2.0 - 1.0)*V_0_Pert_Mag;
        y_vel_pert(i, j, k) += y_vel_prime;
    }
    if (pert_deltaV != 0.0)
    {
        const Real xl = x - prob_lo[0];
        const Real zl = z / pert_ref_height;
        const Real damp = std::exp(-0.5 * zl * zl);
        y_vel_pert(i, j, k) += vfac * damp * z * std::cos(bval * xl);
    }
  });

  // Set the z-velocity
  ParallelForRNG(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept {
    const int dom_lo_z = geomdata.Domain().smallEnd()[2];
    const int dom_hi_z = geomdata.Domain().bigEnd()[2];

    // Set the z-velocity
    if (k == dom_lo_z || k == dom_hi_z+1)
    {
        z_vel_pert(i, j, k) = 0.0;
    }
    else if (W_0_Pert_Mag != 0.0)
    {
        Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
        Real z_vel_prime = (rand_double*2.0 - 1.0)*W_0_Pert_Mag;
        z_vel_pert(i, j, k) = W_0 + z_vel_prime;
    }
  });
}
