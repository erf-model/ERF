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
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);
  pp.query("A_0", parms.A_0);
  pp.query("KE_0", parms.KE_0);

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
  parms.aval = parms.pert_periods_U * Real(2.0) * PI / (probhi[1] - problo[1]);
  parms.bval = parms.pert_periods_V * Real(2.0) * PI / (probhi[0] - problo[0]);
  parms.ufac = parms.pert_deltaU * std::exp(Real(0.5)) / parms.pert_ref_height;
  parms.vfac = parms.pert_deltaV * std::exp(Real(0.5)) / parms.pert_ref_height;

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert (
    const Box&  bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParallelForRNG(bx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
        const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
        const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

        // Define a point (xc,yc,zc) at the center of the domain
        const Real xc = Real(0.5) * (prob_lo[0] + prob_hi[0]);
        const Real yc = Real(0.5) * (prob_lo[1] + prob_hi[1]);
        const Real zc = Real(0.5) * (prob_lo[2] + prob_hi[2]);

        const Real r  = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        // Add temperature perturbations
        if ((z <= parms_d.pert_ref_height) && (parms_d.T_0_Pert_Mag != 0.0)) {
            Real rand_double = Random(engine); // Between 0.0 and 1.0
            state_pert(i, j, k, RhoTheta_comp) = (rand_double*Real(2.0) - Real(1.0))*parms_d.T_0_Pert_Mag;
        }

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state_pert(i, j, k, RhoScalar_comp) = parms_d.A_0 * exp(Real(-10.)*r*r);

        // Set an initial value for KE
        state_pert(i, j, k, RhoKE_comp) = parms_d.KE_0;
    });
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
    // Set the x-velocity
    ParallelForRNG(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
        const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

        // Set the x-velocity
        x_vel_pert(i, j, k) = parms_d.U_0;
        if ((z <= parms_d.pert_ref_height) && (parms_d.U_0_Pert_Mag != 0.0))
        {
            Real rand_double = Random(engine); // Between 0.0 and 1.0
            Real x_vel_prime = (rand_double*Real(2.0) - Real(1.0))*parms_d.U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
        if (parms_d.pert_deltaU != 0.0)
        {
            const Real yl = y - prob_lo[1];
            const Real zl = z / parms_d.pert_ref_height;
            const Real damp = std::exp(Real(-0.5) * zl * zl);
            x_vel_pert(i, j, k) += parms_d.ufac * damp * z * std::cos(parms_d.aval * yl);
        }
    });

    // Set the y-velocity
    ParallelForRNG(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
        const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = parms_d.V_0;
        if ((z <= parms_d.pert_ref_height) && (parms_d.V_0_Pert_Mag != 0.0))
        {
            Real rand_double = Random(engine); // Between 0.0 and 1.0
            Real y_vel_prime = (rand_double*Real(2.0) - Real(1.0))*parms_d.V_0_Pert_Mag;
            y_vel_pert(i, j, k) += y_vel_prime;
        }
        if (parms_d.pert_deltaV != 0.0)
        {
            const Real xl = x - prob_lo[0];
            const Real zl = z / parms_d.pert_ref_height;
            const Real damp = std::exp(Real(-0.5) * zl * zl);
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
            Real rand_double = Random(engine); // Between 0.0 and 1.0
            Real z_vel_prime = (rand_double*Real(2.0) - Real(1.0))*parms_d.W_0_Pert_Mag;
            z_vel_pert(i, j, k) = parms_d.W_0 + z_vel_prime;
        }
    });
}
