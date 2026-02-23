#include <AMReX_ParmParse.H>

#include "ERF_Prob.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* problo,
                const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem (const Real* problo,
                  const Real* probhi)
{
    ParmParse pp("prob");
    Real rho_0 = 1.2; pp.query("rho_0", rho_0);
    Real T_0 = 300.0; pp.query("T_0", T_0);

    init_base_parms(rho_0, T_0);
}

void
Problem::init_custom_pert (
    const Box& /*bx*/,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& /*state_pert*/,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
}

void
Problem::init_custom_pert_vels (
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real const> const& z_nd,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParmParse pp("prob");
    Real U_0 = 0.0; pp.query("U_0", U_0);
    Real V_0 = 0.0; pp.query("V_0", V_0);
    Real U_0_Pert_Mag = 0.0; pp.query("U_0_Pert_Mag", U_0_Pert_Mag);
    Real V_0_Pert_Mag = 0.0; pp.query("V_0_Pert_Mag", V_0_Pert_Mag);
    Real W_0_Pert_Mag = 0.0; pp.query("W_0_Pert_Mag", W_0_Pert_Mag);
    Real pert_deltaU = 0.0; pp.query("pert_deltaU", pert_deltaU);
    Real pert_deltaV = 0.0; pp.query("pert_deltaV", pert_deltaV);
    Real pert_periods_U  = 5.0; pp.query("pert_periods_U", pert_periods_U);
    Real pert_periods_V  = 5.0; pp.query("pert_periods_V", pert_periods_V);;
    Real pert_ref_height = 1.0; pp.query("pert_ref_height", pert_ref_height);

    pp.query("pert_ref_height", pert_ref_height);

    // Real aval = pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
    // Real ufac = pert_deltaU * std::exp(0.5) / pert_ref_height;

    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real bval = pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
    Real vfac = pert_deltaV * std::exp(0.5) / pert_ref_height;

    const bool use_terrain  = (SolverChoice::terrain_type != TerrainType::None);

    // Set the x-velocity
    ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real z = use_terrain ? 0.25*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                                          + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) )
            : prob_lo[2] + (k + 0.5) * dx[2];

        x_vel_pert(i, j, k) = 0.0;
        if ((z <= pert_ref_height) && (U_0_Pert_Mag != 0.0))
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real x_vel_prime = (rand_double*2.0 - 1.0)*U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
    });

    // Set the y-velocity
    ParallelForRNG(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = use_terrain ? 0.25*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                                          + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) )
            : prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = 0.0;
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

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = 1. / dx[0];
    dxInv[1] = 1. / dx[1];
    dxInv[2] = 1. / dx[2];

    // Set the z-velocity from impenetrable condition
    ParallelForRNG(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
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
            z_vel_pert(i, j, k) = z_vel_prime;
        }
    });

    amrex::Gpu::streamSynchronize();
}
