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

Problem::Problem (const amrex::Real* problo,
                  const amrex::Real* probhi)
{
    // Parse params
    ParmParse pp("prob");
    pp.query("rho_0", parms.rho_0);
    pp.query("U_0", parms.U_0);
    pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
    pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
    pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
    pp.query("pert_ref_height", parms.pert_ref_height);
    parms.aval = parms.pert_periods_U * 2.0 * PI / (probhi[1] - problo[1]);
    parms.bval = parms.pert_periods_V * 2.0 * PI / (probhi[0] - problo[0]);
    parms.ufac = parms.pert_deltaU * std::exp(0.5) / parms.pert_ref_height;
    parms.vfac = parms.pert_deltaV * std::exp(0.5) / parms.pert_ref_height;

    init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert (
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& z_nd,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc,
    const int /*lev*/)
{
    const int khi = geomdata.Domain().bigEnd()[2];

    const bool use_moisture = (sc.moisture_type != MoistureType::None);
    const bool use_terrain  = (SolverChoice::terrain_type != TerrainType::None);

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    // Geometry (note we must include these here to get the data on device)
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    // Set the x-velocity
    ParallelForRNG(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real z = use_terrain ? 0.25*( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                                          + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) )
            : prob_lo[2] + (k + 0.5) * dx[2];

        x_vel_pert(i, j, k) = 0.0;
        if ((z <= parms_d.pert_ref_height) && (parms_d.U_0_Pert_Mag != 0.0))
        {
            Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
            Real x_vel_prime = (rand_double*2.0 - 1.0)*parms_d.U_0_Pert_Mag;
            x_vel_pert(i, j, k) += x_vel_prime;
        }
    });

    // Set the y-velocity
    ParallelForRNG(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = use_terrain ? 0.25*( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                                          + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) )
            : prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel_pert(i, j, k) = 0.0;
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

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = 1. / dx[0];
    dxInv[1] = 1. / dx[1];
    dxInv[2] = 1. / dx[2];

    // Set the z-velocity from impenetrable condition
    ParallelForRNG(zbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
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
            z_vel_pert(i, j, k) = z_vel_prime;
        }
    });

    amrex::Gpu::streamSynchronize();
}
