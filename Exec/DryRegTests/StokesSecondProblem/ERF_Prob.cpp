#include "ERF_Prob.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  ParmParse pp("prob");
  Real rho_0 =   1.2; pp.query("rho_0", rho_0);
  Real   T_0 = 300.0; pp.query("T_0", T_0);

  init_base_parms(rho_0, T_0);
}

void
Problem::init_custom_pert (
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // This version perturbs rho but not p -- TODO: CHECK THIS
        state_pert(i, j, k, RhoTheta_comp) = std::pow(1.0,1.0/Gamma) * 101325.0 / 287.0 - p_hse(i,j,k);
    });

    amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_pert_vels (
    const Box& /*xbx*/,
    const Box& /*ybx*/,
    const Box& /*zbx*/,
    Array4<Real      > const& /*x_vel_pert*/,
    Array4<Real      > const& /*y_vel_pert*/,
    Array4<Real      > const& /*z_vel_pert*/,
    Array4<Real const> const& /*z_nd*/,
    GeometryData const& /*geomdata*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& /*sc*/,
    const int /*lev*/)
{ }

#ifdef ERF_USE_TERRAIN_VELOCITY
Real
Problem::compute_terrain_velocity (const Real time)
{
    Real U = 10.0;
    Real omega = 2.0*M_PI*1000.0;
    return U*cos(omega*time);
}
#endif
