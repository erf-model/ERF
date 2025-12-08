#include "ERF_Prob.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* /*problo*/,
                const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem ()
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0"  , parms.T_0);
  pp.query("KE_0" , parms.KE_0);

  pp.query("KE_decay_height", parms.KE_decay_height);
  pp.query("KE_decay_order" , parms.KE_decay_order);

  init_base_parms(parms.rho_0, parms.T_0);
}


void
Problem::init_custom_pert (const amrex::Box&  bx,
                           const amrex::Box& /*xbx*/,
                           const amrex::Box& /*ybx*/,
                           const amrex::Box& /*zbx*/,
                           amrex::Array4<amrex::Real const> const& state,
                           amrex::Array4<amrex::Real      > const& state_pert,
                           amrex::Array4<amrex::Real      > const& /*x_vel_pert*/,
                           amrex::Array4<amrex::Real      > const& /*y_vel_pert*/,
                           amrex::Array4<amrex::Real      > const& /*z_vel_pert*/,
                           amrex::Array4<amrex::Real      > const& /*r_hse*/,
                           amrex::Array4<amrex::Real      > const& /*p_hse*/,
                           amrex::Array4<amrex::Real const> const& /*z_nd*/,
                           amrex::Array4<amrex::Real const> const& z_cc,
                           amrex::GeometryData const& geomdata,
                           amrex::Array4<amrex::Real const> const& /*mf_m*/,
                           amrex::Array4<amrex::Real const> const& /*mf_u*/,
                           amrex::Array4<amrex::Real const> const& /*mf_v*/,
                           const SolverChoice& /*sc*/,
                           const int /*lev*/)
{
    // Initialize KE for PBL model
    ParallelFor(bx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geom info
        const Real* prob_lo = geomdata.ProbLo();
        const Real* dx = geomdata.CellSize();
        const Real z   = (z_cc) ? z_cc(i,j,k) : prob_lo[2] + (k + 0.5) * dx[2];
        const Real dz  = (z_cc) ? z - z_cc(i,j,0) : z - prob_lo[2];

        // Base KE value
        state_pert(i,j,k,RhoKE_comp) = state(i,j,k,Rho_comp) * parms_d.KE_0;

        // Vertical profile
        if (parms_d.KE_decay_height > 0) {
            state_pert(i, j, k, RhoKE_comp) *= max(
                std::pow(1 - min(dz/parms_d.KE_decay_height,1.0), parms_d.KE_decay_order),
                1e-12);
        }
    });
}
