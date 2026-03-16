void init_custom_pert (
    const amrex::Box&  bx,
    amrex::Array4<amrex::Real const> const& state,
    amrex::Array4<amrex::Real      > const& state_pert,
    amrex::Array4<amrex::Real      > const& r_hse,
    amrex::Array4<amrex::Real      > const& p_hse,
    amrex::Array4<amrex::Real const> const& z_nd,
    amrex::Array4<amrex::Real const> const& z_cc,
    amrex::GeometryData const& geomdata,
    amrex::Array4<amrex::Real const> const& mf_m,
    const SolverChoice& sc,
    const int lev) override;

void init_custom_pert_vels (
    const amrex::Box& xbx,
    const amrex::Box& ybx,
    const amrex::Box& zbx,
    amrex::Array4<amrex::Real      > const& x_vel_pert,
    amrex::Array4<amrex::Real      > const& y_vel_pert,
    amrex::Array4<amrex::Real      > const& z_vel_pert,
    amrex::Array4<amrex::Real const> const& z_nd,
    amrex::GeometryData const& geomdata,
    amrex::Array4<amrex::Real const> const& mf_u,
    amrex::Array4<amrex::Real const> const& mf_v,
    const SolverChoice& sc,
    const int lev) override;

AMREX_GPU_DEVICE
static
amrex::Real
erf_vortex_Gaussian(
  amrex::Real x,  amrex::Real y,
  amrex::Real xc, amrex::Real yc,
  amrex::Real R,  amrex::Real beta,
  amrex::Real sigma)
{
  // Evaluate Gaussian function
  const amrex::Real r2 = ((x-xc)*(x-xc) + (y-yc)*(y-yc)) / (R*R);
  return beta * std::exp(-r2/(2.*sigma*sigma));
}

#ifdef ERF_USE_TERRAIN_VELOCITY
amrex::Real
Problem::compute_terrain_velocity (const amrex::Real time)
{
    amrex::Real U = 10.0;
    amrex::Real omega = 2.0*M_PI*1000.0;
    return U*cos(omega*time);
}
#endif
