#include "prob.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
    // Parse params
    amrex::ParmParse pp("prob");
    pp.query("rho_0", parms.rho_0);
    pp.query("T_0", parms.T_0);
    pp.query("A_0", parms.A_0);
    pp.query("B_0", parms.B_0);
    pp.query("u_0", parms.u_0);
    pp.query("v_0", parms.v_0);
    pp.query("rad_0", parms.rad_0);
    pp.query("z0", parms.z0);
    pp.query("zRef", parms.zRef);
    pp.query("uRef", parms.uRef);

    pp.query("xc_frac", parms.xc_frac);
    pp.query("yc_frac", parms.yc_frac);
    pp.query("zc_frac", parms.zc_frac);

    pp.query("prob_type", parms.prob_type);

    init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::erf_init_rayleigh(
    amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
    amrex::Geometry      const& geom,
    std::unique_ptr<MultiFab>& /*z_phys_cc*/)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      rayleigh_ptrs[Rayleigh::tau][k]      = 1.0;
      rayleigh_ptrs[Rayleigh::ubar][k]     = 2.0;
      rayleigh_ptrs[Rayleigh::vbar][k]     = 1.0;
      rayleigh_ptrs[Rayleigh::wbar][k]     = 0.0;
      rayleigh_ptrs[Rayleigh::thetabar][k] = parms.T_0;
  }
}

void
Problem::init_custom_pert(
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
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Define a point (xc,yc,zc) at the center of the domain
    const Real xc = parms.xc_frac * (prob_lo[0] + prob_hi[0]);
    const Real yc = parms.yc_frac * (prob_lo[1] + prob_hi[1]);
    const Real zc = parms.zc_frac * (prob_lo[2] + prob_hi[2]);

    // Define ellipse parameters
    const Real r0   = parms.rad_0 * (prob_hi[0] - prob_lo[0]);

    const Real r3d    = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
    const Real r2d_xy = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
    const Real r2d_xz = std::sqrt((x-xc)*(x-xc) + (z-zc)*(z-zc));

    if (parms.prob_type == 10)
    {
        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain,
        //            + B_0*sin(x)
        state_pert(i, j, k, RhoScalar_comp) = parms.A_0 * exp(-10.*r3d*r3d) + parms.B_0*sin(x);

    } else if (parms.prob_type == 11) {
        state_pert(i, j, k, RhoScalar_comp) = parms.A_0 * 0.25 * (1.0 + std::cos(PI * std::min(r2d_xy, r0) / r0));
    } else if (parms.prob_type == 12) {
        state_pert(i, j, k, RhoScalar_comp) = parms.A_0 * 0.25 * (1.0 + std::cos(PI * std::min(r2d_xz, r0) / r0));
    } else if (parms.prob_type == 13) {
        const Real r0_z = parms.rad_0 * (prob_hi[2] - prob_lo[2]);
        //ellipse for mapfac shear validation
        const Real r2d_xz_ell = std::sqrt((x-xc)*(x-xc)/(r0*r0) + (z-zc)*(z-zc)/(r0_z*r0_z));
        state_pert(i, j, k, RhoScalar_comp) = parms.A_0 * 0.25 * (1.0 + std::cos(PI * std::min(r2d_xz_ell, r0_z) / r0_z));
    } else if (parms.prob_type == 14) {
        state_pert(i, j, k, RhoScalar_comp) = std::cos(PI*x);
    } else {
        // Set scalar = A_0 in a ball of radius r0 and 0 elsewhere
        if (r3d < r0) {
           state_pert(i, j, k, RhoScalar_comp) = parms.A_0;
        } else {
           state_pert(i, j, k, RhoScalar_comp) = 0.0;
        }
    }

    state_pert(i, j, k, RhoScalar_comp) *= parms.rho_0;

    if (use_moisture) {
        state_pert(i, j, k, RhoQ1_comp) = 0.0;
        state_pert(i, j, k, RhoQ2_comp) = 0.0;
    }
  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel_pert(i, j, k) = parms.u_0;

      const Real* prob_lo = geomdata.ProbLo();
      const Real*      dx = geomdata.CellSize();
      const Real        z = prob_lo[2] + (k + 0.5) * dx[2];

      // Set the x-velocity
      x_vel_pert(i, j, k) = parms.u_0 + parms.uRef *
                            std::log((z + parms.z0)/parms.z0)/
                            std::log((parms.zRef +parms.z0)/parms.z0);
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel_pert(i, j, k) = parms.v_0;
  });

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel_pert(i, j, k) = 0.0;
  });
}

    /**
    * Function to perform custom initialization of terrain
    *
    * This version takes a single FArrayBox instead of a MultiFab
    *
    */
    void
    Problem::init_custom_terrain (const amrex::Geometry& /*geom*/,
                                  amrex::FArrayBox& z_phys_nd,
                                  const amrex::Real& /*time*/)
    {
        // Note that this only sets the terrain value at the ground IF k=0 is in the box
        amrex::Print() << "Initializing flat terrain at z=0" << std::endl;

        // Bottom of domain
        int k0 = 0;

        // Grown box with no z range
        amrex::Box bx = z_phys_nd.box();
        amrex::Array4<amrex::Real> const& z_arr = z_phys_nd.array();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            z_arr(i,j,k0) = 2.5;
        });
    }

#if 0
AMREX_GPU_DEVICE
Real
dhdt(int /*i*/, int /*j*/,
     const GpuArray<Real,AMREX_SPACEDIM> /*dx*/,
     const Real /*time_mt*/, const Real /*delta_t*/)
{
    return 0.;
}
#endif
