#include "prob.H"
#include "EOS.H"
#include "TerrainMetrics.H"

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
  pp.query("Ampl", parms.Ampl);
  pp.query("T_0" , parms.T_0);

  init_base_parms(parms.rho_0, parms.T_0);
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
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& z_nd,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
  const bool use_moisture = (sc.moisture_type != MoistureType::None);

  const Real T_sfc    = parms.T_0;
  //const Real thetabar = T_sfc;

  Real H           = geomdata.ProbHi()[2];
  Real Ampl        = parms.Ampl;
  Real wavelength  = 100.;
  Real kp          = 2.0 * PI / wavelength;
  Real g           = CONST_GRAV;
  Real omega       = std::sqrt(g * kp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const auto prob_lo  = geomdata.ProbLo();
      const auto dx       = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.5) * dx[0];

      Real z = z_cc(i,j,k);
      Real z_base = Ampl * std::sin(kp * x);
      z -= z_base;

      Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
      Real p_prime = -(Ampl * omega * omega / kp) * fac * std::sin(kp * x) * r_hse(i,j,k);
      Real p_total = p_hse(i,j,k) + p_prime;

      // Define (rho theta) given pprime
      state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p_total) - getRhoThetagivenP(p_hse(i,j,k));

      // Set scalar = 0 everywhere
      state_pert(i, j, k, RhoScalar_comp) = state_pert(i,j,k,Rho_comp);

      if (use_moisture) {
          state_pert(i, j, k, RhoQ1_comp) = 0.0;
          state_pert(i, j, k, RhoQ2_comp) = 0.0;
      }
  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const auto prob_lo  = geomdata.ProbLo();
      const auto dx       = geomdata.CellSize();

      const Real x = prob_lo[0] + i * dx[0];
      Real z = 0.25 * (z_nd(i,j,k) + z_nd(i,j+1,k) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1));

      Real z_base = Ampl * std::sin(kp * x);
      z -= z_base;

      Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);

      x_vel_pert(i, j, k) = -Ampl * omega * fac * std::sin(kp * x);
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel_pert(i, j, k) = 0.0;
  });

  // Set the z-velocity from impenetrable condition
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const auto dx         = geomdata.CellSize();

      Real x   = (i + 0.5) * dx[0];
      Real z   = 0.25 * ( z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k) );
      Real z_base = Ampl * std::sin(kp * x);

      z -= z_base;
      Real fac = std::sinh( kp * (z - H) ) / std::sinh(kp * H);

      z_vel_pert(i, j, k) = Ampl * omega * fac * std::cos(kp * x);
  });

  amrex::Gpu::streamSynchronize();

}

void
Problem::erf_init_rayleigh (Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
                            Geometry const& geom,
                            std::unique_ptr<MultiFab>& /*z_phys_cc*/)
{
    const int khi = geom.Domain().bigEnd()[2];
    for (int k = 0; k <= khi; k++)
    {
        rayleigh_ptrs[Rayleigh::tau][k]      =   0.0;
        rayleigh_ptrs[Rayleigh::ubar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::vbar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::wbar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::thetabar][k] = 300.0;
    }

    // Damping above k = 60
    for (int k = 60; k <= khi; k++)
    {
        rayleigh_ptrs[Rayleigh::tau][k]      =  10.0;  // Remember that this gets multiplied by dt
        rayleigh_ptrs[Rayleigh::ubar][k]     =   2.0;
        rayleigh_ptrs[Rayleigh::vbar][k]     =   1.0;
        rayleigh_ptrs[Rayleigh::wbar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::thetabar][k] = 300.0;
   }
}

void
Problem::init_custom_terrain (const Geometry& geom,
                              MultiFab& z_phys_nd,
                              const Real& time)
{
    // Check if a valid csv file exists for the terrain
    std::string fname;
    amrex::ParmParse pp("erf");
    auto valid_fname = pp.query("terrain_file_name",fname);
    if (valid_fname) {
        this->read_custom_terrain(fname,geom,z_phys_nd,time);
    } else {

        // Domain cell size and real bounds
        auto dx = geom.CellSizeArray();

        // Domain valid box (z_nd is nodal)
        const amrex::Box& domain = geom.Domain();
        int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
        int domlo_z = domain.smallEnd(2);

        // Number of ghost cells
        int ngrow = z_phys_nd.nGrow();

        // Populate bottom plane
        int k0 = domlo_z;

        Real Ampl        = parms.Ampl;
        Real wavelength  = 100.;
        Real kp          = 2.0 * PI / wavelength;
        Real g           = CONST_GRAV;
        Real omega       = std::sqrt(g * kp);

        for (MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Grown box with no z range
            amrex::Box xybx = mfi.growntilebox(ngrow);
            xybx.setRange(2,0);

            amrex::Array4<Real> const& z_arr = z_phys_nd.array(mfi);

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {

                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

                // Location of nodes
                Real x = ii  * dx[0];

                // Wave height
                Real height = Ampl * std::sin(kp * x - omega * time);

                // Populate terrain height
                z_arr(i,j,k0) = height;
            });
        }
    }
}
