#include "prob.H"
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
    Array4<Real const> const& z_nd,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const int khi = geomdata.Domain().bigEnd()[2];

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

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
    ParallelFor(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = parms_d.u_0;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = 0.0;
    });

    const auto dx = geomdata.CellSize();
    GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = 1. / dx[0];
    dxInv[1] = 1. / dx[1];
    dxInv[2] = 1. / dx[2];

    // Set the z-velocity from impenetrable condition
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
#ifdef ERF_USE_EB
        z_vel_pert(i, j, k) = 0.0;
#else
        z_vel_pert(i, j, k) = WFromOmega(i, j, k, 0.0, x_vel_pert, y_vel_pert, z_nd, dxInv);
#endif
    });

    amrex::Gpu::streamSynchronize();
}

    /**
    * Function to perform custom initialization of terrain
    *
    * This version takes a single FArrayBox instead of a MultiFab
    *
    */
    void
    Problem::init_custom_terrain (const amrex::Geometry& geom,
                                  amrex::FArrayBox& z_phys_nd,
                                  const amrex::Real& /*time*/)
    {
        // Bottom of domain
        int k0 = 0;

        // Domain valid box (z_nd is nodal)
        const amrex::Box& domain = geom.Domain();
        int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;

        const Real* prob_lo = geom.ProbLo();
        const Real* prob_hi = geom.ProbHi();

        const Real* dx = geom.CellSize();

        // User function parameters
        Real a    = 0.5;
        Real num  = 8 * a * a * a;
        Real xcen = 0.5 * (prob_lo[0] + prob_hi[0]);

        // Grown box with no z range
        amrex::Box bx = z_phys_nd.box();
        bx.setRange(2,0);

        amrex::Array4<amrex::Real> const& z_arr = z_phys_nd.array();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Clip indices for ghost-cells
            int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

            // Location of nodes
            Real x = (ii  * dx[0] - xcen);

            // WoA Hill in x-direction
            Real height = num / (x*x + 4 * a * a);

            // Populate terrain height
            z_arr(i,j,k0) = height;

            // z_arr(i,j,k0) = 2.5;
        });
    }
