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
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& state,
    Array4<Real      > const& x_vel,
    Array4<Real      > const& y_vel,
    Array4<Real      > const& z_vel,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& z_nd,
    Array4<Real const> const& z_cc,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice&)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  const Real rho_sfc   = p_0 / (R_d*parms.T_0);
  //const Real thetabar  = parms.T_0;

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = z_cc(i,j,k);

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

    Real L = std::sqrt(
        std::pow((x - parms.x_c)/parms.x_r, 2) +
        std::pow((z - parms.z_c)/parms.z_r, 2)
    );
    if (L <= 1.0) {
        Real dT = parms.T_pert * (std::cos(PI*L) + 1.0)/2.0;

        // Note: dT is a temperature perturbation, theta_perturbed is base state + perturbation in theta
        Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), R_d/parms.C_p);

        // This version perturbs rho but not p
        state(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
    }

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;

    state(i, j, k, RhoQ1_comp) = 0.0;
    state(i, j, k, RhoQ2_comp) = 0.0;
  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel(i, j, k) = parms.U_0;
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel(i, j, k) = 0.0;
  });

  const auto dx = geomdata.CellSize();
  amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
  dxInv[0] = 1. / dx[0];
  dxInv[1] = 1. / dx[1];
  dxInv[2] = 1. / dx[2];

  // Set the z-velocity from impenetrable condition
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = WFromOmega(i, j, k, 0.0, x_vel, y_vel, z_nd, dxInv);
  });

  amrex::Gpu::streamSynchronize();

}

void
Problem::init_custom_terrain(
    const Geometry& geom,
    MultiFab& z_phys_nd,
    const Real& /*time*/)
{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    auto ProbHiArr = geom.ProbHiArray();

    // Domain valid box (z_nd is nodal)
    const amrex::Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    // int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2);

    // User function parameters
    Real a    = 0.5;
    Real num  = 8 * a * a * a;
    Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
    // Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Populate bottom plane
    int k0 = domlo_z;

    for ( amrex::MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        amrex::Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Clip indices for ghost-cells
            int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
            // int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

            // Location of nodes
            Real x = (ii  * dx[0] - xcen);
            // Real y = (jj  * dx[1] - ycen);

            // WoA Hill in x-direction
            Real height = num / (x*x + 4 * a * a);

            // Populate terrain height
            z_arr(i,j,k0) = height;
        });
    }
}
