#include "prob.H"
#include "prob_common.H"

#include "EOS.H"

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
  ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
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
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real const> const&,
    Array4<Real const> const&,
#if defined(ERF_USE_MOISTURE)
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
#elif defined(ERF_USE_WARM_NO_PRECIP)
    Array4<Real      > const&,
    Array4<Real      > const&,
#endif
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice&)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
//const Real z0 = (0.5) * dx[2] + prob_lo[2];
//const Real Tbar = parms.T_0 - z0 * CONST_GRAV / parms.C_p;
//const Real pbar = p_0 * std::pow(Tbar/parms.T_0, R_d/parms.C_p); // from Straka1993
//const Real pbar = p_0 * std::pow(Tbar/parms.T_0, parms.C_p/R_d); // isentropic relation, consistent with exner pressure def
//const Real rhobar = pbar / (R_d*Tbar); // UNUSED

  const Real& rho_sfc   = p_0 / (R_d*parms.T_0);
  const Real& thetabar  = parms.T_0;
  const Real& dz        = geomdata.CellSize()[2];
  const Real& prob_lo_z = geomdata.ProbLo()[2];

  // These are at cell centers (unstaggered)
  Vector<Real> h_r(khi+2);
  Vector<Real> h_p(khi+2);

  amrex::Gpu::DeviceVector<Real> d_r(khi+2);
  amrex::Gpu::DeviceVector<Real> d_p(khi+2);

  init_isentropic_hse(rho_sfc,thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

  amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
  amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

  Real* r = d_r.data();
  Real* p = d_p.data();

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo         = geomdata.ProbLo();
    const auto dx              = geomdata.CellSize();

    const Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const Real z = prob_lo[2] + (k + 0.5) * dx[2];

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const Real Tbar_hse = p[k] / (R_d * r[k]);

    Real L = std::sqrt(
        std::pow((x - parms.x_c)/parms.x_r, 2) +
        std::pow((z - parms.z_c)/parms.z_r, 2)
    );
    Real dT;
    if (L > 1.0) {
        dT = 0.0;
    }
    else {
        dT = parms.T_pert * (std::cos(PI*L) + 1.0)/2.0;
    }

    // Note: dT is a perturbation in temperature, theta_perturbed is theta PLUS perturbation in theta
    Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p[k], R_d/parms.C_p);

    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = std::pow(p[k]/p_0,1.0/Gamma) * p_0 / R_d;
    state(i, j, k, Rho_comp) = state(i, j, k, RhoTheta_comp) / theta_perturbed;

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;

#if defined(ERF_USE_MOISTURE)
    state(i, j, k, RhoQt_comp) = 0.0;
    state(i, j, k, RhoQp_comp) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
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

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_terrain(
    const Geometry& /*geom*/,
    MultiFab& z_phys_nd,
    const Real& /*time*/)
{
    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}
