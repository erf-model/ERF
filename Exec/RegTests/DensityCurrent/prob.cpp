#include "prob.H"
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
    const SolverChoice& sc)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  const Real rho_sfc   = p_0 / (R_d*parms.T_0);
  //const Real thetabar  = parms.T_0;
  const Real dz        = geomdata.CellSize()[2];
  const Real prob_lo_z = geomdata.ProbLo()[2];

  const Real l_x_r = parms.x_r;
  //const Real l_x_r = parms.x_r * mf_u(0,0,0); //used to validate constant msf
  const Real l_z_r = parms.z_r;
  const Real l_x_c = parms.x_c;
  const Real l_z_c = parms.z_c;
  const Real l_Tpt = parms.T_pert;

#if 0
  // These are at cell centers (unstaggered)
  Vector<Real> h_r(khi+2);
  Vector<Real> h_p(khi+2);

  amrex::Gpu::DeviceVector<Real> d_r(khi+2);
  amrex::Gpu::DeviceVector<Real> d_p(khi+2);
#endif

  const Real rdOcp = sc.rdOcp;

  if (z_cc) {

#if 0
    // Create a flat box with same horizontal extent but only one cell in vertical
    Box b2d = surroundingNodes(bx); // Copy constructor
    b2d.setRange(2,0);

    ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
         Array1D<Real,0,255> r;;
         Array1D<Real,0,255> p;;

         init_isentropic_hse_terrain(i,j,rho_sfc,thetabar,&(r(0)),&(p(0)),z_cc,khi);

         for (int k = 0; k <= khi; k++) {
            r_hse(i,j,k) = r(k);
            p_hse(i,j,k) = p(k);
         }
         r_hse(i,j,   -1) = r_hse(i,j,0);
         r_hse(i,j,khi+1) = r_hse(i,j,khi);
      });
#endif

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = z_cc(i,j,k);

#if 0
        // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
        const Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

        Real L = std::sqrt(
            std::pow((x - l_x_c)/l_x_r, 2) +
            std::pow((z - l_z_c)/l_z_r, 2)
        );
        Real dT;
        if (L > 1.0) {
            dT = 0.0;
        }
        else {
            dT = l_Tpt * (std::cos(PI*L) + 1.0)/2.0;
        }

        // Note: dT is a perturbation in temperature, theta_perturbed is theta PLUS perturbation in theta
        Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), rdOcp);

        // This version perturbs rho but not p
        state(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p_hse(i,j,k));
        state(i, j, k, Rho_comp) = state(i, j, k, RhoTheta_comp) / theta_perturbed;
#endif

        Real L = std::sqrt(
            std::pow((x - l_x_c)/l_x_r, 2) +
            std::pow((z - l_z_c)/l_z_r, 2));
        if (L <= 1.0)
        {
            Real dT = l_Tpt * (std::cos(PI*L) + 1.0)/2.0;
            Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

            // Note: dT is a perturbation in temperature, theta_perturbed is base state + perturbation
            Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), rdOcp);

            state(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
        }

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
  } else {

#if 0
      init_isentropic_hse(rho_sfc,thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

      amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
      amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

      Real* r = d_r.data();
      Real* p = d_p.data();
#endif

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

#if 0
        // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
        const Real Tbar_hse = p[k] / (R_d * r[k]);

        Real L = std::sqrt(
            std::pow((x - l_x_c)/l_x_r, 2) +
            std::pow((z - l_z_c)/l_z_r, 2)
        );
        Real dT;
        if (L > 1.0) {
            dT = 0.0;
        }
        else {
            dT = l_Tpt * (std::cos(PI*L) + 1.0)/2.0;
        }

        // Note: dT is a perturbation in temperature, theta_perturbed is theta PLUS perturbation in theta
        Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p[k], rdOcp);

        // This version perturbs rho but not p
        state(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p[k]);
        state(i, j, k, Rho_comp) = state(i, j, k, RhoTheta_comp) / theta_perturbed;

        if ((i==0) && (j==0))
        {
            amrex::Print() << "init_custom_pert("<<i<<","<<j<<","<<k<<")"
                << " r["<<k<<"]=" << r[k] << " r_hse=" << r_hse(i,j,k) << " diff: " << r[k] - r_hse(i,j,k)
                << " p["<<k<<"]=" << p[k] << " p_hse=" << p_hse(i,j,k) << " diff: " << p[k] - p_hse(i,j,k)
                << std::endl;
        }
#endif

        Real L = std::sqrt(
            std::pow((x - l_x_c)/l_x_r, 2) +
            std::pow((z - l_z_c)/l_z_r, 2));
        if (L <= 1.0)
        {
            Real dT = l_Tpt * (std::cos(PI*L) + 1.0)/2.0;
            Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

            // Note: dT is a perturbation in temperature, theta_perturbed is base state + perturbation
            Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), rdOcp);

            state(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
        }

        // Set scalar = 0 everywhere
        state(i, j, k, RhoScalar_comp) = 0.0;

#ifdef ERF_USE_MOISTURE
        state(i, j, k, RhoQt_comp) = 0.0;
        state(i, j, k, RhoQp_comp) = 0.0;
#endif
      });
  }

  const Real u0 = parms.U_0;

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel(i, j, k) = u0;
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

