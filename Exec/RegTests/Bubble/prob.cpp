#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "TileNoZ.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

// TODO: reorder function declarations for consistency

void
Problem::erf_init_dens_hse(
    MultiFab& rho_hse,
    std::unique_ptr<MultiFab>& z_phys_nd,
    std::unique_ptr<MultiFab>& z_phys_cc,
    Geometry const& geom)
{
    // if erf.init_type is specified, then the base density should
    // already have been calculated. In that case, these assumed (dry)
    // surface conditions are probably incorrect.
    AMREX_ALWAYS_ASSERT(parms.T_0 > 0);

    const Real prob_lo_z = geom.ProbLo()[2];
    const Real dz        = geom.CellSize()[2];
    const int khi        = geom.Domain().bigEnd()[2];

    const Real T_sfc    = 300.;
    const Real rho_sfc  = p_0 / (R_d*T_sfc);
    const Real Thetabar = T_sfc;

    // use_terrain = 1
    if (z_phys_nd) {

        if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

        for ( MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > rho_arr  = rho_hse.array(mfi);
            Array4<Real const> z_cc_arr = z_phys_cc->const_array(mfi);

            // Create a flat box with same horizontal extent but only one cell in vertical
            const Box& tbz = mfi.nodaltilebox(2);
            Box b2d = tbz; // Copy constructor
            b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
            b2d.setRange(2,0);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
              Array1D<Real,0,255> r;;
              Array1D<Real,0,255> p;;

              init_isentropic_hse_terrain(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

              for (int k = 0; k <= khi; k++) {
                 rho_arr(i,j,k) = r(k);
              }
              rho_arr(i,j,   -1) = rho_arr(i,j,0);
              rho_arr(i,j,khi+1) = rho_arr(i,j,khi);
            });
        } // mfi

    } else { // use_terrain = 0

        // These are at cell centers (unstaggered)
        Vector<Real> h_r(khi+2);
        Vector<Real> h_p(khi+2);

        amrex::Gpu::DeviceVector<Real> d_r(khi+2);
        amrex::Gpu::DeviceVector<Real> d_p(khi+2);

        init_isentropic_hse(rho_sfc,Thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

        Real* r = d_r.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
          for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              const Box& bx = mfi.growntilebox(1);
              const Array4<Real> rho_hse_arr = rho_hse[mfi].array();
              ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
              {
                  int kk = std::max(k,0);
                  rho_hse_arr(i,j,k) = r[kk];
              });
          } // mfi
    } // no terrain
}

void
Problem::init_custom_prob(
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
    Array4<Real const> const& /*z_nd*/,
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
    const Real thetabar  = parms.T_0;
    const Real dz        = geomdata.CellSize()[2];
    const Real prob_lo_z = geomdata.ProbLo()[2];
    const Real rdOcp     = sc.rdOcp;

    // These are at cell centers (unstaggered)
    Vector<Real> h_r(khi+1);
    Vector<Real> h_p(khi+1);

    amrex::Gpu::DeviceVector<Real> d_r(khi+1);
    amrex::Gpu::DeviceVector<Real> d_p(khi+1);

    amrex::Print() << "Bubble delta T = " << parms.T_pert << " K" << std::endl;
    amrex::Print() << "  centered at ("
        << parms.x_c << " " << parms.y_c << " " << parms.z_c << ")" << std::endl;
    amrex::Print() << "  with extent ("
        << parms.x_r << " " << parms.y_r << " " << parms.z_r << ")" << std::endl;

    if (parms.T_0 <= 0)
    {
        amrex::Print() << "Ignoring parms.T_0 = " << parms.T_0
            << ", background fields should have been initialized with erf.init_type"
            << std::endl;
    }

    if (z_cc) { // nonflat terrain

#if 0
        if (parms.T_0 > 0)
        {
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
        }
#endif

        amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const Real z = z_cc(i,j,k);

            perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k),
                              parms, rdOcp,
                              state(i, j, k, Rho_comp),
                              state(i, j, k, RhoTheta_comp));

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
        if (parms.T_0 > 0)
        {
            init_isentropic_hse(rho_sfc,thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

            Real* r = d_r.data();
            Real* p = d_p.data();

            amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                for (int k = 0; k <= khi; k++) {
                   r_hse(i,j,k) = r[k];
                   p_hse(i,j,k) = p[k];
                }
                r_hse(i,j,   -1) = r_hse(i,j,0);
                r_hse(i,j,khi+1) = r_hse(i,j,khi);
            });
        }
#endif

        amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const Real z = prob_lo[2] + (k + 0.5) * dx[2];

            perturb_rho_theta(x, y, z, p_hse(i,j,k), r_hse(i,j,k),
                              parms, rdOcp,
                              state(i, j, k, Rho_comp),
                              state(i, j, k, RhoTheta_comp));

            state(i, j, k, RhoScalar_comp) = 0.0;

#ifdef ERF_USE_MOISTURE
            state(i, j, k, RhoQt_comp) = 0.0;
            state(i, j, k, RhoQp_comp) = 0.0;
#endif
        });
    }

    const Real u0 = parms.U_0;
    const Real v0 = parms.V_0;
    const Real w0 = parms.W_0;

    // Set the x-velocity
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel(i, j, k) = u0;
    });

    // Set the y-velocity
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel(i, j, k) = v0;
    });

    // Set the z-velocity
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel(i, j, k) = w0;
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

Problem::Problem()
{
  // Parse params
  ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("x_c", parms.x_c);
  pp.query("y_c", parms.y_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("y_r", parms.y_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
  pp.query("T_pert_is_airtemp", parms.T_pert_is_airtemp);
  pp.query("perturb_rho", parms.perturb_rho);
  pp.query("dampcoef", parms.dampcoef);
  pp.query("zdamp", parms.zdamp);
}
