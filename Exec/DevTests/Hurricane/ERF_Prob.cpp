#include "ERF_Prob.H"
#include <ERF_Constants.H>

#include "ERF_Interpolation_Bilinear.H"
#include "ERF_ReadCustomBinaryIC.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem ()
{
  // Parse params
  amrex::ParmParse pp("prob");
}

void
Problem::erf_init_dens_hse_moist (MultiFab& rho_hse,
                                  std::unique_ptr<MultiFab>& z_phys_nd,
                                  Geometry const& geom)
{
    const Real prob_lo_z = geom.ProbLo()[2];
    const Real dz        = geom.CellSize()[2];
    const int khi        = geom.Domain().bigEnd()[2];

    // use_terrain = 1
    if (z_phys_nd) {

        if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

        for ( MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > rho_arr  = rho_hse.array(mfi);
            //Array4<Real const> z_cc_arr = z_phys_cc->const_array(mfi);

            // Create a flat box with same horizontal extent but only one cell in vertical
            const Box& tbz = mfi.nodaltilebox(2);
            Box b2d = tbz; // Copy constructor
            b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
            b2d.setRange(2,0);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
              Array1D<Real,0,255> r;

              //init_isentropic_hse_terrain(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

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
        Vector<Real> h_t(khi+2);
        Vector<Real> h_q_v(khi+2);

        amrex::Gpu::DeviceVector<Real> d_r(khi+2);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());

        Real* r     = d_r.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
          for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              const Box& bx = mfi.growntilebox(1);
              const Array4<Real> rho_hse_arr   = rho_hse[mfi].array();
              ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
              {
                  int kk = std::max(k,0);
                  rho_hse_arr(i,j,k) = 0.0;//r[kk];
              });
          } // mfi
    } // no terrain
}

void
Problem::init_custom_pert (
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
    const SolverChoice& sc,
    const int /*lev*/)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];

  // Call the routine to calculate the 1d background condition

   Vector<Real> h_r(khi+2);
   Vector<Real> h_p(khi+2);
   Vector<Real> h_t(khi+2);
   Vector<Real> h_q_v(khi+2);

   amrex::Gpu::DeviceVector<Real> d_r(khi+2);
   amrex::Gpu::DeviceVector<Real> d_p(khi+2);
   amrex::Gpu::DeviceVector<Real> d_t(khi+2);
   amrex::Gpu::DeviceVector<Real> d_q_v(khi+2);

   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

    Real* t   = d_t.data();
    Real* p   = d_p.data();

    // File to read
    std::string filename;
    ParmParse pp("erf");
    pp.query("hindcast_IC_filename", filename);

    if (filename.empty()) {
        if ( (sc.init_type == InitType::WRFInput) ||
             (sc.init_type == InitType::Metgrid)  ||
             (sc.init_type == InitType::NCFile) ) {
            return;
        } else {
            amrex::Abort("Error: IC_file is not specified in the input file.");
        }
    }

    Vector<Real> latvec_h, lonvec_h;
    Vector<Real> xvec_h, yvec_h, zvec_h;
    Vector<Real> rho_h, uvel_h, vvel_h, wvel_h, theta_h, qv_h, qc_h, qr_h;

    ReadCustomBinaryIC(filename, latvec_h, lonvec_h,
                       xvec_h, yvec_h, zvec_h,
                       rho_h, uvel_h, vvel_h, wvel_h,
                       theta_h, qv_h, qc_h, qr_h);

    int nx = xvec_h.size();
    int ny = yvec_h.size();
    int nz = zvec_h.size();

    amrex::Real dxvec = (xvec_h[nx-1]-xvec_h[0])/(nx-1);
    amrex::Real dyvec = (yvec_h[ny-1]-yvec_h[0])/(ny-1);

    amrex::Gpu::DeviceVector<Real> xvec_d(nx*ny*nz), yvec_d(nx*ny*nz), zvec_d(nx*ny*nz);
    amrex::Gpu::DeviceVector<Real> rho_d(nx*ny*nz), uvel_d(nx*ny*nz), vvel_d(nx*ny*nz), wvel_d(nx*ny*nz),
                                   theta_d(nx*ny*nz), qv_d(nx*ny*nz), qc_d(nx*ny*nz), qr_d(nx*ny*nz);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xvec_h.begin(), xvec_h.end(), xvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yvec_h.begin(), yvec_h.end(), yvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, uvel_h.begin(), uvel_h.end(), uvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, vvel_h.begin(), vvel_h.end(), vvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, wvel_h.begin(), wvel_h.end(), wvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qr_h.begin(), qr_h.end(), qr_d.begin());

    Real* xvec_d_ptr = xvec_d.data();
    Real* yvec_d_ptr = yvec_d.data();
    Real* zvec_d_ptr = zvec_d.data();
    Real* rho_d_ptr   = rho_d.data();
    Real* uvel_d_ptr  = uvel_d.data();
    Real* vvel_d_ptr  = vvel_d.data();
    Real* wvel_d_ptr  = wvel_d.data();
    Real* theta_d_ptr = theta_d.data();
    Real* qv_d_ptr = qv_d.data();
    Real* qc_d_ptr = qc_d.data();
    Real* qr_d_ptr = qr_d.data();

    // Interpolate the data on to the ERF mesh

     ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real x        = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y        = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z        = prob_lo[2] + (k + 0.5) * dx[2];

        // First interpolate where the weather data is available from
        Real tmp_rho, tmp_theta, tmp_qv, tmp_qc, tmp_qr;
        bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                               dxvec, dyvec,
                               nx, ny, nz,
                               x, y, z,
                               rho_d_ptr, tmp_rho);
        bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                               dxvec, dyvec,
                               nx, ny, nz,
                               x, y, z,
                               theta_d_ptr, tmp_theta);

        bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                               dxvec, dyvec,
                               nx, ny, nz,
                               x, y, z,
                               qv_d_ptr, tmp_qv);

        bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                               dxvec, dyvec,
                               nx, ny, nz,
                               x, y, z,
                               qc_d_ptr, tmp_qc);

        bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                               dxvec, dyvec,
                               nx, ny, nz,
                               x, y, z,
                               qr_d_ptr, tmp_qr);

        state_pert(i, j, k, Rho_comp)      = tmp_rho;
        state_pert(i, j, k, RhoTheta_comp) = tmp_rho*tmp_theta - 348.432055749129;

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 0.0;//rho*scalar;

           // mean states
           if (use_moisture) {
               state_pert(i, j, k, RhoQ1_comp) = tmp_rho*tmp_qv;
               state_pert(i, j, k, RhoQ2_comp) = tmp_rho*tmp_qc;
               state_pert(i, j, k, RhoQ3_comp) = tmp_rho*tmp_qr;
           }
    });

  // Set the x-velocity
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x        = prob_lo[0] + i * dx[0];
    const Real y        = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z        = prob_lo[2] + (k + 0.5) * dx[2];

    Real tmp_uvel;
    bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                             dxvec, dyvec,
                             nx, ny, nz,
                             x, y, z,
                             uvel_d_ptr, tmp_uvel);
    x_vel_pert(i,j,k) = tmp_uvel;
  });

  // Set the y-velocity
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x        = prob_lo[0] + (i+0.5) * dx[0];
    const Real y        = prob_lo[1] + j * dx[1];
    const Real z        = prob_lo[2] + (k + 0.5) * dx[2];

    Real tmp_vvel;
    bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                           dxvec, dyvec,
                           nx, ny, nz,
                           x, y, z,
                           vvel_d_ptr, tmp_vvel);

      y_vel_pert(i, j, k) = tmp_vvel;
  });

  // Set the z-velocity
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x        = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y        = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z        = prob_lo[2] + k * dx[2];

    Real tmp_wvel;
    bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                           dxvec, dyvec,
                           nx, ny, nz,
                           x, y, z,
                           wvel_d_ptr, tmp_wvel);

      z_vel_pert(i, j, k) = tmp_wvel;
  });

  Gpu::streamSynchronize();
}

void
Problem::erf_init_rayleigh(
    amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
    amrex::Geometry      const& geom,
    std::unique_ptr<MultiFab>& /*z_phys_nd*/,
    amrex::Real /*zdamp*/)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      rayleigh_ptrs[Rayleigh::ubar][k]     = 2.0;
      rayleigh_ptrs[Rayleigh::vbar][k]     = 1.0;
      rayleigh_ptrs[Rayleigh::wbar][k]     = 0.0;
      rayleigh_ptrs[Rayleigh::thetabar][k] = parms.T_0;
  }
}
