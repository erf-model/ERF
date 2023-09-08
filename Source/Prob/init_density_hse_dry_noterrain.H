/**
 * Initialize hydrostatically balanced density -- no-terrain option only
 *
 * Calls init_isentropic_hse(). Hydrostatic equlibrium (HSE) is satisfied
 * discretely. Note that these routines presume that qv==0 when evaluating the
 * EOS. Both density and pressure in HSE are calculated but only the density is
 * used at this point.
*/
void
erf_init_dens_hse(
    amrex::MultiFab& rho_hse,
    std::unique_ptr<amrex::MultiFab>& z_phys_nd,
    std::unique_ptr<amrex::MultiFab>& z_phys_cc,
    amrex::Geometry const& geom) override
{
    const amrex::Real prob_lo_z = geom.ProbLo()[2];
    const amrex::Real dz        = geom.CellSize()[2];
    const int khi               = geom.Domain().bigEnd()[2];

    const amrex::Real T_sfc    = 300.;
    const amrex::Real rho_sfc  = p_0 / (R_d*T_sfc);
    const amrex::Real Thetabar = T_sfc;

    // These are at cell centers (unstaggered)
    amrex::Vector<amrex::Real> h_r(khi+2);
    amrex::Vector<amrex::Real> h_p(khi+2);

    amrex::Gpu::DeviceVector<amrex::Real> d_r(khi+2);
    amrex::Gpu::DeviceVector<amrex::Real> d_p(khi+2);

    HSEutils::init_isentropic_hse(rho_sfc,Thetabar,h_r.data(),h_p.data(),dz,prob_lo_z,khi);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());

    amrex::Real* r = d_r.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(rho_hse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.growntilebox(1);
      const amrex::Array4<amrex::Real> rho_hse_arr = rho_hse[mfi].array();
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
          int kk = std::max(k,0);
          rho_hse_arr(i,j,k) = r[kk];
      });
    } // mfi
}
