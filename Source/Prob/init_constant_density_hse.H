void
erf_init_dens_hse(amrex::MultiFab& rho_hse,
                  std::unique_ptr<amrex::MultiFab>&,
                  std::unique_ptr<amrex::MultiFab>&,
                  amrex::Geometry const&) override
{
    amrex::Real rho_0 = parms.rho_0;
    for ( amrex::MFIter mfi(rho_hse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
       amrex::Array4<amrex::Real> rho_arr = rho_hse.array(mfi);
       const amrex::Box& gbx = mfi.growntilebox(1);
       ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           rho_arr(i,j,k) = rho_0;
       });
    }
}
