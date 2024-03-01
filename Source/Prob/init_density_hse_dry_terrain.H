/**
 * Initialize hydrostatically balanced density -- assumes non-uniform z
 *
 * Calls init_isentropic_hse_terrain(). Hydrostatic equlibrium (HSE) is satisfied
 * discretely. Note that these routines presume that qv==0 when evaluating the
 * EOS. Both density and pressure in HSE are calculated but only the density is
 * used at this point.
*/
void
erf_init_dens_hse (amrex::MultiFab& rho_hse,
                   std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                   std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                   amrex::Geometry const& geom) override
{
    const int khi               = geom.Domain().bigEnd()[2];

    const amrex::Real T_sfc    = 300.;
    const amrex::Real rho_sfc  = p_0 / (R_d*T_sfc);
    const amrex::Real Thetabar = T_sfc;

    if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

    for ( amrex::MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
    {
        amrex::Array4<amrex::Real      > rho_arr  = rho_hse.array(mfi);
        amrex::Array4<amrex::Real const> z_cc_arr = z_phys_cc->const_array(mfi);

        // Create a flat box with same horizontal extent but only one cell in vertical
        const amrex::Box& tbz = mfi.nodaltilebox(2);
        amrex::Box b2d = tbz; // Copy constructor
        b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
        b2d.setRange(2,0);

        amrex::ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
          amrex::Array1D<amrex::Real,0,255> r;;
          amrex::Array1D<amrex::Real,0,255> p;;

          HSEutils::init_isentropic_hse_terrain(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

          for (int k = 0; k <= khi; k++) {
             rho_arr(i,j,k) = r(k);
          }
          rho_arr(i,j,   -1) = rho_arr(i,j,0);
          rho_arr(i,j,khi+1) = rho_arr(i,j,khi);
        });
    } // mfi
}
