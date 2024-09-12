/**
 * Initialize hydrostatically balanced density
 *
 * Calls init_isentropic_hse_terrain() or init_isentropic_hse() for cases with
 * and without terrain. Hydrostatic equilibrium (HSE) is satisfied discretely.
 * Note that these routines presume that qv==0 when evaluating the EOS. Both
 * density and pressure in HSE are calculated but only the density is used at
 * this point.
*/
void
erf_init_dens_hse (amrex::MultiFab& rho_hse,
                   std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                   std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                   amrex::Geometry const& geom) override
{
    const amrex::Real dz        = geom.CellSize()[2];

    const amrex::Real T_sfc    = 300.;
    const amrex::Real rho_sfc  = p_0 / (R_d*T_sfc);
    const amrex::Real Thetabar = T_sfc;

    const int domlo_z = geom.Domain().smallEnd(2);
    const int domhi_z = geom.Domain().bigEnd(2);

    if (domhi_z > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

    bool use_terrain = (z_phys_nd != nullptr);

    if (use_terrain) {

        for ( amrex::MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
        {
            amrex::Array4<amrex::Real      > rho_arr  = rho_hse.array(mfi);
            amrex::Array4<amrex::Real const> z_cc_arr = (use_terrain) ?
                                                 z_phys_cc->const_array(mfi) : amrex::Array4<amrex::Real const>{};

            // Create a flat box with same horizontal extent but only one cell in vertical
            const amrex::Box& tbz = mfi.nodaltilebox(2);
            amrex::Box b2d = tbz; // Copy constructor
            b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
            b2d.setRange(2,0);

            const int ilo = tbz.smallEnd(0);
            const int jlo = tbz.smallEnd(1);
            const int klo = tbz.smallEnd(2);
            const int khi = tbz.bigEnd(2)-1;

            amrex::Real rho_local_sfc;
            if (klo == 0) {
                rho_local_sfc = rho_sfc;
            } else {
                rho_local_sfc = rho_arr(ilo,jlo,klo-1);
            }

            amrex::ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                amrex::Array1D<amrex::Real,0,255> r;;
                amrex::Array1D<amrex::Real,0,255> p;;
                HSEutils::init_isentropic_hse_terrain(i,j,rho_local_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,klo,khi);

                for (int k = klo; k <= khi; k++) {
                   rho_arr(i,j,k) = r(k);
                }

                // Impose Neumann conditions at bottom and top of domain boundary
                if (klo == domlo_z) {
                    rho_arr(i,j,domlo_z-1) = rho_arr(i,j,domlo_z);
                }
                if (khi == domhi_z) {
                    rho_arr(i,j,domhi_z+1) = rho_arr(i,j,domhi_z);
                }
            });
        } // mfi

    } else {

        // Note that this integrates over the entire domain height,
        //      but only copies into each box as appropriate.
        // There is no assumption that any one box spans the whole vertical height.
        const int klo = geom.Domain().smallEnd(2);
        const int khi = geom.Domain().bigEnd(2);

        // These are at cell centers (unstaggered)
        amrex::Vector<amrex::Real> h_r(khi+2);
        amrex::Vector<amrex::Real> h_p(khi+2);

        amrex::Gpu::DeviceVector<amrex::Real> d_r(khi+2);
        amrex::Gpu::DeviceVector<amrex::Real> d_p(khi+2);

        HSEutils::init_isentropic_hse(rho_sfc,Thetabar,h_r.data(),h_p.data(),dz,klo,khi);

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
    } // no terrain
}
