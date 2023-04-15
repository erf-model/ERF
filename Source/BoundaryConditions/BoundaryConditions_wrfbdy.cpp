#include "ERF.H"
#include "Utils.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
//
// mf is the MultiFab to be filled with data read in from wrfbdy
// time is the time at which the data should be filled
//
void
ERF::fill_from_wrfbdy (const Vector<MultiFab*>& mfs, const Real time)
{
    int lev = 0;

    // Always copy unread vars into relaxation & set region
    int width = wrfbdy_width - 1;

    //
    // *********************************************************************************
    // First we loop over the variables which were NOT read in from wrfbdy
    // *********************************************************************************
    //
    for (int var_idx = Vars::cons; var_idx <= Vars::zvel; var_idx++)
    {
        int icomp = 0;
        int ncomp;
        Box domain = geom[lev].Domain();

        if (var_idx == Vars::zvel) {
           ncomp = 1;
           domain.growHi(2,1); domain.setType(amrex::IndexType(IntVect(0,0,1)));
        } else if (var_idx == Vars::cons) {
            // NOTE: Don't overwrite relaxation zone data!
            ncomp = Cons::NumVars - 2; // Rho & RhoTheta
            icomp = RhoTheta_comp + 1; // Start afer RhoTheta
        }

        if (var_idx == Vars::cons || var_idx == Vars::zvel)
        {
            const auto& dom_lo = amrex::lbound(domain);
            const auto& dom_hi = amrex::ubound(domain);

            MultiFab& mf = *mfs[var_idx];

            IntVect ng_vect = mf.nGrowVect(); ng_vect[2] = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Grown tilebox so we fill exterior ghost cells as well
                Box gbx = mfi.growntilebox(ng_vect);
                const Array4<Real>& dest_arr = mf.array(mfi);

                Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
                compute_interior_ghost_bxs_xy(gbx, domain, width, 0,
                                              bx_xlo, bx_xhi,
                                              bx_ylo, bx_yhi, ng_vect);

                // x-faces (includes y ghost cells)
                ParallelFor(bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    int jj = std::max(j , dom_lo.y);
                        jj = std::min(jj, dom_hi.y);
                    dest_arr(i,j,k,icomp+n) = dest_arr(dom_lo.x+width,jj,k,icomp+n);
                });
                ParallelFor(bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    int jj = std::max(j , dom_lo.y);
                        jj = std::min(jj, dom_hi.y);
                    dest_arr(i,j,k,icomp+n) = dest_arr(dom_hi.x-width,jj,k,icomp+n);
                });

                // y-faces (does not include x ghost cells)
                ParallelFor(bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    dest_arr(i,j,k,icomp+n) = dest_arr(i,dom_lo.y+width,k,icomp+n);
                });
                ParallelFor(bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    dest_arr(i,j,k,icomp+n) = dest_arr(i,dom_hi.y-width,k,icomp+n);
                });
            } // mfi
        } // if (var_idx == zvel or cons)
    } // var_idx

    //
    // *********************************************************************************
    // We now loop over the variables which we did read in from wrfbdy
    // We have read MU and PC but don't use them here so we subtract 2 from NumTypes
    // *********************************************************************************
    //

    // Only populate the set region with read vars
    width = wrfbdy_set_width;

    Real dT = bdy_time_interval;

    int n_time = static_cast<int>(time / dT);
    amrex::Real alpha = (time - n_time * dT) / dT;
    amrex::Real oma   = 1.0 - alpha;

    for (int ivar = 0; ivar < WRFBdyVars::NumTypes-2; ivar++)
    //for (int ivar = 0; ivar <= WRFBdyVars::T; ivar++)
    {
        int icomp   = -1;
        int var_idx = -1;

        if (ivar  == WRFBdyVars::U) {
            var_idx = Vars::xvel;
            icomp   = 0;
        } else if (ivar  == WRFBdyVars::V) {
            var_idx = Vars::yvel;
            icomp   = 0;
        } else if (ivar  == WRFBdyVars::R) {
            var_idx = Vars::cons;
            icomp = Rho_comp;
        } else if (ivar  == WRFBdyVars::T) {
            var_idx = Vars::cons;
            icomp = RhoTheta_comp;
#if defined(ERF_USE_MOISTURE)
        } else if (ivar  == WRFBdyVars::QV) {
            var_idx = Vars::cons;
            icomp = RhoQt_comp;
#elif defined(ERF_USE_WARM_NO_PRECIP)
        } else if (ivar  == WRFBdyVars::QV) {
            var_idx = Vars::cons;
            icomp = RhoQv_comp;
#endif
        } else {
            return;
        }

        MultiFab& mf = *mfs[var_idx];

        IntVect ng_vect = mf.nGrowVect(); ng_vect[2] = 0;

        Box domain = geom[lev].Domain();

        //
        // Note that "domain" is mapped onto the type of box the data is in
        //
        domain.convert(mf.boxArray().ixType());

        const auto& dom_lo = amrex::lbound(domain);
        const auto& dom_hi = amrex::ubound(domain);

        // We have data at fixed time intervals we will call dT
        // Then to interpolate, given time, we can define n = (time/dT)
        // and alpha = (time - n*dT) / dT, then we define the data at time
        // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)
        const auto& bdatxlo_n   = bdy_data_xlo[n_time  ][ivar].const_array();
        const auto& bdatxlo_np1 = bdy_data_xlo[n_time+1][ivar].const_array();
        const auto& bdatxhi_n   = bdy_data_xhi[n_time  ][ivar].const_array();
        const auto& bdatxhi_np1 = bdy_data_xhi[n_time+1][ivar].const_array();
        const auto& bdatylo_n   = bdy_data_ylo[n_time  ][ivar].const_array();
        const auto& bdatylo_np1 = bdy_data_ylo[n_time+1][ivar].const_array();
        const auto& bdatyhi_n   = bdy_data_yhi[n_time  ][ivar].const_array();
        const auto& bdatyhi_np1 = bdy_data_yhi[n_time+1][ivar].const_array();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Grown tilebox so we fill exterior ghost cells as well
            Box gbx = mfi.growntilebox(ng_vect);
            const Array4<Real>& dest_arr = mf.array(mfi);

            // Call w/o interior ghost cells
            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            compute_interior_ghost_bxs_xy(gbx, domain, width, 0,
                                          bx_xlo, bx_xhi,
                                          bx_ylo, bx_yhi, ng_vect);

            // x-faces (includes exterior y ghost cells)
            ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int ii = std::max(i , dom_lo.x);
                int jj = std::max(j , dom_lo.y);
                    jj = std::min(jj, dom_hi.y);
                dest_arr(i,j,k,icomp) = oma   * bdatxlo_n  (ii,jj,k,0)
                                      + alpha * bdatxlo_np1(ii,jj,k,0);
            });
            ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int ii = std::min(i , dom_hi.x);
                int jj = std::max(j , dom_lo.y);
                    jj = std::min(jj, dom_hi.y);
                dest_arr(i,j,k,icomp) = oma   * bdatxhi_n  (ii,jj,k,0)
                                      + alpha * bdatxhi_np1(ii,jj,k,0);
            });

            // y-faces (do not include exterior x ghost cells)
            ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int jj = std::max(j , dom_lo.y);
                dest_arr(i,j,k,icomp) = oma   * bdatylo_n  (i,jj,k,0)
                                      + alpha * bdatylo_np1(i,jj,k,0);
            });
            ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int jj = std::min(j , dom_hi.y);
                dest_arr(i,j,k,icomp) = oma   * bdatyhi_n  (i,jj,k,0)
                                      + alpha * bdatyhi_np1(i,jj,k,0);
            });
          } // mfi
    } // ivar
}
#endif
