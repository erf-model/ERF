#include "ERF.H"

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

    int width = wrfbdy_width;

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
            ncomp = Cons::NumVars;
        }

        if (var_idx == Vars::cons || var_idx == Vars::zvel)
        {

            const auto& dom_lo = amrex::lbound(domain);
            const auto& dom_hi = amrex::ubound(domain);

            MultiFab& mf = *mfs[var_idx];

            IntVect ng_vect = mf.nGrowVect();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const Array4<Real>& dest_arr = mf.array(mfi);
                Box bx  = mfi.tilebox();
                Box gbx = mfi.growntilebox(ng_vect);

                const auto& bx_lo = amrex::lbound(bx);
                const auto& bx_hi = amrex::ubound(bx);

                // Note: "domain" here is the domain for the centering of this variable, i.e.
                //       it is z-face centered for w and cell-centered for the state

                // x-faces
                {
                    if (bx_lo.x == dom_lo.x)
                    {
                        Box bx_xlo(gbx & domain);
                        bx_xlo.setSmall(0,dom_lo.x-ng_vect[0]);
                        bx_xlo.setBig(0,dom_lo.x+width-1);

                        ParallelFor(bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            dest_arr(i,j,k,icomp+n) = dest_arr(dom_lo.x+width,j,k,icomp+n);
                        });
                    } // bx

                    if (bx_hi.x == dom_hi.x)
                    {
                        Box bx_xhi(gbx & domain);
                        bx_xhi.setSmall(0,dom_hi.x-width+1);
                        bx_xhi.setBig(0,dom_hi.x+ng_vect[0]);

                        ParallelFor(bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            dest_arr(i,j,k,icomp+n) = dest_arr(dom_hi.x-width,j,k,icomp+n);
                        });
                    } // bx
                } // x-faces

                // y-faces
                {
                    if (bx_lo.y == dom_lo.y)
                    {
                        Box bx_ylo(gbx & domain);
                        bx_ylo.setSmall(1,dom_lo.y-ng_vect[1]);
                        bx_ylo.setBig(1,dom_lo.y+width-1);

                        ParallelFor(bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            dest_arr(i,j,k,icomp+n) = dest_arr(i,dom_lo.y+width,k,icomp+n);
                        });
                    } // bx

                    if (bx_hi.y == dom_hi.y)
                    {
                        Box bx_yhi(gbx & domain);
                        bx_yhi.setSmall(1,dom_hi.y-width+1);
                        bx_yhi.setBig(1,dom_hi.y+ng_vect[1]);
                        ParallelFor(bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            dest_arr(i,j,k,icomp+n) = dest_arr(i,dom_hi.y-width,k,icomp+n);
                        });
                    } // bx
                } // y-faces

                // corners
                {
                if (bx_lo.x == dom_lo.x && bx_lo.y == dom_lo.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_lo.x - ng_vect[0]);
                    bx_corner.setBig  (0, dom_lo.x - 1);
                    bx_corner.setSmall(1, dom_lo.y - ng_vect[1]);
                    bx_corner.setBig  (1, dom_lo.y - 1);

                    ParallelFor(bx_corner, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        dest_arr(i,j,k,icomp+n) = dest_arr(dom_lo.x,dom_lo.y,k,icomp+n);
                    });
                } // bx

                if (bx_lo.x == dom_lo.x && bx_hi.y == dom_hi.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_lo.x - ng_vect[0]);
                    bx_corner.setBig  (0, dom_lo.x - 1);
                    bx_corner.setSmall(1, dom_hi.y + 1);
                    bx_corner.setBig  (1, dom_hi.y + ng_vect[1]);

                    ParallelFor(bx_corner, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        dest_arr(i,j,k,icomp+n) = dest_arr(dom_lo.x,dom_hi.y,k,icomp+n);
                    });
                } // bx

                if (bx_hi.x == dom_hi.x && bx_lo.y == dom_lo.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_hi.x + 1);
                    bx_corner.setBig  (0, dom_hi.x + ng_vect[0]);
                    bx_corner.setSmall(1, dom_lo.y - ng_vect[1]);
                    bx_corner.setBig  (1, dom_lo.y - 1);

                    ParallelFor(bx_corner, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        dest_arr(i,j,k,icomp+n) = dest_arr(dom_hi.x,dom_lo.y,k,icomp+n);
                    });
                } // bx

                if (bx_hi.x == dom_hi.x && bx_hi.y == dom_hi.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_hi.x + 1);
                    bx_corner.setBig  (0, dom_hi.x + ng_vect[0]);
                    bx_corner.setSmall(1, dom_hi.y + 1);
                    bx_corner.setBig  (1, dom_hi.y + ng_vect[1]);

                    ParallelFor(bx_corner, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        dest_arr(i,j,k,icomp+n) = dest_arr(dom_hi.x,dom_hi.y,k,icomp+n);
                    });
                } // bx

                } // corners
            } // mfi
        } // if (var_idx == zvel or cons)
    } // var_idx

    //
    // *********************************************************************************
    // We now loop over the variables which we did read in from wrfbdy
    // We have read MU and PC but don't use them here so we subtract 2 from NumTypes
    // *********************************************************************************
    //

    Real dT = bdy_time_interval;

    int n_time = static_cast<int>(time / dT);
    amrex::Real alpha = (time - n_time * dT) / dT;
    amrex::Real oma   = 1.0 - alpha;

    for (int ivar = 0; ivar < WRFBdyVars::NumTypes-2; ivar++)
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

        IntVect ng_vect = mf.nGrowVect();

        Box domain = geom[lev].Domain();

        //
        // Note that "domain" is mapped onto the type of box the data is in
        //
        if (var_idx == Vars::xvel) {
           domain.growHi(0,1); domain.setType(amrex::IndexType(IntVect(1,0,0)));
        } else if (var_idx == Vars::yvel) {
           domain.growHi(1,1); domain.setType(amrex::IndexType(IntVect(0,1,0)));
        }

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
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Array4<Real>& dest_arr = mf.array(mfi);
            Box bx  = mfi.tilebox();
            Box gbx = mfi.growntilebox(ng_vect);

            const auto& bx_lo = amrex::lbound(bx);
            const auto& bx_hi = amrex::ubound(bx);

            // Note: "domain" here is the domain for the centering of this variable, i.e.
            //       it is x-face centered for u, y-face-centered for v, and cell-centered for the state

            // x-faces
            {
                if (bx_lo.x == dom_lo.x)
                {
                    Box bx_xlo(gbx & domain);
                    bx_xlo.setSmall(0,dom_lo.x-ng_vect[0]);
                    bx_xlo.setBig(0,dom_lo.x+width-1);

                    ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int ii = std::max(i, dom_lo.x);
                        dest_arr(i,j,k,icomp) = oma   * bdatxlo_n  (ii,j,k,0)
                                              + alpha * bdatxlo_np1(ii,j,k,0);
                    });
                } // bx

                if (bx_hi.x == dom_hi.x)
                {
                    Box bx_xhi(gbx & domain);
                    bx_xhi.setSmall(0,dom_hi.x-width+1);
                    bx_xhi.setBig(0,dom_hi.x+ng_vect[0]);

                    ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int ii = std::min(i, dom_hi.x);
                        dest_arr(i,j,k,icomp) = oma   * bdatxhi_n  (ii,j,k,0)
-                                             + alpha * bdatxhi_np1(ii,j,k,0);
                    });
                } // bx
            } // x-faces

            // y-faces
            {
                if (bx_lo.y == dom_lo.y)
                {
                    Box bx_ylo(gbx & domain);
                    bx_ylo.setSmall(1,dom_lo.y-ng_vect[0]);
                    bx_ylo.setBig(1,dom_lo.y+width-1);

                   ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                   {
                       int jj = std::max(j, dom_lo.y);
                       dest_arr(i,j,k,icomp) = oma   * bdatylo_n  (i,jj,k,0)
                                             + alpha * bdatylo_np1(i,jj,k,0);
                   });
                } // bx

                if (bx_hi.y == dom_hi.y)
                {
                    Box bx_yhi(gbx & domain);
                    bx_yhi.setSmall(1,dom_hi.y-width+1);
                    bx_yhi.setBig(1,dom_hi.y+ng_vect[0]);
                    ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::min(j, dom_hi.y);
                        dest_arr(i,j,k,icomp) = oma   * bdatyhi_n  (i,jj,k,0)
                                              + alpha * bdatyhi_np1(i,jj,k,0);
                    });
                } // bx
            } // y-faces

            // corners
            {
                if (bx_lo.x == dom_lo.x && bx_lo.y == dom_lo.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_lo.x - ng_vect[0]);
                    bx_corner.setBig  (0, dom_lo.x - 1);
                    bx_corner.setSmall(1, dom_lo.y - ng_vect[1]);
                    bx_corner.setBig  (1, dom_lo.y - 1);

                    ParallelFor(bx_corner, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = dest_arr(dom_lo.x,dom_lo.y,k,icomp);
                    });
                } // bx

                if (bx_lo.x == dom_lo.x && bx_hi.y == dom_hi.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_lo.x - ng_vect[0]);
                    bx_corner.setBig  (0, dom_lo.x - 1);
                    bx_corner.setSmall(1, dom_hi.y + 1);
                    bx_corner.setBig  (1, dom_hi.y + ng_vect[1]);

                    ParallelFor(bx_corner, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = dest_arr(dom_lo.x,dom_hi.y,k,icomp);
                    });
                } // bx

                if (bx_hi.x == dom_hi.x && bx_lo.y == dom_lo.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_hi.x + 1);
                    bx_corner.setBig  (0, dom_hi.x + ng_vect[0]);
                    bx_corner.setSmall(1, dom_lo.y - ng_vect[1]);
                    bx_corner.setBig  (1, dom_lo.y - 1);

                    ParallelFor(bx_corner, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = dest_arr(dom_hi.x,dom_lo.y,k,icomp);
                    });
                } // bx

                if (bx_hi.x == dom_hi.x && bx_hi.y == dom_hi.y)
                {
                    Box bx_corner(bx & domain);
                    bx_corner.setSmall(0, dom_hi.x + 1);
                    bx_corner.setBig  (0, dom_hi.x + ng_vect[0]);
                    bx_corner.setSmall(1, dom_hi.y + 1);
                    bx_corner.setBig  (1, dom_hi.y + ng_vect[1]);

                    ParallelFor(bx_corner, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = dest_arr(dom_hi.x,dom_hi.y,k,icomp);
                    });
                } // bx

            } // corners
          } // mfi
    } // ivar
}
#endif
