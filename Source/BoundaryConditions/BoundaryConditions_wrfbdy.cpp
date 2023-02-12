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
    //
    // This uses data read in from WRF real.exe output
    // We assume the data has been read in at the level 0 resolution
    // and that we only fill level 0 MultiFabs
    //

    int lev = 0;

    Real dT = bdy_time_interval;

    int n_time = static_cast<int>(time / dT);
    amrex::Real alpha = (time - n_time * dT) / dT;
    amrex::Real oma   = 1.0 - alpha;

    int ivar;

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
    {
        if (var_idx == Vars::zvel) return;

        MultiFab& mf = *mfs[var_idx];
        const int icomp = 0;
        const int ncomp = 1;

        Box domain = geom[lev].Domain();

        if (var_idx == Vars::xvel) {
           ivar   = WRFBdyVars::U; // U
           domain.growHi(0,1); domain.setType(amrex::IndexType(IntVect(1,0,0)));
        } else if (var_idx == Vars::yvel) {
           ivar   = WRFBdyVars::V; // V
           domain.growHi(1,1); domain.setType(amrex::IndexType(IntVect(0,1,0)));
        } else  if (var_idx == Vars::cons) {
           ivar   = WRFBdyVars::R; // R (we assume that T comes right after R)
        }

        const auto& dom_lo = amrex::lbound(domain);
        const auto& dom_hi = amrex::ubound(domain);

        // We have data at fixed time intervals we will call dT
        // Then to interpolate, given time, we can define n = (time/dT)
        // and alpha = (time - n*dT) / dT, then we define the data at time
        // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)

        for (int nv = 0; nv < ncomp; nv++) {

            const auto& bdatxlo_n   = bdy_data_xlo[n_time  ][ivar+nv].const_array();
            const auto& bdatxlo_np1 = bdy_data_xlo[n_time+1][ivar+nv].const_array();
            const auto& bdatxhi_n   = bdy_data_xhi[n_time  ][ivar+nv].const_array();
            const auto& bdatxhi_np1 = bdy_data_xhi[n_time+1][ivar+nv].const_array();
            const auto& bdatylo_n   = bdy_data_ylo[n_time  ][ivar+nv].const_array();
            const auto& bdatylo_np1 = bdy_data_ylo[n_time+1][ivar+nv].const_array();
            const auto& bdatyhi_n   = bdy_data_yhi[n_time  ][ivar+nv].const_array();
            const auto& bdatyhi_np1 = bdy_data_yhi[n_time+1][ivar+nv].const_array();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(mf); mfi.isValid(); ++mfi)
          {
            const Array4<Real>& dest_arr = mf.array(mfi);
            Box bx = mfi.validbox();

            const auto& bx_lo = amrex::lbound(bx);
            const auto& bx_hi = amrex::ubound(bx);

            // Note: "domain" here is the domain for the centering of this variable, i.e.
            //       it is x-face centered for u, y-face-centered for v, and cell-centered for the state

            // x-faces
            {
                if (bx_lo.x == dom_lo.x)
                {
                    Box bx_xlo(bx & domain);
                    bx_xlo.setBig(0,dom_lo.x);

                    ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp+nv) = oma   * bdatxlo_n  (i,j,k,nv)
                                                 + alpha * bdatxlo_np1(i,j,k,nv);
                    });
                }

                if (bx_hi.x == dom_hi.x)
                {
                    Box bx_xhi(bx & domain);
                    bx_xhi.setSmall(0,dom_hi.x);

                    ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                   {
                       dest_arr(i,j,k,icomp+nv) = oma   * bdatxhi_n  (i,j,k,nv)
                                                + alpha * bdatxhi_np1(i,j,k,nv);
                   });
                }

            } // x-faces

            // y-faces
            {
                if (bx_lo.y == dom_lo.y)
                {
                    Box bx_ylo(bx & domain);
                    bx_ylo.setBig(1,dom_lo.y);

                    ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp+nv) = oma   * bdatylo_n  (i,j,k,nv)
                                                 + alpha * bdatylo_np1(i,j,k,nv);
                    });
                }

                if (bx_hi.y == dom_hi.y)
                {
                    Box bx_yhi(bx & domain);
                    bx_yhi.setSmall(1,dom_hi.y);

                    ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp+nv) = oma   * bdatyhi_n  (i,j,k,nv)
                                                 + alpha * bdatyhi_np1(i,j,k,nv);
                    });
                }
            } // y-faces

          } // mfi
        } // nv
    } // var_idx
}
#endif
