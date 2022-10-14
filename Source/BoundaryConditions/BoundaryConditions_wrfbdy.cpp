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
    const Box& domain = geom[lev].Domain();

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

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

        if (var_idx == Vars::xvel) {
           ivar   = WRFBdyVars::U; // U
        } else if (var_idx == Vars::yvel) {
           ivar   = WRFBdyVars::V; // V
        } else  if (var_idx == Vars::cons) {
           ivar   = WRFBdyVars::R; // R (we assume that T comes right after R)
        }

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

            // x-faces
            {
                Box bx_xlo(bx);
                bx_xlo.setSmall(0,dom_lo.x); bx_xlo.setBig(0,dom_lo.x);
                bx_xlo.setSmall(1,dom_lo.y); bx_xlo.setBig(1,dom_hi.y);
                bx_xlo.setSmall(2,dom_lo.z); bx_xlo.setBig(2,dom_hi.z);

                Box bx_xhi(bx);
                bx_xhi.setSmall(1,dom_lo.y); bx_xhi.setBig(1,dom_hi.y);
                bx_xhi.setSmall(2,dom_lo.z); bx_xhi.setBig(2,dom_hi.z);

                if (var_idx == Vars::xvel) {
                    bx_xhi.setSmall(0,dom_hi.x+1); bx_xhi.setBig(0,dom_hi.x+1);
                } else {
                    bx_xhi.setSmall(0,dom_hi.x); bx_xhi.setBig(0,dom_hi.x);
                }

                ParallelFor(
                  bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int )
                  {
                      dest_arr(i,j,k,icomp+nv) = oma   * bdatxlo_n  (i,j,k,nv)
                                               + alpha * bdatxlo_np1(i,j,k,nv);
                  },
                  bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int )
                  {
                      dest_arr(i,j,k,icomp+nv) = oma   * bdatxhi_n  (i,j,k,nv)
                                               + alpha * bdatxhi_np1(i,j,k,nv);
                  }
                );
            } // x-faces

            // y-faces
            {
                Box bx_ylo(bx);
                bx_ylo.setSmall(0,dom_lo.x); bx_ylo.setBig(0,dom_hi.x);
                bx_ylo.setSmall(1,dom_lo.y); bx_ylo.setBig(1,dom_lo.y);
                bx_ylo.setSmall(2,dom_lo.z); bx_ylo.setBig(2,dom_hi.z);

                Box bx_yhi(bx);
                bx_yhi.setSmall(0,dom_lo.x); bx_yhi.setBig(0,dom_hi.x);
                bx_yhi.setSmall(2,dom_lo.z); bx_yhi.setBig(2,dom_hi.z);

                if (var_idx == Vars::yvel) {
                    bx_yhi.setSmall(1,dom_hi.y+1); bx_yhi.setBig(1,dom_hi.y+1);
                } else {
                    bx_yhi.setSmall(1,dom_hi.y); bx_yhi.setBig(1,dom_hi.y);
                }

                ParallelFor(
                  bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int )
                  {
                      dest_arr(i,j,k,icomp+nv) = oma   * bdatylo_n  (i,j,k,nv)
                                               + alpha * bdatylo_np1(i,j,k,nv);
                  },
                  bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int )
                  {
                      dest_arr(i,j,k,icomp+nv) = oma   * bdatyhi_n  (i,j,k,nv)
                                               + alpha * bdatyhi_np1(i,j,k,nv);
                  }
                );
            } // y-faces
          } // mfi
        } // nv
    } // var_idx
}
#endif
