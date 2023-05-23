#include "ERF.H"
#include "Utils.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
/*
 * Impose boundary conditions using data read in from wrfbdy files
 *
 * @param[out] mfs  Vector of MultiFabs to be filled
 * @param[in] time  time at which the data should be filled
 */

void
ERF::fill_from_wrfbdy (const Vector<MultiFab*>& mfs, const Real time)
{
    int lev = 0;

    // Time interpolation
    Real dT = bdy_time_interval;
    int n_time = static_cast<int>(time / dT);
    amrex::Real alpha = (time - n_time * dT) / dT;
    amrex::Real oma   = 1.0 - alpha;

    // Flags for read vars and index mapping
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
    Vector<int> cons_read = {1, 1, 0, 0, 0, 1, 0};
    Vector<int> cons_map = {WRFBdyVars::R, WRFBdyVars::T, 0, 0, 0, WRFBdyVars::QV, 0};
# else
    Vector<int> cons_read = {1, 1, 0, 0, 0}; // R, RT, RKE, RQKE, RS
    Vector<int> cons_map = {WRFBdyVars::R, WRFBdyVars::T, 0, 0, 0}; // R, RT, RKE, RQKE, RS
#endif

    Vector<Vector<int>> is_read;
    is_read.push_back( cons_read );
    is_read.push_back( {1} ); // xvel
    is_read.push_back( {1} ); // yvel
    is_read.push_back( {0} ); // zvel

    Vector<Vector<int>> ind_map;
    ind_map.push_back( cons_map );
    ind_map.push_back( {WRFBdyVars::U} ); // xvel
    ind_map.push_back( {WRFBdyVars::V} ); // yvel
    ind_map.push_back( {0} );             // zvel

    // Nvars to loop over
    Vector<int> comp_var = {NVAR, 1, 1, 1};

    // Loop over all variable types
    for (int var_idx = Vars::cons; var_idx < Vars::NumTypes; ++var_idx)
    {
        MultiFab& mf = *mfs[var_idx];

        //
        // Note that "domain" is mapped onto the type of box the data is in
        //
        Box domain = geom[lev].Domain();
        domain.convert(mf.boxArray().ixType());
        const auto& dom_lo = amrex::lbound(domain);
        const auto& dom_hi = amrex::ubound(domain);

        // Loop over each component
        for (int comp_idx(0); comp_idx < comp_var[var_idx]; ++comp_idx)
        {
            int width;

            // Variable can be read from wrf bdy
            //------------------------------------
            if (is_read[var_idx][comp_idx])
            {
                width = wrfbdy_set_width;
                int ivar  = ind_map[var_idx][comp_idx];
                IntVect ng_vect = mf.nGrowVect(); ng_vect[2] = 0;

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
                    ParallelFor(bx_xlo, bx_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int ii = std::max(i , dom_lo.x);
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);
                        dest_arr(i,j,k,comp_idx) = oma   * bdatxlo_n  (ii,jj,k,0)
                                                 + alpha * bdatxlo_np1(ii,jj,k,0);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int ii = std::min(i , dom_hi.x);
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);
                        dest_arr(i,j,k,comp_idx) = oma   * bdatxhi_n  (ii,jj,k,0)
                                                 + alpha * bdatxhi_np1(ii,jj,k,0);
                    });

                    // y-faces (do not include exterior x ghost cells)
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y);
                        dest_arr(i,j,k,comp_idx) = oma   * bdatylo_n  (i,jj,k,0)
                                                 + alpha * bdatylo_np1(i,jj,k,0);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::min(j , dom_hi.y);
                        dest_arr(i,j,k,comp_idx) = oma   * bdatyhi_n  (i,jj,k,0)
                                                 + alpha * bdatyhi_np1(i,jj,k,0);
                    });
                } // mfi

            // Variable not read from wrf bdy
            //------------------------------------
            } else {
                width = wrfbdy_width - 1;
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
                    ParallelFor(bx_xlo, bx_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);
                            dest_arr(i,j,k,comp_idx) = dest_arr(dom_lo.x+width,jj,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);
                            dest_arr(i,j,k,comp_idx) = dest_arr(dom_hi.x-width,jj,k,comp_idx);
                    });

                    // y-faces (does not include x ghost cells)
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,comp_idx) = dest_arr(i,dom_lo.y+width,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,comp_idx) = dest_arr(i,dom_hi.y-width,k,comp_idx);
                    });
                } // mfi
            } // is_read
        } // comp
    } // var
}
#endif
