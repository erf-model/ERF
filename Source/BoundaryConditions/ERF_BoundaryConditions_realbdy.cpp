#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
/*
 * Impose boundary conditions using data read in from wrfbdy files
 *
 * @param[out] mfs  Vector of MultiFabs to be filled
 * @param[in] time  time at which the data should be filled
 */

void
ERF::fill_from_realbdy (const Vector<MultiFab*>& mfs,
                        const Real time,
                        bool cons_only,
                        int icomp_cons,
                        int ncomp_cons)
{
    int lev = 0;

    // Time interpolation
    Real dT = bdy_time_interval;
    Real time_since_start = time - start_bdy_time;
    int n_time = static_cast<int>( time_since_start /  dT);
    Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    Real oma   = 1.0 - alpha;

    // Flags for read vars and index mapping
    Vector<int> cons_read = {0, 1, 0,
                             0, 0, 1,
                             0, 0, 0,
                             0, 0};
    Vector<int> cons_map = {Rho_comp,    RealBdyVars::T, RhoKE_comp,
                            RhoQKE_comp, RhoScalar_comp, RealBdyVars::QV,
                            RhoQ2_comp,  RhoQ3_comp,     RhoQ4_comp,
                            RhoQ5_comp,  RhoQ6_comp};

    Vector<Vector<int>> is_read;
    is_read.push_back( cons_read );
    is_read.push_back( {1} ); // xvel
    is_read.push_back( {1} ); // yvel
    is_read.push_back( {0} ); // zvel

    Vector<Vector<int>> ind_map;
    ind_map.push_back( cons_map );
    ind_map.push_back( {RealBdyVars::U} ); // xvel
    ind_map.push_back( {RealBdyVars::V} ); // yvel
    ind_map.push_back( {0} );             // zvel

    // Nvars to loop over
    Vector<int> comp_var = {ncomp_cons, 1, 1, 1};

    // End of vars loop
    int var_idx_end = (cons_only) ? Vars::cons + 1 : Vars::NumTypes;

    // Loop over all variable types
    for (int var_idx = Vars::cons; var_idx < var_idx_end; ++var_idx)
    {
        MultiFab& mf = *mfs[var_idx];

        //
        // Note that "domain" is mapped onto the type of box the data is in
        //
        Box domain = geom[lev].Domain();
        domain.convert(mf.boxArray().ixType());
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);

        // Offset only applies to cons (we may fill a subset of these vars)
        int offset = (var_idx == Vars::cons) ? icomp_cons : 0;

        // Loop over each component
        for (int comp_idx(offset); comp_idx < (comp_var[var_idx]+offset); ++comp_idx)
        {
            int width = real_set_width;

            // Variable can be read from wrf bdy
            //------------------------------------
            if (is_read[var_idx][comp_idx])
            {
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
                        if (var_idx == Vars::cons) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int ii = std::min(i , dom_hi.x);
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);
                        dest_arr(i,j,k,comp_idx) = oma   * bdatxhi_n  (ii,jj,k,0)
                                                 + alpha * bdatxhi_np1(ii,jj,k,0);
                        if (var_idx == Vars::cons) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    });

                    // y-faces (do not include exterior x ghost cells)
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y);
                        dest_arr(i,j,k,comp_idx) = oma   * bdatylo_n  (i,jj,k,0)
                                                 + alpha * bdatylo_np1(i,jj,k,0);
                        if (var_idx == Vars::cons) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::min(j , dom_hi.y);
                        dest_arr(i,j,k,comp_idx) = oma   * bdatyhi_n  (i,jj,k,0)
                                                 + alpha * bdatyhi_np1(i,jj,k,0);
                        if (var_idx == Vars::cons) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    });
                } // mfi

            // Variable not read from wrf bdy
            //------------------------------------
            } else {
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
                        int jj = std::max(j , dom_lo.y+width);
                            jj = std::min(jj, dom_hi.y-width);
                        dest_arr(i,j,k,comp_idx) = dest_arr(dom_lo.x+width,jj,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y+width);
                            jj = std::min(jj, dom_hi.y-width);
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
