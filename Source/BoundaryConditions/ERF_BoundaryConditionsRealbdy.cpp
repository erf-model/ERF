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
                        int ncomp_cons,
                        IntVect ngvect_cons,
                        IntVect ngvect_vels)
{
    int lev = 0;

    // We do not operate on the z ghost cells
    ngvect_cons[2] = 0;
    ngvect_vels[2] = 0;

    // Time interpolation
    Real dT = bdy_time_interval;

    //
    // Note this is because we define "time" to be time since start_bdy_time
    //
    Real time_since_start = time;

    int n_time = static_cast<int>( time_since_start /  dT);
    Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    Real oma   = 1.0 - alpha;

    int n_time_p1 = n_time + 1;
    if ((time == stop_time-start_time) && (alpha==0)) {
        // stop time coincides with final bdy snapshot -- don't try to read in
        // another snapshot
        n_time_p1 = n_time;
    }

    // Data structure for WfromOmega
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom[lev].InvCellSizeArray();

    // Flags for read vars and index mapping
    Vector<int> cons_read = {0, 1, 0, 0,
                             1, 0, 0,
                             0, 0, 0};

    Vector<int> cons_map = {Rho_comp, RealBdyVars::T, RhoKE_comp, RhoScalar_comp,
                            RealBdyVars::QV, RhoQ2_comp, RhoQ3_comp,
                            RhoQ4_comp, RhoQ5_comp, RhoQ6_comp};

    Vector<Vector<int>> is_read;
    is_read.push_back( cons_read );
    is_read.push_back( {1} ); // xvel
    is_read.push_back( {1} ); // yvel
    is_read.push_back( {0} ); // zvel

    Vector<Vector<int>> ind_map;
    ind_map.push_back( cons_map );
    ind_map.push_back( {RealBdyVars::U} ); // xvel
    ind_map.push_back( {RealBdyVars::V} ); // yvel
    ind_map.push_back( {0} );              // zvel NOTE: Loop is forward (need u/v filled first)

    // Nvars to loop over
    Vector<int> comp_var = {ncomp_cons, 1, 1, 1};

    // End of vars loop
    int var_idx_end = (cons_only) ? Vars::cons + 1 : Vars::NumTypes;

    // Extrapolate w from the interior (default) or explicitly set to 0
    bool zero_w = (!cons_only && !real_extrap_w);
    if (zero_w) {
        var_idx_end -= 1;

        mfs[Vars::zvel]->setVal(0.0);
    }

    MultiFab& mf_u = *mfs[Vars::xvel];
    MultiFab& mf_v = *mfs[Vars::yvel];

    // Loop over all variable types
    for (int var_idx = Vars::cons; var_idx < var_idx_end; ++var_idx)
    {
        MultiFab& mf = *mfs[var_idx];

        mf.FillBoundary(geom[lev].periodicity());

        //
        // Note that "domain" is mapped onto the type of box the data is in
        //
        Box domain = geom[lev].Domain();
        domain.convert(mf.boxArray().ixType());
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);

        // Offset only applies to cons (we may fill a subset of these vars)
        int offset = (var_idx == Vars::cons) ? icomp_cons : 0;

        // Ghost cells to be filled
        IntVect ng_vect = (var_idx == Vars::cons) ? ngvect_cons : ngvect_vels;

        // Set region width
        int set_width = real_set_width;

        // Loop over each component
        for (int comp_idx(offset); comp_idx < (comp_var[var_idx]+offset); ++comp_idx)
        {

            // Variable can be read from wrf bdy
            //------------------------------------
            if (is_read[var_idx][comp_idx])
            {
                int ivar  = ind_map[var_idx][comp_idx];

                // We have data at fixed time intervals we will call dT
                // Then to interpolate, given time, we can define n = (time/dT)
                // and alpha = (time - n*dT) / dT, then we define the data at time
                // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)
                const auto& bdatxlo_n   = bdy_data_xlo[n_time   ][ivar].const_array();
                const auto& bdatxlo_np1 = bdy_data_xlo[n_time_p1][ivar].const_array();
                const auto& bdatxhi_n   = bdy_data_xhi[n_time   ][ivar].const_array();
                const auto& bdatxhi_np1 = bdy_data_xhi[n_time_p1][ivar].const_array();
                const auto& bdatylo_n   = bdy_data_ylo[n_time   ][ivar].const_array();
                const auto& bdatylo_np1 = bdy_data_ylo[n_time_p1][ivar].const_array();
                const auto& bdatyhi_n   = bdy_data_yhi[n_time   ][ivar].const_array();
                const auto& bdatyhi_np1 = bdy_data_yhi[n_time_p1][ivar].const_array();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Grown tilebox so we fill exterior ghost cells as well
                    Box gbx = mfi.growntilebox(ng_vect);
                    const Array4<Real>& dest_arr = mf.array(mfi);
                    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
                    realbdy_bc_bxs_xy(gbx, domain, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect);

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
                int width = (var_idx == Vars::zvel) ? real_width : 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Grown tilebox so we fill exterior ghost cells as well
                    Box gbx = mfi.growntilebox(ng_vect);
                    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
                    realbdy_bc_bxs_xy(gbx, domain, width+1,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect);

                    // Bounding
                    int i_xlo = bx_xlo.bigEnd(0)   + 1;
                    int i_xhi = bx_xhi.smallEnd(0) - 1;
                    int j_ylo = bx_ylo.bigEnd(1)   + 1;
                    int j_yhi = bx_yhi.smallEnd(1) - 1;

                    // Destination array
                    const Array4<Real>& dest_arr = mf.array(mfi);

                    // Compute w from Omega
                    const Array4<const Real>& z_nd_arr = z_phys_nd[lev]->const_array(mfi);
                    const Array4<const Real>& u_arr    = mf_u.const_array(mfi);
                    const Array4<const Real>& v_arr    = mf_v.const_array(mfi);
                    const Array4<const Real>& mf_ux    = mapfac[lev][MapFacType::u_x]->const_array(mfi);
                    const Array4<const Real>& mf_vy    = mapfac[lev][MapFacType::v_y]->const_array(mfi);

                    // x-faces (includes y ghost cells)
                    ParallelFor(bx_xlo, bx_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y+width);
                            jj = std::min(jj, dom_hi.y-width);
                        dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,v_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i_xlo,jj,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y+width);
                            jj = std::min(jj, dom_hi.y-width);
                        dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,v_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i_xhi,jj,k,comp_idx);
                    });

                    // y-faces (does not include x ghost cells)
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,u_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i,j_ylo,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                      dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,v_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i,j_yhi,k,comp_idx);
                    });
                } // mfi
            } // is_read
        } // comp
    } // var
}

void
ERF::fill_from_realbdy_upwind (const Vector<MultiFab*>& mfs,
                               const Real time,
                               bool cons_only,
                               int icomp_cons,
                               int ncomp_cons,
                               IntVect ngvect_cons,
                               IntVect ngvect_vels)
{
    int lev = 0;

    // We do not operate on the z ghost cells
    ngvect_cons[2] = 0;
    ngvect_vels[2] = 0;

    // Time interpolation
    Real dT = bdy_time_interval;

    //
    // Note this is because we define "time" to be time since start_bdy_time
    //
    Real time_since_start = time;

    int n_time = static_cast<int>( time_since_start /  dT);
    Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    Real oma   = 1.0 - alpha;

    int n_time_p1 = n_time + 1;
    if ((time == stop_time-start_time) && (alpha==0)) {
        // stop time coincides with final bdy snapshot -- don't try to read in
        // another snapshot
        n_time_p1 = n_time;
    }

    // Data structure for WfromOmega
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom[lev].InvCellSizeArray();

    // Map loop index to variable index
    // NOTE: Loop backwards and do v/u first for WfromOmega
    Vector<int> var_idx_map = {Vars::cons, Vars::zvel, Vars::xvel, Vars::yvel};

    // Variables read from bdy (indexed with var_idx)
    Vector<Vector<int>> is_read;
    is_read.push_back( {0, 1, 0, 0,
                        1, 0, 0,
                        0, 0, 0} );
    is_read.push_back( {1} ); // xvel
    is_read.push_back( {1} ); // yvel
    is_read.push_back( {0} ); // zvel

    // Variable component (indexed with var_idx)
    Vector<Vector<int>> ind_map;
    ind_map.push_back( {Rho_comp, RealBdyVars::T, RhoKE_comp, RhoScalar_comp,
                        RealBdyVars::QV, RhoQ2_comp, RhoQ3_comp,
                        RhoQ4_comp, RhoQ5_comp, RhoQ6_comp} );
    ind_map.push_back( {RealBdyVars::U} ); // xvel
    ind_map.push_back( {RealBdyVars::V} ); // yvel
    ind_map.push_back( {0} );              // zvel

    // Nvars to loop over
    Vector<int> comp_var = {ncomp_cons, 1, 1, 1};

    // End of vars loop
    int var_idx_end = (cons_only) ? Vars::cons + 1 : Vars::NumTypes;

    // Extrapolate w from the interior (default) or explicitly set to 0
    bool zero_w = (!cons_only && !real_extrap_w);
    if (zero_w) {
        var_idx_end -= 1;
        mfs[Vars::zvel]->setVal(0.0);
    }

    // Loop over all variable types -- note we intentionally do
    //      the velocities before the scalars since we may use
    //      the normal velocity on each face for upwinding
    for (int var_lp_idx = var_idx_end-1; var_lp_idx >= Vars::cons; --var_lp_idx)
    {
        int var_idx  = var_idx_map[var_lp_idx];
        MultiFab& mf = *mfs[var_idx];

        mf.FillBoundary(geom[lev].periodicity());

        // CC domain
        Box domain = geom[lev].Domain();
        const auto& dom_cc_lo = lbound(domain);
        const auto& dom_cc_hi = ubound(domain);

        // Map domain onto index type of the data
        domain.convert(mf.boxArray().ixType());
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);
        IntVect iv = domain.type();

        // Offset only applies to cons (we may fill a subset of these vars)
        int offset = (var_idx == Vars::cons) ? icomp_cons : 0;

        // Ghost cells to be filled
        IntVect ng_vect = (var_idx == Vars::cons) ? ngvect_cons : ngvect_vels;

        // Set region width
        int set_width = real_set_width;

        // Loop over each component
        for (int comp_idx(offset); comp_idx < (comp_var[var_idx]+offset); ++comp_idx)
        {
            // Variable can be read from wrf bdy
            //------------------------------------
            if (is_read[var_idx][comp_idx])
            {
                int ivar = ind_map[var_idx][comp_idx];

                // We have data at fixed time intervals we will call dT
                // Then to interpolate, given time, we can define n = (time/dT)
                // and alpha = (time - n*dT) / dT, then we define the data at time
                // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)
                const auto& bdatxlo_n   = bdy_data_xlo[n_time   ][ivar].const_array();
                const auto& bdatxlo_np1 = bdy_data_xlo[n_time_p1][ivar].const_array();
                const auto& bdatxhi_n   = bdy_data_xhi[n_time   ][ivar].const_array();
                const auto& bdatxhi_np1 = bdy_data_xhi[n_time_p1][ivar].const_array();
                const auto& bdatylo_n   = bdy_data_ylo[n_time   ][ivar].const_array();
                const auto& bdatylo_np1 = bdy_data_ylo[n_time_p1][ivar].const_array();
                const auto& bdatyhi_n   = bdy_data_yhi[n_time   ][ivar].const_array();
                const auto& bdatyhi_np1 = bdy_data_yhi[n_time_p1][ivar].const_array();

                // U/V for upwind masking
                const auto& bdatxlo_n_m   = bdy_data_xlo[n_time   ][RealBdyVars::U].const_array();
                const auto& bdatxlo_np1_m = bdy_data_xlo[n_time_p1][RealBdyVars::U].const_array();
                const auto& bdatxhi_n_m   = bdy_data_xhi[n_time   ][RealBdyVars::U].const_array();
                const auto& bdatxhi_np1_m = bdy_data_xhi[n_time_p1][RealBdyVars::U].const_array();
                const auto& bdatylo_n_m   = bdy_data_ylo[n_time   ][RealBdyVars::V].const_array();
                const auto& bdatylo_np1_m = bdy_data_ylo[n_time_p1][RealBdyVars::V].const_array();
                const auto& bdatyhi_n_m   = bdy_data_yhi[n_time   ][RealBdyVars::V].const_array();
                const auto& bdatyhi_np1_m = bdy_data_yhi[n_time_p1][RealBdyVars::V].const_array();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Grown tilebox so we fill exterior ghost cells as well
                    Box gbx = mfi.growntilebox(ng_vect);
                    const Array4<Real>& dest_arr = mf.array(mfi);
                    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
                    realbdy_bc_bxs_xy(gbx, domain, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect);

                    // NOTE: Xlo/hi boxes own corner cells (Ylo/hi)
                    ParallelFor(bx_xlo, bx_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        // Limit for BDY FAB data
                        int ii = std::max(i , dom_lo.x);
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);

                        // Compute bdy velocities
                        int jju = std::min(std::max(j, dom_cc_lo.y), dom_cc_hi.y);
                        int iiv = std::min(std::max(i, dom_cc_lo.x), dom_cc_hi.x);
                        Real u_bdy    = oma   * bdatxlo_n_m  (dom_cc_lo.x,jju,k,0)
                                      + alpha * bdatxlo_np1_m(dom_cc_lo.x,jju,k,0);
                        Real v_bdy_lo = oma   * bdatylo_n_m  (iiv,dom_cc_lo.y,k,0)
                                      + alpha * bdatylo_np1_m(iiv,dom_cc_lo.y,k,0);
                        Real v_bdy_hi = oma   * bdatyhi_n_m  (iiv,dom_cc_hi.y+1,k,0)
                                      + alpha * bdatyhi_np1_m(iiv,dom_cc_hi.y+1,k,0);

                        if ( (u_bdy >= 0.0) ||
                             ((jj == dom_cc_lo.y      ) && (v_bdy_lo >= 0.0)) ||
                             ((jj == dom_cc_hi.y+iv[1]) && (v_bdy_hi <= 0.0)) ) {
                            dest_arr(i,j,k,comp_idx) = oma   * bdatxlo_n  (ii,jj,k,0)
                                                     + alpha * bdatxlo_np1(ii,jj,k,0);
                            if (var_idx == Vars::cons) {
                                dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                            }
                        } else {
                            dest_arr(i,j,k,comp_idx) = dest_arr(dom_lo.x,jj,k,comp_idx);
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        // Limit for BDY FAB data
                        int ii = std::min(i , dom_hi.x);
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);

                        // Compute bdy velocities
                        int jju = std::min(std::max(j, dom_cc_lo.y), dom_cc_hi.y);
                        int iiv = std::min(std::max(i, dom_cc_lo.x), dom_cc_hi.x);
                        Real u_bdy    = oma   * bdatxhi_n_m  (dom_cc_hi.x+1,jju,k,0)
                                      + alpha * bdatxhi_np1_m(dom_cc_hi.x+1,jju,k,0);
                        Real v_bdy_lo = oma   * bdatylo_n_m  (iiv,dom_cc_lo.y,k,0)
                                      + alpha * bdatylo_np1_m(iiv,dom_cc_lo.y,k,0);
                        Real v_bdy_hi = oma   * bdatyhi_n_m  (iiv,dom_cc_hi.y+1,k,0)
                                      + alpha * bdatyhi_np1_m(iiv,dom_cc_hi.y+1,k,0);

                        if ( (u_bdy <= 0.0) ||
                             ((jj == dom_cc_lo.y      ) && (v_bdy_lo >= 0.0)) ||
                             ((jj == dom_cc_hi.y+iv[1]) && (v_bdy_hi <= 0.0)) ) {
                            dest_arr(i,j,k,comp_idx) = oma   * bdatxhi_n  (ii,jj,k,0)
                                                     + alpha * bdatxhi_np1(ii,jj,k,0);
                            if (var_idx == Vars::cons) {
                                dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                            }
                        } else {
                            dest_arr(i,j,k,comp_idx) = dest_arr(dom_hi.x,jj,k,comp_idx);
                        }
                    });

                    // NOTE: Ylo/hi boxes do not own corner cells
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        // Limit for BDY FAB data
                        int jj = std::max(j, dom_lo.y);

                        // Compute bdy velocities
                        int iiv = std::min(std::max(i, dom_cc_lo.x), dom_cc_hi.x);
                        Real v_bdy = oma   * bdatylo_n_m  (iiv,dom_cc_lo.y,k,0)
                                   + alpha * bdatylo_np1_m(iiv,dom_cc_lo.y,k,0);

                        if (v_bdy >= 0.0) {
                            dest_arr(i,j,k,comp_idx) = oma   * bdatylo_n  (i,jj,k,0)
                                                     + alpha * bdatylo_np1(i,jj,k,0);
                            if (var_idx == Vars::cons) {
                                dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                            }
                        } else {
                            dest_arr(i,j,k,comp_idx) = dest_arr(i,dom_lo.y,k,comp_idx);
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        // Limit for BDY FAB data
                        int jj = std::min(j, dom_hi.y);

                        // Compute bdy velocities
                        int iiv = std::min(std::max(i, dom_cc_lo.x), dom_cc_hi.x);
                        Real v_bdy = oma   * bdatyhi_n_m  (iiv,dom_cc_hi.y+1,k,0)
                                   + alpha * bdatyhi_np1_m(iiv,dom_cc_hi.y+1,k,0);

                        if (v_bdy <= 0.0) {
                            dest_arr(i,j,k,comp_idx) = oma   * bdatyhi_n  (i,jj,k,0)
                                                     + alpha * bdatyhi_np1(i,jj,k,0);
                            if (var_idx == Vars::cons) {
                                dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                            }
                        } else {
                            dest_arr(i,j,k,comp_idx) = dest_arr(i,dom_hi.y,k,comp_idx);
                        }
                    });
                } // mfi

            // Variable not read from wrf bdy
            //------------------------------------
            } else {
                MultiFab& mf_u = *mfs[Vars::xvel];
                MultiFab& mf_v = *mfs[Vars::yvel];
                int width = (var_idx == Vars::zvel) ? real_width : 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Grown tilebox so we fill exterior ghost cells as well
                    Box gbx = mfi.growntilebox(ng_vect);
                    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
                    realbdy_bc_bxs_xy(gbx, domain, width+1,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect);

                    // Bounding
                    int i_xlo = bx_xlo.bigEnd(0)   + 1;
                    int i_xhi = bx_xhi.smallEnd(0) - 1;
                    int j_ylo = bx_ylo.bigEnd(1)   + 1;
                    int j_yhi = bx_yhi.smallEnd(1) - 1;

                    // Destination array
                    const Array4<Real>& dest_arr = mf.array(mfi);

                    // Compute w from Omega
                    const Array4<const Real>& z_nd_arr = z_phys_nd[lev]->const_array(mfi);
                    const Array4<const Real>& u_arr    = mf_u.const_array(mfi);
                    const Array4<const Real>& v_arr    = mf_v.const_array(mfi);
                    const Array4<const Real>& mf_ux    = mapfac[lev][MapFacType::u_x]->const_array(mfi);
                    const Array4<const Real>& mf_vy    = mapfac[lev][MapFacType::v_y]->const_array(mfi);

                    // x-faces (includes y ghost cells)
                    ParallelFor(bx_xlo, bx_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y+width);
                            jj = std::min(jj, dom_hi.y-width);
                        dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,v_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i_xlo,jj,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y+width);
                            jj = std::min(jj, dom_hi.y-width);
                        dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,v_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i_xhi,jj,k,comp_idx);
                    });

                    // y-faces (does not include x ghost cells)
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,v_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i,j_ylo,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,comp_idx) = (false) ? //(var_idx == Vars::zvel) ?
                            WFromOmega(i,j,k,0.0,u_arr,v_arr,mf_ux,mf_vy,z_nd_arr,dxInv) :
                            dest_arr(i,j_yhi,k,comp_idx);
                    });
                } // mfi
            } // is_read
        } // comp
    } // var
}
#endif
