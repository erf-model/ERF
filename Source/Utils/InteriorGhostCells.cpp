#include <Utils.H>

using namespace amrex;

PhysBCFunctNoOp void_bc;

/**
 * Get the boxes for looping over interior/exterior ghost cells
 * for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post.
 *
 * @param[in] bx box to intersect with 4 halo regions
 * @param[in] domain box of the whole domain
 * @param[in] width number of cells in (relaxation+specified) zone
 * @param[in] set_width number of cells in (specified) zone
 * @param[out] bx_xlo halo box at x_lo boundary
 * @param[out] bx_xhi halo box at x_hi boundary
 * @param[out] bx_ylo halo box at y_lo boundary
 * @param[out] bx_yhi halo box at y_hi boundary
 * @param[in] ng_vect number of ghost cells in each direction
 * @param[in] get_int_ng flag to get ghost cells inside the domain
 */
void
compute_interior_ghost_bxs_xy (const Box& bx,
                               const Box& domain,
                               const int& width,
                               const int& set_width,
                               Box& bx_xlo,
                               Box& bx_xhi,
                               Box& bx_ylo,
                               Box& bx_yhi,
                               const IntVect& ng_vect,
                               const bool get_int_ng)
{
    AMREX_ALWAYS_ASSERT(bx.ixType() == domain.ixType());

    //
    // NOTE: X-face boxes take ownership of the overlapping region.
    //       With exterior ghost cells (ng_vect != 0), the x-face
    //       boxes will have exterior ghost cells in both x & y.
    //

    // Domain bounds without ghost cells
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // The four boxes surrounding the domain
    Box gdom_xlo(domain); Box gdom_xhi(domain);
    Box gdom_ylo(domain); Box gdom_yhi(domain);

    // Trim the boxes to only include internal ghost cells
    gdom_xlo.setBig(0,dom_lo.x+width-1); gdom_xhi.setSmall(0,dom_hi.x-width+1);
    gdom_ylo.setBig(1,dom_lo.y+width-1); gdom_yhi.setSmall(1,dom_hi.y-width+1);

    // Remove overlapping corners from y-face boxes
    gdom_ylo.setSmall(0,gdom_xlo.bigEnd(0)+1); gdom_ylo.setBig(0,gdom_xhi.smallEnd(0)-1);
    gdom_yhi.setSmall(0,gdom_xlo.bigEnd(0)+1); gdom_yhi.setBig(0,gdom_xhi.smallEnd(0)-1);

    //==================================================================
    // NOTE: 4 boxes now encompass interior cells with thickness 'width'
    //==================================================================

    gdom_xlo.growLo(0,-set_width); gdom_xhi.growHi(0,-set_width);
    gdom_xlo.grow  (1,-set_width); gdom_xhi.grow  (1,-set_width);
    gdom_ylo.growLo(1,-set_width); gdom_yhi.growHi(1,-set_width);

    //==================================================================
    // NOTE: 4 boxes now exclude the regions being set by bndry
    //==================================================================

    // Grow boxes to get external ghost cells only
    gdom_xlo.growLo(0,ng_vect[0]); gdom_xhi.growHi(0,ng_vect[0]);
    gdom_xlo.grow  (1,ng_vect[1]); gdom_xhi.grow  (1,ng_vect[1]);
    gdom_ylo.growLo(1,ng_vect[1]); gdom_yhi.growHi(1,ng_vect[1]);

    // Grow boxes to get internal ghost cells
    if (get_int_ng) {
        gdom_xlo.growHi(0,ng_vect[0]); gdom_xhi.growLo(0,ng_vect[0]);
        gdom_ylo.grow  (0,ng_vect[0]); gdom_yhi.grow  (0,ng_vect[0]);
        gdom_ylo.growHi(1,ng_vect[1]); gdom_yhi.growLo(1,ng_vect[1]);
    }

    // Populate everything
    bx_xlo = (bx & gdom_xlo);
    bx_xhi = (bx & gdom_xhi);
    bx_ylo = (bx & gdom_ylo);
    bx_yhi = (bx & gdom_yhi);
}

/**
 * Compute the RHS in the relaxation zone
 *
 * @param[in] init_type initialization method for this simulation
 * @param[in] bdy_time_interval time interval between boundary condition time stamps
 * @param[in] time    current time
 * @param[in] delta_t timestep
 * @param[in] width   number of cells in (relaxation+specified) zone
 * @param[in] set_width number of cells in (specified) zone
 * @param[in] geom     container for geometric information
 * @param[out] S_rhs   RHS to be computed here
 * @param[in] S_data   current value of the solution
 * @param[in] bdy_data_xlo boundary data on interior of low x-face
 * @param[in] bdy_data_xhi boundary data on interior of high x-face
 * @param[in] bdy_data_ylo boundary data on interior of low y-face
 * @param[in] bdy_data_yhi boundary data on interior of high y-face
 * @param[in] start_bdy_time time of the first boundary data read in
 */
void
wrfbdy_compute_interior_ghost_rhs (const std::string& init_type,
                                   const Real& bdy_time_interval,
                                   const Real& start_bdy_time,
                                   const Real& time,
                                   const Real& delta_t,
                                   int  width,
                                   int  set_width,
                                   const Geometry& geom,
                                   Vector<MultiFab>& S_rhs,
                                   Vector<MultiFab>& S_old_data,
                                   Vector<MultiFab>& S_cur_data,
                                   Vector<Vector<FArrayBox>>& bdy_data_xlo,
                                   Vector<Vector<FArrayBox>>& bdy_data_xhi,
                                   Vector<Vector<FArrayBox>>& bdy_data_ylo,
                                   Vector<Vector<FArrayBox>>& bdy_data_yhi)
{
    BL_PROFILE_REGION("wrfbdy_compute_interior_ghost_RHS()");

    // NOTE: We pass the full width into this routine.
    //       For relaxation, the last cell is a halo
    //       cell for the Laplacian. We remove that
    //       cell here if it is present.

    // The width to do RHS augmentation
    if (width > set_width+1) width -= 1;

    // Relaxation constants
    Real F1 = 1./(10.*delta_t);
    Real F2 = 1./(50.*delta_t);

    // Time interpolation
    Real dT = bdy_time_interval;
    Real time_since_start = time - start_bdy_time;
    int n_time = static_cast<int>( time_since_start /  dT);
    amrex::Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    amrex::Real oma   = 1.0 - alpha;

    /*
    // UNIT TEST DEBUG
    oma = 1.0; alpha = 0.0;
    */

    // Temporary FABs for storage (owned/filled on all ranks)
    FArrayBox U_xlo, U_xhi, U_ylo, U_yhi;
    FArrayBox V_xlo, V_xhi, V_ylo, V_yhi;
    FArrayBox R_xlo, R_xhi, R_ylo, R_yhi;
    FArrayBox T_xlo, T_xhi, T_ylo, T_yhi;

    // Variable index map (WRFBdyVars -> Vars)
    Vector<int> var_map = {Vars::xvel, Vars::yvel, Vars::cons, Vars::cons};
    Vector<int> ivar_map = {IntVar::xmom, IntVar::ymom, IntVar::cons, IntVar::cons};

    // Variable icomp map
    Vector<int> comp_map = {0, 0, Rho_comp, RhoTheta_comp};

    // Indices
    int BdyEnd = RealBdyVars::NumTypes-1;
    int  ivarU = RealBdyVars::U;
    int  ivarV = RealBdyVars::V;
    int  ivarR = RealBdyVars::R;
    int  ivarT = RealBdyVars::T;


    // Size the FABs
    //==========================================================
    for (int ivar(ivarU); ivar < BdyEnd; ivar++) {
        int var_idx = var_map[ivar];
        Box domain  = geom.Domain();
        domain.convert(S_cur_data[var_idx].boxArray().ixType());

        // Grown domain to get the 4 halo boxes w/ ghost cells
        // NOTE: 2 ghost cells needed for U -> rho*U
        IntVect ng_vect{2,2,0};
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, 0,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        // Size the FABs
        if (ivar  == ivarU) {
            U_xlo.resize(bx_xlo,1,The_Async_Arena()); U_xhi.resize(bx_xhi,1,The_Async_Arena());
            U_ylo.resize(bx_ylo,1,The_Async_Arena()); U_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else if (ivar  == ivarV) {
            V_xlo.resize(bx_xlo,1,The_Async_Arena()); V_xhi.resize(bx_xhi,1,The_Async_Arena());
            V_ylo.resize(bx_ylo,1,The_Async_Arena()); V_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else if (ivar  == ivarR) {
            R_xlo.resize(bx_xlo,1,The_Async_Arena()); R_xhi.resize(bx_xhi,1,The_Async_Arena());
            R_ylo.resize(bx_ylo,1,The_Async_Arena()); R_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else if (ivar  == ivarT){
            T_xlo.resize(bx_xlo,1,The_Async_Arena()); T_xhi.resize(bx_xhi,1,The_Async_Arena());
            T_ylo.resize(bx_ylo,1,The_Async_Arena()); T_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else {
            continue;
        }
    } // ivar

    // Populate FABs from bdy interpolation
    //==========================================================
    for (int ivar(ivarU); ivar < BdyEnd; ivar++) {
        int var_idx = var_map[ivar];
        Box domain  = geom.Domain();
        domain.convert(S_cur_data[var_idx].boxArray().ixType());
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);

        // NOTE: 2 ghost cells needed here. The first
        //       ghost cell is to access the Laplacian
        //       halo cell. The second ghost cell is
        //       for averaging u -> rho*u.
        IntVect ng_vect{2,2,0};
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, 0,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
        Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
        if (ivar  == ivarU) {
            arr_xlo = U_xlo.array(); arr_xhi = U_xhi.array();
            arr_ylo = U_ylo.array(); arr_yhi = U_yhi.array();
        } else if (ivar  == ivarV) {
            arr_xlo = V_xlo.array(); arr_xhi = V_xhi.array();
            arr_ylo = V_ylo.array(); arr_yhi = V_yhi.array();
        } else if (ivar  == ivarR) {
            arr_xlo = R_xlo.array(); arr_xhi = R_xhi.array();
            arr_ylo = R_ylo.array(); arr_yhi = R_yhi.array();
        } else if (ivar  == ivarT){
            arr_xlo = T_xlo.array(); arr_xhi = T_xhi.array();
            arr_ylo = T_ylo.array(); arr_yhi = T_yhi.array();
        } else {
            continue;
        }

        // Boundary data at fixed time intervals
        const auto& bdatxlo_n   = bdy_data_xlo[n_time  ][ivar].const_array();
        const auto& bdatxlo_np1 = bdy_data_xlo[n_time+1][ivar].const_array();
        const auto& bdatxhi_n   = bdy_data_xhi[n_time  ][ivar].const_array();
        const auto& bdatxhi_np1 = bdy_data_xhi[n_time+1][ivar].const_array();
        const auto& bdatylo_n   = bdy_data_ylo[n_time  ][ivar].const_array();
        const auto& bdatylo_np1 = bdy_data_ylo[n_time+1][ivar].const_array();
        const auto& bdatyhi_n   = bdy_data_yhi[n_time  ][ivar].const_array();
        const auto& bdatyhi_np1 = bdy_data_yhi[n_time+1][ivar].const_array();

        // NOTE: width is now one less than the total bndy width
        //       if we have a relaxation zone; so we can access
        //       dom_lo/hi +- width. If we do not have a relax
        //       zone, this offset is set_width - 1.
        int offset = set_width - 1;
        if (width > set_width) offset = width;

        // Populate with interpolation (protect from ghost cells)
        ParallelFor(bx_xlo, bx_xhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_lo.x+offset);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            arr_xlo(i,j,k) = oma   * bdatxlo_n  (ii,jj,k,0)
                           + alpha * bdatxlo_np1(ii,jj,k,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_hi.x-offset);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            arr_xhi(i,j,k) = oma   * bdatxhi_n  (ii,jj,k,0)
                           + alpha * bdatxhi_np1(ii,jj,k,0);
        });

        ParallelFor(bx_ylo, bx_yhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_lo.y+offset);
            arr_ylo(i,j,k) = oma   * bdatylo_n  (ii,jj,k,0)
                           + alpha * bdatylo_np1(ii,jj,k,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_hi.y-offset);
                jj = std::min(jj, dom_hi.y);
            arr_yhi(i,j,k) = oma   * bdatyhi_n  (ii,jj,k,0)
                           + alpha * bdatyhi_np1(ii,jj,k,0);
        });
    } // ivar


    // Velocity to momentum
    //==========================================================
    for (int ivar(ivarU); ivar <= ivarV; ivar++) {
        int ivar_idx = ivar_map[ivar];
        Box domain   = geom.Domain();
        domain.convert(S_cur_data[ivar_idx].boxArray().ixType());

        // NOTE: 1 ghost cell needed here. This first
        //       ghost cell is to access the Laplacian
        //       halo cell. We will touch the second ghost
        //       cell when averaging u -> rho*u.
        IntVect ng_vect{1,1,0};
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, 0,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        Array4<Real> rarr_xlo = R_xlo.array();  Array4<Real> rarr_xhi = R_xhi.array();
        Array4<Real> rarr_ylo = R_ylo.array();  Array4<Real> rarr_yhi = R_yhi.array();

        Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
        Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
        if (ivar  == ivarU) {
            arr_xlo = U_xlo.array(); arr_xhi = U_xhi.array();
            arr_ylo = U_ylo.array(); arr_yhi = U_yhi.array();

            ParallelFor(bx_xlo, bx_xhi,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xlo(i-1,j,k) + rarr_xlo(i,j,k) );
                arr_xlo(i,j,k) *= rho_interp;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xhi(i-1,j,k) + rarr_xhi(i,j,k) );
                arr_xhi(i,j,k) *= rho_interp;
            });

            ParallelFor(bx_ylo, bx_yhi,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_ylo(i-1,j,k) + rarr_ylo(i,j,k) );
                arr_ylo(i,j,k) *= rho_interp;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_yhi(i-1,j,k) + rarr_yhi(i,j,k) );
                arr_yhi(i,j,k) *= rho_interp;
            });
        } else {
            arr_xlo = V_xlo.array(); arr_xhi = V_xhi.array();
            arr_ylo = V_ylo.array(); arr_yhi = V_yhi.array();

            ParallelFor(bx_xlo, bx_xhi,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xlo(i,j-1,k) + rarr_xlo(i,j,k) );
                arr_xlo(i,j,k) *= rho_interp;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xhi(i,j-1,k) + rarr_xhi(i,j,k) );
                arr_xhi(i,j,k) *= rho_interp;
            });

            ParallelFor(bx_ylo, bx_yhi,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_ylo(i,j-1,k) + rarr_ylo(i,j,k) );
                arr_ylo(i,j,k) *= rho_interp;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_yhi(i,j-1,k) + rarr_yhi(i,j,k) );
                arr_yhi(i,j,k) *= rho_interp;
            });
        }
    } // ivar


    // Compute RHS in specified region
    //==========================================================
    if (set_width > 0 ) {
        for (int ivar(ivarU); ivar < BdyEnd; ivar++) {
            int ivar_idx = ivar_map[ivar];
            int icomp    = comp_map[ivar];

            Box domain = geom.Domain();
            domain.convert(S_old_data[ivar_idx].boxArray().ixType());
            const auto& dom_hi = ubound(domain);
            const auto& dom_lo = lbound(domain);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(S_old_data[ivar_idx],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box tbx = mfi.tilebox();
                Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
                compute_interior_ghost_bxs_xy(tbx, domain, width, 0,
                                              tbx_xlo, tbx_xhi,
                                              tbx_ylo, tbx_yhi);

                Array4<Real> rhs_arr; Array4<Real> data_arr;
                Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
                Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
                if (ivar  == ivarU) {
                    arr_xlo  = U_xlo.array(); arr_xhi = U_xhi.array();
                    arr_ylo  = U_ylo.array(); arr_yhi = U_yhi.array();
                    rhs_arr  = S_rhs[IntVar::xmom].array(mfi);
                    data_arr = S_old_data[IntVar::xmom].array(mfi);
                } else if (ivar  == ivarV) {
                    arr_xlo  = V_xlo.array(); arr_xhi = V_xhi.array();
                    arr_ylo  = V_ylo.array(); arr_yhi = V_yhi.array();
                    rhs_arr  = S_rhs[IntVar::ymom].array(mfi);
                    data_arr = S_old_data[IntVar::ymom].array(mfi);
                } else if (ivar  == ivarR) {
                    arr_xlo  = R_xlo.array(); arr_xhi = R_xhi.array();
                    arr_ylo  = R_ylo.array(); arr_yhi = R_yhi.array();
                    rhs_arr  = S_rhs[IntVar::cons].array(mfi);
                    data_arr = S_old_data[IntVar::cons].array(mfi);
                } else if (ivar  == ivarT){
                    arr_xlo  = T_xlo.array(); arr_xhi = T_xhi.array();
                    arr_ylo  = T_ylo.array(); arr_yhi = T_yhi.array();
                    rhs_arr  = S_rhs[IntVar::cons].array(mfi);
                    data_arr = S_old_data[IntVar::cons].array(mfi);
                } else {
                    continue;
                }

                wrfbdy_set_rhs_in_spec_region(delta_t, icomp, 1,
                                              width, set_width, dom_lo, dom_hi,
                                              tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi,
                                              arr_xlo, arr_xhi, arr_ylo, arr_yhi,
                                              data_arr, rhs_arr);
            } // mfi
        } // ivar
    } // set_width


    // Compute RHS in relaxation region
    //==========================================================
    if (width > set_width) {
        for (int ivar(ivarU); ivar < BdyEnd; ivar++) {
            int ivar_idx = ivar_map[ivar];
            int icomp    = comp_map[ivar];

            Box domain = geom.Domain();
            domain.convert(S_cur_data[ivar_idx].boxArray().ixType());
            const auto& dom_hi = ubound(domain);
            const auto& dom_lo = lbound(domain);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(S_cur_data[ivar_idx],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box tbx = mfi.tilebox();
                Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
                compute_interior_ghost_bxs_xy(tbx, domain, width, set_width,
                                              tbx_xlo, tbx_xhi,
                                              tbx_ylo, tbx_yhi);

                Array4<Real> rhs_arr; Array4<Real> data_arr;
                Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
                Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
                if (ivar  == ivarU) {
                    arr_xlo  = U_xlo.array(); arr_xhi = U_xhi.array();
                    arr_ylo  = U_ylo.array(); arr_yhi = U_yhi.array();
                    rhs_arr  = S_rhs[IntVar::xmom].array(mfi);
                    data_arr = S_cur_data[IntVar::xmom].array(mfi);
                } else if (ivar  == ivarV) {
                    arr_xlo  = V_xlo.array(); arr_xhi = V_xhi.array();
                    arr_ylo  = V_ylo.array(); arr_yhi = V_yhi.array();
                    rhs_arr  = S_rhs[IntVar::ymom].array(mfi);
                    data_arr = S_cur_data[IntVar::ymom].array(mfi);
                } else if (ivar  == ivarR) {
                    arr_xlo  = R_xlo.array(); arr_xhi = R_xhi.array();
                    arr_ylo  = R_ylo.array(); arr_yhi = R_yhi.array();
                    rhs_arr  = S_rhs[IntVar::cons].array(mfi);
                    data_arr = S_cur_data[IntVar::cons].array(mfi);
                } else if (ivar  == ivarT){
                    arr_xlo  = T_xlo.array(); arr_xhi = T_xhi.array();
                    arr_ylo  = T_ylo.array(); arr_yhi = T_yhi.array();
                    rhs_arr  = S_rhs[IntVar::cons].array(mfi);
                    data_arr = S_cur_data[IntVar::cons].array(mfi);
                } else {
                    continue;
                }

                wrfbdy_compute_laplacian_relaxation(icomp, 1,
                                                    width, set_width, dom_lo, dom_hi, F1, F2,
                                                    tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi,
                                                    arr_xlo, arr_xhi, arr_ylo, arr_yhi,
                                                    data_arr, rhs_arr);

                /*
                // UNIT TEST DEBUG
                compute_interior_ghost_bxs_xy(tbx, domain, width+1, 0,
                                              tbx_xlo, tbx_xhi,
                                              tbx_ylo, tbx_yhi);
                ParallelFor(tbx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (arr_xlo(i,j,k) != data_arr(i,j,k,icomp)) {
                        amrex::Print() << "ERROR XLO: " << ivar << ' ' << icomp << ' ' << IntVect(i,j,k) << "\n";
                        exit(0);
                    }
                });
                ParallelFor(tbx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (arr_xhi(i,j,k) != data_arr(i,j,k,icomp)) {
                        amrex::Print() << "ERROR XHI: " << ivar << ' ' << icomp << ' ' << IntVect(i,j,k) << "\n";
                        exit(0);
                    }
                });
                ParallelFor(tbx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (arr_ylo(i,j,k) != data_arr(i,j,k,icomp)) {
                        amrex::Print() << "ERROR YLO: " << ivar << ' ' << icomp << ' ' << IntVect(i,j,k) << "\n";
                        exit(0);
                    }
                });
                ParallelFor(tbx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (arr_yhi(i,j,k) != data_arr(i,j,k,icomp)) {
                        amrex::Print() << "ERROR YHI: " << ivar << ' ' << icomp << ' ' << IntVect(i,j,k) << "\n";
                        exit(0);
                    }
                });
                */
            } // mfi
        } // ivar
    } // width
}

/**
 * Compute the RHS in the fine relaxation zone
 *
 * @param[in]  time      current time
 * @param[in]  delta_t   timestep
 * @param[in]  width     number of cells in (relaxation+specified) zone
 * @param[in]  set_width number of cells in (specified) zone
 * @param[in]  FPr_c     cons fine patch container
 * @param[in]  FPr_u     uvel fine patch container
 * @param[in]  FPr_v     vvel fine patch container
 * @param[in]  FPr_w     wvel fine patch container
 * @param[in]  boxes_at_level boxes at current level
 * @param[in]  domain_bcs_type boundary condition types
 * @param[out] S_rhs     RHS to be computed here
 * @param[in]  S_data    current value of the solution
 */
void
fine_compute_interior_ghost_rhs (const Real& time,
                                 const Real& delta_t,
                                 const int& width,
                                 const int& set_width,
                                 const Geometry& geom,
                                 ERFFillPatcher* FPr_c,
                                 ERFFillPatcher* FPr_u,
                                 ERFFillPatcher* FPr_v,
                                 ERFFillPatcher* FPr_w,
                                 Vector<BCRec>& domain_bcs_type,
                                 Vector<MultiFab>& S_rhs_f,
                                 Vector<MultiFab>& S_data_f)
{
    BL_PROFILE_REGION("fine_compute_interior_ghost_RHS()");

    // Relaxation constants
    Real F1 = 1./(10.*delta_t);
    Real F2 = 1./(50.*delta_t);

    // Vector of MFs to hold data (dm differs w/ fine patch)
    Vector<MultiFab> fmf_p_v;

    // Loop over the variables
    for (int ivar_idx = 0; ivar_idx < IntVar::NumVars; ++ivar_idx)
    {
        // Fine mfs
        MultiFab& fmf = S_data_f[ivar_idx];
        MultiFab& rhs = S_rhs_f [ivar_idx];

        // NOTE: These temporary MFs and copy operations are horrible
        //       for memory usage and efficiency. However, we need to
        //       have access to ghost cells in the cons array to convert
        //       from primitive u/v/w to momentum. Furthermore, the BA
        //       for the fine patches in ERFFillPatcher don't match the
        //       BA for the data/RHS. For this reason, the data is copied
        //       to a vector of MFs (with ghost cells) so the BAs match
        //       the BA of data/RHS and we have access to rho to convert
        //       prim to conserved.

        // Temp MF on box (distribution map differs w/ fine patch)
        int num_var = fmf.nComp();
        fmf_p_v.emplace_back(fmf.boxArray(), fmf.DistributionMap(), num_var, fmf.nGrowVect());
        MultiFab& fmf_p = fmf_p_v[ivar_idx];
        amrex::MultiFab::Copy(fmf_p,fmf, 0, 0, num_var, fmf.nGrowVect());

        // Integer mask MF
        int set_mask_val;
        int relax_mask_val;
        iMultiFab* mask;

        // Fill fine patch on interior halo region
        //==========================================================
        if (ivar_idx == IntVar::cons)
        {
            FPr_c->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_c->GetMask();
            set_mask_val   = FPr_c->GetSetMaskVal();
            relax_mask_val = FPr_c->GetRelaxMaskVal();
        }
        else if (ivar_idx == IntVar::xmom)
        {
            FPr_u->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_u->GetMask();
            set_mask_val   = FPr_u->GetSetMaskVal();
            relax_mask_val = FPr_u->GetRelaxMaskVal();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(fmf_p,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();

                const Array4<Real>& prim_arr = fmf_p.array(mfi);
                const Array4<const Real>& rho_arr  = fmf_p_v[0].const_array(mfi);
                const Array4<const int>&  mask_arr = mask->const_array(mfi);

                ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (mask_arr(i,j,k) == relax_mask_val) {
                        Real rho_interp = 0.5 * ( rho_arr(i-1,j,k) + rho_arr(i,j,k) );
                        prim_arr(i,j,k) *= rho_interp;
                    }
                });
            } // mfi
        }
        else if (ivar_idx == IntVar::ymom)
        {
            FPr_v->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_v->GetMask();
            set_mask_val   = FPr_v->GetSetMaskVal();
            relax_mask_val = FPr_v->GetRelaxMaskVal();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(fmf_p,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();

                const Array4<Real>& prim_arr = fmf_p.array(mfi);
                const Array4<const Real>& rho_arr  = fmf_p_v[0].const_array(mfi);
                const Array4<const int>&  mask_arr = mask->const_array(mfi);

                ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (mask_arr(i,j,k) == relax_mask_val) {
                        Real rho_interp = 0.5 * ( rho_arr(i,j-1,k) + rho_arr(i,j,k) );
                        prim_arr(i,j,k) *= rho_interp;
                    }
                });
            } // mfi
        }
        else if (ivar_idx == IntVar::zmom)
        {
            FPr_w->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_w->GetMask();
            set_mask_val   = FPr_w->GetSetMaskVal();
            relax_mask_val = FPr_w->GetRelaxMaskVal();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(fmf_p,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();

                const Array4<Real>& prim_arr = fmf_p.array(mfi);
                const Array4<const Real>& rho_arr  = fmf_p_v[0].const_array(mfi);
                const Array4<const int>&  mask_arr = mask->const_array(mfi);

                ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (mask_arr(i,j,k) == relax_mask_val) {
                        Real rho_interp = 0.5 * ( rho_arr(i,j,k-1) + rho_arr(i,j,k) );
                        prim_arr(i,j,k) *= rho_interp;
                    }
                });
            } // mfi
        } else {
            amrex::Abort("Dont recognize this variable type in fine_compute_interior_ghost_RHS");
        }


        // Zero RHS in set region
        //==========================================================
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(rhs,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box vbx = mfi.validbox();
            const Array4<Real>& rhs_arr  = rhs.array(mfi);
            const Array4<const int>& mask_arr = mask->const_array(mfi);

            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (mask_arr(i,j,k) == set_mask_val) {
                    rhs_arr(i,j,k) = 0.0;
                }
            });
        } // mfi

        // For Laplacian stencil
        rhs.FillBoundary(geom.periodicity());


        // Compute RHS in relaxation region
        //==========================================================
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(fmf_p,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box vbx = mfi.validbox();
            const Array4<Real>&        rhs_arr = rhs.array(mfi);
            const Array4<const Real>& fine_arr = fmf_p.const_array(mfi);
            const Array4<const Real>& data_arr = fmf.const_array(mfi);
            const Array4<const int>&  mask_arr = mask->const_array(mfi);

            const auto& vbx_lo = lbound(vbx);
            const auto& vbx_hi = ubound(vbx);

            int icomp = 0;

            int Spec_z  = set_width;
            int Relax_z = width - Spec_z;
            amrex::Real num   = amrex::Real(Spec_z + Relax_z);
            amrex::Real denom = amrex::Real(Relax_z - 1);
            ParallelFor(vbx, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
               if (mask_arr(i,j,k) == relax_mask_val) {

                   // Indices
                   Real n_ind;
                   int ii{width-1}; int jj{width-1};
                   bool near_x_lo_wall{false}; bool near_x_hi_wall{false};
                   bool near_y_lo_wall{false}; bool near_y_hi_wall{false};
                   bool mask_x_found{false};   bool mask_y_found{false};

                   // Near x-wall
                   if ((i-vbx_lo.x) < width) {
                       near_x_lo_wall = true;
                       ii = i-vbx_lo.x;
                       if (mask_arr(vbx_lo.x,j,k) == 2) mask_x_found = true;
                   } else if ((vbx_hi.x-i) < width) {
                       near_x_hi_wall = true;
                       ii = vbx_hi.x-i;
                       if (mask_arr(vbx_hi.x,j,k) == 2) mask_x_found = true;
                   }

                   // Near y-wall
                   if ((j-vbx_lo.y) < width) {
                       near_y_lo_wall = true;
                       jj = j-vbx_lo.y;
                       if (mask_arr(i,vbx_lo.y,k) == 2) mask_y_found = true;
                   } else if ((vbx_hi.y-j) < width) {
                       near_y_hi_wall = true;
                       jj = vbx_hi.y-j;
                       if (mask_arr(i,vbx_hi.y,k) == 2) mask_y_found = true;
                   }

                   // Found a nearby masked cell (valid n_ind)
                   if (mask_x_found && mask_y_found) {
                       n_ind = std::min(ii,jj) + 1.0;
                   } else if (mask_x_found) {
                       n_ind = ii + 1.0;
                   } else if (mask_y_found) {
                       n_ind = jj + 1.0;
                   // Pesky corner cell
                   } else {
                       if (near_x_lo_wall || near_x_hi_wall) {
                           Real dj_min{width-1.0};
                           int j_lb = std::max(vbx_lo.y,j-width);
                           int j_ub = std::min(vbx_hi.y,j+width);
                           int li   = (near_x_lo_wall) ? vbx_lo.x : vbx_hi.x;
                           for (int lj(j_lb); lj<=j_ub; ++lj) {
                               if (mask_arr(li,lj,k) == 2) {
                                   mask_y_found = true;
                                   dj_min = std::min(dj_min,(Real) std::abs(lj-j));
                               }
                           }
                           if (mask_y_found) {
                               Real mag = sqrt( amrex::Real(dj_min*dj_min + ii*ii) );
                               n_ind = std::min(mag,width-1.0) + 1.0;
                           } else {
                               amrex::Abort("Mask not found near x wall!");
                           }
                       } else if (near_y_lo_wall || near_y_hi_wall) {
                           Real di_min{width-1.0};
                           int i_lb = std::max(vbx_lo.x,i-width);
                           int i_ub = std::min(vbx_hi.x,i+width);
                           int lj   = (near_y_lo_wall) ? vbx_lo.y : vbx_hi.y;
                           for (int li(i_lb); li<=i_ub; ++li) {
                               if (mask_arr(li,lj,k) == 2) {
                                   mask_x_found = true;
                                   di_min = std::min(di_min,(Real) std::abs(li-i));
                               }
                           }
                           if (mask_x_found) {
                               Real mag = sqrt( amrex::Real(di_min*di_min + jj*jj) );
                               n_ind = std::min(mag,width-1.0) + 1.0;
                           } else {
                               amrex::Abort("Mask not found near y wall!");
                           }
                       } else {
                           amrex::Abort("Relaxation cell must be near a wall!");
                       }
                   }

                   amrex::Real Factor   = (num - n_ind)/denom;
                   amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp) + delta_t*rhs_arr(i  , j  , k  ,n+icomp);
                   amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i+1, j  , k  ,n+icomp);
                   amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i-1, j  , k  ,n+icomp);
                   amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp) + delta_t*rhs_arr(i  , j+1, k  ,n+icomp);
                   amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp) + delta_t*rhs_arr(i  , j-1, k  ,n+icomp);
                   amrex::Real delta    = fine_arr(i  ,j  ,k,n) - d;
                   amrex::Real delta_xp = fine_arr(i+1,j  ,k,n) - d_ip1;
                   amrex::Real delta_xm = fine_arr(i-1,j  ,k,n) - d_im1;
                   amrex::Real delta_yp = fine_arr(i  ,j+1,k,n) - d_jp1;
                   amrex::Real delta_ym = fine_arr(i  ,j-1,k,n) - d_jm1;
                   amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
                   rhs_arr(i,j,k,n) += (F1*delta - F2*Laplacian) * Factor;
               }
            });
        } // mfi
    } // ivar_idx
}

