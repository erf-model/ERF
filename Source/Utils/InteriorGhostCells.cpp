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
compute_interior_ghost_bxs_xy(const Box& bx,
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
wrfbdy_compute_interior_ghost_RHS(const Real& bdy_time_interval,
                                  const Real& start_bdy_time,
                                  const Real& time,
                                  const Real& delta_t,
                                  const int&  width,
                                  const int&  set_width,
                                  const Geometry& geom,
                                  Vector<MultiFab>& S_rhs,
                                  Vector<MultiFab>& S_data,
                                  Vector<Vector<FArrayBox>>& bdy_data_xlo,
                                  Vector<Vector<FArrayBox>>& bdy_data_xhi,
                                  Vector<Vector<FArrayBox>>& bdy_data_ylo,
                                  Vector<Vector<FArrayBox>>& bdy_data_yhi)
{
    BL_PROFILE_REGION("wrfbdy_compute_interior_ghost_RHS()");

    // Relaxation constants
    Real F1 = 1./(10.*delta_t);
    Real F2 = 1./(50.*delta_t);

    // Time interpolation
    Real dT = bdy_time_interval;
    Real time_since_start = (time - start_bdy_time) / 1.e10;
    int n_time = static_cast<int>( time_since_start / dT);
    amrex::Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    amrex::Real oma   = 1.0 - alpha;

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

    // End of loop for WRFBdyVars
    int WRFBdyEnd = WRFBdyVars::NumTypes-3;


    // Size the FABs
    //==========================================================
    for (int ivar(WRFBdyVars::U); ivar < WRFBdyEnd; ivar++)
    {
        int var_idx = var_map[ivar];
        Box domain  = geom.Domain();
        domain.convert(S_data[var_idx].boxArray().ixType());

        // Grown domain to get the 4 halo boxes w/ ghost cells
        // NOTE: 2 ghost cells needed for U -> rho*U
        IntVect ng_vect{2,2,0};
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        if (ivar  == WRFBdyVars::U) {
            U_xlo.resize(bx_xlo,1); U_xhi.resize(bx_xhi,1);
            U_ylo.resize(bx_ylo,1); U_yhi.resize(bx_yhi,1);
        } else if (ivar  == WRFBdyVars::V) {
            V_xlo.resize(bx_xlo,1); V_xhi.resize(bx_xhi,1);
            V_ylo.resize(bx_ylo,1); V_yhi.resize(bx_yhi,1);
        } else if (ivar  == WRFBdyVars::R) {
            R_xlo.resize(bx_xlo,1); R_xhi.resize(bx_xhi,1);
            R_ylo.resize(bx_ylo,1); R_yhi.resize(bx_yhi,1);
        } else if (ivar  == WRFBdyVars::T){
            T_xlo.resize(bx_xlo,1); T_xhi.resize(bx_xhi,1);
            T_ylo.resize(bx_ylo,1); T_yhi.resize(bx_yhi,1);
        } else {
            continue;
        }
    } // ivar

    // Elixir to avoid destruction before kernel end
    Elixir U_xlo_eli = U_xlo.elixir(); Elixir U_xhi_eli = U_xhi.elixir();
    Elixir U_ylo_eli = U_ylo.elixir(); Elixir U_yhi_eli = U_yhi.elixir();

    Elixir V_xlo_eli = V_xlo.elixir(); Elixir V_xhi_eli = V_xhi.elixir();
    Elixir V_ylo_eli = V_ylo.elixir(); Elixir V_yhi_eli = V_yhi.elixir();

    Elixir R_xlo_eli = R_xlo.elixir(); Elixir R_xhi_eli = R_xhi.elixir();
    Elixir R_ylo_eli = R_ylo.elixir(); Elixir R_yhi_eli = R_yhi.elixir();

    Elixir T_xlo_eli = T_xlo.elixir(); Elixir T_xhi_eli = T_xhi.elixir();
    Elixir T_ylo_eli = T_ylo.elixir(); Elixir T_yhi_eli = T_yhi.elixir();


    // Populate FABs from boundary interpolation
    //==========================================================
    for (int ivar(WRFBdyVars::U); ivar < WRFBdyEnd; ivar++)
    {
        int var_idx = var_map[ivar];
        Box domain  = geom.Domain();
        domain.convert(S_data[var_idx].boxArray().ixType());
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);

        // Grown domain to get the 4 halo boxes w/ ghost cells
        // NOTE: 2 ghost cells needed for U -> rho*U
        IntVect ng_vect{2,2,0};
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
        Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
        if (ivar  == WRFBdyVars::U) {
            arr_xlo = U_xlo.array(); arr_xhi = U_xhi.array();
            arr_ylo = U_ylo.array(); arr_yhi = U_yhi.array();
        } else if (ivar  == WRFBdyVars::V) {
            arr_xlo = V_xlo.array(); arr_xhi = V_xhi.array();
            arr_ylo = V_ylo.array(); arr_yhi = V_yhi.array();
        } else if (ivar  == WRFBdyVars::R) {
            arr_xlo = R_xlo.array(); arr_xhi = R_xhi.array();
            arr_ylo = R_ylo.array(); arr_yhi = R_yhi.array();
        } else if (ivar  == WRFBdyVars::T){
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

        // Populate with interpolation (protect from ghost cells)
        amrex::ParallelFor(bx_xlo, bx_xhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_lo.x+width-1);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            arr_xlo(i,j,k) = oma   * bdatxlo_n  (ii,jj,k,0)
                           + alpha * bdatxlo_np1(ii,jj,k,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_hi.x-width+1);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            arr_xhi(i,j,k) = oma   * bdatxhi_n  (ii,jj,k,0)
                           + alpha * bdatxhi_np1(ii,jj,k,0);
        });

        amrex::ParallelFor(bx_ylo, bx_yhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_lo.y+width-1);
            arr_ylo(i,j,k) = oma   * bdatylo_n  (ii,jj,k,0)
                           + alpha * bdatylo_np1(ii,jj,k,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_hi.y-width+1);
                jj = std::min(jj, dom_hi.y);
            arr_yhi(i,j,k) = oma   * bdatyhi_n  (ii,jj,k,0)
                           + alpha * bdatyhi_np1(ii,jj,k,0);
        });
    } // ivar


    // Velocity to momentum
    //==========================================================
    for (int ivar(WRFBdyVars::U); ivar <= WRFBdyVars::V; ivar++)
    {
        int ivar_idx = ivar_map[ivar];
        Box domain   = geom.Domain();
        domain.convert(S_data[ivar_idx].boxArray().ixType());

        // Grown domain to get the 4 halo boxes w/ ghost cells
        // NOTE: 1 ghost cell needed for Laplacian
        IntVect ng_vect{1,1,0};
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        Array4<Real> rarr_xlo = R_xlo.array();  Array4<Real> rarr_xhi = R_xhi.array();
        Array4<Real> rarr_ylo = R_ylo.array();  Array4<Real> rarr_yhi = R_yhi.array();

        Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
        Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
        if (ivar  == WRFBdyVars::U) {
            arr_xlo = U_xlo.array(); arr_xhi = U_xhi.array();
            arr_ylo = U_ylo.array(); arr_yhi = U_yhi.array();

            amrex::ParallelFor(bx_xlo, bx_xhi,
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

            amrex::ParallelFor(bx_ylo, bx_yhi,
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

            amrex::ParallelFor(bx_xlo, bx_xhi,
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

            amrex::ParallelFor(bx_ylo, bx_yhi,
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


    // Zero RHS in set region
    //==========================================================
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_rhs[IntVar::cons],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<Real> rhs_cons = S_rhs[IntVar::cons].array(mfi);
        const Array4<Real> rhs_xmom = S_rhs[IntVar::xmom].array(mfi);
        const Array4<Real> rhs_ymom = S_rhs[IntVar::ymom].array(mfi);
        const Array4<Real> rhs_zmom = S_rhs[IntVar::zmom].array(mfi);
        {
            Box tbx = mfi.tilebox();
            Box domain = geom.Domain();
            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            compute_interior_ghost_bxs_xy(tbx, domain, set_width, 0,
                                          bx_xlo, bx_xhi,
                                          bx_ylo, bx_yhi);
            zero_RHS_in_set_region(Rho_comp, 2, bx_xlo, bx_xhi, bx_ylo, bx_yhi, rhs_cons);
        }

        {
            Box tbx = mfi.tilebox(IntVect(1,0,0));
            Box domain = geom.Domain();
            domain.convert(S_rhs[IntVar::xmom].boxArray().ixType());
            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            compute_interior_ghost_bxs_xy(tbx, domain, set_width, 0,
                                          bx_xlo, bx_xhi,
                                          bx_ylo, bx_yhi);
            zero_RHS_in_set_region(0, 1, bx_xlo, bx_xhi, bx_ylo, bx_yhi, rhs_xmom);
        }

        {
            Box tbx = mfi.tilebox(IntVect(0,1,0));
            Box domain = geom.Domain();
            domain.convert(S_rhs[IntVar::ymom].boxArray().ixType());
            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            compute_interior_ghost_bxs_xy(tbx, domain, set_width, 0,
                                          bx_xlo, bx_xhi,
                                          bx_ylo, bx_yhi);
            zero_RHS_in_set_region(0, 1, bx_xlo, bx_xhi, bx_ylo, bx_yhi, rhs_ymom);
        }

        {
            Box tbx = mfi.tilebox(IntVect(0,0,1));
            Box domain = geom.Domain();
            domain.convert(S_rhs[IntVar::zmom].boxArray().ixType());
            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            compute_interior_ghost_bxs_xy(tbx, domain, set_width, 0,
                                          bx_xlo, bx_xhi,
                                          bx_ylo, bx_yhi);
            zero_RHS_in_set_region(0, 1, bx_xlo, bx_xhi, bx_ylo, bx_yhi, rhs_zmom);
        }
    } // mfi


    // Compute RHS in relaxation region
    //==========================================================
    for (int ivar(WRFBdyVars::U); ivar < WRFBdyEnd; ivar++)
    {
        int ivar_idx = ivar_map[ivar];
        int icomp    = comp_map[ivar];

        Box domain = geom.Domain();
        domain.convert(S_data[ivar_idx].boxArray().ixType());
        const auto& dom_hi = ubound(domain);
        const auto& dom_lo = lbound(domain);

        // For Laplacian stencil
        S_rhs[ivar_idx].FillBoundary(geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(S_data[ivar_idx],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();
            Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
            compute_interior_ghost_bxs_xy(tbx, domain, width, 0,
                                          tbx_xlo, tbx_xhi,
                                          tbx_ylo, tbx_yhi);

            Array4<Real> rhs_arr; Array4<Real> data_arr;
            Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
            Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
            if (ivar  == WRFBdyVars::U) {
                arr_xlo  = U_xlo.array(); arr_xhi = U_xhi.array();
                arr_ylo  = U_ylo.array(); arr_yhi = U_yhi.array();
                rhs_arr  = S_rhs[IntVar::xmom].array(mfi);
                data_arr = S_data[IntVar::xmom].array(mfi);
            } else if (ivar  == WRFBdyVars::V) {
                arr_xlo  = V_xlo.array(); arr_xhi = V_xhi.array();
                arr_ylo  = V_ylo.array(); arr_yhi = V_yhi.array();
                rhs_arr  = S_rhs[IntVar::ymom].array(mfi);
                data_arr = S_data[IntVar::ymom].array(mfi);
            } else if (ivar  == WRFBdyVars::R) {
                arr_xlo  = R_xlo.array(); arr_xhi = R_xhi.array();
                arr_ylo  = R_ylo.array(); arr_yhi = R_yhi.array();
                rhs_arr  = S_rhs[IntVar::cons].array(mfi);
                data_arr = S_data[IntVar::cons].array(mfi);
            } else if (ivar  == WRFBdyVars::T){
                arr_xlo  = T_xlo.array(); arr_xhi = T_xhi.array();
                arr_ylo  = T_ylo.array(); arr_yhi = T_yhi.array();
                rhs_arr  = S_rhs[IntVar::cons].array(mfi);
                data_arr = S_data[IntVar::cons].array(mfi);
            } else {
                continue;
            }

            compute_Laplacian_relaxation(delta_t, icomp, 1, width, set_width, dom_lo, dom_hi, F1, F2,
                                         tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi,
                                         arr_xlo, arr_xhi, arr_ylo, arr_yhi,
                                         data_arr, rhs_arr);
        } // mfi
    } // ivar
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
fine_compute_interior_ghost_RHS(const Real& time,
                                const Real& delta_t,
                                const int& width,
                                const int& set_width,
                                const Geometry& geom,
                                ERFFillPatcher* FPr_c,
                                ERFFillPatcher* FPr_u,
                                ERFFillPatcher* FPr_v,
                                ERFFillPatcher* FPr_w,
                                Vector<Box>& boxes_at_level,
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
        int num_var = 1;
        if (ivar_idx == IntVar::cons) num_var = 2;
        fmf_p_v.emplace_back(fmf.boxArray(), fmf.DistributionMap(), num_var, fmf.nGrowVect());
        MultiFab& fmf_p = fmf_p_v[ivar_idx];
        amrex::MultiFab::Copy(fmf_p,fmf, 0, 0, num_var, fmf.nGrowVect());

        // Integer mask MF
        iMultiFab* mask;

        // Fill fine patch on interior halo region
        //==========================================================
        if (ivar_idx == IntVar::cons)
        {
            FPr_c->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask = FPr_c->GetMask();
        }
        else if (ivar_idx == IntVar::xmom)
        {
            FPr_u->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask = FPr_u->GetMask();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(fmf_p,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();
                const Array4<Real>& prim_arr = fmf_p.array(mfi);
                const Array4<const Real>& rho_arr  = fmf_p_v[0].const_array(mfi);
                const Array4<const int>&  mask_arr = mask->const_array(mfi);

                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (mask_arr(i,j,k) == 1) {
                        Real rho_interp = 0.5 * ( rho_arr(i-1,j,k) + rho_arr(i,j,k) );
                        prim_arr(i,j,k) *= rho_interp;
                    }
                });
            } // mfi
        }
        else if (ivar_idx == IntVar::ymom)
        {
            FPr_v->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask = FPr_v->GetMask();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(fmf_p,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();
                const Array4<Real>& prim_arr = fmf_p.array(mfi);
                const Array4<const Real>& rho_arr  = fmf_p_v[0].const_array(mfi);
                const Array4<const int>&  mask_arr = mask->const_array(mfi);

                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (mask_arr(i,j,k) == 1) {
                        Real rho_interp = 0.5 * ( rho_arr(i,j-1,k) + rho_arr(i,j,k) );
                        prim_arr(i,j,k) *= rho_interp;
                    }
                });
            } // mfi
        }
        else if (ivar_idx == IntVar::zmom)
        {
            FPr_w->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask = FPr_w->GetMask();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(fmf_p,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();
                const Array4<Real>& prim_arr = fmf_p.array(mfi);
                const Array4<const Real>& rho_arr  = fmf_p_v[0].const_array(mfi);
                const Array4<const int>&  mask_arr = mask->const_array(mfi);

                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (mask_arr(i,j,k) == 1) {
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

            amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (mask_arr(i,j,k) == 2) {
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
            amrex::ParallelFor(vbx, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
               if (mask_arr(i,j,k) == 1) {

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
                       if ((near_x_lo_wall || near_x_hi_wall) && (near_y_lo_wall || near_y_hi_wall))
                           amrex::Abort("Mask must be found if you are near a corner!");

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
                               /*
                               amrex::Print() << "X CORNER FOUND: "
                                              << IntVect(i,j,k) << ' '
                                              << vbx << ' '
                                              << dj_min << ' '
                                              << ii << ' '
                                              << n_ind << "\n";
                               */
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
                               /*
                               amrex::Print() << "Y CORNER FOUND: "
                                              << IntVect(i,j,k) << ' '
                                              << vbx << ' '
                                              << di_min << ' '
                                              << jj << ' '
                                              << n_ind << "\n";
                               */
                           } else {
                               amrex::Abort("Mask not found near y wall!");
                           }
                       } else {
                           amrex::Abort("Relaxation cell must be near a wall!");
                       }
                   }

                   /*
                   if(n_ind<2 || n_ind>4) {
                       amrex::Print() << "SHIT: " << ivar_idx << ' ' << IntVect(i,j,k) << ' '
                                      << vbx << ' ' << n_ind << "\n";
                       exit(0);
                   }
                   */

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
                   amrex::Real HOLD = rhs_arr(i,j,k,n);
                   rhs_arr(i,j,k,n) += (F1*delta - F2*Laplacian) * Factor;

                   /*
                   Real eps = std::numeric_limits<Real>::epsilon();
                   if (delta   /(fine_arr(i  ,j  ,k,n)+eps)    > 0.1 ||
                       delta_xp/(fine_arr(i+1,j  ,k,n)+eps) > 0.1 ||
                       delta_xm/(fine_arr(i-1,j  ,k,n)+eps) > 0.1 ||
                       delta_yp/(fine_arr(i  ,j+1,k,n)+eps) > 0.1 ||
                       delta_ym/(fine_arr(i  ,j-1,k,n)+eps) > 0.1 ) {
                       amrex::Print() << "ERROR: " << ivar_idx << ' ' << n << ' ' << vbx << ' ' << IntVect(i,j,k) << ' ' << delta << ' ' << delta_xp << ' '
                                      << delta_xm << ' ' << delta_yp << ' ' << delta_ym << ' ' << Factor << ' ' << n_ind << ' ' << (F1*delta - F2*Laplacian) * Factor << ' ' << HOLD << "\n";
                       amrex::Print() << fine_arr(i  ,j  ,k,n) << ' ' << data_arr(i  ,j  ,k  ,n+icomp) << ' ' << d     << ' ' << rhs_arr(i  , j  , k  ,n+icomp) << "\n";
                       amrex::Print() << fine_arr(i+1,j  ,k,n) << ' ' << data_arr(i+1,j  ,k  ,n+icomp) << ' ' << d_ip1 << ' ' << rhs_arr(i+1, j  , k  ,n+icomp) << "\n";
                       amrex::Print() << fine_arr(i-1,j  ,k,n) << ' ' << data_arr(i-1,j  ,k  ,n+icomp) << ' ' << d_im1 << ' ' << rhs_arr(i-1, j  , k  ,n+icomp) << "\n";
                       amrex::Print() << fine_arr(i  ,j+1,k,n) << ' ' << data_arr(i  ,j+1,k  ,n+icomp) << ' ' << d_jp1 << ' ' << rhs_arr(i  , j+1, k  ,n+icomp) << "\n";
                       amrex::Print() << fine_arr(i  ,j-1,k,n) << ' ' << data_arr(i  ,j-1,k  ,n+icomp) << ' ' << d_jm1 << ' ' << rhs_arr(i  , j-1, k  ,n+icomp) << "\n";
                       amrex::Print() << "\n";
                       //exit(0);
                   }
                   */

               }
           });
        } // mfi
    } // ivar_idx
}

/**
 * Update the solution in the relaxation zone
 *
 * @param[in]  delta_t timestep
 * @param[in]  width   number of cells in (relaxation+specified) zone
 * @param[in]  geom    container for geometric information
 * @param[in]  S_rhs  RHS to be added here
 * @param[in]  S_old  previous value of the solution
 * @param[out] S_data new value of the solution defined here
 */

void
update_interior_ghost(const int& lev,
                      const Real& delta_t,
                      const int&  width,
                      Vector<Box>& boxes_at_level,
                      Vector<MultiFab>& S_rhs,
                      const Vector<MultiFab>& S_old,
                      Vector<MultiFab>& S_data,
                      ERFFillPatcher* FPr_c,
                      ERFFillPatcher* FPr_u,
                      ERFFillPatcher* FPr_v,
                      ERFFillPatcher* FPr_w)
{
    BL_PROFILE_REGION("update_interior_ghost()");

    Vector<int> ncomp_map = {2, 1, 1, 1};
    for (int ivar(IntVar::cons); ivar < IntVar::NumVars; ivar++)
    {
        int ncomp = ncomp_map[ivar];
        IndexType m_ixt = S_data[ivar].boxArray().ixType();

        // WRF BDY set region
        if (lev==0) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(S_data[ivar],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box tbx = mfi.tilebox();
                Array4<Real> rhs_arr       = S_rhs[ivar].array(mfi);
                Array4<const Real> old_arr = S_old[ivar].array(mfi);
                Array4<Real> new_arr       = S_data[ivar].array(mfi);

                for (int g_ind(0); g_ind<boxes_at_level.size(); ++g_ind)
                {
                    Box domain = boxes_at_level[g_ind];
                    domain.convert(m_ixt);
                    Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
                    compute_interior_ghost_bxs_xy(tbx, domain, width, 0,
                                                  tbx_xlo, tbx_xhi,
                                                  tbx_ylo, tbx_yhi);

                    amrex::ParallelFor(tbx_xlo, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        new_arr(i,j,k,n) = old_arr(i,j,k,n) + delta_t * rhs_arr(i,j,k,n);
                    },
                    tbx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        new_arr(i,j,k,n) = old_arr(i,j,k,n) + delta_t * rhs_arr(i,j,k,n);
                    });

                    amrex::ParallelFor(tbx_ylo, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        new_arr(i,j,k,n) = old_arr(i,j,k,n) + delta_t * rhs_arr(i,j,k,n);
                    },
                    tbx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        new_arr(i,j,k,n) = old_arr(i,j,k,n) + delta_t * rhs_arr(i,j,k,n);
                    });
                } // boxes_at_level
            } // mfi

        // CRSE-FINE set region
        } else {
            iMultiFab* mask;

            if (ivar == IntVar::cons) {
                mask = FPr_c->GetMask();
            } else if (ivar == IntVar::xmom) {
                mask = FPr_u->GetMask();
            } else if (ivar == IntVar::ymom) {
                mask = FPr_v->GetMask();
            } else if (ivar == IntVar::zmom) {
                mask = FPr_w->GetMask();
            } else {
                amrex::Abort("Dont recognize this variable type in update_interior_ghost");
            }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(S_data[ivar],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();
                Array4<Real> rhs_arr       = S_rhs[ivar].array(mfi);
                Array4<const Real> old_arr = S_old[ivar].array(mfi);
                Array4<Real> new_arr       = S_data[ivar].array(mfi);
                const Array4<const int> mask_arr = mask->const_array(mfi);

                amrex::ParallelFor(vbx, ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (mask_arr(i,j,k) == 2) {
                        rhs_arr(i,j,k,n) = 0.0;
                        new_arr(i,j,k,n) = old_arr(i,j,k,n);
                    }
                });
            } // mfi
        } // lev
    } // ivar
}
