#include <Utils.H>

using namespace amrex;

//
// Get the boxes for looping over interior/exterior ghost cells
// for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post.
//
// NOTE: X-face boxes take ownership of the overlapping region.
//       With exterior ghost cells (ng_vect != 0), the x-face
//       boxes will have exterior ghost cells in both x & y.
//
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

void
compute_interior_ghost_RHS(const Real& bdy_time_interval,
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
    BL_PROFILE_REGION("compute_interior_ghost_RHS()");

    // Relaxation constants
    Real F1 = 1./(10.*delta_t);
    Real F2 = 1./(50.*delta_t);

    // Compute interpolation factors
    Real dT = bdy_time_interval;
    int n_time = static_cast<int>(time / dT);
    amrex::Real alpha = (time - n_time * dT) / dT;
    amrex::Real oma   = 1.0 - alpha;

    // Temporary FABs for storage (owned/filled on all ranks)
    FArrayBox U_xlo, U_xhi, U_ylo, U_yhi;
    FArrayBox V_xlo, V_xhi, V_ylo, V_yhi;
    FArrayBox R_xlo, R_xhi, R_ylo, R_yhi;
    FArrayBox T_xlo, T_xhi, T_ylo, T_yhi;

    // Variable index map
    Vector<int> var_map = {Vars::xvel, Vars::yvel, Vars::cons, Vars::cons};

    // Variable icomp map
    Vector<int> comp_map = {0, 0, Rho_comp, RhoTheta_comp};

    // Size the FABs for each variable
    //==========================================================
    for (int ivar(WRFBdyVars::U); ivar <= WRFBdyVars::T; ivar++)
    {
        // Convert the domain to the ixtype of the variable
        int var_idx = var_map[ivar];
        Box domain  = geom.Domain();
        domain.convert(S_data[var_idx].boxArray().ixType());

        // Grown domain to get the 4 halo boxes w/ ghost cells
        // NOTE: 2 ghost cells needed for U -> rho*U
        IntVect ng_vect{2,2,0};
        Box gdom(domain); gdom.grow(ng_vect);

        // 4 halo boxes w/ interior ghost cells
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        // Size the FABs
        if (ivar  == WRFBdyVars::U) {
            U_xlo.resize(bx_xlo,1); U_xhi.resize(bx_xhi,1);
            U_ylo.resize(bx_ylo,1); U_yhi.resize(bx_yhi,1);
        } else if (ivar  == WRFBdyVars::V) {
            V_xlo.resize(bx_xlo,1); V_xhi.resize(bx_xhi,1);
            V_ylo.resize(bx_ylo,1); V_yhi.resize(bx_yhi,1);
        } else if (ivar  == WRFBdyVars::R) {
            R_xlo.resize(bx_xlo,1); R_xhi.resize(bx_xhi,1);
            R_ylo.resize(bx_ylo,1); R_yhi.resize(bx_yhi,1);
        } else {
            T_xlo.resize(bx_xlo,1); T_xhi.resize(bx_xhi,1);
            T_ylo.resize(bx_ylo,1); T_yhi.resize(bx_yhi,1);
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
    for (int ivar(WRFBdyVars::U); ivar <= WRFBdyVars::T; ivar++)
    {
        // Convert the domain to the ixtype of the variable
        int var_idx = var_map[ivar];
        Box domain  = geom.Domain();
        domain.convert(S_data[var_idx].boxArray().ixType());
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);

        // Grown domain to get the 4 halo boxes w/ ghost cells
        // NOTE: 2 ghost cells needed for U -> rho*U
        IntVect ng_vect{2,2,0};
        Box gdom(domain); gdom.grow(ng_vect);

        // 4 halo boxes w/ interior ghost cells
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        // Get the FAB array4s
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
        } else {
            arr_xlo = T_xlo.array(); arr_xhi = T_xhi.array();
            arr_ylo = T_ylo.array(); arr_yhi = T_yhi.array();
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
        amrex::ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_lo.x+width-1);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            arr_xlo(i,j,k) = oma   * bdatxlo_n  (ii,jj,k,0)
                           + alpha * bdatxlo_np1(ii,jj,k,0);
        });
        amrex::ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_hi.x-width+1);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            arr_xhi(i,j,k) = oma   * bdatxhi_n  (ii,jj,k,0)
                           + alpha * bdatxhi_np1(ii,jj,k,0);
        });
        amrex::ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_lo.y+width-1);
            arr_ylo(i,j,k) = oma   * bdatylo_n  (ii,jj,k,0)
                           + alpha * bdatylo_np1(ii,jj,k,0);
        });
        amrex::ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
        // Convert the domain to the ixtype of the variable
        int var_idx = var_map[ivar];
        Box domain  = geom.Domain();
        domain.convert(S_data[var_idx].boxArray().ixType());

        // Grown domain to get the 4 halo boxes w/ ghost cells
        // NOTE: 1 ghost cell needed for Laplacian
        IntVect ng_vect{1,1,0};
        Box gdom(domain); gdom.grow(ng_vect);

        // 4 halo boxes w/ interior ghost cells
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        compute_interior_ghost_bxs_xy(gdom, domain, width, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect, true);

        // Density FAB array4s
        Array4<Real> rarr_xlo = R_xlo.array();  Array4<Real> rarr_xhi = R_xhi.array();
        Array4<Real> rarr_ylo = R_ylo.array();  Array4<Real> rarr_yhi = R_yhi.array();

        // Get the FAB array4s
        Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
        Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
        if (ivar  == WRFBdyVars::U) {
            arr_xlo = U_xlo.array(); arr_xhi = U_xhi.array();
            arr_ylo = U_ylo.array(); arr_yhi = U_yhi.array();

            amrex::ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xlo(i-1,j,k) + rarr_xlo(i,j,k) );
                arr_xlo(i,j,k) *= rho_interp;
            });
            amrex::ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xhi(i-1,j,k) + rarr_xhi(i,j,k) );
                arr_xhi(i,j,k) *= rho_interp;
            });
            amrex::ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_ylo(i-1,j,k) + rarr_ylo(i,j,k) );
                arr_ylo(i,j,k) *= rho_interp;
            });
            amrex::ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_yhi(i-1,j,k) + rarr_yhi(i,j,k) );
                arr_yhi(i,j,k) *= rho_interp;
            });
        } else {
            arr_xlo = V_xlo.array(); arr_xhi = V_xhi.array();
            arr_ylo = V_ylo.array(); arr_yhi = V_yhi.array();

            amrex::ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xlo(i,j-1,k) + rarr_xlo(i,j,k) );
                arr_xlo(i,j,k) *= rho_interp;
            });
            amrex::ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_xhi(i,j-1,k) + rarr_xhi(i,j,k) );
                arr_xhi(i,j,k) *= rho_interp;
            });
            amrex::ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_ylo(i,j-1,k) + rarr_ylo(i,j,k) );
                arr_ylo(i,j,k) *= rho_interp;
            });
            amrex::ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rho_interp = 0.5 * ( rarr_yhi(i,j-1,k) + rarr_yhi(i,j,k) );
                arr_yhi(i,j,k) *= rho_interp;
            });
        }
    } // ivar


    // Compute the RHS (fill all CONS & W first)
    //==========================================================
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_rhs[IntVar::cons],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        {
            // 4 halo boxes w/ no ghost cells for CC vars
            Box domain = geom.Domain();
            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            compute_interior_ghost_bxs_xy(domain, domain, width, 0,
                                          bx_xlo, bx_xhi,
                                          bx_ylo, bx_yhi);
            S_rhs[IntVar::cons][mfi].template setVal<RunOn::Device>(0.,bx_xlo,Rho_comp,NVAR);
            S_rhs[IntVar::cons][mfi].template setVal<RunOn::Device>(0.,bx_xhi,Rho_comp,NVAR);
            S_rhs[IntVar::cons][mfi].template setVal<RunOn::Device>(0.,bx_ylo,Rho_comp,NVAR);
            S_rhs[IntVar::cons][mfi].template setVal<RunOn::Device>(0.,bx_yhi,Rho_comp,NVAR);
        }

        {
            // 4 halo boxes w/ no ghost cells for CC vars
            Box domain = geom.Domain();
            domain.convert(S_rhs[IntVar::zmom].boxArray().ixType());
            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            compute_interior_ghost_bxs_xy(domain, domain, width, 0,
                                          bx_xlo, bx_xhi,
                                          bx_ylo, bx_yhi);
            S_rhs[IntVar::zmom][mfi].template setVal<RunOn::Device>(0.,bx_xlo);
            S_rhs[IntVar::zmom][mfi].template setVal<RunOn::Device>(0.,bx_xhi);
            S_rhs[IntVar::zmom][mfi].template setVal<RunOn::Device>(0.,bx_ylo);
            S_rhs[IntVar::zmom][mfi].template setVal<RunOn::Device>(0.,bx_yhi);
        }
    } // mfi

    for (int ivar(WRFBdyVars::U); ivar <= WRFBdyVars::T; ivar++)
    {
        // Variable and comp maps
        int var_idx =  var_map[ivar];
        int icomp   = comp_map[ivar];

        // Convert the domain to the ixtype of the variable
        Box domain = geom.Domain();
        domain.convert(S_data[var_idx].boxArray().ixType());
        const auto& dom_hi = ubound(domain);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(S_data[var_idx],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box tbx = mfi.tilebox();

            // Intersection of tilebox and 4 boxes w/o ANY ghost cells
            // NOTE: 0 ghost cells needed
            Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
            compute_interior_ghost_bxs_xy(tbx, domain, width, 0,
                                          tbx_xlo, tbx_xhi,
                                          tbx_ylo, tbx_yhi);

            // Get the array4s
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
            } else {
                arr_xlo  = T_xlo.array(); arr_xhi = T_xhi.array();
                arr_ylo  = T_ylo.array(); arr_yhi = T_yhi.array();
                rhs_arr  = S_rhs[IntVar::cons].array(mfi);
                data_arr = S_data[IntVar::cons].array(mfi);
            }

            // RHS computation
            int Spec_z  = set_width;
            int Relax_z = width - Spec_z + 1;
            Real num    = Real(Spec_z + Relax_z);
            Real denom  = Real(Relax_z - 1);
            amrex::ParallelFor(tbx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Corners with x boxes
                int j_lo = std::min(j,width-1);
                int j_hi = std::min(dom_hi.y-j,width-1);
                int jj   = std::min(j_lo,j_hi);
                int n    = std::min(i,jj) + 1;
                if (n <= Spec_z) {
                    rhs_arr(i,j,k,icomp) = 0.0;
                } else {
                    Real Factor   = (num - Real(n))/denom;
                    Real delta    = arr_xlo(i  ,j  ,k) - data_arr(i  ,j  ,k,icomp);
                    Real delta_xp = arr_xlo(i+1,j  ,k) - data_arr(i+1,j  ,k,icomp);
                    Real delta_xm = arr_xlo(i-1,j  ,k) - data_arr(i-1,j  ,k,icomp);
                    Real delta_yp = arr_xlo(i  ,j+1,k) - data_arr(i  ,j+1,k,icomp);
                    Real delta_ym = arr_xlo(i  ,j-1,k) - data_arr(i  ,j-1,k,icomp);
                    Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
                    rhs_arr(i,j,k,icomp) = (F1*delta - F2*Laplacian) * Factor;
                }
            });
            amrex::ParallelFor(tbx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Corners with x boxes
                int j_lo = std::min(j,width-1);
                int j_hi = std::min(dom_hi.y-j,width-1);
                int jj   = std::min(j_lo,j_hi);
                int n    = std::min(dom_hi.x-i,jj) + 1;
                if (n <= Spec_z) {
                    rhs_arr(i,j,k,icomp) = 0.0;
                } else {
                    Real Factor   = (num - Real(n))/denom;
                    Real delta    = arr_xhi(i  ,j  ,k) - data_arr(i  ,j  ,k,icomp);
                    Real delta_xp = arr_xhi(i+1,j  ,k) - data_arr(i+1,j  ,k,icomp);
                    Real delta_xm = arr_xhi(i-1,j  ,k) - data_arr(i-1,j  ,k,icomp);
                    Real delta_yp = arr_xhi(i  ,j+1,k) - data_arr(i  ,j+1,k,icomp);
                    Real delta_ym = arr_xhi(i  ,j-1,k) - data_arr(i  ,j-1,k,icomp);
                    Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
                    rhs_arr(i,j,k,icomp) = (F1*delta - F2*Laplacian) * Factor;
                }
            });
            amrex::ParallelFor(tbx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // No corners for y boxes
                int n = j + 1;
                if (n <= Spec_z) {
                    rhs_arr(i,j,k,icomp) = 0.0;
                } else {
                    Real Factor   = (num - Real(n))/denom;
                    Real delta    = arr_ylo(i  ,j  ,k) - data_arr(i  ,j  ,k,icomp);
                    Real delta_xp = arr_ylo(i+1,j  ,k) - data_arr(i+1,j  ,k,icomp);
                    Real delta_xm = arr_ylo(i-1,j  ,k) - data_arr(i-1,j  ,k,icomp);
                    Real delta_yp = arr_ylo(i  ,j+1,k) - data_arr(i  ,j+1,k,icomp);
                    Real delta_ym = arr_ylo(i  ,j-1,k) - data_arr(i  ,j-1,k,icomp);
                    Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
                    rhs_arr(i,j,k,icomp) = (F1*delta - F2*Laplacian) * Factor;
                }
            });
            amrex::ParallelFor(tbx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // No corners for y boxes
                int n = dom_hi.y - j + 1;
                if (n <= Spec_z) {
                    rhs_arr(i,j,k,icomp) = 0.0;
                } else {
                    Real Factor   = (num - Real(n))/denom;
                    Real delta    = arr_yhi(i  ,j  ,k) - data_arr(i  ,j  ,k,icomp);
                    Real delta_xp = arr_yhi(i+1,j  ,k) - data_arr(i+1,j  ,k,icomp);
                    Real delta_xm = arr_yhi(i-1,j  ,k) - data_arr(i-1,j  ,k,icomp);
                    Real delta_yp = arr_yhi(i  ,j+1,k) - data_arr(i  ,j+1,k,icomp);
                    Real delta_ym = arr_yhi(i  ,j-1,k) - data_arr(i  ,j-1,k,icomp);
                    Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
                    rhs_arr(i,j,k,icomp) = (F1*delta - F2*Laplacian) * Factor;
                }
            });
        } // mfi
    } // ivar
}

// TODO: Figure out if this is needed for fast integrator
//       and how to incorporate if it is needed.
void
update_interior_ghost(const Real& delta_t,
                      const int&  width,
                      const Geometry& geom,
                      Vector<MultiFab>& S_rhs,
                      const Vector<MultiFab>& S_old,
                      Vector<MultiFab>& S_data)
{
    BL_PROFILE_REGION("update_interior_ghost()");

    // Copy of the Domain
    Box domain = geom.Domain();

    // Variable index map
    Vector<int> var_map = {Vars::xvel, Vars::yvel, Vars::cons, Vars::cons};

    // Variable icomp map
    Vector<int> comp_map = {0, 0, Rho_comp, RhoTheta_comp};

    // Update the vars on the 4 halo regions
    //==========================================================
    for (int ivar(WRFBdyVars::U); ivar <= WRFBdyVars::T; ivar++)
    {
        // Variable and comp maps
        int var_idx =  var_map[ivar];
        int icomp   = comp_map[ivar];

        // Convert the domain to the ixtype of the variable
        domain.convert(S_data[var_idx].boxArray().ixType());

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(S_data[var_idx],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box tbx = mfi.tilebox();

            // Intersection of tilebox and 4 boxes w/o ANY ghost cells
            Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
            compute_interior_ghost_bxs_xy(tbx, domain, width, 0,
                                          tbx_xlo, tbx_xhi,
                                          tbx_ylo, tbx_yhi);

            // Get the array4s
            Array4<Real> rhs_arr; Array4<const Real> old_arr; Array4<Real> new_arr;
            if (ivar  == WRFBdyVars::U) {
                rhs_arr = S_rhs[IntVar::xmom].array(mfi);
                old_arr = S_old[IntVar::xmom].array(mfi);
                new_arr = S_data[IntVar::xmom].array(mfi);
            } else if (ivar  == WRFBdyVars::V) {
                rhs_arr = S_rhs[IntVar::ymom].array(mfi);
                old_arr = S_old[IntVar::ymom].array(mfi);
                new_arr = S_data[IntVar::ymom].array(mfi);
            } else if (ivar  == WRFBdyVars::R) {
                rhs_arr = S_rhs[IntVar::cons].array(mfi);
                old_arr = S_old[IntVar::cons].array(mfi);
                new_arr = S_data[IntVar::cons].array(mfi);
            } else {
                rhs_arr = S_rhs[IntVar::cons].array(mfi);
                old_arr = S_old[IntVar::cons].array(mfi);
                new_arr = S_data[IntVar::cons].array(mfi);
            }

            // Update computation
            amrex::ParallelFor(tbx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                new_arr(i,j,k,icomp) = old_arr(i,j,k,icomp) + delta_t * rhs_arr(i,j,k,icomp);
            });
            amrex::ParallelFor(tbx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                new_arr(i,j,k,icomp) = old_arr(i,j,k,icomp) + delta_t * rhs_arr(i,j,k,icomp);
            });
            amrex::ParallelFor(tbx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                new_arr(i,j,k,icomp) = old_arr(i,j,k,icomp) + delta_t * rhs_arr(i,j,k,icomp);
            });
            amrex::ParallelFor(tbx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                new_arr(i,j,k,icomp) = old_arr(i,j,k,icomp) + delta_t * rhs_arr(i,j,k,icomp);
            });
        } // mfi
    } // ivar
}
