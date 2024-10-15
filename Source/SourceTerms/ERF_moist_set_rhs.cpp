#if defined(ERF_USE_NETCDF)

#include <ERF_Src_headers.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Function for setting the slow variables in the "specified" zones at the domain boundary
 *
*/

void
moist_set_rhs (const Box& tbx,
               const Box& gtbx,
               const Array4<Real const>& old_cons,
               const Array4<Real const>& new_cons,
               const Array4<Real      >& cell_rhs,
               const Real& bdy_time_interval,
               const Real& start_bdy_time,
               const Real& new_stage_time,
               const Real& dt,
               int  width,
               int  set_width,
               const Box& domain,
               Vector<Vector<FArrayBox>>& bdy_data_xlo,
               Vector<Vector<FArrayBox>>& bdy_data_xhi,
               Vector<Vector<FArrayBox>>& bdy_data_ylo,
               Vector<Vector<FArrayBox>>& bdy_data_yhi)
{
    // NOTE: We pass the full width into this routine. For relaxation, the last cell is a halo
    //       cell for the Laplacian. We remove that cell here if it is present.

    // The width to do RHS augmentation
    if (width > set_width+1) width -= 2;

    // Relaxation constants
    Real F1 = 1./(10.*dt);
    Real F2 = 1./(50.*dt);

    // Domain bounds
    const auto& dom_hi = ubound(domain);
    const auto& dom_lo = lbound(domain);

    // Time interpolation
    Real dT = bdy_time_interval;
    Real time_since_start = new_stage_time - start_bdy_time;
    int n_time = static_cast<int>( time_since_start /  dT);
    Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    Real oma   = 1.0 - alpha;

    /*
    // UNIT TEST DEBUG
    oma = 1.0; alpha = 0.0;
    */

    // NOTE: These sizing of the temporary BDY FABS is
    //       GLOBAL and occurs over the entire BDY region.

    // Size the FABs
    //==========================================================
    IntVect ng_vect{2,2,0};
    Box gdom(domain); gdom.grow(ng_vect);
    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
    compute_interior_ghost_bxs_xy(gdom, domain, width, 0,
                                  bx_xlo, bx_xhi,
                                  bx_ylo, bx_yhi,
                                  ng_vect, true);

    // Temporary FABs for storage (owned/filled on all ranks)
    FArrayBox QV_xlo, QV_xhi, QV_ylo, QV_yhi;
    QV_xlo.resize(bx_xlo,1,The_Async_Arena()); QV_xhi.resize(bx_xhi,1,The_Async_Arena());
    QV_ylo.resize(bx_ylo,1,The_Async_Arena()); QV_yhi.resize(bx_yhi,1,The_Async_Arena());


    // NOTE: These operations use the BDY FABS and RHO. The
    //       use of RHO to go from PRIM -> CONS requires that
    //       these operations be LOCAL. So we have allocated
    //       enough space to do global operations (1 rank) but
    //       will fill a subset of that data that the rank owns.

    // Populate FABs from bdy interpolation (primitive vars)
    //==========================================================
    const auto& bdatxlo_n   = bdy_data_xlo[n_time  ][WRFBdyVars::QV].const_array();
    const auto& bdatxlo_np1 = bdy_data_xlo[n_time+1][WRFBdyVars::QV].const_array();
    const auto& bdatxhi_n   = bdy_data_xhi[n_time  ][WRFBdyVars::QV].const_array();
    const auto& bdatxhi_np1 = bdy_data_xhi[n_time+1][WRFBdyVars::QV].const_array();
    const auto& bdatylo_n   = bdy_data_ylo[n_time  ][WRFBdyVars::QV].const_array();
    const auto& bdatylo_np1 = bdy_data_ylo[n_time+1][WRFBdyVars::QV].const_array();
    const auto& bdatyhi_n   = bdy_data_yhi[n_time  ][WRFBdyVars::QV].const_array();
    const auto& bdatyhi_np1 = bdy_data_yhi[n_time+1][WRFBdyVars::QV].const_array();

    // Get Array4 of interpolated values
    Array4<Real> arr_xlo = QV_xlo.array();  Array4<Real> arr_xhi = QV_xhi.array();
    Array4<Real> arr_ylo = QV_ylo.array();  Array4<Real> arr_yhi = QV_yhi.array();

    Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
    compute_interior_ghost_bxs_xy(gtbx, domain, width, 0,
                                  tbx_xlo, tbx_xhi,
                                  tbx_ylo, tbx_yhi,
                                  ng_vect, true);

    // NOTE: width is now one less than the total bndy width
    //       if we have a relaxation zone; so we can access
    //       dom_lo/hi +- width. If we do not have a relax
    //       zone, this offset is set_width - 1.
    int offset = set_width - 1;
    if (width > set_width) offset = width;

    // Populate with interpolation (protect from ghost cells)
    ParallelFor(tbx_xlo, tbx_xhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::max(i , dom_lo.x);
            ii = std::min(ii, dom_lo.x+offset);
        int jj = std::max(j , dom_lo.y);
            jj = std::min(jj, dom_hi.y);
        arr_xlo(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatxlo_n  (ii,jj,k)
                                                    + alpha * bdatxlo_np1(ii,jj,k) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::max(i , dom_hi.x-offset);
            ii = std::min(ii, dom_hi.x);
        int jj = std::max(j , dom_lo.y);
            jj = std::min(jj, dom_hi.y);
        arr_xhi(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatxhi_n  (ii,jj,k)
                                                    + alpha * bdatxhi_np1(ii,jj,k) );
    });

    ParallelFor(tbx_ylo, tbx_yhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::max(i , dom_lo.x);
            ii = std::min(ii, dom_hi.x);
        int jj = std::max(j , dom_lo.y);
            jj = std::min(jj, dom_lo.y+offset);
        arr_ylo(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatylo_n  (ii,jj,k)
                                                    + alpha * bdatylo_np1(ii,jj,k) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::max(i , dom_lo.x);
            ii = std::min(ii, dom_hi.x);
        int jj = std::max(j , dom_hi.y-offset);
            jj = std::min(jj, dom_hi.y);
        arr_yhi(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatyhi_n  (ii,jj,k)
                                                    + alpha * bdatyhi_np1(ii,jj,k) );
    });


    // NOTE: We pass 'old_cons' here since the tendencies are with
    //       respect to the start of the RK integration.

    // Compute RHS in specified region
    //==========================================================
    if (set_width > 0) {
        compute_interior_ghost_bxs_xy(tbx, domain, width, 0,
                                      tbx_xlo, tbx_xhi,
                                      tbx_ylo, tbx_yhi,
                                      ng_vect, true);
        wrfbdy_set_rhs_in_spec_region(dt, RhoQ1_comp, 1,
                                      width, set_width, dom_lo, dom_hi,
                                      tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi,
                                      arr_xlo, arr_xhi, arr_ylo, arr_yhi,
                                      old_cons, cell_rhs);
    }


    // NOTE: We pass 'new_cons' here since it has its ghost cells
    //       populated and we are only operating on RhoQv; thus,
    //       we do not need the updated fast quantities.

    // Compute RHS in relaxation region
    //==========================================================
    if (width > set_width) {
        compute_interior_ghost_bxs_xy(tbx, domain, width, set_width,
                                      tbx_xlo, tbx_xhi,
                                      tbx_ylo, tbx_yhi);
        wrfbdy_compute_laplacian_relaxation(RhoQ1_comp, 1,
                                            width, set_width, dom_lo, dom_hi, F1, F2,
                                            tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi,
                                            arr_xlo, arr_xhi, arr_ylo, arr_yhi,
                                            new_cons, cell_rhs);
    }

    /*
    // UNIT TEST DEBUG
    compute_interior_ghost_bxs_xy(tbx, domain, width+1, 0,
                                  bx_xlo, bx_xhi,
                                  bx_ylo, bx_yhi);
    ParallelFor(bx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (arr_xlo(i,j,k) != new_cons(i,j,k,RhoQ1_comp)) {
            Print() << "ERROR XLO: " <<  RhoQ1_comp << ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    ParallelFor(bx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (arr_xhi(i,j,k) != new_cons(i,j,k,RhoQ1_comp)) {
            Print() << "ERROR XHI: " << RhoQ1_comp<< ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    ParallelFor(bx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (arr_ylo(i,j,k) != new_cons(i,j,k,RhoQ1_comp)) {
            Print() << "ERROR YLO: " << RhoQ1_comp << ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    ParallelFor(bx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (arr_yhi(i,j,k) != new_cons(i,j,k,RhoQ1_comp)) {
            Print() << "ERROR YHI: " << RhoQ1_comp << ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    */
} // moist_set_rhs
#endif
