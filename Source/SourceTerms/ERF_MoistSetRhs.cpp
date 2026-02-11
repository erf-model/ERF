#if defined(ERF_USE_NETCDF)

#include <ERF_SrcHeaders.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Function for setting the slow variables in the "specified" zones at the domain boundary
 *
*/

void
moist_set_rhs (const Geometry& geom,
               const Box& tbx,
               const Array4<Real const>& old_cons,
               const Array4<Real const>& new_cons,
               const Array4<Real      >& cell_rhs,
               const Real& bdy_time_interval,
               const Real& new_stage_time,
               const Real& dt,
               const Real & stop_time_elapsed,
               int  width,
               int  set_width,
               bool do_upwind,
               const Box& domain,
               Vector<Vector<FArrayBox>>& bdy_data_xlo,
               Vector<Vector<FArrayBox>>& bdy_data_xhi,
               Vector<Vector<FArrayBox>>& bdy_data_ylo,
               Vector<Vector<FArrayBox>>& bdy_data_yhi)
{
    // Relaxation constants
    Real F1 = 1./dt;

    // Domain bounds
    const auto& dom_hi = ubound(domain);
    const auto& dom_lo = lbound(domain);
    auto dx = geom.CellSizeArray();
    auto ProbHi = geom.ProbHiArray();
    auto ProbLo = geom.ProbLoArray();

    // Time interpolation
    Real dT = bdy_time_interval;

    // NOTE: This is because we define "time" to be time since start_bdy_time
    Real time_since_start = new_stage_time;

    int n_time = static_cast<int>( time_since_start /  dT);
    Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    Real oma   = 1.0 - alpha;

    int n_time_p1 = n_time + 1;
    if ((new_stage_time == stop_time_elapsed) && (alpha==0)) {
        // stop time coincides with final bdy snapshot -- don't try to read in
        // another snapshot
        n_time_p1 = n_time;
    }

    /*
    // UNIT TEST DEBUG
    oma = 1.0; alpha = 0.0;
    */

    // NOTE: The sizing of the temporary BDY FABS is
    //       GLOBAL and occurs over the entire BDY region.

    // Size the FABs
    //==========================================================
    // NOTE: No ghost cells, we force mask to be idx type (0,0,0)
    IntVect ng_vect(0);
    Box gdom(domain); gdom.grow(ng_vect);
    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
    realbdy_interior_bxs_xy(gdom, domain, width,
                            bx_xlo, bx_xhi,
                            bx_ylo, bx_yhi,
                            0, ng_vect, true);

    // Temporary FABs for storage (owned/filled on all ranks)
    FArrayBox QV_xlo, QV_xhi, QV_ylo, QV_yhi;
    QV_xlo.resize(bx_xlo,1,The_Async_Arena()); QV_xhi.resize(bx_xhi,1,The_Async_Arena());
    QV_ylo.resize(bx_ylo,1,The_Async_Arena()); QV_yhi.resize(bx_yhi,1,The_Async_Arena());

    // Masks for upwinding
    FArrayBox U_xlo, U_xhi, V_xlo, V_xhi, V_ylo, V_yhi;
    U_xlo.resize(convert(bx_xlo,IntVect(1,0,0)),1,The_Async_Arena());
    U_xhi.resize(convert(bx_xhi,IntVect(1,0,0)),1,The_Async_Arena());
    V_xlo.resize(convert(bx_xlo,IntVect(0,1,0)),1,The_Async_Arena());
    V_xhi.resize(convert(bx_xhi,IntVect(0,1,0)),1,The_Async_Arena());
    V_ylo.resize(convert(bx_ylo,IntVect(0,1,0)),1,The_Async_Arena());
    V_yhi.resize(convert(bx_yhi,IntVect(0,1,0)),1,The_Async_Arena());

    // Populate FABs from bdy interpolation (primitive vars)
    //==========================================================
    const auto& bdatxlo_n   = bdy_data_xlo[n_time   ][WRFBdyVars::QV].const_array();
    const auto& bdatxlo_np1 = bdy_data_xlo[n_time_p1][WRFBdyVars::QV].const_array();
    const auto& bdatxhi_n   = bdy_data_xhi[n_time   ][WRFBdyVars::QV].const_array();
    const auto& bdatxhi_np1 = bdy_data_xhi[n_time_p1][WRFBdyVars::QV].const_array();
    const auto& bdatylo_n   = bdy_data_ylo[n_time   ][WRFBdyVars::QV].const_array();
    const auto& bdatylo_np1 = bdy_data_ylo[n_time_p1][WRFBdyVars::QV].const_array();
    const auto& bdatyhi_n   = bdy_data_yhi[n_time   ][WRFBdyVars::QV].const_array();
    const auto& bdatyhi_np1 = bdy_data_yhi[n_time_p1][WRFBdyVars::QV].const_array();

    const auto& bdatxlo_n_u   = bdy_data_xlo[n_time   ][WRFBdyVars::U].const_array();
    const auto& bdatxlo_np1_u = bdy_data_xlo[n_time_p1][WRFBdyVars::U].const_array();
    const auto& bdatxhi_n_u   = bdy_data_xhi[n_time   ][WRFBdyVars::U].const_array();
    const auto& bdatxhi_np1_u = bdy_data_xhi[n_time_p1][WRFBdyVars::U].const_array();

    const auto& bdatxlo_n_v   = bdy_data_xlo[n_time   ][WRFBdyVars::V].const_array();
    const auto& bdatxlo_np1_v = bdy_data_xlo[n_time_p1][WRFBdyVars::V].const_array();
    const auto& bdatxhi_n_v   = bdy_data_xhi[n_time   ][WRFBdyVars::V].const_array();
    const auto& bdatxhi_np1_v = bdy_data_xhi[n_time_p1][WRFBdyVars::V].const_array();

    const auto& bdatylo_n_v   = bdy_data_ylo[n_time   ][WRFBdyVars::V].const_array();
    const auto& bdatylo_np1_v = bdy_data_ylo[n_time_p1][WRFBdyVars::V].const_array();
    const auto& bdatyhi_n_v   = bdy_data_yhi[n_time   ][WRFBdyVars::V].const_array();
    const auto& bdatyhi_np1_v = bdy_data_yhi[n_time_p1][WRFBdyVars::V].const_array();

    // Get Array4 of interpolated values
    Array4<Real> arr_xlo = QV_xlo.array();  Array4<Real> arr_xhi = QV_xhi.array();
    Array4<Real> arr_ylo = QV_ylo.array();  Array4<Real> arr_yhi = QV_yhi.array();

    Array4<Real> u_xlo = U_xlo.array();  Array4<Real> u_xhi = U_xhi.array();
    Array4<Real> v_xlo = V_xlo.array();  Array4<Real> v_xhi = V_xhi.array();
    Array4<Real> v_ylo = V_ylo.array();  Array4<Real> v_yhi = V_yhi.array();

    Box gtbx = grow(tbx,ng_vect);
    Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
    realbdy_interior_bxs_xy(gtbx, domain, width,
                            tbx_xlo, tbx_xhi,
                            tbx_ylo, tbx_yhi,
                            0, ng_vect, true);

    // Limiting offset
    int offset = set_width - 1;
    if (width > set_width) offset = width - 1;

    // Populate with interpolation (protect from ghost cells)
    ParallelFor(tbx_xlo, tbx_xhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_lo.x), dom_lo.x+offset);
        int jj = std::min(std::max(j , dom_lo.y), dom_hi.y       );
        arr_xlo(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatxlo_n  (ii,jj,k)
                                                    + alpha * bdatxlo_np1(ii,jj,k) );
        u_xlo(i,j,k) = ( oma * bdatxlo_n_u(ii,jj,k) + alpha * bdatxlo_np1_u(ii,jj,k) );
        v_xlo(i,j,k) = ( oma * bdatxlo_n_v(ii,jj,k) + alpha * bdatxlo_np1_v(ii,jj,k) );
        if (j == dom_hi.y) {
            v_xlo(i,dom_hi.y+1,k) = ( oma   * bdatxlo_n_v  (ii,dom_hi.y+1,k)
                                    + alpha * bdatxlo_np1_v(ii,dom_hi.y+1,k) );
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_hi.x-offset), dom_hi.x);
        int jj = std::min(std::max(j , dom_lo.y       ), dom_hi.y);
        arr_xhi(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatxhi_n  (ii,jj,k)
                                                    + alpha * bdatxhi_np1(ii,jj,k) );
        // NOTE: correct for idx type mismatch with u bdy data
        u_xhi(i+1,j,k) = ( oma * bdatxhi_n_u(ii+1,jj,k) + alpha * bdatxhi_np1_u(ii+1,jj,k) );
        v_xhi(i  ,j,k) = ( oma * bdatxhi_n_v(ii  ,jj,k) + alpha * bdatxhi_np1_v(ii,jj,k) );
        if (j == dom_hi.y) {
            v_xhi(i,dom_hi.y+1,k) = ( oma   * bdatxhi_n_v  (ii,dom_hi.y+1,k)
                                    + alpha * bdatxhi_np1_v(ii,dom_hi.y+1,k) );
        }
    });

    ParallelFor(tbx_ylo, tbx_yhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_lo.x), dom_hi.x       );
        int jj = std::min(std::max(j , dom_lo.y), dom_lo.y+offset);
        arr_ylo(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatylo_n  (ii,jj,k)
                                                    + alpha * bdatylo_np1(ii,jj,k) );
        v_ylo(i,j,k) = ( oma * bdatylo_n_v(ii,jj,k) + alpha * bdatylo_np1_v(ii,jj,k) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_lo.x       ), dom_hi.x);
        int jj = std::min(std::max(j , dom_hi.y-offset), dom_hi.y);
            jj = std::min(jj, dom_hi.y);
        arr_yhi(i,j,k) = new_cons(i,j,k,Rho_comp) * ( oma   * bdatyhi_n  (ii,jj,k)
                                                    + alpha * bdatyhi_np1(ii,jj,k) );
        // NOTE: correct for idx type mismatch with v bdy data
        v_yhi(i,j+1,k) = ( oma * bdatyhi_n_v(ii,jj+1,k) + alpha * bdatyhi_np1_v(ii,jj+1,k) );
    });


    // NOTE: We pass 'old_cons' here since the tendencies are with
    //       respect to the start of the RK integration.

    // Compute RHS in specified region
    //==========================================================
    if (set_width > 0) {
        realbdy_interior_bxs_xy(tbx, domain, width,
                                tbx_xlo, tbx_xhi,
                                tbx_ylo, tbx_yhi);
        realbdy_set_rhs_in_spec_region(dt, RhoQ1_comp, 1,
                                       width, set_width-1, set_width-1,
                                       domain, domain,
                                       tbx_xlo , tbx_xhi , tbx_ylo , tbx_yhi ,
                                       arr_xlo , arr_xhi , arr_ylo , arr_yhi ,
                                       u_xlo, u_xhi, v_xlo, v_xhi, v_ylo, v_yhi,
                                       old_cons, cell_rhs, do_upwind);
    }


    // NOTE: We pass 'new_cons' here since it has its ghost cells
    //       populated and we are only operating on RhoQv; thus,
    //       we do not need the updated fast quantities.

    // Compute RHS in relaxation region
    //==========================================================
    if (width > set_width) {
        realbdy_interior_bxs_xy(tbx, domain, width,
                                tbx_xlo, tbx_xhi,
                                tbx_ylo, tbx_yhi,
                                set_width, ng_vect);
        realbdy_compute_relaxation(RhoQ1_comp, 1,
                                   width, dx, ProbLo, ProbHi, F1, domain,
                                   tbx_xlo , tbx_xhi , tbx_ylo , tbx_yhi ,
                                   arr_xlo , arr_xhi , arr_ylo , arr_yhi ,
                                   u_xlo, u_xhi, v_xlo, v_xhi, v_ylo, v_yhi,
                                   new_cons, cell_rhs, do_upwind);
    }

    /*
    // UNIT TEST DEBUG
    realbdy_interior_bxs_xy(tbx, domain, width,
                            tbx_xlo, tbx_xhi,
                            tbx_ylo, tbx_yhi);
    ParallelFor(tbx_xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (std::fabs(arr_xlo(i,j,k) - new_cons(i,j,k,RhoQ1_comp)) > 1.0e-7) {
            Print() << "ERROR XLO: " <<  RhoQ1_comp << ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    ParallelFor(tbx_xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (std::fabs(arr_xhi(i,j,k) - new_cons(i,j,k,RhoQ1_comp)) > 1.0e-7) {
            Print() << "ERROR XHI: " << RhoQ1_comp<< ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    ParallelFor(tbx_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (std::fabs(arr_ylo(i,j,k) - new_cons(i,j,k,RhoQ1_comp))> 1.0e-7) {
            Print() << "ERROR YLO: " << RhoQ1_comp << ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    ParallelFor(tbx_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (std::fabs(arr_yhi(i,j,k) - new_cons(i,j,k,RhoQ1_comp))> 1.0e-7) {
            Print() << "ERROR YHI: " << RhoQ1_comp << ' ' << IntVect(i,j,k) << "\n";
            exit(0);
        }
    });
    exit(0);
    */
} // moist_set_rhs
#endif
