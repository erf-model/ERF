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
               const Array4<Real const>& new_cons,
               const Array4<Real      >& cell_rhs,
               const Real& time,
               const Real& dt,
               const Real& start_bdy_time,
               const Real& final_bdy_time,
               const Real& bdy_time_interval,
               const Real& nudge_factor,
               int  width,
               const Box& domain,
               Vector<Vector<FArrayBox>>& bdy_data_xlo,
               Vector<Vector<FArrayBox>>& bdy_data_xhi,
               Vector<Vector<FArrayBox>>& bdy_data_ylo,
               Vector<Vector<FArrayBox>>& bdy_data_yhi,
               std::unique_ptr<ReadBndryPlanes>& m_r2d)
{
    // Get bndry data
    int bdy_comp = BCVars::RhoQ1_bc_comp;
    Array4<Real> bdatxlo, bdatxhi, bdatylo, bdatyhi;
    if (m_r2d) {
        Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time);
        bdatxlo = (*bndry_data[0])[0].array();
        bdatylo = (*bndry_data[1])[0].array();
        bdatxhi = (*bndry_data[3])[0].array();
        bdatyhi = (*bndry_data[4])[0].array();
    }

    //
    // Note that time (= start_time+old_stage_time)  is measured as total time
    //           start_bdy_time and final_bdy_time are also measured as total time
    //

    // Relaxation constants
    Real F1 = one/(nudge_factor*dt);

    // Domain bounds
    const auto& dom_hi = ubound(domain);
    const auto& dom_lo = lbound(domain);
    auto dx = geom.CellSizeArray();
    auto ProbHi = geom.ProbHiArray();
    auto ProbLo = geom.ProbLoArray();

    // Time interpolation
    Real dT = bdy_time_interval;

    int n_time    = static_cast<int>( (time-start_bdy_time) /  dT);
    int n_time_p1 = n_time + 1;
    Real alpha    = ((time-start_bdy_time) - n_time * dT) / dT;

    // Do not over run the last bdy file
    if (time >= final_bdy_time) {
      n_time    = static_cast<int>( (final_bdy_time - start_bdy_time)/ dT);
      n_time_p1 = n_time;
      alpha     = zero;
    }

    AMREX_ALWAYS_ASSERT( alpha >= zero && alpha <= one);
    Real oma   = one - alpha;

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
                            ng_vect, true);

    // Temporary FABs for storage (owned/filled on all ranks)
    FArrayBox QV_xlo, QV_xhi, QV_ylo, QV_yhi;
    QV_xlo.resize(bx_xlo,1,The_Async_Arena()); QV_xhi.resize(bx_xhi,1,The_Async_Arena());
    QV_ylo.resize(bx_ylo,1,The_Async_Arena()); QV_yhi.resize(bx_yhi,1,The_Async_Arena());

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

    // Get Array4 of interpolated values
    Array4<Real> arr_xlo = QV_xlo.array();  Array4<Real> arr_xhi = QV_xhi.array();
    Array4<Real> arr_ylo = QV_ylo.array();  Array4<Real> arr_yhi = QV_yhi.array();

    Box gtbx = grow(tbx,ng_vect);
    Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
    realbdy_interior_bxs_xy(gtbx, domain, width,
                            tbx_xlo, tbx_xhi,
                            tbx_ylo, tbx_yhi,
                            ng_vect, true);

    // Limiting offset
    int offset = width - 1;

    // Populate with interpolation (protect from ghost cells)
    ParallelFor(tbx_xlo, tbx_xhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_lo.x), dom_lo.x+offset);
        int jj = std::min(std::max(j , dom_lo.y), dom_hi.y       );
        arr_xlo(i,j,k) = (bdatxlo) ? new_cons(i,j,k,Rho_comp) * bdatxlo(ii,jj,k,bdy_comp) :
            new_cons(i,j,k,Rho_comp) * ( oma   * bdatxlo_n  (ii,jj,k)
                                       + alpha * bdatxlo_np1(ii,jj,k) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_hi.x-offset), dom_hi.x);
        int jj = std::min(std::max(j , dom_lo.y       ), dom_hi.y);
        arr_xhi(i,j,k) = (bdatxhi) ? new_cons(i,j,k,Rho_comp) * bdatxhi(ii,jj,k,bdy_comp) :
            new_cons(i,j,k,Rho_comp) * ( oma   * bdatxhi_n  (ii,jj,k)
                                       + alpha * bdatxhi_np1(ii,jj,k) );
    });

    ParallelFor(tbx_ylo, tbx_yhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_lo.x), dom_hi.x       );
        int jj = std::min(std::max(j , dom_lo.y), dom_lo.y+offset);
        arr_ylo(i,j,k) = (bdatylo) ? new_cons(i,j,k,Rho_comp) * bdatylo(ii,jj,k,bdy_comp) :
            new_cons(i,j,k,Rho_comp) * ( oma   * bdatylo_n  (ii,jj,k)
                                       + alpha * bdatylo_np1(ii,jj,k) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int ii = std::min(std::max(i , dom_lo.x       ), dom_hi.x);
        int jj = std::min(std::max(j , dom_hi.y-offset), dom_hi.y);
            jj = std::min(jj, dom_hi.y);
            arr_yhi(i,j,k) = (bdatyhi) ? new_cons(i,j,k,Rho_comp) * bdatyhi(ii,jj,k,bdy_comp) :
                new_cons(i,j,k,Rho_comp) * ( oma   * bdatyhi_n  (ii,jj,k)
                                           + alpha * bdatyhi_np1(ii,jj,k) );
    });

    // Compute RHS in relaxation region
    //==========================================================
    realbdy_interior_bxs_xy(tbx, domain, width,
                            tbx_xlo, tbx_xhi,
                            tbx_ylo, tbx_yhi,
                            ng_vect);
    realbdy_compute_relaxation(RhoQ1_comp, 1,
                               width, dx, ProbLo, ProbHi, F1,
                               tbx_xlo , tbx_xhi , tbx_ylo , tbx_yhi ,
                               arr_xlo , arr_xhi , arr_ylo , arr_yhi ,
                               new_cons, cell_rhs);
} // moist_set_rhs
#endif
