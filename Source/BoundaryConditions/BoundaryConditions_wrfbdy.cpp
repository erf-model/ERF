#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

#ifdef ERF_USE_NETCDF
//
// dest_arr is the Array4 to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals
//     so this follows the BCVars enum
//
void ERFPhysBCFunct::fill_from_wrfbdy (int lev, const Box& bx, const Array4<Real>& dest_arr,
                                       const int icomp, const int bccomp, const int ncomp,
                                       const Box& domain, const BCRec* bc_ptr, const Real time)
{
    //
    // This uses data read in from WRF real.exe output
    // We assume the data has been read in at the level 0 resolution
    //

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    // This is a total HACK
    amrex::Real dT = 1.0;

    int n = time / dT;
    amrex::Real alpha = (time - n * dT) / dT;
    amrex::Real oma   = 1.0 - alpha;

    int ifab;

    int ingested;

    if (bccomp == BCVars::xvel_bc) {
       amrex::Print() << "FILLING FROM WRFBDY: XVEL BC " << icomp << " " << bccomp << " " << ncomp << std::endl;
       ifab = 0; // U
       AMREX_ALWAYS_ASSERT(ncomp == 1);
    } else if (bccomp == BCVars::yvel_bc) {
       amrex::Print() << "FILLING FROM WRFBDY: YVEL BC " << icomp << " " << bccomp << " " << ncomp << std::endl;
       ifab = 1; // V
       AMREX_ALWAYS_ASSERT(ncomp == 1);
    } else if (bccomp == BCVars::zvel_bc) {
       amrex::Print() << "FILLING FROM WRFBDY: ZVEL BC " << icomp << " " << bccomp << " " << ncomp << std::endl;
       ifab = 2; // W
       AMREX_ALWAYS_ASSERT(ncomp == 1);
    } else  if (bccomp == BCVars::RhoTheta_bc_comp) {
       amrex::Print() << "FILLING FROM WRFBDY: CONS BC " << bccomp << " " << bccomp << " " << icomp << " " << ncomp << std::endl;
       ifab = 3; // T
       AMREX_ALWAYS_ASSERT(ncomp == 1);
    } else {
       amrex::Print() << "In fill_from_wrfbdy: icomp = " << icomp << " , bccomp = " << bccomp << " , ncomp = " << ncomp << std::endl;
       amrex::Abort("Don't know this bccomp in fill_from_wrfbdy");
    }

    //
    // We have data at fixed time intervals we will call dT
    // Then to interpolate, given time, we can define n = (time/dT)
    // and alpha = (time - n*dT) / dT, then we define the data at time
    // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)

    //
    const auto& bdatxlo_n   = m_bdy_data_xlo[n  ][ifab].const_array();
    const auto& bdatxlo_np1 = m_bdy_data_xlo[n+1][ifab].const_array();
    const auto& bdatxhi_n   = m_bdy_data_xhi[n  ][ifab].const_array();
    const auto& bdatxhi_np1 = m_bdy_data_xhi[n+1][ifab].const_array();
    const auto& bdatylo_n   = m_bdy_data_ylo[n  ][ifab].const_array();
    const auto& bdatylo_np1 = m_bdy_data_ylo[n+1][ifab].const_array();
    const auto& bdatyhi_n   = m_bdy_data_yhi[n  ][ifab].const_array();
    const auto& bdatyhi_np1 = m_bdy_data_yhi[n+1][ifab].const_array();

    // Fill here all the boundary conditions which are supplied by
    // planes we have read in and are interpolating in time
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
            int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = oma   * bdatxlo_n  (dom_lo.x-1,jb,kb,bccomp+n)
                                    + alpha * bdatxlo_np1(dom_lo.x-1,jb,kb,bccomp+n);
        }
        if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
            int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = oma   * bdatylo_n  (ib,dom_lo.y-1,kb,bccomp+n)
                                    + alpha * bdatylo_np1(ib,dom_lo.y-1,kb,bccomp+n);
        }
        if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
            // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            // dest_arr(i,j,k,icomp+n) = oma   * bdatzlo_n  (ib,jb,dom_lo.z-1,bccomp+n)
            //                         + alpha * bdatzlo_np1(ib,jb,dom_lo.z-1,bccomp+n);
        }
        if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir_ingested) {
            int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = oma   * bdatxhi_n  (dom_hi.x+1,jb,kb,bccomp+n)
                                    + alpha * bdatxhi_np1(dom_hi.x+1,jb,kb,bccomp+n);
        }
        if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir_ingested) {
            int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = oma   * bdatyhi_n  (ib,dom_hi.y+1,kb,bccomp+n)
                                    + alpha * bdatyhi_np1(ib,dom_hi.y+1,kb,bccomp+n);
        }
        if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir_ingested) {
            // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            // dest_arr(i,j,k,icomp+n) = oma   * bdatzhi_n  (ib,jb,dom_hi.z+1,bccomp+n)
            //                         + alpha * bdatzhi_np1(ib,jb,dom_hi.z+1,bccomp+n);
        }

        if (bccomp == BCVars::xvel_bc)
        {
            if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
                int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                dest_arr(i,j,k,icomp+n) = oma * bdatxlo_n  (dom_lo.x-1,jb,kb,bccomp+n)
                                        + alpha * bdatxlo_np1(dom_lo.x-1,jb,kb,bccomp+n);
            }
        }
        if (bccomp == BCVars::yvel_bc)
        {
            if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
                int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                dest_arr(i,j,k,icomp+n) = oma * bdatylo_n  (ib,dom_lo.y-1,kb,bccomp+n)
                                        + alpha * bdatylo_np1(ib,dom_lo.y-1,kb,bccomp+n);
            }
        }
        if (bccomp == BCVars::zvel_bc)
        {
            if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
                // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                // dest_arr(i,j,k,icomp+n) = oma * bdatzlo_n  (ib,jb,dom_lo.z-1,bccomp+n)
                //                         + alpha * bdatzlo_np1(ib,jb,dom_lo.z-1,bccomp+n);
            }
        }

     }); // ParallelFor
}
#endif
