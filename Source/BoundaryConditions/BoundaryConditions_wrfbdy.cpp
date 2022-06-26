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
                                       const int icomp, const int bccomp_in, const int ncomp,
                                       const Box& domain, const BCRec* bc_ptr,
                                       const Real time, const Real dT)
{
    //
    // This uses data read in from WRF real.exe output
    // We assume the data has been read in at the level 0 resolution
    //

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    int n = static_cast<int>(time / dT);
    amrex::Real alpha = (time - n * dT) / dT;
    amrex::Real oma   = 1.0 - alpha;

    int ifab;
    // NOTE: we currently only store 1 component for each of bdy_data_xlo etc of each type, so for now we set
    //       bccomp = 0 to use in bdy_data_xlo etc ... but if we start to fill more than one component of cons
    //       we will need to fix this
    int bccomp = 0;

    AMREX_ALWAYS_ASSERT(ncomp == 1);

    if (bccomp_in == BCVars::xvel_bc) {
       amrex::Print() << "FILLING XVEL BC FROM WRFBDY " << icomp << " " << bccomp_in << " " << ncomp << std::endl;
       ifab = 0; // U
    } else if (bccomp_in == BCVars::yvel_bc) {
       amrex::Print() << "FILLING XVEL BC FROM WRFBDY " << icomp << " " << bccomp_in << " " << ncomp << std::endl;
       ifab = 1; // V
    } else if (bccomp_in == BCVars::zvel_bc) {
       amrex::Print() << "FILLING XVEL BC FROM WRFBDY " << icomp << " " << bccomp_in << " " << ncomp << std::endl;
       ifab = 2; // W
    } else  if (bccomp_in == BCVars::RhoTheta_bc_comp) {
       amrex::Print() << "FILLING RHO THETA BC FROM WRFBDY " << icomp << " " << bccomp_in << " " << ncomp << std::endl;
       ifab = 3; // T
    } else {
       amrex::Print() << "In fill_from_wrfbdy: icomp = " << icomp << " , bccomp = " << bccomp_in << " , ncomp = " << ncomp << std::endl;
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
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (i < dom_lo.x) {
            int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp) = oma   * bdatxlo_n  (dom_lo.x-1,jb,kb,bccomp)
                                  + alpha * bdatxlo_np1(dom_lo.x-1,jb,kb,bccomp);
        }
        if (j < dom_lo.y) {
            int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp) = oma   * bdatylo_n  (ib,dom_lo.y-1,kb,bccomp)
                                  + alpha * bdatylo_np1(ib,dom_lo.y-1,kb,bccomp);
        }
        if (k < dom_lo.z) {
            // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            // dest_arr(i,j,k,icomp) = oma   * bdatzlo_n  (ib,jb,dom_lo.z-1,bccomp)
            //                       + alpha * bdatzlo_np1(ib,jb,dom_lo.z-1,bccomp);
        }
        if (i > dom_hi.x) {
            int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp) = oma   * bdatxhi_n  (dom_hi.x+1,jb,kb,bccomp)
                                  + alpha * bdatxhi_np1(dom_hi.x+1,jb,kb,bccomp);
        }
        if (j > dom_hi.y) {
            int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp) = oma   * bdatyhi_n  (ib,dom_hi.y+1,kb,bccomp)
                                  + alpha * bdatyhi_np1(ib,dom_hi.y+1,kb,bccomp);
        }
        if (k > dom_hi.z) {
            // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            // dest_arr(i,j,k,icomp) = oma   * bdatzhi_n  (ib,jb,dom_hi.z+1,bccomp)
            //                       + alpha * bdatzhi_np1(ib,jb,dom_hi.z+1,bccomp);
        }

        if (bccomp == BCVars::xvel_bc)
        {
            if (i == dom_lo.x) {
                int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                dest_arr(i,j,k,icomp) = oma * bdatxlo_n  (dom_lo.x-1,jb,kb,bccomp)
                                    + alpha * bdatxlo_np1(dom_lo.x-1,jb,kb,bccomp);
            }
        }
        if (bccomp == BCVars::yvel_bc)
        {
            if (j == dom_lo.y) {
                int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                dest_arr(i,j,k,icomp) = oma * bdatylo_n  (ib,dom_lo.y-1,kb,bccomp)
                                    + alpha * bdatylo_np1(ib,dom_lo.y-1,kb,bccomp);
            }
        }
        if (bccomp == BCVars::zvel_bc)
        {
            if (k == dom_lo.z) {
                // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                // dest_arr(i,j,k,icomp) = oma * bdatzlo_n  (ib,jb,dom_lo.z-1,bccomp)
                //                         + alpha * bdatzlo_np1(ib,jb,dom_lo.z-1,bccomp);
            }
        }

     }); // ParallelFor
}
#endif
