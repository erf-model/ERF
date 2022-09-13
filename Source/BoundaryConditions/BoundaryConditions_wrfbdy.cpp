#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

using namespace amrex;

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


    int ivar;
    // NOTE: we currently only store 1 component for each of bdy_data_xlo etc of each type, so for now we set
    //       bccomp = 0 to use in bdy_data_xlo etc ... but if we start to fill more than one component of cons
    //       we will need to fix this
    int bccomp = 0;

    AMREX_ALWAYS_ASSERT(ncomp == 1);

    if (bccomp_in == BCVars::xvel_bc) {
       ivar = WRFBdyVars::U; // U
    } else if (bccomp_in == BCVars::yvel_bc) {
       ivar = WRFBdyVars::V; // V
    } else  if (bccomp_in == BCVars::Rho_bc_comp) {
       ivar = WRFBdyVars::R; // R
    } else  if (bccomp_in == BCVars::RhoTheta_bc_comp) {
       ivar = WRFBdyVars::T; // T
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
    const auto& bdatxlo_n   = m_bdy_data_xlo[n  ][ivar].const_array();
    const auto& bdatxlo_np1 = m_bdy_data_xlo[n+1][ivar].const_array();
    const auto& bdatxhi_n   = m_bdy_data_xhi[n  ][ivar].const_array();
    const auto& bdatxhi_np1 = m_bdy_data_xhi[n+1][ivar].const_array();
    const auto& bdatylo_n   = m_bdy_data_ylo[n  ][ivar].const_array();
    const auto& bdatylo_np1 = m_bdy_data_ylo[n+1][ivar].const_array();
    const auto& bdatyhi_n   = m_bdy_data_yhi[n  ][ivar].const_array();
    const auto& bdatyhi_np1 = m_bdy_data_yhi[n+1][ivar].const_array();

    // x-faces
    {
        Box bx_xlo(bx);
        bx_xlo.setSmall(0,dom_lo.x); bx_xlo.setBig(0,dom_lo.x);
        bx_xlo.setSmall(1,dom_lo.y); bx_xlo.setBig(1,dom_hi.y);
        bx_xlo.setSmall(2,dom_lo.z); bx_xlo.setBig(2,dom_hi.z);

        Box bx_xhi(bx);
        bx_xhi.setSmall(1,dom_lo.y); bx_xhi.setBig(1,dom_hi.y);
        bx_xhi.setSmall(2,dom_lo.z); bx_xhi.setBig(2,dom_hi.z);
        if (bccomp == BCVars::xvel_bc) {
            bx_xhi.setSmall(0,dom_hi.x+1); bx_xhi.setBig(0,dom_hi.x+1);
        } else {
            bx_xhi.setSmall(0,dom_hi.x); bx_xhi.setBig(0,dom_hi.x);
        }

        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                dest_arr(i,j,k,icomp) = oma   * bdatxlo_n  (i,j,k,bccomp)
                                      + alpha * bdatxlo_np1(i,j,k,bccomp);
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                dest_arr(i,j,k,icomp) = oma   * bdatxhi_n  (i,j,k,bccomp)
                                      + alpha * bdatxhi_np1(i,j,k,bccomp);
            }
        );
    }

    // y-faces
    {
        Box bx_ylo(bx);
        bx_ylo.setSmall(0,dom_lo.x); bx_ylo.setBig(0,dom_hi.x);
        bx_ylo.setSmall(1,dom_lo.y); bx_ylo.setBig(1,dom_lo.y);
        bx_ylo.setSmall(2,dom_lo.z); bx_ylo.setBig(2,dom_hi.z);

        Box bx_yhi(bx);
        bx_yhi.setSmall(0,dom_lo.x); bx_yhi.setBig(0,dom_hi.x);
        bx_yhi.setSmall(2,dom_lo.z); bx_yhi.setBig(2,dom_hi.z);

        if (bccomp == BCVars::yvel_bc) {
            bx_yhi.setSmall(1,dom_hi.y+1); bx_yhi.setBig(1,dom_hi.y+1);
        } else {
            bx_yhi.setSmall(1,dom_hi.y); bx_yhi.setBig(1,dom_hi.y);
        }

        ParallelFor(
           bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                dest_arr(i,j,k,icomp) = oma   * bdatylo_n  (i,j,k,bccomp)
                                      + alpha * bdatylo_np1(i,j,k,bccomp);
            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                dest_arr(i,j,k,icomp) = oma   * bdatyhi_n  (i,j,k,bccomp)
                                      + alpha * bdatyhi_np1(i,j,k,bccomp);
            }
        );
    }
}
#endif
