#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

using namespace amrex;

//
// dest_arr is the Array4 to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals
//     so this follows the BCVars enum
//
void ERFPhysBCFunct::fill_from_bndryregs (int lev, const Box& bx, const Array4<Real>& dest_arr,
                                          const int icomp, const int bccomp, const int ncomp,
                                          const Box& domain, const BCRec* bc_ptr, const Real time)
{
    // This uses data read in as BndryRegisters from a previous ERF run
    AMREX_ALWAYS_ASSERT(m_r2d);

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    amrex::Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time);

    const auto& bdatxlo = (*bndry_data[0])[lev].const_array();
    const auto& bdatylo = (*bndry_data[1])[lev].const_array();
    const auto& bdatxhi = (*bndry_data[3])[lev].const_array();
    const auto& bdatyhi = (*bndry_data[4])[lev].const_array();

    // const auto& bdatzlo = (*bndry_data[2])[lev].const_array();
    // const auto& bdatzhi = (*bndry_data[5])[lev].const_array();

    // Fill here all the boundary conditions which are supplied by
    // planes we have read in and are interpolating in time
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
            int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = bdatxlo(dom_lo.x-1,jb,kb,bccomp+n);
        }
        if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
            int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = bdatylo(ib,dom_lo.y-1,kb,bccomp+n);
        }
        if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
            // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            // dest_arr(i,j,k,icomp+n) = bdatzlo(ib,jb,dom_lo.z-1,bccomp+n);
        }
        if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir_ingested) {
            int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = bdatxhi(dom_hi.x+1,jb,kb,bccomp+n);
        }
        if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir_ingested) {
            int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
            dest_arr(i,j,k,icomp+n) = bdatyhi(ib,dom_hi.y+1,kb,bccomp+n);
        }
        if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir_ingested) {
            // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
            // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
            // dest_arr(i,j,k,icomp+n) = bdatzhi(ib,jb,dom_hi.z+1,bccomp+n);
        }

        if (bccomp == BCVars::xvel_bc)
        {
            if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
                int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                dest_arr(i,j,k,icomp+n) = bdatxlo(dom_lo.x-1,jb,kb,bccomp+n);
            }
        } // xvel
        if (bccomp == BCVars::yvel_bc)
        {
            if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
                int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                dest_arr(i,j,k,icomp+n) = bdatylo(ib,dom_lo.y-1,kb,bccomp+n);
            }
        } // yvel
        if (bccomp == BCVars::zvel_bc)
        {
            if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
                // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                // dest_arr(i,j,k,icomp+n) = bdatzlo(ib,jb,dom_lo.z-1,bccomp+n);
            }
        } // zvel

    }); // ParallelFor
}
