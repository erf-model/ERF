#include "ERF.H"

using namespace amrex;

//
// This routine uses data read in as BndryRegisters from a previous ERF run
//
// mf is the MultiFab to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals
//     so this follows the BCVars enum
//
void
ERF::fill_from_bndryregs (MultiFab& mf, const int icomp, const int bccomp, const int ncomp,
                          const Real time)
{
    //
    // We now assume that if we read in on one face, we read in on all faces
    //
    AMREX_ALWAYS_ASSERT(m_r2d);

    int lev = 0;
    const Box& domain = geom[lev].Domain();

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    amrex::Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time);

    const auto& bdatxlo = (*bndry_data[0])[lev].const_array();
    const auto& bdatylo = (*bndry_data[1])[lev].const_array();
    const auto& bdatxhi = (*bndry_data[3])[lev].const_array();
    const auto& bdatyhi = (*bndry_data[4])[lev].const_array();

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Array4<Real>& dest_arr = mf.array(mfi);
        Box bx = mfi.validbox();

        // x-faces
        {
            Box bx_xlo(bx);
            bx_xlo.setSmall(0,dom_lo.x-1); bx_xlo.setBig(0,dom_lo.x-1);
            bx_xlo.setSmall(1,dom_lo.y); bx_xlo.setBig(1,dom_hi.y);
            bx_xlo.setSmall(2,dom_lo.z); bx_xlo.setBig(2,dom_hi.z);
            if (bccomp == BCVars::xvel_bc) {
                bx_xlo.setSmall(0,dom_lo.x); bx_xlo.setBig(0,dom_lo.x);
            } else {
                bx_xlo.setSmall(0,dom_lo.x-1); bx_xlo.setBig(0,dom_lo.x-1);
            }

            Box bx_xhi(bx);
            bx_xhi.setSmall(1,dom_lo.y  ); bx_xhi.setBig(1,dom_hi.y  );
            bx_xhi.setSmall(2,dom_lo.z  ); bx_xhi.setBig(2,dom_hi.z  );
            bx_xhi.setSmall(0,dom_hi.x+1); bx_xhi.setBig(0,dom_hi.x+1);

            ParallelFor(
                bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                    dest_arr(i,j,k,icomp+n) = bdatxlo(dom_lo.x-1,jb,kb,bccomp+n);
                },
                bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                    dest_arr(i,j,k,icomp+n) = bdatxhi(dom_hi.x+1,jb,kb,bccomp+n);
                }
            );
        } // x-faces

        // y-faces
        {
            Box bx_ylo(bx);
            bx_ylo.setSmall(0,dom_lo.x); bx_ylo.setBig(0,dom_hi.x);
            bx_ylo.setSmall(1,dom_lo.y); bx_ylo.setBig(1,dom_lo.y);
            if (bccomp == BCVars::yvel_bc) {
                bx_ylo.setSmall(1,dom_lo.y); bx_ylo.setBig(1,dom_lo.y);
            } else {
                bx_ylo.setSmall(1,dom_lo.y-1); bx_ylo.setBig(1,dom_lo.y-1);
            }

            Box bx_yhi(bx);
            bx_yhi.setSmall(0,dom_lo.x  ); bx_yhi.setBig(0,dom_hi.x);
            bx_yhi.setSmall(2,dom_lo.z  ); bx_yhi.setBig(2,dom_hi.z);
            bx_yhi.setSmall(1,dom_hi.y+1); bx_yhi.setBig(1,dom_hi.y+1);


            ParallelFor(
               bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                    dest_arr(i,j,k,icomp+n) = bdatylo(ib,dom_lo.y-1,kb,bccomp+n);
                },
                bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                    dest_arr(i,j,k,icomp+n) = bdatyhi(ib,dom_hi.y+1,kb,bccomp+n);
                }
            );
        } // y-faces

    } // mf
}
