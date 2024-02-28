#include "ERF.H"

using namespace amrex;

/*
 * Impose boundary conditions using data read in as BndryRegisters from a previous ERF run
 *
 * @param[out] mfs  Vector of MultiFabs to be filled
 * @param[in]  time time at which the data should be filled
 */

void
ERF::fill_from_bndryregs (const Vector<MultiFab*>& mfs, const Real time)
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

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5
    const auto& bdatxlo = (*bndry_data[0])[lev].const_array();
    const auto& bdatylo = (*bndry_data[1])[lev].const_array();
    const auto& bdatxhi = (*bndry_data[3])[lev].const_array();
    const auto& bdatyhi = (*bndry_data[4])[lev].const_array();

    int bccomp;

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
    {
        MultiFab& mf = *mfs[var_idx];
        const int icomp = 0;
        const int ncomp = mf.nComp();

        if (var_idx == Vars::xvel) {
           bccomp = BCVars::xvel_bc;
        } else if (var_idx == Vars::yvel) {
           bccomp = BCVars::yvel_bc;
        } else if (var_idx == Vars::zvel) {
           bccomp = BCVars::zvel_bc;
        } else  if (var_idx == Vars::cons) {
           bccomp = BCVars::cons_bc;
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Array4<Real>& dest_arr = mf.array(mfi);
            Box bx = mfi.growntilebox();

            // x-faces
            {
            Box bx_xlo(bx); bx_xlo.setBig(0,dom_lo.x-1);
            if (var_idx == Vars::xvel) bx_xlo.setBig(0,dom_lo.x);

            Box bx_xhi(bx); bx_xhi.setSmall(0,dom_hi.x+1);
            if (var_idx == Vars::xvel) bx_xhi.setSmall(0,dom_hi.x);

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
            Box bx_ylo(bx); bx_ylo.setBig  (1,dom_lo.y-1);
            if (var_idx == Vars::yvel) bx_ylo.setBig(1,dom_lo.y);

            Box bx_yhi(bx); bx_yhi.setSmall(1,dom_hi.y+1);
            if (var_idx == Vars::yvel) bx_yhi.setSmall(1,dom_hi.y);

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
    } // var_idx
}
