#ifndef ERF_LSM_INTERP_H_
#define ERF_LSM_INTERP_H_

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <ERF_Constants.H>

/**
 * Sets lsm_undefined in every fine cell whose coarse parent is not a valid LSM value.
 *
 * A land-surface field holds lsm_undefined where the LSM has no value, and its consumers
 * fall back to MOST there.  cell_cons_interp keeps the valid cells smooth, and its MC
 * limiter takes the smaller one-sided difference, so a valid cell next to lsm_undefined
 * gets a slope from its valid side only and moves less than that difference.  An
 * lsm_undefined parent is a local maximum, so it gets zero slope and its fine cells stay
 * exactly lsm_undefined, as long as no neighbour exceeds it.  A neighbour that does
 * (+inf, or memory that was never set) gives a slope of order lsm_undefined, and its fine
 * cells come out as large values below lsm_undefined (about 0.5*lsm_undefined at ratio 2)
 * that pass is_valid_lsm_value.  This marks those cells from a piecewise-constant copy of
 * the coarse field, so a fine cell is either a linearly interpolated value from a valid
 * parent or exactly lsm_undefined.
 *
 * Only the cells that InterpFromCoarseLevel fills with the same nghost and
 * nghost_outside_domain are tested; the others are left as they are.
 */
inline void
lsm_mask_invalid_parents (amrex::MultiFab& fine,
                          const amrex::IntVect& nghost,
                          const amrex::IntVect& nghost_outside_domain,
                          const amrex::MultiFab& crse, int ncomp,
                          const amrex::Geometry& cgeom, const amrex::Geometry& fgeom,
                          const amrex::IntVect& ratio,
                          const amrex::Vector<amrex::BCRec>& bcs, int bcscomp)
{
    // Cells this leaves unfilled keep zero, which counts as valid, so they are not masked
    amrex::MultiFab parent(fine.boxArray(), fine.DistributionMap(), ncomp, nghost);
    parent.setVal(amrex::Real(0));
    amrex::InterpFromCoarseLevel(parent, nghost, nghost_outside_domain,
                                 crse, 0, 0, ncomp, cgeom, fgeom, ratio,
                                 &amrex::pc_interp, bcs, bcscomp);

    for (amrex::MFIter mfi(fine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box gbx = mfi.growntilebox(nghost);
        const auto fine_arr   = fine.array(mfi);
        const auto parent_arr = parent.const_array(mfi);
        amrex::ParallelFor(gbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (!is_valid_lsm_value(parent_arr(i,j,k,n)) ||
                !is_valid_lsm_value(fine_arr(i,j,k,n))) {
                fine_arr(i,j,k,n) = lsm_undefined;
            }
        });
    }
}

/**
 * Interpolates a land-surface field from a coarser level without letting lsm_undefined
 * leak into the valid cells: linear interpolation, then lsm_mask_invalid_parents.
 * Components 0 to ncomp-1 of crse go to components 0 to ncomp-1 of fine.
 */
inline void
lsm_interp_from_coarse (amrex::MultiFab& fine,
                        const amrex::IntVect& nghost,
                        const amrex::IntVect& nghost_outside_domain,
                        const amrex::MultiFab& crse, int ncomp,
                        const amrex::Geometry& cgeom, const amrex::Geometry& fgeom,
                        const amrex::IntVect& ratio,
                        const amrex::Vector<amrex::BCRec>& bcs, int bcscomp)
{
    amrex::InterpFromCoarseLevel(fine, nghost, nghost_outside_domain,
                                 crse, 0, 0, ncomp, cgeom, fgeom, ratio,
                                 &amrex::cell_cons_interp, bcs, bcscomp);
    lsm_mask_invalid_parents(fine, nghost, nghost_outside_domain, crse, ncomp,
                             cgeom, fgeom, ratio, bcs, bcscomp);
}

#endif
