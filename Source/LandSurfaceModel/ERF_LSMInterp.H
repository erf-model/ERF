#ifndef ERF_LSM_INTERP_H_
#define ERF_LSM_INTERP_H_

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <ERF_Constants.H>

/**
 * Corrects a land-surface field that cell_cons_interp has filled from a coarser level so
 * that no lsm_undefined coarse cell affects a valid fine cell.
 *
 * A land-surface field holds lsm_undefined where the LSM has no value, and its consumers
 * fall back to MOST there.  cell_cons_interp builds each coarse cell's slopes, and limits
 * them, from its 3x3(x3) neighbourhood in the refined directions.  An lsm_undefined
 * neighbour makes both of these meaningless: the slope from the valid side alone can
 * still push a fine value past every valid coarse value, and a neighbour above
 * lsm_undefined (+inf, or memory that was never set) can turn an lsm_undefined parent
 * into large values that pass is_valid_lsm_value.  So each fine cell becomes
 *   - lsm_undefined, if its parent or its interpolated value is not valid;
 *   - its parent's value (piecewise constant), if any coarse cell in that stencil is not
 *     valid;
 *   - the linearly interpolated value, which is bounded by the valid stencil, otherwise.
 * Stencil cells outside a non-periodic domain boundary are not tested, since
 * cell_cons_interp takes one-sided slopes there.
 *
 * The parents come from a single ParallelCopy of crse onto the coarsened fine layout.
 * Only the cells that InterpFromCoarseLevel fills with the same nghost and
 * nghost_outside_domain are tested; the others are left as they are.
 */
inline void
lsm_mask_invalid_parents (amrex::MultiFab& fine,
                          const amrex::IntVect& nghost,
                          const amrex::IntVect& nghost_outside_domain,
                          const amrex::MultiFab& crse, int ncomp,
                          const amrex::Geometry& cgeom, const amrex::Geometry& fgeom,
                          const amrex::IntVect& ratio)
{
    // Coarse cells under the fine cells to be tested, plus the interpolation stencil
    amrex::IntVect stencil(0), cng(0);
    for (int d = 0; d < 3; ++d) {
        stencil[d] = (ratio[d] > 1) ? 1 : 0;
        cng[d] = (nghost[d] + ratio[d] - 1) / ratio[d] + stencil[d];
    }

    // Coarse cells this leaves unfilled inside the domain are lsm_undefined, so the fine
    // cells next to them are not trusted
    amrex::MultiFab parent(amrex::coarsen(fine.boxArray(), ratio), fine.DistributionMap(),
                           ncomp, cng);
    parent.setVal(lsm_undefined);
    parent.ParallelCopy(crse, 0, 0, ncomp, amrex::IntVect(0), cng, cgeom.periodicity());

    // These fields are planar, so only their own z extent can hold stencil cells
    const amrex::Box crse_extent = crse.boxArray().minimalBox();
    amrex::Box cdomain = cgeom.Domain();
    cdomain.setRange(2, crse_extent.smallEnd(2), crse_extent.length(2));
    amrex::Box fdomain = fgeom.Domain();
    for (int d = 0; d < 3; ++d) {
        if (cgeom.isPeriodic(d)) { cdomain.grow(d, cng[d]); }
        fdomain.grow(d, fgeom.isPeriodic(d) ? nghost[d] : nghost_outside_domain[d]);
    }

    for (amrex::MFIter mfi(fine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box gbx = mfi.growntilebox(nghost) & fdomain;
        if (gbx.isEmpty()) { continue; }
        const auto fine_arr   = fine.array(mfi);
        const auto parent_arr = parent.const_array(mfi);
        amrex::ParallelFor(gbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const amrex::IntVect ic = amrex::coarsen(amrex::IntVect(i,j,k), ratio);
            const amrex::Real p = parent_arr(ic, n);
            if (!is_valid_lsm_value(p) || !is_valid_lsm_value(fine_arr(i,j,k,n))) {
                fine_arr(i,j,k,n) = lsm_undefined;
                return;
            }
            bool stencil_valid = true;
            for (int koff = -stencil[2]; koff <= stencil[2]; ++koff) {
            for (int joff = -stencil[1]; joff <= stencil[1]; ++joff) {
            for (int ioff = -stencil[0]; ioff <= stencil[0]; ++ioff) {
                const amrex::IntVect c = ic + amrex::IntVect(ioff, joff, koff);
                if (cdomain.contains(c) && !is_valid_lsm_value(parent_arr(c, n))) {
                    stencil_valid = false;
                }
            }}}
            if (!stencil_valid) {
                fine_arr(i,j,k,n) = p;
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
                             cgeom, fgeom, ratio);
}

#endif
