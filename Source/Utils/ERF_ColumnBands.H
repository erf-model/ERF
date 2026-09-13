#ifndef ERF_COLUMNBANDS_H_
#define ERF_COLUMNBANDS_H_

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

//
// A column integration (e.g. of the hydrostatic base state) runs upward from the bottom of each
//    box.  A box stacked on another box of the same level -- the BoxArray is split in z, e.g.
//    amr.max_grid_size below the number of cells in z -- has to continue from the values the
//    box below it reached.  Such an integration visits the boxes in bands of equal lowest z
//    index, bottom up, and calls fill_below_band before each band but the first.
//

/**
 * The distinct lowest z indices of the boxes of a BoxArray, bottom up.
 *
 * @param[in] ba BoxArray of the level
 * @return lowest z index of every band of boxes, in increasing order
 */
amrex::Vector<int> column_bands (const amrex::BoxArray& ba);

/**
 * Fill the cells just below the boxes of mf whose lowest z index is klo_band, widened laterally
 * by lateral_ng, from the boxes of mf below that band, which must already have been integrated.
 * A cell inside a box below takes that box's value; otherwise a cell in the lateral ghost region
 * of a box below takes the value that box holds there.  A cell with no box of mf below it (e.g.
 * under a refined patch that does not reach the ground) keeps the value it has.
 *
 * @param[in,out] mf         field being integrated; needs lateral_ng ghost cells plus one in z
 * @param[in]     icomp      first component to fill
 * @param[in]     ncomp      number of components to fill
 * @param[in]     klo_band   lowest z index of the band about to be integrated
 * @param[in]     lateral_ng lateral width of the filled slabs; its z entry must be zero
 * @param[in]     geom       geometry of the level (lateral periodicity is honoured)
 */
void fill_below_band (amrex::MultiFab& mf, int icomp, int ncomp, int klo_band,
                      const amrex::IntVect& lateral_ng, const amrex::Geometry& geom);

#endif
