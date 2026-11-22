#ifndef ERF_LENGTH_SCALE_CELL_SIZE_H_
#define ERF_LENGTH_SCALE_CELL_SIZE_H_

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

/**
 * Cell sizes for formulas that use a lateral cell size as a length scale rather than to
 * divide a difference: LES filter widths, scale-aware parameterizations, drag rates and
 * per-area densities.
 *
 * On a domain with one cell in y and more than one in x (an x-z run), dy is arbitrary, so
 * dx is used in its place; likewise dy replaces dx on a domain with one cell in x.  The
 * formula then keeps its 3-D form, as if the unresolved direction had the spacing of the
 * resolved one.  A single column (one cell in both) keeps the dx and dy of the inputs as
 * its nominal column size.
 *
 * The level's own domain decides: a lateral direction with one cell is not refined
 * (amr.ref_ratio_vect = r 1 r), so it has one cell on every level, and the spacing copied
 * is that of the level.
 *
 * @param[in] geom geometry of the level
 * @return dx, dy and dz with the substitution above
 */
inline amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
LengthScaleCellSize (const amrex::Geometry& geom)
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const amrex::Box& domain = geom.Domain();
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    if (ny == 1 && nx > 1) {
        dx[1] = dx[0];
    } else if (nx == 1 && ny > 1) {
        dx[0] = dx[1];
    }
    return dx;
}

/**
 * Inverse of LengthScaleCellSize.
 *
 * @param[in] geom geometry of the level
 * @return 1/dx, 1/dy and 1/dz with the substitution of LengthScaleCellSize
 */
inline amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
LengthScaleInvCellSize (const amrex::Geometry& geom)
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const amrex::Box& domain = geom.Domain();
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    if (ny == 1 && nx > 1) {
        dxinv[1] = dxinv[0];
    } else if (nx == 1 && ny > 1) {
        dxinv[0] = dxinv[1];
    }
    return dxinv;
}

#endif
