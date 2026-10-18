#ifndef ERF_TERRAIN_FILE_CHECKS_H_
#define ERF_TERRAIN_FILE_CHECKS_H_

/**
 * \file ERF_TerrainFileChecks.H
 *
 * Checks of the height tables read from erf.terrain_file_name, erf.terrain_file_name_USGS
 * and erf.buildings_file_name (ProblemBase::read_custom_terrain), a start-up summary of each
 * table, and the warning for steep terrain next to an inflow or outflow face.
 *
 * Table layouts (nx, ny >= 2):
 *  - ERF text format: x has nx values, y has ny values, and height n sits at
 *    (x[n / ny], y[n % ny]), i.e. y varies fastest.
 *  - USGS format: x has one value per point, y one value per row of nx points, and
 *    height n sits at (x[n], y[n / nx]), i.e. x varies fastest.
 */

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <string>

namespace erf_terrain_file {

/**
 * The first problem found in a height table, or an empty string if there is none.
 * Problems: fewer than two points in x or y, list sizes that do not match nx and ny,
 * a value that is not finite, x or y values that do not increase strictly (in the
 * USGS format: the x values of each row, and the y values from row to row), and a
 * table that does not cover the domain: each end of the x and y lists may stop at
 * most one table spacing inside prob_lo or prob_hi (the layout of a file with one
 * height per cell starting at prob_lo), with an allowance of 1e-3 of that spacing for
 * rounded coordinates.  Beyond the table the reader only extends the heights from its
 * edge: linearly before the first point, and (ERF format) by holding the last value
 * after the last point; USGS extends linearly at both ends.
 *
 * @param[in] is_usgs  USGS layout (x varies fastest, one x per point)
 * @param[in] nx, ny   table size
 * @param[in] x, y, z  table coordinates and heights
 * @param[in] prob_lo, prob_hi  domain extent
 */
std::string
TableProblem (bool is_usgs, int nx, int ny,
              const amrex::Vector<amrex::Real>& x,
              const amrex::Vector<amrex::Real>& y,
              const amrex::Vector<amrex::Real>& z,
              const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& prob_lo,
              const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& prob_hi);

/**
 * Text describing a valid height table: its size and extent, the order its heights are
 * read in, where the lowest and highest heights sit, and the largest slopes |dh/dx| and
 * |dh/dy| between neighbouring table points with where they sit.  Meant to be printed once
 * per file, so that a transposed or mis-ordered file can be spotted at start-up.
 *
 * @param[in] fname    file name, for the first line
 * @param[in] is_usgs  USGS layout
 * @param[in] nx, ny   table size
 * @param[in] x, y, z  table coordinates and heights (must pass TableProblem)
 */
std::string
TableSummary (const std::string& fname, bool is_usgs, int nx, int ny,
              const amrex::Vector<amrex::Real>& x,
              const amrex::Vector<amrex::Real>& y,
              const amrex::Vector<amrex::Real>& z);

/** Terrain slope |grad h| above which a band next to an inflow or outflow face is reported. */
constexpr amrex::Real open_face_slope_warn = amrex::Real(0.3);

/**
 * Physical slope |grad h| of the terrain surface over cell (i,j) from the four surface
 * nodes at k: each of dh/dx and dh/dy is the mean over the cell's two edges in that
 * direction.
 *
 * @param[in] i, j, k   cell indices (k = the surface node level)
 * @param[in] z_nd      nodal heights
 * @param[in] dx, dy    horizontal cell sizes
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
SurfaceSlope (int i, int j, int k, const amrex::Array4<const amrex::Real>& z_nd,
              amrex::Real dx, amrex::Real dy) noexcept
{
    using amrex::Math::abs;
    const amrex::Real hx = amrex::Real(0.5) * ( abs(z_nd(i+1,j  ,k) - z_nd(i,j  ,k))
                                              + abs(z_nd(i+1,j+1,k) - z_nd(i,j+1,k)) ) / dx;
    const amrex::Real hy = amrex::Real(0.5) * ( abs(z_nd(i  ,j+1,k) - z_nd(i  ,j,k))
                                              + abs(z_nd(i+1,j+1,k) - z_nd(i+1,j,k)) ) / dy;
    return std::sqrt(hx*hx + hy*hy);
}

/**
 * Largest SurfaceSlope, on one level, over the cells within band cells of each lateral
 * domain face, in the order xlo, xhi, ylo, yhi; -1 for a face whose band holds no cell
 * of this level.  Reduced over all ranks.
 *
 * @param[in] z_phys_nd nodal heights of the level (surface at the domain's lowest k)
 * @param[in] geom      geometry of the level
 * @param[in] band      width of the band in cells (>= 1)
 */
amrex::Array<amrex::Real,4>
FaceBandSlopes (const amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom, int band);

/**
 * Warning text for the inflow and outflow faces of one level whose band slope (from
 * FaceBandSlopes) exceeds open_face_slope_warn, or an empty string if there is none.
 * The caller prints it only when no steep-terrain option is on at that level.
 *
 * @param[in] lev       level, for the text
 * @param[in] band      band width in cells, for the text
 * @param[in] slopes    band slopes in the order xlo, xhi, ylo, yhi (-1 = no cells)
 * @param[in] open      whether each face is an inflow or outflow face
 * @param[in] bc_names  boundary type of each face, for the text
 */
std::string
OpenFaceSlopeWarning (int lev, int band,
                      const amrex::Array<amrex::Real,4>& slopes,
                      const amrex::Array<bool,4>& open,
                      const amrex::Array<std::string,4>& bc_names);

} // namespace erf_terrain_file

#endif
