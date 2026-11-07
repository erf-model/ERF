#ifndef ERF_BNDRY_PLANE_INDEX_H
#define ERF_BNDRY_PLANE_INDEX_H

#include <cmath>

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

/**
 * The index along a boundary plane that fills the ghost cell at index i: inside the domain, i
 * itself; outside it, the periodic image when the direction is periodic, so that a ghost cell
 * beside the inflow face holds what the plane holds at the periodic image (as the interior cells
 * across the seam do); else the nearest cell of the plane. A plane holds the domain's cells only
 * along its face, so its ghost cells in the other directions come from here.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int bndry_plane_index (int i, int lo, int hi, bool periodic) noexcept
{
    if (i >= lo && i <= hi) { return i; }
    if (periodic) {
        const int n = hi - lo + 1;
        const int r = (i - lo) % n;
        return lo + ((r < 0) ? r + n : r);
    }
    return (i < lo) ? lo : hi;
}

/**
 * Where a fine index sits among the level-0 plane indices along one tangential direction:
 * the two plane indices to blend and the weight on the second.  The plane holds cell-centred
 * values in every direction (the writer averages the face velocities to cell centres), so a
 * cell-centred index and a nodal (face) index are both located by their positions relative
 * to the plane's cell centres.  At level 0 (ratio 1) a cell-centred index is the plane cell
 * itself with zero weight on its neighbour, and node j is the average of plane cells j - 1
 * and j, the two cells it separates.  Both indices are brought into the plane's extent, which
 * is the level-0 domain: along a periodic direction they take the periodic image, as the
 * interior cells across the seam do (bndry_plane_index); otherwise they are clamped.
 */
struct PlaneIndex { int lo; int hi; amrex::Real w; };

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PlaneIndex plane_index (int i, int ratio, bool nodal, int plo, int phi, bool periodic) noexcept
{
    using amrex::Real;
    const Real x  = (nodal ? Real(i) : Real(i) + Real(0.5)) / Real(ratio) - Real(0.5);
    const int  i0 = static_cast<int>(std::floor(x));
    const Real w  = x - Real(i0);
    return { bndry_plane_index(i0  , plo, phi, periodic),
             bndry_plane_index(i0+1, plo, phi, periodic), w };
}

#endif
