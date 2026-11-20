#ifndef ERF_ACTUATOR_GEOMETRY_H
#define ERF_ACTUATOR_GEOMETRY_H

/** Face heights and face volumes of the staggered grid, shared by the actuator sampler and spreader. */

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace erf_actuator {

/**
 * Return the height (m) of the centre of face (i,j,k) of the staggered grid dir (0: x face,
 * 1: y face, 2: z face): the mean of the four mesh nodes of that face on a stretched or
 * terrain-following mesh (znd, has_znd true), or the nominal height on a uniform-dz mesh:
 * zlo + k dz for a z face (z face k is the bottom face of cell k), zlo + (k + 1/2) dz otherwise.
 * Shared by the velocity sampling and the force spreading so that both see the same face positions.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real face_height (int dir, int i, int j, int k,
                         amrex::Array4<amrex::Real const> const& znd, bool has_znd,
                         amrex::Real zlo, amrex::Real dz)
{
    using amrex::Real;
    if (!has_znd) {
        return (dir == 2) ? zlo + k * dz : zlo + (k + Real(0.5)) * dz;
    }
    if (dir == 0) {
        return Real(0.25) * (znd(i,j,k) + znd(i,j+1,k) + znd(i,j,k+1) + znd(i,j+1,k+1));
    }
    if (dir == 1) {
        return Real(0.25) * (znd(i,j,k) + znd(i+1,j,k) + znd(i,j,k+1) + znd(i+1,j,k+1));
    }
    return Real(0.25) * (znd(i,j,k) + znd(i+1,j,k) + znd(i,j+1,k) + znd(i+1,j+1,k));
}

/**
 * Return the volume (m^3) attached to face (i,j,k) of the staggered grid dir: half of each
 * neighbouring cell, with the cell volumes dx dy dz detJ on a stretched or terrain-following mesh
 * (detJ = 1 otherwise). Reads detj at i-1 (x), j-1 (y) or k-1 (z): detJ needs one ghost cell at a
 * box's low faces.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real face_volume (int dir, int i, int j, int k,
                         amrex::Array4<amrex::Real const> const& detj, bool has_detj,
                         amrex::Real dxdydz)
{
    using amrex::Real;
    if (!has_detj) { return dxdydz; }
    if (dir == 0) { return Real(0.5) * (detj(i-1,j,k) + detj(i,j,k)) * dxdydz; }
    if (dir == 1) { return Real(0.5) * (detj(i,j-1,k) + detj(i,j,k)) * dxdydz; }
    return Real(0.5) * (detj(i,j,k-1) + detj(i,j,k)) * dxdydz;
}

} // namespace erf_actuator

#endif
