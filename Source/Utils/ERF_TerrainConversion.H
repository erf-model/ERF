#ifndef ERF_TERRAIN_CONVERSION_H_
#define ERF_TERRAIN_CONVERSION_H_

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

// Particle position storage convention:
//  * pos(0), pos(1): physical x, y (no horizontal mapping in ERF).
//  * pos(AMREX_SPACEDIM-1): computational vertical coordinate zeta.
//
// For uniform-z grids the computational zeta equals the physical z, so the
// helpers below are identity transforms.
//
// For terrain-fitted grids the physical z at a column (x, y) is bilinearly
// interpolated from z_phys_nd(i, j, k_face).

namespace ERF {
namespace ParticlePos {

// Bilinear interpolation of the height array at (x, y) for a given vertical
// node index k_face.  Returns the physical z of node k_face at column (x, y).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
z_face_at_xy (amrex::Real x, amrex::Real y, int k_face,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& plo,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxi,
              amrex::Array4<amrex::Real const> const& height_arr) noexcept
{
    const amrex::Real lx = (x - plo[0]) * dxi[0];
    const amrex::Real ly = (y - plo[1]) * dxi[1];
    const int  ix = static_cast<int>(amrex::Math::floor(lx));
    const int  iy = static_cast<int>(amrex::Math::floor(ly));
    const amrex::Real fx = lx - static_cast<amrex::Real>(ix);
    const amrex::Real fy = ly - static_cast<amrex::Real>(iy);
    return height_arr(ix  , iy  , k_face) * (amrex::Real(1) - fx) * (amrex::Real(1) - fy)
         + height_arr(ix+1, iy  , k_face) *                fx     * (amrex::Real(1) - fy)
         + height_arr(ix  , iy+1, k_face) * (amrex::Real(1) - fx) *                fy
         + height_arr(ix+1, iy+1, k_face) *                fx     *                fy;
}

// Computational zeta -> physical z at column (x, y).
//   zeta is in the same units as the cell-centred z coordinate (i.e. it ranges
//   over the level's z domain box scaled by dxi[2]).  The integer cell index
//   k = floor((zeta - plo_z) * dxi_z) is the same as for a uniform-z grid;
//   the sub-cell fraction is interpolated linearly between the two terrain
//   nodes that bound the cell.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
z_from_zeta (amrex::Real x, amrex::Real y, amrex::Real zeta,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& plo,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxi,
             amrex::Array4<amrex::Real const> const& height_arr) noexcept
{
    if (!height_arr) {
        return zeta;  // no terrain: identity
    }
    const amrex::Real lz = (zeta - plo[AMREX_SPACEDIM-1]) * dxi[AMREX_SPACEDIM-1];
    const int  k  = static_cast<int>(amrex::Math::floor(lz));
    const amrex::Real fz = lz - static_cast<amrex::Real>(k);
    const amrex::Real z_lo = z_face_at_xy(x, y, k    , plo, dxi, height_arr);
    const amrex::Real z_hi = z_face_at_xy(x, y, k + 1, plo, dxi, height_arr);
    return (amrex::Real(1) - fz) * z_lo + fz * z_hi;
}

// Physical z -> computational zeta at column (x, y).
//   Walks the terrain nodes at (x, y) to find the cell whose physical z range
//   straddles the requested z, then linearly interpolates the sub-cell
//   fraction.  k_max is the largest valid cell index in the level.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
zeta_from_z (amrex::Real x, amrex::Real y, amrex::Real z,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& plo,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxi,
             amrex::Array4<amrex::Real const> const& height_arr,
             int k_max) noexcept
{
    if (!height_arr) {
        return z;  // no terrain: identity
    }
    // Linear walk to find k where z_face[k] <= z < z_face[k+1]
    int k = 0;
    while (k < k_max) {
        const amrex::Real z_hi = z_face_at_xy(x, y, k + 1, plo, dxi, height_arr);
        if (z < z_hi) { break; }
        ++k;
    }
    const amrex::Real z_lo = z_face_at_xy(x, y, k    , plo, dxi, height_arr);
    const amrex::Real z_hi = z_face_at_xy(x, y, k + 1, plo, dxi, height_arr);
    const amrex::Real fz = (z_hi > z_lo)
                         ? (z - z_lo) / (z_hi - z_lo)
                         : amrex::Real(0);
    return plo[AMREX_SPACEDIM-1]
         + (static_cast<amrex::Real>(k) + fz) / dxi[AMREX_SPACEDIM-1];
}

} // namespace ParticlePos
} // namespace ERF

#endif
