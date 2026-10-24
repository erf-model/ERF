#ifndef ERF_PROB_HEIGHTS_H_
#define ERF_PROB_HEIGHTS_H_

/**
 * Physical heights at the face centres of a cell, from the nodal heights z_nd.
 *
 * On a mesh other than ConstantDz, geometry's prob_lo[2], prob_hi[2] and CellSize(2) only
 * bound and step the computational coordinate, so problem setups take their heights from
 * z_nd (and z_cc) there.
 */

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace erf_prob_heights {

//! Height of the centre of x-face (i,j,k)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
z_xface (int i, int j, int k, amrex::Array4<amrex::Real const> const& z_nd) noexcept
{
    return amrex::Real(0.25) * ( z_nd(i,j  ,k) + z_nd(i,j  ,k+1)
                               + z_nd(i,j+1,k) + z_nd(i,j+1,k+1) );
}

//! Height of the centre of y-face (i,j,k)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
z_yface (int i, int j, int k, amrex::Array4<amrex::Real const> const& z_nd) noexcept
{
    return amrex::Real(0.25) * ( z_nd(i  ,j,k) + z_nd(i  ,j,k+1)
                               + z_nd(i+1,j,k) + z_nd(i+1,j,k+1) );
}

//! Height of the centre of z-face (i,j,k)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
z_zface (int i, int j, int k, amrex::Array4<amrex::Real const> const& z_nd) noexcept
{
    return amrex::Real(0.25) * ( z_nd(i,j  ,k) + z_nd(i+1,j  ,k)
                               + z_nd(i,j+1,k) + z_nd(i+1,j+1,k) );
}

} // namespace erf_prob_heights

#endif
