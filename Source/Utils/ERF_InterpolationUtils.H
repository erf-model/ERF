#ifndef ERF_INTERPOLATION_UTILS_H
#define ERF_INTERPOLATION_UTILS_H

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include "ERF_Constants.H"
#include "ERF_TerrainConversion.H"

/**
 * @file ERF_InterpolationUtils.H
 * @brief Utility functions for field interpolation at particle positions
 */

namespace ERF {
namespace Interpolation {

/**
 * @brief Interpolate multiple field values at a particle position using CIC.
 */
template <typename ParticleType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void interpolateFields(
    const ParticleType& p,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
    const amrex::Array4<const amrex::Real>* field_arrays,
    amrex::ParticleReal* result_values,
    int num_fields
) {
    for (int i = 0; i < num_fields; ++i) {
        cic_interpolate(p, plo, dxi, field_arrays[i], &result_values[i], 1);
    }
}

/**
 * @brief Test whether a particle's vertical CIC stencil reads past the local
 *        fab.  cic_interpolate reads cell-centered data at k0+1 and k0+2 with
 *        k0 in {pk-1, pk}, so the maximum read offset is pk+2.  Returns true
 *        if the stencil is out of bounds and the caller should skip the
 *        particle.
 */
template <typename ParticleType, typename Array4Type>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool stencilOutOfBoundsZ(
    const ParticleType& p,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
    const Array4Type& zheight
) noexcept
{
    const int pk = int(amrex::Math::floor((p.pos(AMREX_SPACEDIM-1) - plo[AMREX_SPACEDIM-1])
                                          * dxi[AMREX_SPACEDIM-1]));
    return (pk < zheight.begin[2] || pk + 2 >= zheight.end[2]);
}

} // namespace Interpolation
} // namespace ERF

#endif // ERF_INTERPOLATION_UTILS_H
