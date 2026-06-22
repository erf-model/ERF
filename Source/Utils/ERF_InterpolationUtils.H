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

} // namespace Interpolation
} // namespace ERF

#endif // ERF_INTERPOLATION_UTILS_H
