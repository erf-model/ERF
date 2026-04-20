#ifndef ERF_INTERPOLATION_UTILS_H
#define ERF_INTERPOLATION_UTILS_H

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include "ERF_Constants.H"

/**
 * @file ERF_InterpolationUtils.H
 * @brief Utility functions for field interpolation at particle positions
 */

namespace ERF {
namespace Interpolation {

/**
 * @brief Helper function to interpolate multiple field values at a particle position
 *
 * This function interpolates values from multiple field arrays at a particle's position
 * using the cloud-in-cell (CIC) interpolation method. It supports both periodic and
 * non-periodic z-boundaries.
 *
 * @tparam ParticleType The type of particle
 * @param p The particle
 * @param plo The domain lower bounds
 * @param dxi The inverse of cell spacing
 * @param field_arrays Array of field arrays to interpolate from
 * @param result_values Output array to store interpolated values
 * @param num_fields Number of fields to interpolate
 * @param is_periodic_z Whether the z-boundary is periodic
 * @param zheight Height field array (only used for non-periodic z-boundary)
 */
template <typename ParticleType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void interpolateFields(
    const ParticleType& p,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
    const amrex::Array4<const amrex::Real>* field_arrays,
    amrex::ParticleReal* result_values,
    int num_fields,
    int is_periodic_z = 1,
    const amrex::Array4<amrex::Real>* zheight = nullptr
) {
    // Iterate over each field and interpolate its value at the particle position
    for (int i = 0; i < num_fields; ++i) {
        if (is_periodic_z != 0) {
            cic_interpolate(p, plo, dxi, field_arrays[i], &result_values[i], 1);
        } else {
            cic_interpolate_mapped_z(p, plo, dxi, field_arrays[i], *zheight, &result_values[i], 1);
        }
    }
}

} // namespace Interpolation
} // namespace ERF

#endif // ERF_INTERPOLATION_UTILS_H
