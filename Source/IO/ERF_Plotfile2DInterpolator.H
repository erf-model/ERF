/**
 * \file ERF_Plotfile2DInterpolator.H
 */
#ifndef ERF_PLOTFILE2DINTERPOLATOR_H_
#define ERF_PLOTFILE2DINTERPOLATOR_H_

#include <AMReX_MultiFab.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ERF_Plotfile2DSampledLevel.H"

// Vertical sampling helpers for ERF 2D plotfiles. These helpers copy model
// levels or linearly interpolate cell-centered fields onto one target level.
// The 2D writer calls this layer after it allocates output components.

namespace plotfile2d
{

struct SampledBracket
{
    int klo = -1;
    int khi = -1;
    bool found = false;
};

struct SampledWindSources
{
    const amrex::MultiFab* xvel = nullptr;
    const amrex::MultiFab* yvel = nullptr;
    const amrex::MultiFab* zvel = nullptr;
    const amrex::MultiFab* cos_alpha = nullptr;
    const amrex::MultiFab* sin_alpha = nullptr;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sampled_target_is_bracketed (amrex::Real target,
                                  amrex::Real c0,
                                  amrex::Real c1) noexcept
{
    return (c0 <= target && target <= c1) ||
           (c1 <= target && target <= c0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real linear_interpolate (amrex::Real lo_value,
                                amrex::Real hi_value,
                                amrex::Real lo_coord,
                                amrex::Real hi_coord,
                                amrex::Real target) noexcept
{
    if (hi_coord == lo_coord) {
        return lo_value;
    }

    const amrex::Real weight = (target - lo_coord) / (hi_coord - lo_coord);
    return lo_value + weight * (hi_value - lo_value);
}

void fill_sampled_level_component (amrex::MultiFab& dst,
                                   int dst_comp,
                                   const Plotfile2DOutputDescriptor& descriptor,
                                   const amrex::MultiFab& cons,
                                   const amrex::MultiFab* z_phys_cc,
                                   const amrex::MultiFab& z_phys_nd,
                                   bool have_z_phys_cc,
                                   const MoistureComponentIndices& moisture_indices,
                                   int klo,
                                   int khi,
                                   const SampledWindSources& wind_sources = {},
                                   // Base state pressure -- pass this only for an anelastic
                                   //    run, in which case it is used in place of the
                                   //    compressible EOS for the pressure field
                                   const amrex::MultiFab* p_hse = nullptr);

} // namespace plotfile2d

#endif
