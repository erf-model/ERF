#ifndef ERF_PLOTFILE2DINTERPOLATOR_H_
#define ERF_PLOTFILE2DINTERPOLATOR_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

// Vertical sampling helpers for ERF 2D plotfiles. These helpers copy model
// levels or linearly interpolate cell-centered fields onto one target level.

namespace plotfile2d
{

struct SampledBracket
{
    int klo = -1;
    int khi = -1;
    bool found = false;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SampledBracket find_sampled_bracket (amrex::Real target,
                                     amrex::Real lo_value,
                                     amrex::Real hi_value) noexcept
{
    SampledBracket bracket;
    if ((lo_value <= target && target <= hi_value) ||
        (hi_value <= target && target <= lo_value)) {
        bracket.klo = 0;
        bracket.khi = 1;
        bracket.found = true;
    }
    return bracket;
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

} // namespace plotfile2d

#endif
