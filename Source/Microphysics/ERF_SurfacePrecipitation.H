#ifndef ERF_SURFACE_PRECIPITATION_H_
#define ERF_SURFACE_PRECIPITATION_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

// Shared microphysics-to-IO contract for surface precipitation accumulation.
// Scheme-native accumulators may be stored in depth-like restart units, so the
// output layer applies the conversion factor to a common kg/m^2 basis.

struct SurfacePrecipAccumulationSource
{
    const amrex::MultiFab* accum = nullptr;
    amrex::Real native_to_kg_m2 = amrex::Real(1.0);
};

struct SurfacePrecipAccumulationSources
{
    SurfacePrecipAccumulationSource rain{};
    SurfacePrecipAccumulationSource snow{};
    SurfacePrecipAccumulationSource graupel{};
    SurfacePrecipAccumulationSource hail{};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool surface_precip_has_source (const SurfacePrecipAccumulationSource& src) noexcept
{
    return src.accum != nullptr;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool surface_precip_has_any_source (const SurfacePrecipAccumulationSources& src) noexcept
{
    return surface_precip_has_source(src.rain)
        || surface_precip_has_source(src.snow)
        || surface_precip_has_source(src.graupel)
        || surface_precip_has_source(src.hail);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool surface_precip_has_any_frozen_source (const SurfacePrecipAccumulationSources& src) noexcept
{
    return surface_precip_has_source(src.snow)
        || surface_precip_has_source(src.graupel)
        || surface_precip_has_source(src.hail);
}

#endif
