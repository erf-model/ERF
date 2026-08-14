#ifndef ERF_SURFACE_PRECIPITATION_H_
#define ERF_SURFACE_PRECIPITATION_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

// Shared microphysics-to-IO contract for surface precipitation accumulation.
// Scheme-native precipitation accumulators do not all have the same semantics.
// Some schemes store explicit species accumulations; others store total
// precipitation plus frozen-species subset accumulations. The IO layer uses
// these typed slots to build public kg/m^2 fields without double-counting.

/**
 * @brief Precipitation accumulation field and its conversion factor.
 * @var SurfacePrecipAccumulationSource::accum Pointer to the accumulation MultiFab.
 * @var SurfacePrecipAccumulationSource::native_to_kg_m2 Conversion factor from native units to kg/m^2.
 */
struct SurfacePrecipAccumulationSource
{
    const amrex::MultiFab* accum = nullptr;
    amrex::Real native_to_kg_m2 = amrex::Real(1.0);
};

/**
 * @brief Collection of precipitation accumulation sources for different species.
 * @var SurfacePrecipAccumulationSources::total Total precipitation source.
 * @var SurfacePrecipAccumulationSources::rain Rain precipitation source.
 * @var SurfacePrecipAccumulationSources::snow Snow precipitation source.
 * @var SurfacePrecipAccumulationSources::graupel Graupel precipitation source.
 * @var SurfacePrecipAccumulationSources::hail Hail precipitation source.
 */
struct SurfacePrecipAccumulationSources
{
    SurfacePrecipAccumulationSource total{};
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
bool surface_precip_has_total_source (const SurfacePrecipAccumulationSources& src) noexcept
{
    return surface_precip_has_source(src.total);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool surface_precip_has_explicit_rain_source (const SurfacePrecipAccumulationSources& src) noexcept
{
    return surface_precip_has_source(src.rain);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool surface_precip_has_any_source (const SurfacePrecipAccumulationSources& src) noexcept
{
    return surface_precip_has_source(src.total)
        || surface_precip_has_source(src.rain)
        || surface_precip_has_source(src.snow)
        || surface_precip_has_source(src.graupel)
        || surface_precip_has_source(src.hail);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool surface_precip_has_any_frozen_source (const SurfacePrecipAccumulationSources& src) noexcept
{
    // Public frozen precipitation is available whenever the scheme has any
    // precipitation source, because the IO layer can derive the frozen field
    // from the total or species accumulators and fill zero for warm-rain cases.
    return surface_precip_has_any_source(src);
}

#endif
