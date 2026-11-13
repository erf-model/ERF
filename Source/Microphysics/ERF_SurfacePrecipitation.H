#ifndef ERF_SURFACE_PRECIPITATION_H_
#define ERF_SURFACE_PRECIPITATION_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
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

/**
 * @brief Total liquid-water-equivalent surface precipitation accumulation [kg/m2]
 *        of every column, from the scheme-native accumulators.
 *
 * Fills accum_2d, a one-component MultiFab on the level's boxes coarsened to one
 * cell in z (box n of accum_2d is box n of the accumulators, with the same
 * owner), from the k = klo plane of the sources. Schemes that provide a total
 * give it directly; schemes that provide the species give their sum, so
 * nothing is counted twice (the same rule as the 2D plotfile precipitation
 * diagnostics). Without any source the field is zero. A box whose 3D box does
 * not hold the k = klo plane (a fine level that does not span z, or the upper
 * box of a column split in z) has no ground to read and stays zero, so no
 * accumulator is read outside its box. Used by the fire layer's
 * fuel moisture (ERF_FirePrecip.H) and the force-restore hydrology
 * (ERF_ForceRestoreHydrology.H).
 */
inline void surface_precip_accum_k0 (amrex::MultiFab& accum_2d,
                                     const SurfacePrecipAccumulationSources& sources,
                                     int klo)
{
    using namespace amrex;
    accum_2d.setVal(Real(0.0));
    if (!surface_precip_has_any_source(sources)) { return; }

    const bool has_total   = surface_precip_has_source(sources.total);
    const bool has_rain    = surface_precip_has_source(sources.rain);
    const bool has_snow    = surface_precip_has_source(sources.snow);
    const bool has_graupel = surface_precip_has_source(sources.graupel);
    const bool has_hail    = surface_precip_has_source(sources.hail);
    const Real f_total   = sources.total.native_to_kg_m2;
    const Real f_rain    = sources.rain.native_to_kg_m2;
    const Real f_snow    = sources.snow.native_to_kg_m2;
    const Real f_graupel = sources.graupel.native_to_kg_m2;
    const Real f_hail    = sources.hail.native_to_kg_m2;

    // Whether box n of a source holds the plane k = klo.
    auto holds_klo = [klo] (const SurfacePrecipAccumulationSource& src, int n) {
        if (!surface_precip_has_source(src)) { return true; }
        const Box& b = src.accum->boxArray()[n];
        return b.smallEnd(2) <= klo && klo <= b.bigEnd(2);
    };

    for (MFIter mfi(accum_2d, false); mfi.isValid(); ++mfi) {
        const int n = mfi.index();
        if (!(holds_klo(sources.total, n) && holds_klo(sources.rain, n) &&
              holds_klo(sources.snow, n) && holds_klo(sources.graupel, n) &&
              holds_klo(sources.hail, n))) { continue; }
        Array4<Real> a = accum_2d.array(mfi);
        Array4<const Real> total   = has_total   ? sources.total.accum->const_array(mfi)   : Array4<const Real>{};
        Array4<const Real> rain    = has_rain    ? sources.rain.accum->const_array(mfi)    : Array4<const Real>{};
        Array4<const Real> snow    = has_snow    ? sources.snow.accum->const_array(mfi)    : Array4<const Real>{};
        Array4<const Real> graupel = has_graupel ? sources.graupel.accum->const_array(mfi) : Array4<const Real>{};
        Array4<const Real> hail    = has_hail    ? sources.hail.accum->const_array(mfi)    : Array4<const Real>{};
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) noexcept
        {
            const int i = iv[0], j = iv[1];
            Real v;
            if (has_total) {
                v = total(i, j, klo) * f_total;
            } else {
                v = Real(0.0);
                if (has_rain)    { v += rain(i, j, klo)    * f_rain; }
                if (has_snow)    { v += snow(i, j, klo)    * f_snow; }
                if (has_graupel) { v += graupel(i, j, klo) * f_graupel; }
                if (has_hail)    { v += hail(i, j, klo)    * f_hail; }
            }
            a(i, j, 0) = v;
        });
    }
}

#endif
