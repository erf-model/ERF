#ifndef ERF_SURFACE_DIAGNOSTIC_SOURCE_H_
#define ERF_SURFACE_DIAGNOSTIC_SOURCE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace surface_diagnostics
{

// Diagnostic-only provenance codes written to 2D plotfiles.
// These numeric values are a user-visible plotfile convention. Renumbering
// them is a breaking output change. The mask is cell-centered and must not
// affect flux calculations, boundary conditions, state updates, or staggered
// stress-face source mixing.
enum class SurfaceDiagnosticSource : int
{
    Missing              = 0,
    SurfaceLayerLand     = 1,
    LSMLand              = 2,
    SurfaceLayerFallback = 3,
    SurfaceLayerSea      = 4,
    Custom               = 5,
    RICO                 = 6
};

/**
 * Convert a diagnostic source enum to its integer code.
 *
 * @param[in] source diagnostic source enum value
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
to_int (SurfaceDiagnosticSource source) noexcept
{
    return static_cast<int>(source);
}

/**
 * Convert a diagnostic source enum to the Real value written to plotfiles.
 *
 * @param[in] source diagnostic source enum value
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
to_plot_value (SurfaceDiagnosticSource source) noexcept
{
    return static_cast<amrex::Real>(to_int(source));
}

/**
 * Classify the provenance of a scalar surface flux.
 *
 * @param[in] is_custom whether the custom flux mode is active
 * @param[in] is_rico whether the RICO flux mode is active
 * @param[in] is_land whether the cell is classified as land
 * @param[in] has_lsm_flux whether an LSM flux array is available
 * @param[in] lsm_flux_is_valid whether the LSM flux value is valid
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SurfaceDiagnosticSource
classify_scalar_source (bool is_custom,
                        bool is_rico,
                        bool is_land,
                        bool has_lsm_flux,
                        bool lsm_flux_is_valid) noexcept
{
    if (is_land && has_lsm_flux && lsm_flux_is_valid) {
        return SurfaceDiagnosticSource::LSMLand;
    }
    if (is_custom) {
        return SurfaceDiagnosticSource::Custom;
    }
    if (is_rico) {
        return SurfaceDiagnosticSource::RICO;
    }
    if (is_land && has_lsm_flux && !lsm_flux_is_valid) {
        return SurfaceDiagnosticSource::SurfaceLayerFallback;
    }
    if (is_land) {
        return SurfaceDiagnosticSource::SurfaceLayerLand;
    }
    return SurfaceDiagnosticSource::SurfaceLayerSea;
}

/**
 * Return the string name for a diagnostic source enum value.
 *
 * @param[in] source diagnostic source enum value
 */
inline const char*
source_name (SurfaceDiagnosticSource source) noexcept
{
    switch (source) {
    case SurfaceDiagnosticSource::Missing:              return "missing";
    case SurfaceDiagnosticSource::SurfaceLayerLand:     return "surface_layer_land";
    case SurfaceDiagnosticSource::LSMLand:              return "lsm_land";
    case SurfaceDiagnosticSource::SurfaceLayerFallback: return "surface_layer_fallback";
    case SurfaceDiagnosticSource::SurfaceLayerSea:      return "surface_layer_sea";
    case SurfaceDiagnosticSource::Custom:               return "custom";
    case SurfaceDiagnosticSource::RICO:                 return "rico";
    }
    return "unknown";
}

} // namespace surface_diagnostics

#endif
