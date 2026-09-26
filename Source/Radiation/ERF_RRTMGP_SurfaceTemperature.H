#ifndef ERF_RRTMGP_SURFACE_TEMPERATURE_H_
#define ERF_RRTMGP_SURFACE_TEMPERATURE_H_

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_SurfaceTemperature.H>

namespace rrtmgp {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
inverse_exner (const amrex::Real pressure, const amrex::Real rdOcp)
{
    return one / getExnergivenP(pressure, rdOcp);
}

/**
 * Resolve the RRTMGP surface-temperature input from ERF's surface sources.
 *
 * RRTMGP and LSM surface-temperature fields contain absolute temperature,
 * while the SurfaceLayer field contains potential temperature.  The latter
 * must therefore be converted with the physical surface pressure diagnosed
 * from the lowest atmospheric cell before it is handed to RRTMGP.
 *
 * @param[in]  is_land                 true when the column is over land
 * @param[in]  has_lsm_t_sfc           true when an LSM t_sfc field exists
 * @param[in]  valid_lsm_t_sfc         true when the LSM value is not a sentinel
 * @param[in]  lsm_t_sfc               LSM t_sfc value, when present
 * @param[in]  has_surface_layer       true when a SurfaceLayer field exists
 * @param[in]  surface_layer_theta     SurfaceLayer potential temperature
 * @param[in]  surface_pressure        physical surface pressure diagnosed from
 *                                     the lowest atmospheric cell
 * @param[in]  rdOcp                   configured Rd/cp exponent
 * @param[in]  default_t_sfc           configured absolute-temperature fallback
 * @param[out] t_sfc                   resolved absolute temperature for RRTMGP
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
resolve_surface_temperature (const bool is_land,
                              const bool has_lsm_t_sfc,
                              const bool valid_lsm_t_sfc,
                              const amrex::Real lsm_t_sfc,
                              const bool has_surface_layer,
                              const amrex::Real surface_layer_theta,
                              const amrex::Real surface_pressure,
                              const amrex::Real rdOcp,
                              const amrex::Real default_t_sfc,
                              amrex::Real& t_sfc)
{
    // A valid land-model value is already absolute temperature.  A valid LSM
    // value over water follows the existing fallback chain instead.
    const bool use_lsm_t_sfc = is_land && has_lsm_t_sfc && valid_lsm_t_sfc;

    if (use_lsm_t_sfc) {
        t_sfc = lsm_t_sfc;
    } else if (has_surface_layer &&
               erf_surface_temperature::theta_to_temperature(
                   surface_layer_theta, surface_pressure, rdOcp, t_sfc)) {
        // SurfaceLayer::get_t_surf() returns theta.  RRTMGP requires T.
        // The helper also rejects an invalid pressure or nonphysical theta.
    } else {
        t_sfc = default_t_sfc;
    }

}

} // namespace rrtmgp

#endif
