#ifndef ERF_RRTMGP_SURFACE_TEMPERATURE_H_
#define ERF_RRTMGP_SURFACE_TEMPERATURE_H_

#include <ERF_Constants.H>
#include <ERF_EOS.H>

namespace rrtmgp {

/**
 * Resolve the RRTMGP surface-temperature input from ERF's surface sources.
 *
 * RRTMGP and LSM surface-temperature fields contain absolute temperature,
 * while the SurfaceLayer field contains potential temperature.  The latter
 * must therefore be converted with the pressure of the lowest atmospheric
 * cell before it is handed to RRTMGP.
 *
 * @param[in]  is_land                 true when the column is over land
 * @param[in]  has_lsm_t_sfc           true when an LSM t_sfc field exists
 * @param[in]  valid_lsm_t_sfc         true when the LSM value is not a sentinel
 * @param[in]  lsm_t_sfc               LSM t_sfc value, when present
 * @param[in]  has_surface_layer       true when a SurfaceLayer field exists
 * @param[in]  surface_layer_theta     SurfaceLayer potential temperature
 * @param[in]  surface_pressure        pressure in the lowest atmospheric cell
 * @param[in]  default_t_sfc           configured absolute-temperature fallback
 * @param[out] t_sfc                   resolved absolute temperature for RRTMGP
 * @param[in,out] lsm_t_sfc_out        optional LSM fallback/writeback field
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
                              const amrex::Real default_t_sfc,
                              amrex::Real& t_sfc,
                              amrex::Real* lsm_t_sfc_out)
{
    // A valid land-model value is already absolute temperature.  A valid LSM
    // value over water follows the existing fallback chain instead.
    const bool use_lsm_t_sfc = is_land && has_lsm_t_sfc && valid_lsm_t_sfc;

    if (use_lsm_t_sfc) {
        t_sfc = lsm_t_sfc;
    } else if (has_surface_layer) {
        // SurfaceLayer::get_t_surf() returns theta.  RRTMGP requires T.
        t_sfc = getTgivenPandTh(surface_pressure, surface_layer_theta, RdoCp);
    } else {
        t_sfc = default_t_sfc;
    }

    // Preserve the existing fallback/writeback behavior, but write the value
    // in the LSM contract: absolute temperature.
    if (lsm_t_sfc_out != nullptr && !use_lsm_t_sfc) {
        *lsm_t_sfc_out = t_sfc;
    }
}

} // namespace rrtmgp

#endif
