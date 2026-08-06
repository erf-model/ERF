#ifndef ERF_TWO_STREAM_SW_H_
#define ERF_TWO_STREAM_SW_H_

#include <AMReX_GpuControl.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

/**
 * @file ERF_TwoStreamSW.H
 * @brief Shortwave (solar) radiation using Beer-Lambert direct-beam model.
 *
 * Implements a simplified, clear-sky shortwave radiation model for Phase 1
 * using the Beer-Lambert law:
 *
 *   F_dir(z) = S0 * cos(zenith) * exp(-tau_cumulative(z) / cos(zenith))
 *
 * where tau_cumulative is the cumulative optical depth from TOA down to level z.
 * For Phase 1 (uniform optical depth per layer):
 *
 *   tau_cumulative(k) = tau_per_layer * k
 *
 * Heating rates are computed from flux divergence:
 *
 *   dF/dz = -d(F_dir)/dz = F_dir(z) * (tau_per_layer / (rho * cp * cos(zenith)))
 *
 * References:
 * -----------
 * - Beer, A., 1852: Bestimmung der Absorption des rothen Lichts in farbigen
 *   Flüssigkeiten. Annalen der Physik und Chemie, 86(78), 78-88.
 *
 * - Bird, R. E., R. L. Hulstrom, and J. A. Lewis, 1984: Terrestrial solar
 *   spectral irradiance: Comparison of proposed standard spectra. J. Climate
 *   Appl. Meteor., 23, 539-545.
 */

/**
 * @brief Compute Beer-Lambert direct-beam flux at a given altitude.
 *
 * Computes the downwelling direct shortwave flux at a single vertical level
 * (or column) using the Beer-Lambert law:
 *
 *   F_dir(z) = S0 * cos(zenith) * exp(-tau_cumulative(z) / cos(zenith))
 *
 * where tau_cumulative is the optical depth from TOA to level z.
 *
 * This is a GPU-safe inline function intended for use within device-side
 * parallel kernels.
 *
 * @tparam D Direction of calculation: -1 for downward (TOA → surface),
 *           +1 for upward (surface → TOA). Phase 1 only uses downward.
 *
 * @param[in] tau_cumulative Cumulative optical depth from TOA to this level [unitless].
 * @param[in] S0 Solar constant at TOA [W/m^2].
 * @param[in] cos_zenith Cosine of solar zenith angle [unitless, 0 ≤ cos_zenith ≤ 1].
 *
 * @return Downwelling direct-beam flux [W/m^2], or 0 if cos_zenith ≤ 0 (night).
 *
 * @note If cos_zenith ≤ 0 (sun below horizon), returns 0.
 * @note Unstable if cos_zenith is exactly 0 (grazing incidence).
 *       At zenith angle ≈ 89°, optical depth effectively approaches infinity.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_sw_direct_flux(amrex::Real tau_cumulative, amrex::Real S0,
                       amrex::Real cos_zenith)
{
  if (cos_zenith <= 0.0) {
    return 0.0;
  }
  return S0 * cos_zenith * std::exp(-tau_cumulative / cos_zenith);
}

/**
 * @brief Compute shortwave heating rate from flux divergence.
 *
 * Given the downwelling direct-beam fluxes at the top and bottom of a layer,
 * computes the radiative heating rate within that layer as:
 *
 *   dF/dz = -[F_dir(k) - F_dir(k+1)] / dz
 *   Q_rad = -dF/dz / (rho * cp)
 *
 * The factor (rho * cp)^{-1} converts radiative flux divergence [W/m^3] into
 * a temperature tendency [K/s].
 *
 * @param[in] flux_top Downwelling flux at the top of the layer [W/m^2].
 * @param[in] flux_bot Downwelling flux at the bottom of the layer [W/m^2].
 * @param[in] dz Vertical thickness of the layer [m]. Must be positive.
 * @param[in] rho Density [kg/m^3]. Must be positive.
 * @param[in] cp Specific heat at constant pressure [J/(kg·K)]. Must be positive.
 *
 * @return Heating rate [K/s]. Positive values indicate warming; negative indicate cooling.
 *
 * @note This function assumes the layer is the control volume with uniform
 *       properties. No sub-layer variation is computed.
 * @note The heating rate formula follows from energy conservation:
 *       dT/dt = -(1/rho/cp) * (∂F_rad/∂z)
 *            = -(1/rho/cp) * (F_out - F_in) / dz
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_sw_heating_rate(amrex::Real flux_top, amrex::Real flux_bot,
                        amrex::Real dz, amrex::Real rho, amrex::Real cp)
{
  if (dz <= 0.0 || rho <= 0.0 || cp <= 0.0) {
    return 0.0;
  }
  // Flux divergence: (W/m^2) / (m) = W/m^3
  amrex::Real flux_divergence = (flux_top - flux_bot) / dz;
  // Heating rate: (W/m^3) / (kg/m^3 * J/(kg*K)) = K/s
  return flux_divergence / (rho * cp);
}

#endif // ERF_TWO_STREAM_SW_H_
