#ifndef ERF_RADIATION_PARAMS_H_
#define ERF_RADIATION_PARAMS_H_

#include <AMReX_GpuControl.H>

/**
 * @file ERF_RadiationParams.H
 * @brief Parameters for two-stream radiation solvers.
 *
 * This file contains optical properties and parameters for the two-stream
 * radiation model. Implements constant, uniform optical properties
 * (Beer-Lambert SW and gray-gas LW). Future work will extend this to include
 * vertical profiles, clouds, scattering, etc.
 *
 * References:
 * -----------
 * - Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam, 1989:
 *   Rapid calculation of radiative heating rates and photodissociation rates
 *   in inhomogeneous multiple scattering atmospheres. J. Geophys. Res., 94,
 *   16387-16405. https://doi.org/10.1029/JD094iD13p16387
 *
 * - Beer, A., 1852: Bestimmung der Absorption des rothen Lichts in farbigen
 *   Flüssigkeiten. Annalen der Physik und Chemie, 86(78), 78-88.
 *
 * - Kirchhoff, G., 1860: Ueber den Zusammenhang zwischen den
 *   Emissionsvermögen und den Absorptionsvermögen der Körper für Wärmestrahlung.
 *   Monatsberichte der Akademie der Wissenschaften zu Berlin, 783-787.
 */

/**
 * @class RadiationParams
 * @brief Container for optical properties in the two-stream model.
 *
 * Provides constant, uniform optical properties across space and time.
 * Future work will add:
 * - Per-layer vertical profiles of optical depth
 * - Cloud optical properties and coverage
 * - Scattering effects
 * - Gaseous absorption bands
 *
 * All quantities are marked as placeholders to guide future development.
 */
struct RadiationParams
{
  /**
   * @brief Shortwave optical depth per layer (unitless).
   *
   * Used in Beer-Lambert direct-beam attenuation:
   *   F_dir(k) = S0 * cos(zenith) * exp(-tau_sw_cumulative(k) / cos(zenith))
   *
   * Uniform constant value for all layers.
   * Future work: vertically varying profile based on clouds/aerosols.
   */
  amrex::Real tau_sw = 0.05;

  /**
   * @brief Longwave optical depth per layer (unitless).
   *
   * Used in gray-gas LW two-stream sweeps:
   *   F_up(k+1) = F_up(k) * exp(-tau_lw_k) + sigma * T(k)^4 * (1 - exp(-tau_lw_k))
   *
   * Uniform constant value for all layers.
   * Future work: vertically varying profile based on humidity/clouds.
   */
  amrex::Real tau_lw = 1.0;

  /**
   * @brief Solar zenith angle in radians.
   *
   * Relates to solar zenith angle input (in degrees) via:
   *   cos_zenith = cos(zenith_rad)
   *
   * This is used to scale the Beer-Lambert optical depth and to weight
   * the direct-beam flux.
   */
  amrex::Real cos_zenith = 0.5;  // Default to 60° zenith angle

  /**
   * @brief Solar constant (top-of-atmosphere irradiance) in W/m^2.
   *
   * Typical value: 1361 W/m^2 (TSI).
   * Constant global value.
   * Future work: account for solar cycle variability.
   */
  amrex::Real S0 = 1361.0;

  /**
   * @brief Stefan-Boltzmann constant (W/(m^2·K^4)).
   *
   * Fundamental constant for thermal radiation emission:
   *   F_rad = sigma * T^4
   *
   * Value: 5.670374419e-8 W/(m^2·K^4).
   * This is typically provided by ERF_Constants.H.
   */
  amrex::Real sigma = 5.670374419e-8;

  /**
   * @brief Constructor with default parameters.
   *
   * Initializes optical properties to sensible defaults.
   * Parameters can be overridden via the RadChoice initialization.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  RadiationParams()
    : tau_sw(0.05), tau_lw(1.0), cos_zenith(0.5),
      S0(1361.0), sigma(5.670374419e-8)
  {}

  /**
   * @brief Constructor with explicit parameters.
   *
   * @param[in] tau_sw_in Shortwave optical depth per layer.
   * @param[in] tau_lw_in Longwave optical depth per layer.
   * @param[in] cos_zenith_in Cosine of solar zenith angle.
   * @param[in] S0_in Solar constant in W/m^2.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  RadiationParams(amrex::Real tau_sw_in, amrex::Real tau_lw_in,
                  amrex::Real cos_zenith_in, amrex::Real S0_in)
    : tau_sw(tau_sw_in), tau_lw(tau_lw_in), cos_zenith(cos_zenith_in),
      S0(S0_in), sigma(5.670374419e-8)
  {}
};

#endif // ERF_RADIATION_PARAMS_H_
