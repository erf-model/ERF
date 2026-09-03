#ifndef ERF_TWO_STREAM_LW_H_
#define ERF_TWO_STREAM_LW_H_

#include <AMReX_GpuControl.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

/**
 * @file ERF_TwoStreamLW.H
 * @brief Longwave (thermal) radiation using gray-gas two-stream model.
 *
 * Implements a simplified, clear-sky longwave radiation model
 * using the gray-gas two-stream approximation (Toon et al. 1989):
 *
 * Vertical orientation follows ERF: the vertical index k increases upward,
 * so layer k = kmin touches the surface and k = kmax the top of the domain.
 * Fluxes live on layer interfaces m = k - kmin (bottom of layer k) and
 * m + 1 (top of layer k); m = 0 is the surface, m = nlev the TOA.
 *
 * UPWARD SWEEP (from surface to TOA):
 * -----------
 *   F_up(m+1) = F_up(m) * exp(-tau_lw) + sigma * T(k)^4 * (1 - exp(-tau_lw))
 *
 *   Initial condition: F_up(0) = emissivity * sigma * T_surface^4
 *                            + (1 - emissivity) * F_down(0)
 *   (gray surface: emitted plus reflected downwelling flux, so the
 *   downward sweep must be completed before the upward sweep starts)
 *
 * DOWNWARD SWEEP (from TOA to surface):
 * ---------------
 *   F_down(m) = F_down(m+1) * exp(-tau_lw) + sigma * T(k)^4 * (1 - exp(-tau_lw))
 *
 *   Initial condition: F_down(nlev) = 0 (no incoming from space)
 *
 * NET FLUX AND HEATING:
 * ---------------------
 *   F_net(m) = F_up(m) - F_down(m)          (positive upward)
 *   Q_lw(k)  = -(1 / (rho * cp)) * [F_net(m+1) - F_net(m)] / dz(k)   [K/s]
 *
 *   A layer that emits more than it absorbs has F_net increasing with
 *   height, so Q_lw < 0 (radiative cooling).
 *
 * Layer temperature T(k) is the absolute temperature obtained from
 * rho*theta through the Exner function (see get_temperature_from_rhotheta()
 * in ERF_AdvanceTwoStreamRadiation.cpp), not the potential temperature.
 *
 * MODEL SIMPLIFICATIONS:
 * -------------------------
 * - Uniform optical depth (tau_lw) per layer
 * - No scattering
 * - No clouds
 * - Isothermal test mode available: if isothermal_test=true, force all T = T_iso_K
 *   This gives F_up = F_down = sigma * T_iso_K^4 everywhere, zero net flux.
 *
 * References:
 * -----------
 * - Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam, 1989:
 *   Rapid calculation of radiative heating rates and photodissociation rates
 *   in inhomogeneous multiple scattering atmospheres. J. Geophys. Res., 94,
 *   16387-16405. https://doi.org/10.1029/JD094iD13p16387
 *
 * - Kirchhoff, G., 1860: Ueber den Zusammenhang zwischen den
 *   Emissionsvermögen und den Absorptionsvermögen der Körper für Wärmestrahlung.
 *   Monatsberichte der Akademie der Wissenschaften zu Berlin, 783-787.
 */

/**
 * @brief Compute thermal (LW) intensity for a given temperature.
 *
 * Using the Stefan-Boltzmann law, compute the upwelling or downwelling
 * thermal radiation intensity:
 *
 *   I_rad = sigma * T^4
 *
 * where sigma is the Stefan-Boltzmann constant and T is absolute temperature.
 *
 * This is a GPU-safe inline function.
 *
 * @param[in] T Absolute temperature [K]. Must be positive.
 * @param[in] sigma Stefan-Boltzmann constant [W/(m^2·K^4)].
 *
 * @return Radiative intensity [W/m^2]. Always non-negative.
 *
 * @note If T ≤ 0, returns 0 (unphysical, but defensive).
 * @note For typical terrestrial temperatures (200-400 K), this ranges
 *       from ~1 W/m^2 to ~1500 W/m^2.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_thermal_intensity(amrex::Real T, amrex::Real sigma)
{
  if (T <= 0.0) {
    return 0.0;
  }
  return sigma * T * T * T * T;
}

/**
 * @brief Compute transmission coefficient through an absorbing layer.
 *
 * For a layer with optical depth tau_lw, the fraction of incident radiation
 * that is transmitted (not absorbed or scattered) is:
 *
 *   transmit = exp(-tau_lw)
 *
 * The fraction absorbed/scattered is:
 *
 *   absorb = 1 - exp(-tau_lw)
 *
 * This is used in the two-stream upward/downward sweep formulas.
 *
 * @param[in] tau_lw Optical depth of the layer [unitless]. Must be ≥ 0.
 *
 * @return Transmission fraction [0, 1]. Satisfies transmit + absorb = 1.
 *
 * @note If tau_lw = 0, returns 1 (fully transparent).
 * @note If tau_lw >> 1, returns ≈ 0 (fully opaque).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_lw_transmit(amrex::Real tau_lw)
{
  if (tau_lw < 0.0) {
    return 1.0;  // Unphysical, but defensive
  }
  return std::exp(-tau_lw);
}

/**
 * @brief Compute upwelling LW flux in one layer of a two-stream sweep.
 *
 * Implements the upward flux formula in the gray-gas two-stream model:
 *
 *   F_up_current = F_up_below * exp(-tau_lw) + sigma * T^4 * (1 - exp(-tau_lw))
 *
 * This is called iteratively from the surface (where
 * F_up_surface = emissivity * sigma * T_surface^4 + (1 - emissivity) * F_down_surface)
 * upward to the TOA.
 *
 * @param[in] F_up_below Upwelling flux from the layer below [W/m^2].
 * @param[in] T_layer Temperature of the current layer [K].
 * @param[in] sigma Stefan-Boltzmann constant [W/(m^2·K^4)].
 * @param[in] tau_lw Optical depth of the current layer [unitless].
 *
 * @return Upwelling flux at the top of the current layer [W/m^2].
 *
 * @note This function assumes each layer acts as a semi-infinite slab
 *       with uniform temperature and optical depth.
 * @note Called within a device-side kernel; must be GPU-safe.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_lw_flux_up(amrex::Real F_up_below, amrex::Real T_layer,
                   amrex::Real sigma, amrex::Real tau_lw)
{
  amrex::Real transmit = compute_lw_transmit(tau_lw);
  amrex::Real absorb = 1.0 - transmit;
  amrex::Real I_layer = compute_thermal_intensity(T_layer, sigma);
  return F_up_below * transmit + I_layer * absorb;
}

/**
 * @brief Compute downwelling LW flux in one layer of a two-stream sweep.
 *
 * Implements the downward flux formula in the gray-gas two-stream model:
 *
 *   F_down_current = F_down_above * exp(-tau_lw) + sigma * T^4 * (1 - exp(-tau_lw))
 *
 * This is called iteratively from the TOA (where F_down_toa = 0)
 * downward to the surface.
 *
 * @param[in] F_down_above Downwelling flux from the layer above [W/m^2].
 * @param[in] T_layer Temperature of the current layer [K].
 * @param[in] sigma Stefan-Boltzmann constant [W/(m^2·K^4)].
 * @param[in] tau_lw Optical depth of the current layer [unitless].
 *
 * @return Downwelling flux at the bottom of the current layer [W/m^2].
 *
 * @note This function assumes each layer acts as a semi-infinite slab
 *       with uniform temperature and optical depth.
 * @note Called within a device-side kernel; must be GPU-safe.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_lw_flux_down(amrex::Real F_down_above, amrex::Real T_layer,
                     amrex::Real sigma, amrex::Real tau_lw)
{
  amrex::Real transmit = compute_lw_transmit(tau_lw);
  amrex::Real absorb = 1.0 - transmit;
  amrex::Real I_layer = compute_thermal_intensity(T_layer, sigma);
  return F_down_above * transmit + I_layer * absorb;
}

/**
 * @brief Compute LW heating rate from net flux divergence.
 *
 * Given the net LW flux (F_up - F_down, positive upward) at the top and
 * bottom interfaces of a layer, compute the radiative heating rate:
 *
 *   dF_net/dz = [F_net_top - F_net_bot] / dz
 *   Q_lw = -dF_net/dz / (rho * cp)  [K/s]
 *
 * Net upward flux increasing with height means the layer loses more
 * energy through its top than it gains through its bottom, i.e. the layer
 * cools (Q_lw < 0). This is the same sign convention as the SW heating
 * rate: energy converging into the layer warms it.
 *
 * **Nonuniform dz Support & Sanity Checks**
 *
 * This function accepts per-level dz values to support terrain-aware or
 * nonuniform vertical grids. When dz <= 0, or rho <= 0, or cp <= 0,
 * returns 0 (heating rate is undefined or unphysical).
 *
 * The computed heating rate is guarded against NaN/Inf via the same
 * parameter validation.
 *
 * @param[in] F_net_top Net flux at the top of the layer [W/m^2].
 * @param[in] F_net_bot Net flux at the bottom of the layer [W/m^2].
 * @param[in] dz Vertical thickness of the layer [m]. Must be positive.
 *              For nonuniform grids, use per-level spacing.
 * @param[in] rho Density [kg/m^3]. Must be positive.
 * @param[in] cp Specific heat at constant pressure [J/(kg·K)]. Must be positive.
 *
 * @return Heating rate [K/s]. Positive values indicate warming.
 *         Returns 0 if inputs are unphysical (dz, rho, cp <= 0).
 *
 * @note The sign convention is: positive dF_net/dz (net upward flux
 *       increasing with height) corresponds to cooling (dT/dt < 0).
 * @note For terrain-aware grids, pass dz = z_cc(k+1) - z_cc(k) for each level.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_lw_heating_rate(amrex::Real F_net_top, amrex::Real F_net_bot,
                        amrex::Real dz, amrex::Real rho, amrex::Real cp)
{
  if (dz <= 0.0 || rho <= 0.0 || cp <= 0.0) {
    return 0.0;
  }
  // Net flux divergence: (W/m^2) / (m) = W/m^3
  amrex::Real flux_divergence = (F_net_top - F_net_bot) / dz;
  // Heating rate: -(W/m^3) / (kg/m^3 * J/(kg*K)) = K/s
  // Check for NaN/Inf as defensive sanity check
  amrex::Real heating = -flux_divergence / (rho * cp);
  if (!std::isfinite(heating)) {
    return 0.0;
  }
  return heating;
}

#endif // ERF_TWO_STREAM_LW_H_
