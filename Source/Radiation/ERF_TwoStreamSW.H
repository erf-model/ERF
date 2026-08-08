#ifndef ERF_TWO_STREAM_SW_H_
#define ERF_TWO_STREAM_SW_H_

#include <AMReX_GpuControl.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

/**
 * @file ERF_TwoStreamSW.H
 * @brief Shortwave (solar) radiation using Beer-Lambert direct-beam model
 * plus (Phase 4) a two-stream diffuse (scattering) approximation.
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
 * Phase 4 adds a diffuse SW flux contribution using the Meador-Weaver
 * (1980) two-stream scattering approximation, parameterized by the
 * single-scattering albedo (omega) and asymmetry factor (g) of each
 * layer. When omega == 0 (the Phase 1-3 default), the diffuse
 * contribution below is identically zero and the total SW flux reduces
 * EXACTLY to the pre-Phase-4 direct-beam-only result.
 *
 * References:
 * -----------
 * - Beer, A., 1852: Bestimmung der Absorption des rothen Lichts in farbigen
 *   Flüssigkeiten. Annalen der Physik und Chemie, 86(78), 78-88.
 *
 * - Bird, R. E., R. L. Hulstrom, and J. A. Lewis, 1984: Terrestrial solar
 *   spectral irradiance: Comparison of proposed standard spectra. J. Climate
 *   Appl. Meteor., 23, 539-545.
 *
 * - Meador, W. E., and W. R. Weaver, 1980: Two-stream approximations to
 *   radiative transfer in planetary atmospheres: A unified description
 *   of existing methods and a new improvement. J. Atmos. Sci., 37,
 *   630-643.
 *
 * - Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam, 1989:
 *   Rapid calculation of radiative heating rates and photodissociation
 *   rates in inhomogeneous multiple scattering atmospheres. J. Geophys.
 *   Res., 94, 16287-16301.
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
 * **Phase 9: Nonuniform dz Support & Sanity Checks**
 *
 * This function accepts per-level dz values to support terrain-aware or
 * nonuniform vertical grids. When dz <= 0, or rho <= 0, or cp <= 0,
 * returns 0 (heating rate is undefined or unphysical).
 *
 * The computed heating rate is guarded against NaN/Inf via the same
 * parameter validation.
 *
 * @param[in] flux_top Downwelling flux at the top of the layer [W/m^2].
 * @param[in] flux_bot Downwelling flux at the bottom of the layer [W/m^2].
 * @param[in] dz Vertical thickness of the layer [m]. Must be positive.
 *              For nonuniform grids (Phase 9+), use per-level spacing.
 * @param[in] rho Density [kg/m^3]. Must be positive.
 * @param[in] cp Specific heat at constant pressure [J/(kg·K)]. Must be positive.
 *
 * @return Heating rate [K/s]. Positive values indicate warming; negative indicate cooling.
 *         Returns 0 if inputs are unphysical (dz, rho, cp <= 0).
 *
 * @note This function assumes the layer is the control volume with uniform
 *       properties. No sub-layer variation is computed.
 * @note The heating rate formula follows from energy conservation:
 *       dT/dt = -(1/rho/cp) * (∂F_rad/∂z)
 *            = -(1/rho/cp) * (F_out - F_in) / dz
 * @note For terrain-aware grids, pass dz = z_cc(k) - z_cc(k+1) for each level.
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
  // Check for NaN/Inf as defensive sanity check
  amrex::Real heating = flux_divergence / (rho * cp);
  if (!std::isfinite(heating)) {
    return 0.0;
  }
  return heating;
}

/**
 * @brief (Phase 4) Compute the diffuse (scattered) downwelling SW flux at
 * the base of a single homogeneous layer using the Meador-Weaver (1980)
 * two-stream scattering approximation, driven by the direct-beam flux
 * incident at the top of that layer.
 *
 * This models the single-scattering source of diffuse radiation
 * generated within one layer by the attenuating direct beam, plus
 * multiple-scattering amplification via the layer's diffuse
 * reflectance/transmittance. It is intended to be accumulated layer by
 * layer in a vertical sweep (added to the transmitted diffuse flux from
 * the layer above) to build up the total diffuse downwelling flux
 * profile, analogous to how tau_cumulative is accumulated for the
 * direct beam in compute_sw_direct_flux().
 *
 * Meador-Weaver two-stream coefficients (delta-scaled quadrature closure,
 * gamma set for the "quadrature" two-stream method):
 *   gamma1 = (7 - omega*(4 + 3*g)) / 4
 *   gamma2 = -(1 - omega*(4 - 3*g)) / 4
 *   gamma3 = (2 - 3*g*cos_zenith) / 4
 *   gamma4 = 1 - gamma3
 *
 *   lambda = sqrt(gamma1^2 - gamma2^2)
 *   Gamma  = gamma2 / (gamma1 + lambda)      (equivalent to (gamma1-lambda)/gamma2)
 *
 *   Tdir   = exp(-tau / cos_zenith)          (direct-beam transmittance of layer)
 *   Tdiff  = exp(-lambda * tau)              (diffuse-beam transmittance of layer)
 *
 * Layer diffuse reflectance/transmittance (no surface albedo, isolated layer):
 *   R_diff = Gamma * (1 - Tdiff^2) / (1 - Gamma^2 * Tdiff^2)
 *   T_diff = Tdiff * (1 - Gamma^2) / (1 - Gamma^2 * Tdiff^2)
 *
 * Particular (single-scattering) solution source terms driven by the
 * attenuating direct beam (mu0 = cos_zenith), following Toon et al.
 * (1989) Eq. 27-28 form, simplified for downwelling diffuse flux at the
 * layer base assuming no diffuse flux incident from above this layer
 * (i.e., this function returns ONLY the diffuse flux generated within
 * this layer; the caller is responsible for adding the diffuse flux
 * transmitted through from layers above, scaled by T_diff, if a full
 * multi-layer sweep is desired):
 *
 *   C_down = omega * F_dir_top * cos_zenith / (lambda^2 - 1/cos_zenith^2)
 *            * [ (gamma1 - 1/cos_zenith) * Tdir  ... see Toon et al. 1989 ]
 *
 * For simplicity and robustness (avoiding singularities when
 * lambda ≈ 1/cos_zenith), this implementation uses the standard
 * single-scattering-source approximation:
 *
 *   F_diff_base = omega * (F_dir_top - F_dir_top * Tdir) * (gamma3 + gamma4) / 2
 *                 + F_dir_top * Tdir * R_diff * gamma3
 *
 * which reduces to the physically-required limits:
 *   - omega == 0  => F_diff_base == 0 exactly (pure direct-beam, Phase 1-3)
 *   - omega -> 1, g -> 1 (strong forward scattering) => most of the
 *     scattered energy remains in the forward (downward) direction
 *
 * This is a GPU-safe inline function intended for use within device-side
 * parallel kernels.
 *
 * @param[in] tau Optical depth of this single layer (NOT cumulative) [unitless].
 * @param[in] F_dir_top Direct-beam downwelling flux incident at the TOP
 *   of this layer [W/m^2]. Typically the output of
 *   compute_sw_direct_flux() evaluated at the cumulative optical depth
 *   at the top of this layer.
 * @param[in] cos_zenith Cosine of solar zenith angle [unitless, 0 < cos_zenith <= 1].
 * @param[in] omega Single-scattering albedo of this layer, in [0, 1].
 *   omega == 0.0 means a purely absorbing layer; the function returns
 *   exactly 0.0 in that case, with NO further calculation performed.
 * @param[in] g Asymmetry factor of this layer, in [-1, 1]. 0 = isotropic
 *   scattering; close to 1 = strongly forward-scattering (typical for
 *   cloud droplets, g ~ 0.85).
 *
 * @return Diffuse downwelling SW flux generated within this layer,
 *   evaluated at the layer base [W/m^2]. Always >= 0.
 *
 * @note If cos_zenith <= 0 (night) or omega <= 0 (purely absorbing),
 *   returns exactly 0.0 — this is the Phase 1-3 direct-beam-only limit.
 * @note This single-layer formulation does not itself accumulate
 *   diffuse transmission from layers above; the calling vertical sweep
 *   is expected to add contributions layer-by-layer (see
 *   ERF_AdvanceTwoStreamRadiation.cpp Phase 4 section).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_sw_diffuse_flux(amrex::Real tau, amrex::Real F_dir_top,
                        amrex::Real cos_zenith, amrex::Real omega,
                        amrex::Real g)
{
  // Exact zero-scattering limit: byte-identical to Phase 1-3 behavior.
  if (omega <= 0.0 || cos_zenith <= 0.0 || F_dir_top <= 0.0 || tau <= 0.0) {
    return 0.0;
  }

  // Clip inputs defensively to valid physical ranges.
  amrex::Real w0 = omega;
  if (w0 > 1.0) w0 = 1.0;
  amrex::Real asym = g;
  if (asym > 1.0) asym = 1.0;
  if (asym < -1.0) asym = -1.0;

  // Meador-Weaver (1980) quadrature two-stream gamma coefficients.
  amrex::Real gamma1 = (7.0 - w0 * (4.0 + 3.0 * asym)) / 4.0;
  amrex::Real gamma2 = -(1.0 - w0 * (4.0 - 3.0 * asym)) / 4.0;
  amrex::Real gamma3 = (2.0 - 3.0 * asym * cos_zenith) / 4.0;
  amrex::Real gamma4 = 1.0 - gamma3;

  amrex::Real lambda_sq = gamma1 * gamma1 - gamma2 * gamma2;
  if (lambda_sq < 1.0e-12) lambda_sq = 1.0e-12;  // Defensive: avoid sqrt of ~0/negative
  amrex::Real lambda = std::sqrt(lambda_sq);

  // Gamma (reflection amplitude ratio); guard against gamma1+lambda == 0.
  amrex::Real denom_Gamma = gamma1 + lambda;
  amrex::Real Gamma = (std::abs(denom_Gamma) > 1.0e-12)
                         ? (gamma2 / denom_Gamma)
                         : 0.0;

  amrex::Real Tdir  = std::exp(-tau / cos_zenith);
  amrex::Real Tdiff = std::exp(-lambda * tau);

  amrex::Real Gamma_sq = Gamma * Gamma;
  amrex::Real Tdiff_sq = Tdiff * Tdiff;
  amrex::Real denom = 1.0 - Gamma_sq * Tdiff_sq;
  if (std::abs(denom) < 1.0e-12) denom = 1.0e-12;  // Defensive

  amrex::Real R_diff = Gamma * (1.0 - Tdiff_sq) / denom;

  // Direct-beam flux absorbed/scattered within this layer.
  amrex::Real F_dir_bot = F_dir_top * Tdir;
  amrex::Real F_dir_attenuated = F_dir_top - F_dir_bot;

  // Single-scattering source of diffuse flux at the layer base: energy
  // scattered from the direct beam as it passes through the layer,
  // partitioned into the downward hemisphere via (gamma3+gamma4)/2 = 0.5
  // (by construction gamma3+gamma4=1, so this reduces to a factor of 0.5,
  // i.e. symmetric single-scattering source), plus additional downward
  // diffuse flux from the direct beam that reached the layer base being
  // partially backscattered by the layer's own diffuse reflectance.
  amrex::Real F_diff_base = w0 * F_dir_attenuated * 0.5
                             + F_dir_bot * R_diff * gamma3;

  if (F_diff_base < 0.0) F_diff_base = 0.0;  // Defensive: non-negative flux

  return F_diff_base;
}

#endif // ERF_TWO_STREAM_SW_H_
