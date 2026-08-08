#ifndef ERF_RAD_STRUCT_H_
#define ERF_RAD_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <ERF_TurbStruct.H>

/**
 * @file ERF_RadStruct.H
 * @brief Radiation model type and control parameters.
 */

/**
 * @brief Radiation model type.
 *
 * Selects which radiative transfer solver to use, if any.
 */
AMREX_ENUM(RadType, None, TwoStream, RRTMGP);

/**
 * @brief Optical depth profile type for the two-stream solver.
 *
 * - Constant: tau_per_layer / tau_lw_per_layer applied uniformly at every
 *   level (Phase 1/2 behavior, byte-identical).
 * - CloudLayer: an additional cloud optical depth is added within
 *   [cloud_base_height_m, cloud_top_height_m], on top of the clear-sky
 *   background tau_per_layer / tau_lw_per_layer.
 */
AMREX_ENUM(TauProfileType, Constant, CloudLayer);

/**
 * @brief Container holding radiation-related choices and parameters.
 *
 * Manages the runtime selection and configuration of radiation solvers.
 * Phase 1/2: clear-sky, no-scattering physics, uniform optical depth.
 * Phase 3: adds height-varying optical depth (cloud layer) and cloud
 * fraction masking (blend of clear-sky and cloudy-column fluxes).
 * Phase 4: adds a diffuse (scattering) SW component via a two-stream
 * scattering approximation, parameterized by single-scattering albedo
 * and asymmetry factor, with independent clear-sky and cloud values.
 * Later phases will extend with RRTMGP support.
 */
struct RadChoice
{
public:
  /**
   * @brief Radiation model type (None, TwoStream, RRTMGP).
   */
  RadType rad_type = RadType::None;

  /**
   * @brief Enable shortwave (solar) radiation computation.
   * Only used when rad_type == RadType::TwoStream.
   */
  bool sw_enabled = true;

  /**
   * @brief Enable longwave (thermal) radiation computation.
   * Only used when rad_type == RadType::TwoStream.
   */
  bool lw_enabled = true;

  /**
   * @brief Shortwave optical depth per layer (constant for all layers).
   * Used in Beer-Lambert direct-beam formula.
   * Only used when rad_type == RadType::TwoStream.
   */
  amrex::Real tau_per_layer = 0.05;

  /**
   * @brief Longwave optical depth per layer (constant for all layers).
   * Used in gray-gas LW two-stream solver.
   * Only used when rad_type == RadType::TwoStream.
   */
  amrex::Real tau_lw_per_layer = 1.0;

  /**
   * @brief Solar zenith angle in degrees.
   * Only used when rad_type == RadType::TwoStream and sw_enabled == true.
   */
  amrex::Real solar_zenith_deg = 45.0;

  /**
   * @brief Solar constant (top-of-atmosphere irradiance) in W/m^2.
   * Only used when rad_type == RadType::TwoStream and sw_enabled == true.
   */
  amrex::Real S0 = 1361.0;

  /**
   * @brief Enable isothermal LW test mode.
   * Forces all temperatures to T_iso_K for analytical verification.
   * Only used when rad_type == RadType::TwoStream and lw_enabled == true.
   */
  bool isothermal_test = false;

  /**
   * @brief Isothermal temperature for LW test (Kelvin).
   * Only used when rad_type == RadType::TwoStream and isothermal_test == true.
   */
  amrex::Real T_iso_K = 288.15;

  /**
   * @brief Debug verbosity level (0=off, 1+=debug output).
   * Gated debug prints only emit when this is >= 1.
   */
  int verbosity = 0;

  /**
   * @brief Output file path for radiation diagnostics CSV.
   * Paths are relative to the simulation working directory.
   */
  std::string diag_file = "radiation_diag.dat";

  /**
   * @brief (Phase 7) Master switch for radiation diagnostics emission.
   * When false, no diagnostics output (stdout/tagged/regtest/csv).
   * Default true preserves Phase 6 behavior.
   */
  bool diag_enable = true;

  /**
   * @brief (Phase 7) Enable human-readable stdout diagnostics block.
   * Controls tagged [RAD][...] debug line emission (when verbosity >= 1).
   * Default true preserves Phase 6 behavior.
   */
  bool diag_stdout_enable = true;

  /**
   * @brief (Phase 7) Enable tagged [RAD][...] debug lines.
   * Controls emission of bracketed-tag debug output.
   * Default true preserves Phase 6 behavior.
   */
  bool diag_tagged_enable = true;

  /**
   * @brief (Phase 7) Enable RADIATION_DIAG: line emission.
   * Controls regtest-parsing diagnostics line.
   * Default true preserves Phase 6 behavior.
   */
  bool diag_regtest_line_enable = true;

  /**
   * @brief (Phase 7) Enable CSV file append behavior.
   * When false, no file writes; other enabled streams may still print.
   * Default true preserves Phase 6 behavior.
   */
  bool diag_csv_enable = true;

  /**
   * @brief (Phase 7) Call-site filtering mode: "both" (pre+post rows),
   * "pre_only", or "post_only". Default "both" preserves Phase 6 behavior.
   */
  std::string diag_callsite_mode = "both";

  /**
   * @brief (Phase 7) Tolerance for time equality in duplicate guard [s].
   * Times within this tolerance are considered equal.
   * Default 1e-12 preserves Phase 6 behavior.
   */
  amrex::Real diag_dedup_tol = 1.0e-12;

  /**
   * @brief (Phase 3) Optical depth profile type: "constant" (default,
   * byte-identical to Phase 1/2) or "cloud_layer" (adds height-varying
   * cloud optical depth between cloud_base_height_m and cloud_top_height_m).
   */
  TauProfileType tau_profile_type = TauProfileType::Constant;

  /**
   * @brief (Phase 3) Cloud layer base height [m]. Only used when
   * tau_profile_type == CloudLayer.
   */
  amrex::Real cloud_base_height_m = 500.0;

  /**
   * @brief (Phase 3) Cloud layer top height [m]. Only used when
   * tau_profile_type == CloudLayer. Must be >= cloud_base_height_m.
   */
  amrex::Real cloud_top_height_m = 1000.0;

  /**
   * @brief (Phase 3) Additional optical depth per layer contributed by the
   * cloud, added on top of the clear-sky tau_per_layer/tau_lw_per_layer
   * background within [cloud_base_height_m, cloud_top_height_m]. Only used
   * when tau_profile_type == CloudLayer.
   */
  amrex::Real cloud_tau_per_layer = 0.5;

  /**
   * @brief (Phase 3) Cloud fraction in [0, 1] used to blend clear-sky and
   * cloudy-column fluxes: F = (1 - cloud_fraction) * F_clear +
   * cloud_fraction * F_cloudy. Default 0.0 (fully clear sky), which
   * reduces exactly to Phase 2d behavior regardless of tau_profile_type.
   */
  amrex::Real cloud_fraction = 0.0;

  /**
   * @brief (Phase 4) Clear-sky (background gas) single-scattering albedo
   * in [0, 1] for the SW diffuse (scattering) two-stream calculation.
   * Default 0.0 means the clear-sky column is purely absorbing, and the
   * diffuse SW contribution from the clear-sky column is exactly zero —
   * reducing EXACTLY to the Phase 1-3 direct-beam-only result.
   * Only used when rad_type == RadType::TwoStream and sw_enabled == true.
   */
  amrex::Real single_scattering_albedo = 0.0;

  /**
   * @brief (Phase 4) Clear-sky (background gas) scattering asymmetry
   * factor in [-1, 1] (0 = isotropic scattering, close to 1 = strongly
   * forward-scattering). Only meaningful when single_scattering_albedo
   * > 0. Default 0.0.
   */
  amrex::Real asymmetry_factor = 0.0;

  /**
   * @brief (Phase 4) Cloud single-scattering albedo in [0, 1], used
   * instead of single_scattering_albedo for levels where the Phase 3
   * cloud-layer optical depth enhancement applies (i.e., within
   * [cloud_base_height_m, cloud_top_height_m] on the "cloudy" column
   * evaluation). Default 0.0 preserves byte-identical Phase 3 behavior
   * for SW_Cloud_Layer unless explicitly configured otherwise. Typical
   * liquid water cloud values are close to 1 (e.g., 0.9999).
   * Only used when tau_profile_type == CloudLayer.
   */
  amrex::Real cloud_single_scattering_albedo = 0.0;

  /**
   * @brief (Phase 4) Cloud scattering asymmetry factor in [-1, 1], used
   * instead of asymmetry_factor for levels where the Phase 3 cloud-layer
   * optical depth enhancement applies. Default 0.0. Typical liquid water
   * cloud droplet values are around 0.85 (strongly forward-scattering).
   * Only meaningful when cloud_single_scattering_albedo > 0.
   */
  amrex::Real cloud_asymmetry_factor = 0.0;

  /**
   * @brief (Phase 11) Shortwave surface albedo fallback [0, 1] for per-column
   * heterogeneous-surface-property resolution. Used when per-column LSM/radiation
   * interface albedo field is unavailable or invalid. Default 0.3 represents
   * typical land/water reflectivity. Only used when hetero surface fields not
   * present or disabled.
   */
  amrex::Real surface_albedo_sw = 0.3;

  /**
   * @brief (Phase 11) Longwave surface emissivity fallback [0, 1] for per-column
   * heterogeneous-surface-property resolution. Used when per-column LSM/radiation
   * interface emissivity field is unavailable or invalid. Default 0.99 represents
   * typical land/water/ice emissivity (nearly blackbody). Only used when hetero
   * surface fields not present or disabled.
   */
  amrex::Real surface_emissivity_lw = 0.99;

  /**
   * @brief (Phase 11) Surface temperature fallback [K] for per-column heterogeneous
   * -surface-property resolution. Used as LW boundary condition when per-column
   * LSM/radiation interface surface-temperature field is unavailable or invalid.
   * Default 300.0 K (typical surface temp). Only used when hetero t_sfc field
   * not present or disabled.
   */
  amrex::Real surface_temp_k = 300.0;

  /**
   * @brief (Phase 12) Enable dynamic shortwave optical depth diagnosis from moisture/clouds.
   * When true, SW optical depth is computed per-level from qv and qc.
   * When false (default), uses static tau_per_layer (Phase 11 behavior).
   * Only used when rad_type == RadType::TwoStream and sw_enabled == true.
   */
  bool tau_sw_dynamic_enable = false;

  /**
   * @brief (Phase 12) Enable dynamic longwave optical depth diagnosis from moisture/clouds.
   * When true, LW optical depth is computed per-level from qv and qc.
   * When false (default), uses static tau_lw_per_layer (Phase 11 behavior).
   * Only used when rad_type == RadType::TwoStream and lw_enabled == true.
   */
  bool tau_lw_dynamic_enable = false;

  /**
   * @brief (Phase 12) Shortwave optical depth coefficient for water vapor [dimensionless].
   * Dynamic SW tau = tau_sw_coeff_qv * qv + tau_sw_coeff_qc * qc + tau_per_layer.
   * Default 0.0 means tau is purely from static tau_per_layer (no dynamic component).
   * Only used when tau_sw_dynamic_enable == true.
   */
  amrex::Real tau_sw_coeff_qv = 0.0;

  /**
   * @brief (Phase 12) Shortwave optical depth coefficient for cloud liquid water [dimensionless].
   * Default 0.0. When both coefficients are 0, dynamic path still invoked but yields
   * tau_per_layer unchanged, preserving Phase 11 behavior.
   * Only used when tau_sw_dynamic_enable == true.
   */
  amrex::Real tau_sw_coeff_qc = 0.0;

  /**
   * @brief (Phase 12) Longwave optical depth coefficient for water vapor [dimensionless].
   * Dynamic LW tau = tau_lw_coeff_qv * qv + tau_lw_coeff_qc * qc + tau_lw_per_layer.
   * Default 0.0 means tau is purely from static tau_lw_per_layer (no dynamic component).
   * Only used when tau_lw_dynamic_enable == true.
   */
  amrex::Real tau_lw_coeff_qv = 0.0;

  /**
   * @brief (Phase 12) Longwave optical depth coefficient for cloud liquid water [dimensionless].
   * Default 0.0. When both coefficients are 0, dynamic path still invoked but yields
   * tau_lw_per_layer unchanged, preserving Phase 11 behavior.
   * Only used when tau_lw_dynamic_enable == true.
   */
  amrex::Real tau_lw_coeff_qc = 0.0;

  /**
   * @brief (Phase 14) Enable prognostic cloud fraction diagnosis from RH/qc.
   * When true, cloud fraction is computed per-level from relative humidity and cloud water.
   * When false (default), uses static cloud_fraction (Phase 13 behavior).
   * Only used when rad_type == RadType::TwoStream.
   */
  bool cloud_fraction_prog_enable = false;

  /**
   * @brief (Phase 14) Minimum relative humidity threshold for cloud fraction diagnosis [0, 1].
   * Cloud fraction ramps from 0 at rh < cloud_fraction_rh_min to 1 at rh >= cloud_fraction_rh_max.
   * Default 0.0 means cloud fraction depends primarily on qc content.
   * Only used when cloud_fraction_prog_enable == true.
   */
  amrex::Real cloud_fraction_rh_min = 0.0;

  /**
   * @brief (Phase 14) Maximum relative humidity threshold for cloud fraction diagnosis [0, 1].
   * Must be >= cloud_fraction_rh_min. Cloud fraction saturates to 1 at this RH.
   * Default 1.0 means cloud fraction increases linearly from rh_min to 1.0 RH.
   * Only used when cloud_fraction_prog_enable == true.
   */
  amrex::Real cloud_fraction_rh_max = 1.0;

  /**
   * @brief (Phase 14) Scaling coefficient for cloud liquid water contribution to cloud fraction [dimensionless].
   * Cloud fraction has additional contribution from qc: cf += cloud_fraction_qc_scale * qc.
   * Default 1.0e-3 converts qc [kg/kg] to dimensionless cloud fraction contribution.
   * Only used when cloud_fraction_prog_enable == true.
   */
  amrex::Real cloud_fraction_qc_scale = 1.0e-3;

  /**
   * @brief (Phase 14) Enable temporal smoothing of prognostic cloud fraction.
   * When true, diagnosed cf(k) is smoothed with exponential moving average (EMA).
   * When false (default), diagnosed cf(k) is applied directly without smoothing.
   * Only used when cloud_fraction_prog_enable == true.
   */
  bool cloud_fraction_smooth_enable = false;

  /**
   * @brief (Phase 14) Temporal smoothing parameter for cloud fraction [0, 1].
   * Controls EMA blending: cf_smooth = alpha * cf_new + (1 - alpha) * cf_old.
   * alpha = 0.0 means no smoothing (cf_old retained); alpha = 1.0 means instant update.
   * Default 0.0 preserves feature-off behavior when smoothing enabled but alpha=0.
   * Only used when cloud_fraction_smooth_enable == true.
   */
  amrex::Real cloud_fraction_smooth_alpha = 0.0;

  /**
   * @brief Read radiation options from the input parameter database.
   *
   * @param[in] lev AMR level index.
   * @param[in] max_level Maximum AMR level configured for the run.
   * @param[in] pp_prefix ParmParse prefix for the ERF input namespace (e.g., "erf").
   *
   * Queries:
   * - erf.radiation_type (enum, required if TwoStream is desired)
   * - erf.radiation.sw_enabled (bool, optional, default true)
   * - erf.radiation.lw_enabled (bool, optional, default true)
   * - erf.radiation.tau_per_layer (real, optional)
   * - erf.radiation.tau_lw_per_layer (real, optional)
   * - erf.radiation.solar_zenith (real, optional)
   * - erf.radiation.S0 (real, optional)
   * - erf.radiation.isothermal_test (bool, optional)
   * - erf.radiation.T_iso_K (real, optional)
   * - erf.radiation.v (int, optional, debug verbosity)
   * - erf.radiation.diag_file (string, optional)
   * - erf.radiation.tau_profile_type (string "constant"|"cloud_layer", optional, Phase 3)
   * - erf.radiation.cloud_base_height_m (real, optional, Phase 3)
   * - erf.radiation.cloud_top_height_m (real, optional, Phase 3)
   * - erf.radiation.cloud_tau_per_layer (real, optional, Phase 3)
   * - erf.radiation.cloud_fraction (real, optional, Phase 3)
   * - erf.radiation.single_scattering_albedo (real, optional, Phase 4)
   * - erf.radiation.asymmetry_factor (real, optional, Phase 4)
   * - erf.radiation.cloud_single_scattering_albedo (real, optional, Phase 4)
   * - erf.radiation.cloud_asymmetry_factor (real, optional, Phase 4)
   * - erf.radiation.diag_enable (bool, optional, Phase 7, default true)
   * - erf.radiation.diag_stdout_enable (bool, optional, Phase 7, default true)
   * - erf.radiation.diag_tagged_enable (bool, optional, Phase 7, default true)
   * - erf.radiation.diag_regtest_line_enable (bool, optional, Phase 7, default true)
   * - erf.radiation.diag_csv_enable (bool, optional, Phase 7, default true)
   * - erf.radiation.diag_callsite_mode (string, optional, Phase 7, default "both")
   * - erf.radiation.diag_dedup_tol (real, optional, Phase 7, default 1e-12)
   * - erf.radiation.tau_sw_dynamic_enable (bool, optional, Phase 12, default false)
   * - erf.radiation.tau_lw_dynamic_enable (bool, optional, Phase 12, default false)
   * - erf.radiation.tau_sw_coeff_qv (real, optional, Phase 12, default 0.0)
   * - erf.radiation.tau_sw_coeff_qc (real, optional, Phase 12, default 0.0)
   * - erf.radiation.tau_lw_coeff_qv (real, optional, Phase 12, default 0.0)
   * - erf.radiation.tau_lw_coeff_qc (real, optional, Phase 12, default 0.0)
   * - erf.radiation.cloud_fraction_prog_enable (bool, optional, Phase 14, default false)
   * - erf.radiation.cloud_fraction_rh_min (real, optional, Phase 14, default 0.0)
   * - erf.radiation.cloud_fraction_rh_max (real, optional, Phase 14, default 1.0)
   * - erf.radiation.cloud_fraction_qc_scale (real, optional, Phase 14, default 1.0e-3)
   * - erf.radiation.cloud_fraction_smooth_enable (bool, optional, Phase 14, default false)
   * - erf.radiation.cloud_fraction_smooth_alpha (real, optional, Phase 14, default 0.0)
   *
   * Sub-option queries are only performed if rad_type == RadType::TwoStream.
   */
  void init_params(int lev, int max_level, std::string pp_prefix)
  {
    amrex::ParmParse pp(pp_prefix);

    // Query radiation type, using helper to support per-level specification
    // Reuse the same helper used in TurbStruct for enum handling
    query_one_or_per_level_enum_case_insensitive(
      pp, "radiation_type", rad_type, lev, max_level);

    // Only query sub-options if TwoStream radiation is selected
    if (rad_type == RadType::TwoStream) {
      pp.query("radiation.sw_enabled", sw_enabled);
      pp.query("radiation.lw_enabled", lw_enabled);
      pp.query("radiation.tau_per_layer", tau_per_layer);
      pp.query("radiation.tau_lw_per_layer", tau_lw_per_layer);
      pp.query("radiation.solar_zenith", solar_zenith_deg);
      pp.query("radiation.S0", S0);
      pp.query("radiation.isothermal_test", isothermal_test);
      pp.query("radiation.T_iso_K", T_iso_K);
      pp.query("radiation.v", verbosity);
      pp.query("radiation.diag_file", diag_file);

      // Phase 3: cloud optical properties (all optional; defaults preserve
      // exact Phase 2d behavior)
      std::string tau_profile_type_str = "constant";
      pp.query("radiation.tau_profile_type", tau_profile_type_str);
      if (tau_profile_type_str == "cloud_layer") {
        tau_profile_type = TauProfileType::CloudLayer;
      } else {
        tau_profile_type = TauProfileType::Constant;
      }
      pp.query("radiation.cloud_base_height_m", cloud_base_height_m);
      pp.query("radiation.cloud_top_height_m", cloud_top_height_m);
      pp.query("radiation.cloud_tau_per_layer", cloud_tau_per_layer);
      pp.query("radiation.cloud_fraction", cloud_fraction);

      // Parameter validation (R5/R5'): clip cloud_fraction to [0, 1]
      if (cloud_fraction < 0.0) cloud_fraction = 0.0;
      if (cloud_fraction > 1.0) cloud_fraction = 1.0;
      if (cloud_tau_per_layer < 0.0) cloud_tau_per_layer = 0.0;
      if (cloud_top_height_m < cloud_base_height_m) {
        cloud_top_height_m = cloud_base_height_m;
      }

      // Phase 4: scattering parameters (all optional; defaults preserve
      // exact Phase 3 direct-beam-only behavior)
      pp.query("radiation.single_scattering_albedo", single_scattering_albedo);
      pp.query("radiation.asymmetry_factor", asymmetry_factor);
      pp.query("radiation.cloud_single_scattering_albedo", cloud_single_scattering_albedo);
      pp.query("radiation.cloud_asymmetry_factor", cloud_asymmetry_factor);

      // Parameter validation (R5/R5'): clip albedos to [0, 1] and
      // asymmetry factors to [-1, 1]
      if (single_scattering_albedo < 0.0) single_scattering_albedo = 0.0;
      if (single_scattering_albedo > 1.0) single_scattering_albedo = 1.0;
      if (asymmetry_factor < -1.0) asymmetry_factor = -1.0;
      if (asymmetry_factor > 1.0) asymmetry_factor = 1.0;
      if (cloud_single_scattering_albedo < 0.0) cloud_single_scattering_albedo = 0.0;
      if (cloud_single_scattering_albedo > 1.0) cloud_single_scattering_albedo = 1.0;
      if (cloud_asymmetry_factor < -1.0) cloud_asymmetry_factor = -1.0;
      if (cloud_asymmetry_factor > 1.0) cloud_asymmetry_factor = 1.0;

      // Phase 7: runtime diagnostics controls (all optional; defaults preserve
      // exact Phase 6 behavior)
      pp.query("radiation.diag_enable", diag_enable);
      pp.query("radiation.diag_stdout_enable", diag_stdout_enable);
      pp.query("radiation.diag_tagged_enable", diag_tagged_enable);
      pp.query("radiation.diag_regtest_line_enable", diag_regtest_line_enable);
      pp.query("radiation.diag_csv_enable", diag_csv_enable);
      pp.query("radiation.diag_callsite_mode", diag_callsite_mode);
      pp.query("radiation.diag_dedup_tol", diag_dedup_tol);

      // Parameter validation (Phase 7): validate callsite mode
      if (diag_callsite_mode != "both" &&
          diag_callsite_mode != "pre_only" &&
          diag_callsite_mode != "post_only") {
        diag_callsite_mode = "both";  // Default on invalid input
      }
      // Ensure dedup tolerance is non-negative
      if (diag_dedup_tol < 0.0) diag_dedup_tol = 0.0;

      // Phase 11: surface heterogeneity fallback parameters (all optional;
      // defaults preserve Phase 1-10 scalar behavior; only used when hetero
      // LSM/radiation interface fields are absent or invalid)
      pp.query("radiation.surface_albedo_sw", surface_albedo_sw);
      pp.query("radiation.surface_emissivity_lw", surface_emissivity_lw);
      pp.query("radiation.surface_temp_k", surface_temp_k);

      // Parameter validation (Phase 11): clamp albedo/emissivity to [0, 1],
      // ensure surface temp is positive and finite
      if (surface_albedo_sw < 0.0) surface_albedo_sw = 0.0;
      if (surface_albedo_sw > 1.0) surface_albedo_sw = 1.0;
      if (surface_emissivity_lw < 0.0) surface_emissivity_lw = 0.0;
      if (surface_emissivity_lw > 1.0) surface_emissivity_lw = 1.0;
      if (surface_temp_k <= 0.0) surface_temp_k = 300.0;  // Ensure positive
      if (!std::isfinite(surface_temp_k)) surface_temp_k = 300.0;  // Ensure finite

      // Phase 12: dynamic optical depth diagnosis from moisture/clouds
      // (all optional; defaults disable dynamic tau, preserving Phase 11 behavior)
      pp.query("radiation.tau_sw_dynamic_enable", tau_sw_dynamic_enable);
      pp.query("radiation.tau_lw_dynamic_enable", tau_lw_dynamic_enable);
      pp.query("radiation.tau_sw_coeff_qv", tau_sw_coeff_qv);
      pp.query("radiation.tau_sw_coeff_qc", tau_sw_coeff_qc);
      pp.query("radiation.tau_lw_coeff_qv", tau_lw_coeff_qv);
      pp.query("radiation.tau_lw_coeff_qc", tau_lw_coeff_qc);

      // Parameter validation (Phase 12): clamp dynamic tau coefficients
      // to nonnegative (negative would give unphysical behavior)
      if (tau_sw_coeff_qv < 0.0) tau_sw_coeff_qv = 0.0;
      if (tau_sw_coeff_qc < 0.0) tau_sw_coeff_qc = 0.0;
      if (tau_lw_coeff_qv < 0.0) tau_lw_coeff_qv = 0.0;
      if (tau_lw_coeff_qc < 0.0) tau_lw_coeff_qc = 0.0;

      // Phase 14: prognostic cloud fraction diagnosis from RH/qc
      // (all optional; defaults disable prognostic cf, preserving Phase 13 behavior)
      pp.query("radiation.cloud_fraction_prog_enable", cloud_fraction_prog_enable);
      pp.query("radiation.cloud_fraction_rh_min", cloud_fraction_rh_min);
      pp.query("radiation.cloud_fraction_rh_max", cloud_fraction_rh_max);
      pp.query("radiation.cloud_fraction_qc_scale", cloud_fraction_qc_scale);
      pp.query("radiation.cloud_fraction_smooth_enable", cloud_fraction_smooth_enable);
      pp.query("radiation.cloud_fraction_smooth_alpha", cloud_fraction_smooth_alpha);

      // Parameter validation (Phase 14): clamp cloud fraction parameters
      // to physically reasonable ranges
      if (cloud_fraction_rh_min < 0.0) cloud_fraction_rh_min = 0.0;
      if (cloud_fraction_rh_min > 1.0) cloud_fraction_rh_min = 1.0;
      if (cloud_fraction_rh_max < 0.0) cloud_fraction_rh_max = 0.0;
      if (cloud_fraction_rh_max > 1.0) cloud_fraction_rh_max = 1.0;
      if (cloud_fraction_rh_max < cloud_fraction_rh_min) {
        cloud_fraction_rh_max = cloud_fraction_rh_min;
      }
      if (cloud_fraction_qc_scale < 0.0) cloud_fraction_qc_scale = 0.0;
      if (cloud_fraction_smooth_alpha < 0.0) cloud_fraction_smooth_alpha = 0.0;
      if (cloud_fraction_smooth_alpha > 1.0) cloud_fraction_smooth_alpha = 1.0;
  }
}

  /**
   * @brief Check if two-stream radiation is active.
   *
   * @return true if rad_type == RadType::TwoStream, false otherwise.
   */
  [[nodiscard]] bool uses_two_stream_radiation() const noexcept
  {
    return rad_type == RadType::TwoStream;
  }
};

#endif // ERF_RAD_STRUCT_H_
