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
 * Later phases will extend with scattering and RRTMGP support.
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
