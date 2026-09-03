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
 * @brief Aerosol optical depth profile type for the two-stream solver (Phase 15).
 *
 * - Constant: uniform aerosol tau_per_layer applied at every level.
 * - Exponential: aerosol tau decays exponentially with height via
 *   tau_aerosol(k) = tau_surface * exp(-z(k) / scale_height).
 * - Table: aerosol tau prescribed per level via lookup (future extension).
 */
AMREX_ENUM(AerosolProfileType, Constant, Exponential, Table);

/**
 * @brief How the clear-sky optical depth of a layer is obtained.
 *
 * - PerLayer: the fixed tau_per_layer / tau_lw_per_layer for every layer
 *   regardless of its thickness (the original model; the column optical
 *   depth then scales with the number of vertical cells).
 * - Mass: from the layer mass path, tau = rho dz (k_dry + k_vapor qv + ...),
 *   so the column optical depth is set by the atmosphere, not by the grid.
 *   Shortwave uses sw_kabs_dry, sw_kscat_dry (Rayleigh, scattering),
 *   sw_kabs_vapor and sw_kext_cloud with sw_cloud_omega / sw_cloud_g;
 *   longwave uses lw_kabs_dry, lw_kabs_vapor and lw_kabs_cloud.
 */
AMREX_ENUM(TauModel, PerLayer, Mass);

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
   * @brief Scale the solar constant by the Earth-Sun distance factor
   * (d0/d)^2 for day_of_year (Spencer 1971). About +3.4% in early January
   * and -3.3% in early July. Default false keeps S0 fixed.
   */
  bool earth_sun_distance_enable = false;

  /**
   * @brief Optical depth model: "per_layer" (default, fixed tau per layer)
   * or "mass" (from the layer mass path in both bands; see TauModel).
   */
  TauModel tau_model = TauModel::PerLayer;

  /** @brief Mass model, SW: dry-air gray absorption coefficient [m^2/kg]. */
  amrex::Real sw_kabs_dry = 4.0e-6;

  /** @brief Mass model, SW: dry-air (Rayleigh) scattering coefficient [m^2/kg], omega = 1, g = 0. */
  amrex::Real sw_kscat_dry = 3.0e-6;

  /** @brief Mass model, SW: water-vapor gray absorption coefficient [m^2/kg]. */
  amrex::Real sw_kabs_vapor = 4.0e-3;

  /** @brief Mass model, SW: cloud-water extinction coefficient [m^2/kg] (~1.5/r_e with r_e = 10 um). */
  amrex::Real sw_kext_cloud = 150.0;

  /** @brief Mass model, SW: cloud-water single-scattering albedo. */
  amrex::Real sw_cloud_omega = 0.9999;

  /** @brief Mass model, SW: cloud-water asymmetry factor. */
  amrex::Real sw_cloud_g = 0.85;

  /**
   * @brief Gray longwave optical depth from the layer mass path instead of
   * the fixed tau_lw_per_layer (also selected by tau_model = mass):
   *   tau_lw(k) = rho dz (lw_kabs_dry + lw_kabs_vapor qv + lw_kabs_cloud qc)
   * The cloud term is the Stephens (1978) emissivity, eps = 1 - exp(-a LWP)
   * with a = 0.158 m^2/g, written as an optical depth. The column optical
   * depth is then independent of the vertical resolution. The cloud-band,
   * moisture-coefficient and aerosol additions still apply on top.
   * Default false keeps tau_lw_per_layer.
   */
  bool lw_mass_absorption_enable = false;

  /** @brief Dry-air gray LW mass absorption coefficient [m^2/kg]. */
  amrex::Real lw_kabs_dry = 1.0e-4;

  /** @brief Water-vapor gray LW mass absorption coefficient [m^2/kg]. */
  amrex::Real lw_kabs_vapor = 0.1;

  /** @brief Cloud-water LW mass absorption coefficient [m^2/kg] (Stephens 1978: 0.158 m^2/g). */
  amrex::Real lw_kabs_cloud = 158.0;

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
   * @brief Shortwave surface albedo for diffuse light [0,1]. A negative
   * value (default) uses surface_albedo_sw for both direct and diffuse light.
   */
  amrex::Real surface_albedo_sw_diffuse = -1.0;

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
   * @brief (Phase 17) Master switch for simplified SEB infrastructure.
   * When false (default), no SEB MultiFabs are allocated and behavior is unchanged.
   */
  bool seb_enable = false;

  /**
   * @brief (Phase 17) Fallback SEB surface net shortwave flux [W/m^2].
   */
  amrex::Real seb_sw_flux_default = 0.0;

  /**
   * @brief (Phase 17) Fallback SEB surface net longwave flux [W/m^2].
   */
  amrex::Real seb_lw_flux_default = 0.0;

  /**
   * @brief (Phase 17) Fallback SEB sensible heat flux [W/m^2].
   */
  amrex::Real seb_hfx_default = 0.0;

  /**
   * @brief (Phase 17) Fallback SEB latent heat flux [W/m^2].
   */
  amrex::Real seb_lh_default = 0.0;

  /**
   * @brief (Phase 17) Fallback SEB ground heat flux [W/m^2].
   */
  amrex::Real seb_grdflx_default = 0.0;

  /**
   * @brief (Phase 17) Fallback SEB surface moisture [kg/kg].
   */
  amrex::Real seb_q_sfc_default = 0.0;

  /**
   * @brief (Phase 17) Fallback SEB deep soil temperature [K].
   */
  amrex::Real seb_t_deep_default = 300.0;

  /**
   * @brief (Phase 17) Fallback SEB deep soil moisture [kg/kg].
   */
  amrex::Real seb_q_deep_default = 0.0;

  /**
   * @brief (Phase 18) Enable diagnostic SEB residual computation and reporting.
   * When true, computes and reports SEB residual = (SW_net + LW_net) - H - LE - G.
   * When false (default), no residual computation or new diagnostic output.
   * Requires seb_enable=true to be meaningful; if false but seb_enable=true, auto-enables internally.
   * Note: This is diagnostic-only; no prognostic T_s or flux update (Phase 19).
   */
  bool seb_diagnostic_enable = false;

  /**
   * @brief (Phase 19b) Enable prognostic SEB surface temperature and moisture evolution.
   * When true, integrates SEB residual forward in time to evolve T_s and q_s using
   * force-restore formulation. When false (default), no prognostic update occurs.
   * Requires seb_enable=true and seb_diagnostic_enable=true; auto-enables internally if needed.
   * Only active when Noah-MP is NOT driving LSM fields at that level.
   */
  bool seb_prognostic_enable = false;

  /**
   * @brief (Phase 19b) Effective surface heat capacity [J/(m^2*K)].
   * Used in force-restore prognostic T_s equation: dT_s/dt = SEB_residual / C_s - (2*pi/tau)*(T_s - T_deep).
   * Default 2.0e4 J/(m^2*K) represents a thin surface-active layer.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_surface_heat_capacity = 2.0e4;

  /**
   * @brief (Phase 19b) Force-restore timescale for surface temperature [s].
   * Relaxation time constant toward T_deep in prognostic T_s equation.
   * Default 86400.0 s (1 day) provides weak damping on diurnal timescales.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_restore_timescale_s = 86400.0;

  /**
   * @brief (Phase 19b) Effective surface moisture layer depth [m].
   * Used in force-restore prognostic q_s equation: dq_s/dt = -(LE/(L_v*rho_w*d_s)) - (1/tau_q)*(q_s - q_deep).
   * Default 0.1 m represents a thin surface-active moisture layer.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_moisture_layer_depth_m = 0.1;

  /**
   * @brief (Phase 19b) Force-restore timescale for surface moisture [s].
   * Relaxation time constant toward q_deep in prognostic q_s equation.
   * Default 86400.0 s (1 day) provides weak damping on diurnal timescales.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_moisture_restore_timescale_s = 86400.0;

  /**
   * @brief (Phase 19b) Minimum clamping bound for prognostic surface temperature [K].
   * After each Euler update, T_s is clamped to [seb_prognostic_t_min_k, seb_prognostic_t_max_k]
   * to prevent instability if dt is large relative to C_s/tau.
   * Default 200.0 K (far below realistic surface temperatures).
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_t_min_k = 200.0;

  /**
   * @brief (Phase 19b) Maximum clamping bound for prognostic surface temperature [K].
   * After each Euler update, T_s is clamped to [seb_prognostic_t_min_k, seb_prognostic_t_max_k].
   * Default 340.0 K (above typical surface temperatures, but prevents unrealistic excursions).
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_t_max_k = 340.0;

  /**
   * @brief (Phase 19b) Minimum clamping bound for prognostic surface moisture [kg/kg].
   * After each Euler update, q_s is clamped to [seb_prognostic_q_min, seb_prognostic_q_max].
   * Default 0.0 kg/kg prevents negative moisture.
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_q_min = 0.0;

  /**
   * @brief (Phase 19b) Maximum clamping bound for prognostic surface moisture [kg/kg].
   * After each Euler update, q_s is clamped to [seb_prognostic_q_min, seb_prognostic_q_max].
   * Default 1.0 kg/kg ensures physically meaningful moisture.
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_q_max = 1.0;


  /**
   * @brief (Phase 12) Enable dynamic longwave optical depth diagnosis from moisture/clouds.
   * When true, LW optical depth is computed per-level from qv and qc.
   * When false (default), uses static tau_lw_per_layer (Phase 11 behavior).
   * Only used when rad_type == RadType::TwoStream and lw_enabled == true.
   */
  bool tau_lw_dynamic_enable = false;
  bool tau_sw_dynamic_enable = false;

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
   * @brief (Phase 15) Enable prescribed bulk aerosol/turbidity optical depth.
   * When true, aerosol optical depth is added on top of existing tau_base + cloud + dynamic-moisture contributions.
   * When false (default), no aerosol contribution (Phase 14B behavior preserved).
   * Only used when rad_type == RadType::TwoStream and sw_enabled == true (LW optional, Phase 16+).
   */
  bool aerosol_enable = false;

  /**
   * @brief (Phase 15) Aerosol optical depth profile type: Constant, Exponential, or Table.
   * - Constant: uniform aerosol tau per layer (aerosol_tau_per_layer).
   * - Exponential: tau decays with height via tau_aerosol(k) = aerosol_tau_surface * exp(-z(k) / aerosol_scale_height_m).
   * - Table: per-level tau prescribed via lookup (future extension, currently not supported).
   * Default Constant. Only used when aerosol_enable == true.
   */
  AerosolProfileType aerosol_profile_type = AerosolProfileType::Constant;

  /**
   * @brief (Phase 15) Constant aerosol optical depth per layer [dimensionless].
   * Used for AerosolProfileType::Constant profile. Default 0.0 (no aerosol contribution).
   * Only used when aerosol_enable == true and aerosol_profile_type == Constant.
   */
  amrex::Real aerosol_tau_per_layer = 0.0;

  /**
   * @brief (Phase 15) Scale height for exponential aerosol profile [m].
   * Controls decay rate: tau_aerosol(k) = aerosol_tau_surface * exp(-z(k) / scale_height).
   * Default 2000.0 m (typical scale height for tropospheric aerosols).
   * Only used when aerosol_enable == true and aerosol_profile_type == Exponential.
   */
  amrex::Real aerosol_scale_height_m = 2000.0;

  /**
   * @brief (Phase 15) Total-column aerosol optical depth at surface [dimensionless].
   * Used for AerosolProfileType::Exponential profile as the surface boundary condition.
   * Default 0.0 (no aerosol contribution).
   * Only used when aerosol_enable == true and aerosol_profile_type == Exponential.
   */
  amrex::Real aerosol_tau_surface = 0.0;

  /**
   * @brief (Phase 16) Enable time-varying solar geometry from diurnal cycle.
   * When true, solar zenith angle is computed dynamically from simulation time,
   * latitude, longitude, and day-of-year using astronomical formulas.
   * When false (default), uses static solar_zenith_deg (Phase 15 behavior).
   * Only used when rad_type == RadType::TwoStream and sw_enabled == true.
   */
  bool solar_geometry_dynamic_enable = false;

  /**
   * @brief (Phase 16) Site latitude [degrees, -90 to 90, positive north].
   * Used to compute solar zenith angle when solar_geometry_dynamic_enable == true.
   * Clamped to [-90, 90] in init_params().
   * Default 0.0 (equator). Only used when solar_geometry_dynamic_enable == true.
   */
  amrex::Real latitude_deg = 0.0;

  /**
   * @brief (Phase 16) Site longitude [degrees, -180 to 180, positive east].
   * Used to compute solar hour angle when solar_geometry_dynamic_enable == true.
   * Clamped to [-180, 180] in init_params().
   * Default 0.0 (prime meridian). Only used when solar_geometry_dynamic_enable == true.
   */
  amrex::Real longitude_deg = 0.0;

  /**
   * @brief (Phase 16) Reference day-of-year at simulation start [1-366].
   * Used to compute solar declination and equation-of-time correction.
   * Default 172 (approx. June 21, summer solstice in Northern Hemisphere).
   * Clamped to [1, 366] in init_params().
   * Only used when solar_geometry_dynamic_enable == true.
   */
  amrex::Real day_of_year = 172.0;

  /**
   * @brief (Phase 16) Time zone offset from UTC [hours].
   * Used to convert UTC simulation time to local solar time.
   * Positive east of Greenwich, negative west. Default 0.0 (UTC).
   * No explicit bounds; values outside ±12 hours are allowed (e.g., for
   * dateline regions). Only used when solar_geometry_dynamic_enable == true.
   */
  amrex::Real time_zone_offset_hours = 0.0;

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
   * - erf.radiation.earth_sun_distance_enable (bool, optional, default false)
   * - erf.radiation.tau_model (string "per_layer"|"mass", optional, default "per_layer")
   * - erf.radiation.sw_kabs_dry, sw_kscat_dry, sw_kabs_vapor, sw_kext_cloud (real, optional, mass model)
   * - erf.radiation.sw_cloud_omega, sw_cloud_g (real, optional, mass model)
   * - erf.radiation.lw_mass_absorption_enable (bool, optional, default false)
   * - erf.radiation.lw_kabs_dry (real, optional, default 1.0e-4 m^2/kg)
   * - erf.radiation.lw_kabs_vapor (real, optional, default 0.1 m^2/kg)
   * - erf.radiation.lw_kabs_cloud (real, optional, default 158.0 m^2/kg)
   * - erf.radiation.surface_albedo_sw_diffuse (real, optional, default -1 = same as direct)
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
   * - erf.radiation.aerosol_enable (bool, optional, Phase 15, default false)
   * - erf.radiation.aerosol_profile_type (string "constant"|"exponential"|"table", optional, Phase 15, default "constant")
   * - erf.radiation.aerosol_tau_per_layer (real, optional, Phase 15, default 0.0)
   * - erf.radiation.aerosol_scale_height_m (real, optional, Phase 15, default 2000.0)
   * - erf.radiation.aerosol_tau_surface (real, optional, Phase 15, default 0.0)
   * - erf.radiation.solar_geometry_dynamic_enable (bool, optional, Phase 16, default false)
   * - erf.radiation.latitude_deg (real, optional, Phase 16, default 0.0, clamped to [-90, 90])
   * - erf.radiation.longitude_deg (real, optional, Phase 16, default 0.0, clamped to [-180, 180])
   * - erf.radiation.day_of_year (real, optional, Phase 16, default 172.0, clamped to [1, 366])
   * - erf.radiation.time_zone_offset_hours (real, optional, Phase 16, default 0.0)
   * - erf.radiation.seb_enable (bool, optional, Phase 17, default false)
   * - erf.radiation.seb_sw_flux_default (real, optional, Phase 17, default 0.0)
   * - erf.radiation.seb_lw_flux_default (real, optional, Phase 17, default 0.0)
   * - erf.radiation.seb_hfx_default (real, optional, Phase 17, default 0.0)
   * - erf.radiation.seb_lh_default (real, optional, Phase 17, default 0.0)
   * - erf.radiation.seb_grdflx_default (real, optional, Phase 17, default 0.0)
   * - erf.radiation.seb_q_sfc_default (real, optional, Phase 17, default 0.0)
   * - erf.radiation.seb_t_deep_default (real, optional, Phase 17, default 300.0)
   * - erf.radiation.seb_q_deep_default (real, optional, Phase 17, default 0.0)
   * - erf.radiation.seb_diagnostic_enable (bool, optional, Phase 18, default false)
   * - erf.radiation.seb_prognostic_enable (bool, optional, Phase 19b, default false)
   * - erf.radiation.seb_surface_heat_capacity (real, optional, Phase 19b, default 2.0e4, validated > 0)
   * - erf.radiation.seb_restore_timescale_s (real, optional, Phase 19b, default 86400.0, validated > 0)
   * - erf.radiation.seb_moisture_layer_depth_m (real, optional, Phase 19b, default 0.1, validated > 0)
   * - erf.radiation.seb_moisture_restore_timescale_s (real, optional, Phase 19b, default 86400.0, validated > 0)
   * - erf.radiation.seb_prognostic_t_min_k (real, optional, Phase 19b, default 200.0)
   * - erf.radiation.seb_prognostic_t_max_k (real, optional, Phase 19b, default 340.0)
   * - erf.radiation.seb_prognostic_q_min (real, optional, Phase 19b, default 0.0, clamped to [0, 1])
   * - erf.radiation.seb_prognostic_q_max (real, optional, Phase 19b, default 1.0, clamped to [0, 1])
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
      pp.query("radiation.earth_sun_distance_enable", earth_sun_distance_enable);
      std::string tau_model_str = "per_layer";
      pp.query("radiation.tau_model", tau_model_str);
      if (tau_model_str == "mass") {
        tau_model = TauModel::Mass;
      } else {
        tau_model = TauModel::PerLayer;
      }
      pp.query("radiation.sw_kabs_dry", sw_kabs_dry);
      pp.query("radiation.sw_kscat_dry", sw_kscat_dry);
      pp.query("radiation.sw_kabs_vapor", sw_kabs_vapor);
      pp.query("radiation.sw_kext_cloud", sw_kext_cloud);
      pp.query("radiation.sw_cloud_omega", sw_cloud_omega);
      pp.query("radiation.sw_cloud_g", sw_cloud_g);
      if (sw_kabs_dry < 0.0) sw_kabs_dry = 0.0;
      if (sw_kscat_dry < 0.0) sw_kscat_dry = 0.0;
      if (sw_kabs_vapor < 0.0) sw_kabs_vapor = 0.0;
      if (sw_kext_cloud < 0.0) sw_kext_cloud = 0.0;
      if (sw_cloud_omega < 0.0) sw_cloud_omega = 0.0;
      if (sw_cloud_omega > 1.0) sw_cloud_omega = 1.0;
      if (sw_cloud_g < -1.0) sw_cloud_g = -1.0;
      if (sw_cloud_g > 1.0) sw_cloud_g = 1.0;
      pp.query("radiation.lw_mass_absorption_enable", lw_mass_absorption_enable);
      pp.query("radiation.lw_kabs_dry", lw_kabs_dry);
      pp.query("radiation.lw_kabs_vapor", lw_kabs_vapor);
      pp.query("radiation.lw_kabs_cloud", lw_kabs_cloud);
      if (lw_kabs_dry < 0.0) lw_kabs_dry = 0.0;
      if (lw_kabs_vapor < 0.0) lw_kabs_vapor = 0.0;
      if (lw_kabs_cloud < 0.0) lw_kabs_cloud = 0.0;
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
      pp.query("radiation.surface_albedo_sw_diffuse", surface_albedo_sw_diffuse);
      if (surface_albedo_sw_diffuse > 1.0) surface_albedo_sw_diffuse = 1.0;

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

      // Phase 15: prescribed bulk aerosol/turbidity optical depth
      // (all optional; defaults disable aerosol, preserving Phase 14B behavior)
      pp.query("radiation.aerosol_enable", aerosol_enable);
       
      // Parse aerosol profile type as string
      std::string aerosol_profile_type_str = "constant";
      pp.query("radiation.aerosol_profile_type", aerosol_profile_type_str);
      if (aerosol_profile_type_str == "exponential") {
        aerosol_profile_type = AerosolProfileType::Exponential;
      } else if (aerosol_profile_type_str == "table") {
        aerosol_profile_type = AerosolProfileType::Table;
      } else {
        aerosol_profile_type = AerosolProfileType::Constant;
      }
       
      pp.query("radiation.aerosol_tau_per_layer", aerosol_tau_per_layer);
      pp.query("radiation.aerosol_scale_height_m", aerosol_scale_height_m);
      pp.query("radiation.aerosol_tau_surface", aerosol_tau_surface);

      // Parameter validation (Phase 15): clamp aerosol parameters
      // to physically reasonable ranges
      if (aerosol_tau_per_layer < 0.0) aerosol_tau_per_layer = 0.0;
      if (aerosol_tau_surface < 0.0) aerosol_tau_surface = 0.0;
      if (aerosol_scale_height_m <= 0.0) aerosol_scale_height_m = 2000.0;  // Ensure positive

      // Phase 16: time-varying solar geometry from diurnal cycle
      // (all optional; defaults disable dynamic solar geometry, preserving Phase 15 behavior)
      pp.query("radiation.solar_geometry_dynamic_enable", solar_geometry_dynamic_enable);
      pp.query("radiation.latitude_deg", latitude_deg);
      pp.query("radiation.longitude_deg", longitude_deg);
      pp.query("radiation.day_of_year", day_of_year);
      pp.query("radiation.time_zone_offset_hours", time_zone_offset_hours);

      // Parameter validation (Phase 16): clamp latitude/longitude/day-of-year to valid ranges
      // Latitude must be in [-90, 90]
      if (!std::isfinite(latitude_deg)) latitude_deg = 0.0;
      if (latitude_deg < -90.0) latitude_deg = -90.0;
      if (latitude_deg > 90.0) latitude_deg = 90.0;

      // Longitude must be in [-180, 180]; wrap if outside
      if (!std::isfinite(longitude_deg)) longitude_deg = 0.0;
      while (longitude_deg > 180.0) longitude_deg -= 360.0;
      while (longitude_deg < -180.0) longitude_deg += 360.0;

      // Day-of-year must be in [1, 366]
      if (!std::isfinite(day_of_year)) day_of_year = 172.0;
      if (day_of_year < 1.0) day_of_year = 1.0;
      if (day_of_year > 366.0) day_of_year = 366.0;

      // Time zone offset: no strict bounds, but ensure finite
      if (!std::isfinite(time_zone_offset_hours)) time_zone_offset_hours = 0.0;

      // Phase 17: simplified SEB infrastructure fallback parameters
      pp.query("radiation.seb_enable", seb_enable);
      pp.query("radiation.seb_sw_flux_default", seb_sw_flux_default);
      pp.query("radiation.seb_lw_flux_default", seb_lw_flux_default);
      pp.query("radiation.seb_hfx_default", seb_hfx_default);
      pp.query("radiation.seb_lh_default", seb_lh_default);
      pp.query("radiation.seb_grdflx_default", seb_grdflx_default);
      pp.query("radiation.seb_q_sfc_default", seb_q_sfc_default);
      pp.query("radiation.seb_t_deep_default", seb_t_deep_default);
      pp.query("radiation.seb_q_deep_default", seb_q_deep_default);

      // Phase 18: SEB diagnostic mode
      pp.query("radiation.seb_diagnostic_enable", seb_diagnostic_enable);

      // Auto-enable seb_enable if seb_diagnostic_enable is requested but seb_enable is false
      if (seb_diagnostic_enable && !seb_enable) {
          seb_enable = true;  // Internal auto-enable for safety/simplicity
      }

      // Phase 19b: SEB prognostic mode parameters
      pp.query("radiation.seb_prognostic_enable", seb_prognostic_enable);
      pp.query("radiation.seb_surface_heat_capacity", seb_surface_heat_capacity);
      pp.query("radiation.seb_restore_timescale_s", seb_restore_timescale_s);
      pp.query("radiation.seb_moisture_layer_depth_m", seb_moisture_layer_depth_m);
      pp.query("radiation.seb_moisture_restore_timescale_s", seb_moisture_restore_timescale_s);
      pp.query("radiation.seb_prognostic_t_min_k", seb_prognostic_t_min_k);
      pp.query("radiation.seb_prognostic_t_max_k", seb_prognostic_t_max_k);
      pp.query("radiation.seb_prognostic_q_min", seb_prognostic_q_min);
      pp.query("radiation.seb_prognostic_q_max", seb_prognostic_q_max);

      // Auto-enable seb_enable and seb_diagnostic_enable if seb_prognostic_enable is requested
      if (seb_prognostic_enable && !seb_enable) {
          seb_enable = true;
      }
      if (seb_prognostic_enable && !seb_diagnostic_enable) {
          seb_diagnostic_enable = true;
      }

      if (!std::isfinite(seb_sw_flux_default)) seb_sw_flux_default = 0.0;
      if (!std::isfinite(seb_lw_flux_default)) seb_lw_flux_default = 0.0;
      if (!std::isfinite(seb_hfx_default)) seb_hfx_default = 0.0;
      if (!std::isfinite(seb_lh_default)) seb_lh_default = 0.0;
      if (!std::isfinite(seb_grdflx_default)) seb_grdflx_default = 0.0;
      if (!std::isfinite(seb_q_sfc_default)) seb_q_sfc_default = 0.0;
      if (!std::isfinite(seb_t_deep_default) || seb_t_deep_default <= 0.0) seb_t_deep_default = 300.0;
      if (!std::isfinite(seb_q_deep_default)) seb_q_deep_default = 0.0;

      seb_sw_flux_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_sw_flux_default));
      seb_lw_flux_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_lw_flux_default));
      seb_hfx_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_hfx_default));
      seb_lh_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_lh_default));
      seb_grdflx_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_grdflx_default));
      seb_q_sfc_default = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), seb_q_sfc_default));
      seb_q_deep_default = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), seb_q_deep_default));

      // Parameter validation (Phase 19b): validate prognostic parameters
      if (!std::isfinite(seb_surface_heat_capacity) || seb_surface_heat_capacity <= 0.0) {
          seb_surface_heat_capacity = 2.0e4;  // Reset to default
      }
      if (!std::isfinite(seb_restore_timescale_s) || seb_restore_timescale_s <= 0.0) {
          seb_restore_timescale_s = 86400.0;  // Reset to default
      }
      if (!std::isfinite(seb_moisture_layer_depth_m) || seb_moisture_layer_depth_m <= 0.0) {
          seb_moisture_layer_depth_m = 0.1;  // Reset to default
      }
      if (!std::isfinite(seb_moisture_restore_timescale_s) || seb_moisture_restore_timescale_s <= 0.0) {
          seb_moisture_restore_timescale_s = 86400.0;  // Reset to default
      }

      // Clamp temperature bounds
      if (!std::isfinite(seb_prognostic_t_min_k)) seb_prognostic_t_min_k = 200.0;
      if (!std::isfinite(seb_prognostic_t_max_k)) seb_prognostic_t_max_k = 340.0;
      if (seb_prognostic_t_min_k > seb_prognostic_t_max_k) {
          amrex::Real tmp = seb_prognostic_t_min_k;
          seb_prognostic_t_min_k = seb_prognostic_t_max_k;
          seb_prognostic_t_max_k = tmp;
      }

      // Clamp moisture bounds
      if (!std::isfinite(seb_prognostic_q_min)) seb_prognostic_q_min = 0.0;
      if (!std::isfinite(seb_prognostic_q_max)) seb_prognostic_q_max = 1.0;
      seb_prognostic_q_min = amrex::max(0.0, amrex::min(1.0, seb_prognostic_q_min));
      seb_prognostic_q_max = amrex::max(0.0, amrex::min(1.0, seb_prognostic_q_max));
      if (seb_prognostic_q_min > seb_prognostic_q_max) {
          amrex::Real tmp = seb_prognostic_q_min;
          seb_prognostic_q_min = seb_prognostic_q_max;
          seb_prognostic_q_max = tmp;
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
