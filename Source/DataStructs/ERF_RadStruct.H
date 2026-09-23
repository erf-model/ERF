#ifndef ERF_RAD_STRUCT_H_
#define ERF_RAD_STRUCT_H_

#include <cmath>
#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <ERF_IndexDefines.H>
#include <ERF_TurbStruct.H>

/**
 * @file ERF_RadStruct.H
 * @brief Radiation model type and control parameters.
 */

/**
 * @brief Optical depth profile type for the two-stream solver.
 *
 * - Constant: tau_per_layer / tau_lw_per_layer applied uniformly at every
 *   level, independent of layer thickness.
 * - CloudLayer: an additional cloud optical depth is added within
 *   [cloud_base_height_m, cloud_top_height_m], on top of the clear-sky
 *   background tau_per_layer / tau_lw_per_layer.
 */
AMREX_ENUM(TauProfileType, Constant, CloudLayer);

/**
 * @brief Aerosol optical depth profile type for the two-stream solver.
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
 *
 * The two-stream solver builds up from a clear-sky, no-scattering column of
 * uniform optical depth. On top of that it offers height-varying optical
 * depth (a cloud layer) with cloud-fraction masking, which blends clear-sky
 * and cloudy-column fluxes, and a diffuse shortwave component from a
 * two-stream scattering approximation parameterized by single-scattering
 * albedo and asymmetry factor, with independent clear-sky and cloud values.
 * Every addition is off by default, so the defaults reproduce the plain
 * clear-sky column.
 */
struct RadChoice
{
public:
  /**
   * @brief True when erf.radiation_model = TwoStream. Set by init_params from
   * the model SolverChoice read; the erf.radiation.* sub-options are read only then.
   */
  bool enabled = false;

  /**
   * @brief Enable shortwave (solar) radiation computation.
   * Only used by the two-stream model.
   */
  bool sw_enabled = true;

  /**
   * @brief Enable longwave (thermal) radiation computation.
   * Only used by the two-stream model.
   */
  bool lw_enabled = true;

  /**
   * @brief Shortwave optical depth per layer (constant for all layers).
   * Used in Beer-Lambert direct-beam formula.
   * Only used by the two-stream model.
   */
  amrex::Real tau_per_layer = 0.05;

  /**
   * @brief Longwave optical depth per layer (constant for all layers).
   * Used in gray-gas LW two-stream solver.
   * Only used by the two-stream model.
   */
  amrex::Real tau_lw_per_layer = 1.0;

  /**
   * @name Inputs shared with the RRTMGP interface
   * The two-stream model reads the sun, site and surface-temperature inputs
   * under the same keys, with the same meaning, as the RRTMGP interface
   * (Source/Radiation/RRTMGP/ERF_Radiation.H), so a deck switches
   * between the two models by changing erf.radiation_model alone.
   * @{
   */

  /**
   * @brief erf.fixed_solar_zenith_angle: the COSINE of the solar zenith
   * angle (mu0), as RRTMGP takes it, applied to every column when > 0.
   * Otherwise (default) the sun follows the calendar: start_datetime sets
   * the date and time, the position is evaluated at rad_cons_lat/lon or at
   * the per-column latitude and longitude when the grid carries them.
   */
  amrex::Real fixed_solar_zenith_angle = -9999.0;

  /**
   * @brief erf.fixed_total_solar_irradiance: top-of-atmosphere irradiance
   * [W/m^2] when >= 0. Otherwise (default) 1360.9 W/m^2 scaled by the
   * Earth-Sun distance factor of the date, as RRTMGP does.
   */
  amrex::Real fixed_total_solar_irradiance = -9999.0;

  /**
   * @brief erf.rad_t_sfc: surface temperature [K] used as the longwave
   * boundary condition wherever no land-surface model, surface layer or
   * prognostic surface energy balance supplies one; also the initial value
   * of the prognostic surface temperature. Required, as for RRTMGP.
   */
  amrex::Real rad_t_sfc = -1.0;

  /** @brief erf.rad_cons_lat / erf.rad_cons_lon: site latitude and longitude [degrees] when the grid has no lat/lon fields (RRTMGP defaults). */
  amrex::Real rad_cons_lat = 39.809860;
  amrex::Real rad_cons_lon = -98.555183;

  /**
   * @brief erf.rad_orbital_year and the orbital overrides: negative means
   * the year of the timestamp and the Berger (1978) values for that year.
   */
  int rad_orbital_year = -9999;
  amrex::Real rad_orbital_eccentricity = -9999.0;
  amrex::Real rad_orbital_obliquity = -9999.0;
  amrex::Real rad_orbital_mvelp = -9999.0;
  /** @} */

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
   * @brief Master switch for radiation diagnostics emission.
   * When false, no diagnostics output (stdout/tagged/regtest/csv).
   * Default false: a run that did not ask for the diagnostics should not
   * write a CSV row and a stdout block on every step. The per-stream
   * switches below stay on by default, so setting this one input to true
   * turns the whole set on.
   */
  bool diag_enable = false;

  /**
   * @brief Enable human-readable stdout diagnostics block.
   * Controls tagged [RAD][...] debug line emission (when verbosity >= 1).
   * Default true.
   */
  bool diag_stdout_enable = true;

  /**
   * @brief Enable tagged [RAD][...] debug lines.
   * Controls emission of bracketed-tag debug output.
   * Default true.
   */
  bool diag_tagged_enable = true;

  /**
   * @brief Enable RADIATION_DIAG: line emission.
   * Controls regtest-parsing diagnostics line.
   * Default true.
   */
  bool diag_regtest_line_enable = true;

  /**
   * @brief Enable CSV file append behavior.
   * When false, no file writes; other enabled streams may still print.
   * Default true.
   */
  bool diag_csv_enable = true;

  /**
   * @brief Call-site filtering mode: "both" (pre+post rows),
   * "pre_only", or "post_only".
   */
  std::string diag_callsite_mode = "both";

  /**
   * @brief Tolerance for time equality in duplicate guard [s].
   * Times within this tolerance are considered equal.
   * Default 1e-12.
   */
  amrex::Real diag_dedup_tol = 1.0e-12;

  /**
   * @brief Optical depth profile type: "constant" (default,
   * the thickness-independent tau) or "cloud_layer" (adds height-varying
   * cloud optical depth between cloud_base_height_m and cloud_top_height_m).
   */
  TauProfileType tau_profile_type = TauProfileType::Constant;

  /**
   * @brief Cloud layer base height [m]. Only used when
   * tau_profile_type == CloudLayer.
   */
  amrex::Real cloud_base_height_m = 500.0;

  /**
   * @brief Cloud layer top height [m]. Only used when
   * tau_profile_type == CloudLayer. Must be >= cloud_base_height_m.
   */
  amrex::Real cloud_top_height_m = 1000.0;

  /**
   * @brief Additional optical depth per layer contributed by the
   * cloud, added on top of the clear-sky tau_per_layer/tau_lw_per_layer
   * background within [cloud_base_height_m, cloud_top_height_m]. Only used
   * when tau_profile_type == CloudLayer.
   */
  amrex::Real cloud_tau_per_layer = 0.5;

  /**
   * @brief Cloud fraction in [0, 1] used to blend clear-sky and
   * cloudy-column fluxes: F = (1 - cloud_fraction) * F_clear +
   * cloud_fraction * F_cloudy. Default 0.0 (fully clear sky), which
   * reduces exactly to the clear-sky column regardless of tau_profile_type.
   */
  amrex::Real cloud_fraction = 0.0;

  /**
   * @brief Clear-sky (background gas) single-scattering albedo
   * in [0, 1] for the SW diffuse (scattering) two-stream calculation.
   * Default 0.0 means the clear-sky column is purely absorbing, and the
   * diffuse SW contribution from the clear-sky column is exactly zero —
   * reducing exactly to the direct-beam-only result.
   * Only used by the two-stream model and sw_enabled == true.
   */
  amrex::Real single_scattering_albedo = 0.0;

  /**
   * @brief Clear-sky (background gas) scattering asymmetry
   * factor in [-1, 1] (0 = isotropic scattering, close to 1 = strongly
   * forward-scattering). Only meaningful when single_scattering_albedo
   * > 0. Default 0.0.
   */
  amrex::Real asymmetry_factor = 0.0;

  /**
   * @brief Cloud single-scattering albedo in [0, 1], used
   * instead of single_scattering_albedo for levels where the
   * cloud-layer optical depth enhancement applies (i.e., within
   * [cloud_base_height_m, cloud_top_height_m] on the "cloudy" column
   * evaluation). Default 0.0 leaves the cloud band non-scattering
   * for SW_Cloud_Layer unless explicitly configured otherwise. Typical
   * liquid water cloud values are close to 1 (e.g., 0.9999).
   * Only used when tau_profile_type == CloudLayer.
   */
  amrex::Real cloud_single_scattering_albedo = 0.0;

  /**
   * @brief Cloud scattering asymmetry factor in [-1, 1], used
   * instead of asymmetry_factor for levels where the cloud-layer
   * optical depth enhancement applies. Default 0.0. Typical liquid water
   * cloud droplet values are around 0.85 (strongly forward-scattering).
   * Only meaningful when cloud_single_scattering_albedo > 0.
   */
  amrex::Real cloud_asymmetry_factor = 0.0;

  /**
   * @brief Shortwave surface albedo fallback [0, 1] for per-column
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
   * @brief Longwave surface emissivity fallback [0, 1] for per-column
   * heterogeneous-surface-property resolution. Used when per-column LSM/radiation
   * interface emissivity field is unavailable or invalid. Default 0.99 represents
   * typical land/water/ice emissivity (nearly blackbody). Only used when hetero
   * surface fields not present or disabled.
   */
  amrex::Real surface_emissivity_lw = 0.99;

  /**
   * @brief Master switch for simplified SEB infrastructure.
   * When false (default), no SEB MultiFabs are allocated and behavior is unchanged.
   */
  bool seb_enable = false;

  /**
   * @brief Take the SEB's net surface shortwave and longwave fluxes from the
   * two-stream sweep itself (per column, absorbed fluxes positive) wherever
   * the LSM does not expose them. When false (default) those fluxes are the
   * LSM fields or, without an LSM, the seb_sw_flux_default and
   * seb_lw_flux_default constants. Sensible, latent and ground heat fluxes
   * are unaffected.
   */
  bool seb_use_radiation_fluxes = false;

  /**
   * @brief Fallback SEB surface net shortwave flux [W/m^2].
   */
  amrex::Real seb_sw_flux_default = 0.0;

  /**
   * @brief Fallback SEB surface net longwave flux [W/m^2].
   */
  amrex::Real seb_lw_flux_default = 0.0;

  /**
   * @brief Fallback SEB sensible heat flux [W/m^2].
   */
  amrex::Real seb_hfx_default = 0.0;

  /**
   * @brief Fallback SEB latent heat flux [W/m^2].
   */
  amrex::Real seb_lh_default = 0.0;

  /**
   * @brief Fallback SEB ground heat flux [W/m^2].
   */
  amrex::Real seb_grdflx_default = 0.0;

  /**
   * @brief Fallback SEB surface moisture [kg/kg].
   */
  amrex::Real seb_q_sfc_default = 0.0;

  /**
   * @brief Fallback SEB deep soil temperature [K].
   */
  amrex::Real seb_t_deep_default = 300.0;

  /**
   * @brief Fallback SEB deep soil moisture [kg/kg].
   */
  amrex::Real seb_q_deep_default = 0.0;

  /**
   * @brief Enable diagnostic SEB residual computation and reporting.
   * When true, computes and reports SEB residual = (SW_net + LW_net) - H - LE - G.
   * When false (default), no residual computation or new diagnostic output.
   * Requires seb_enable=true to be meaningful; if false but seb_enable=true, auto-enables internally.
   * Note: this is diagnostic-only; it makes no prognostic T_s or flux update.
   * Set seb_prognostic_enable for that.
   */
  bool seb_diagnostic_enable = false;

  /**
   * @brief Enable prognostic SEB surface temperature and moisture evolution.
   * When true, integrates SEB residual forward in time to evolve T_s and q_s using
   * force-restore formulation. When false (default), no prognostic update occurs.
   * Requires seb_enable=true and seb_diagnostic_enable=true; auto-enables internally if needed.
   * Only active when Noah-MP is NOT driving LSM fields at that level.
   */
  bool seb_prognostic_enable = false;

  /**
   * @brief Effective surface heat capacity [J/(m^2*K)].
   * Used in force-restore prognostic T_s equation: dT_s/dt = SEB_residual / C_s - (2*pi/tau)*(T_s - T_deep).
   * Default 2.0e4 J/(m^2*K) represents a thin surface-active layer.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_surface_heat_capacity = 2.0e4;

  /**
   * @brief Force-restore timescale for surface temperature [s].
   * Relaxation time constant toward T_deep in prognostic T_s equation.
   * Default 86400.0 s (1 day) provides weak damping on diurnal timescales.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_restore_timescale_s = 86400.0;

  /**
   * @brief Effective surface moisture layer depth [m].
   * Used in force-restore prognostic q_s equation: dq_s/dt = -(LE/(L_v*rho_w*d_s)) - (1/tau_q)*(q_s - q_deep).
   * Default 0.1 m represents a thin surface-active moisture layer.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_moisture_layer_depth_m = 0.1;

  /**
   * @brief Force-restore timescale for surface moisture [s].
   * Relaxation time constant toward q_deep in prognostic q_s equation.
   * Default 86400.0 s (1 day) provides weak damping on diurnal timescales.
   * Only used when seb_prognostic_enable=true. Validated > 0 in init_params().
   */
  amrex::Real seb_moisture_restore_timescale_s = 86400.0;

  /**
   * @brief Minimum clamping bound for prognostic surface temperature [K].
   * After each Euler update, T_s is clamped to [seb_prognostic_t_min_k, seb_prognostic_t_max_k]
   * to prevent instability if dt is large relative to C_s/tau.
   * Default 200.0 K (far below realistic surface temperatures).
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_t_min_k = 200.0;

  /**
   * @brief Maximum clamping bound for prognostic surface temperature [K].
   * After each Euler update, T_s is clamped to [seb_prognostic_t_min_k, seb_prognostic_t_max_k].
   * Default 340.0 K (above typical surface temperatures, but prevents unrealistic excursions).
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_t_max_k = 340.0;

  /**
   * @brief Minimum clamping bound for prognostic surface moisture [kg/kg].
   * After each Euler update, q_s is clamped to [seb_prognostic_q_min, seb_prognostic_q_max].
   * Default 0.0 kg/kg prevents negative moisture.
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_q_min = 0.0;

  /**
   * @brief Maximum clamping bound for prognostic surface moisture [kg/kg].
   * After each Euler update, q_s is clamped to [seb_prognostic_q_min, seb_prognostic_q_max].
   * Default 1.0 kg/kg ensures physically meaningful moisture.
   * Only used when seb_prognostic_enable=true.
   */
  amrex::Real seb_prognostic_q_max = 1.0;


  /**
   * @brief Enable dynamic longwave optical depth diagnosis from moisture/clouds.
   * When true, LW optical depth is computed per-level from qv and qc.
   * When false (default), uses the static tau_lw_per_layer.
   * Only used by the two-stream model and lw_enabled == true.
   */
  bool tau_lw_dynamic_enable = false;
  bool tau_sw_dynamic_enable = false;

  /**
   * @brief Shortwave optical depth coefficient for water vapor [dimensionless].
   * Dynamic SW tau = tau_sw_coeff_qv * qv + tau_sw_coeff_qc * qc + tau_per_layer.
   * Default 0.0 means tau is purely from static tau_per_layer (no dynamic component).
   * Only used when tau_sw_dynamic_enable == true.
   */
  amrex::Real tau_sw_coeff_qv = 0.0;

  /**
   * @brief Shortwave optical depth coefficient for cloud liquid water [dimensionless].
   * Default 0.0. When both coefficients are 0, dynamic path still invoked but yields
   * tau_per_layer unchanged.
   * Only used when tau_sw_dynamic_enable == true.
   */
  amrex::Real tau_sw_coeff_qc = 0.0;

  /**
   * @brief Longwave optical depth coefficient for water vapor [dimensionless].
   * Dynamic LW tau = tau_lw_coeff_qv * qv + tau_lw_coeff_qc * qc + tau_lw_per_layer.
   * Default 0.0 means tau is purely from static tau_lw_per_layer (no dynamic component).
   * Only used when tau_lw_dynamic_enable == true.
   */
  amrex::Real tau_lw_coeff_qv = 0.0;

  /**
   * @brief Longwave optical depth coefficient for cloud liquid water [dimensionless].
   * Default 0.0. When both coefficients are 0, dynamic path still invoked but yields
   * tau_lw_per_layer unchanged.
   * Only used when tau_lw_dynamic_enable == true.
   */
  amrex::Real tau_lw_coeff_qc = 0.0;

  /**
   * @brief Enable prognostic cloud fraction diagnosis from RH/qc.
   * When true, cloud fraction is computed per-level from relative humidity and cloud water.
   * When false (default), uses the static cloud_fraction.
   * Only used by the two-stream model.
   */
  bool cloud_fraction_prog_enable = false;

  /**
   * @brief Minimum relative humidity threshold for cloud fraction diagnosis [0, 1].
   * Cloud fraction ramps from 0 at rh < cloud_fraction_rh_min to 1 at rh >= cloud_fraction_rh_max.
   * Default 0.0 means cloud fraction depends primarily on qc content.
   * Only used when cloud_fraction_prog_enable == true.
   */
  amrex::Real cloud_fraction_rh_min = 0.0;

  /**
   * @brief Maximum relative humidity threshold for cloud fraction diagnosis [0, 1].
   * Must be >= cloud_fraction_rh_min. Cloud fraction saturates to 1 at this RH.
   * Default 1.0 means cloud fraction increases linearly from rh_min to 1.0 RH.
   * Only used when cloud_fraction_prog_enable == true.
   */
  amrex::Real cloud_fraction_rh_max = 1.0;

  /**
   * @brief Cloud water [kg/kg] at which the liquid-water term alone gives a
   * cloud fraction of 1: cf += min(1, qc / cloud_fraction_qc_scale).
   * Default 1.0e-3 kg/kg, so a typical stratus cloud water saturates the term.
   * Only used when cloud_fraction_prog_enable == true.
   */
  amrex::Real cloud_fraction_qc_scale = 1.0e-3;

  /**
   * @brief Enable prescribed bulk aerosol/turbidity optical depth.
   * When true, aerosol optical depth is added on top of existing tau_base + cloud + dynamic-moisture contributions.
   * When false (default), there is no aerosol contribution.
   * Only used by the two-stream model and sw_enabled == true.
   */
  bool aerosol_enable = false;

  /**
   * @brief Aerosol optical depth profile type: Constant, Exponential, or Table.
   * - Constant: uniform aerosol tau per layer (aerosol_tau_per_layer).
   * - Exponential: tau decays with height via tau_aerosol(k) = aerosol_tau_surface * exp(-z(k) / aerosol_scale_height_m).
   * - Table: per-level tau prescribed via lookup (future extension, currently not supported).
   * Default Constant. Only used when aerosol_enable == true.
   */
  AerosolProfileType aerosol_profile_type = AerosolProfileType::Constant;

  /**
   * @brief Constant aerosol optical depth per layer [dimensionless].
   * Used for AerosolProfileType::Constant profile. Default 0.0 (no aerosol contribution).
   * Only used when aerosol_enable == true and aerosol_profile_type == Constant.
   */
  amrex::Real aerosol_tau_per_layer = 0.0;

  /**
   * @brief Scale height for exponential aerosol profile [m].
   * Controls decay rate: tau_aerosol(k) = aerosol_tau_surface * exp(-z(k) / scale_height).
   * Default 2000.0 m (typical scale height for tropospheric aerosols).
   * Only used when aerosol_enable == true and aerosol_profile_type == Exponential.
   */
  amrex::Real aerosol_scale_height_m = 2000.0;

  /**
   * @brief Total-column aerosol optical depth at surface [dimensionless].
   * Used for AerosolProfileType::Exponential profile as the surface boundary condition.
   * Default 0.0 (no aerosol contribution).
   * Only used when aerosol_enable == true and aerosol_profile_type == Exponential.
   */
  amrex::Real aerosol_tau_surface = 0.0;

  /**
   * @brief Read radiation options from the input parameter database.
   *
   * @param[in] lev AMR level index.
   * @param[in] max_level Maximum AMR level configured for the run.
   * @param[in] pp_prefix ParmParse prefix for the ERF input namespace (e.g., "erf").
   *
   * Queries:
   * - erf.radiation_model = TwoStream (read by SolverChoice) selects the model
   * - erf.radiation.sw_enabled (bool, optional, default true)
   * - erf.radiation.lw_enabled (bool, optional, default true)
   * - erf.radiation.tau_per_layer (real, optional)
   * - erf.radiation.tau_lw_per_layer (real, optional)
   * - erf.fixed_solar_zenith_angle, erf.fixed_total_solar_irradiance, erf.rad_t_sfc
   *   (required), erf.rad_cons_lat, erf.rad_cons_lon, erf.rad_orbital_year,
   *   erf.rad_orbital_eccentricity, erf.rad_orbital_obliquity, erf.rad_orbital_mvelp:
   *   the RRTMGP interface's keys, with its meanings (see the members above)
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
   * - erf.radiation.tau_profile_type (string "constant"|"cloud_layer", optional)
   * - erf.radiation.cloud_base_height_m (real, optional)
   * - erf.radiation.cloud_top_height_m (real, optional)
   * - erf.radiation.cloud_tau_per_layer (real, optional)
   * - erf.radiation.cloud_fraction (real, optional)
   * - erf.radiation.single_scattering_albedo (real, optional)
   * - erf.radiation.asymmetry_factor (real, optional)
   * - erf.radiation.cloud_single_scattering_albedo (real, optional)
   * - erf.radiation.cloud_asymmetry_factor (real, optional)
   * - erf.radiation.diag_enable (bool, optional, default false)
   * - erf.radiation.diag_stdout_enable (bool, optional, default true)
   * - erf.radiation.diag_tagged_enable (bool, optional, default true)
   * - erf.radiation.diag_regtest_line_enable (bool, optional, default true)
   * - erf.radiation.diag_csv_enable (bool, optional, default true)
   * - erf.radiation.diag_callsite_mode (string, optional, default "both")
   * - erf.radiation.diag_dedup_tol (real, optional, default 1e-12)
   * - erf.radiation.tau_sw_dynamic_enable (bool, optional, default false)
   * - erf.radiation.tau_lw_dynamic_enable (bool, optional, default false)
   * - erf.radiation.tau_sw_coeff_qv (real, optional, default 0.0)
   * - erf.radiation.tau_sw_coeff_qc (real, optional, default 0.0)
   * - erf.radiation.tau_lw_coeff_qv (real, optional, default 0.0)
   * - erf.radiation.tau_lw_coeff_qc (real, optional, default 0.0)
   * - erf.radiation.cloud_fraction_prog_enable (bool, optional, default false)
   * - erf.radiation.cloud_fraction_rh_min (real, optional, default 0.0)
   * - erf.radiation.cloud_fraction_rh_max (real, optional, default 1.0)
   * - erf.radiation.cloud_fraction_qc_scale (real, optional, default 1.0e-3)
   * - erf.radiation.aerosol_enable (bool, optional, default false)
   * - erf.radiation.aerosol_profile_type (string "constant"|"exponential"|"table", optional, default "constant")
   * - erf.radiation.aerosol_tau_per_layer (real, optional, default 0.0)
   * - erf.radiation.aerosol_scale_height_m (real, optional, default 2000.0)
   * - erf.radiation.aerosol_tau_surface (real, optional, default 0.0)
   * - erf.radiation.seb_enable (bool, optional, default false)
   * - erf.radiation.seb_sw_flux_default (real, optional, default 0.0)
   * - erf.radiation.seb_lw_flux_default (real, optional, default 0.0)
   * - erf.radiation.seb_hfx_default (real, optional, default 0.0)
   * - erf.radiation.seb_lh_default (real, optional, default 0.0)
   * - erf.radiation.seb_grdflx_default (real, optional, default 0.0)
   * - erf.radiation.seb_q_sfc_default (real, optional, default 0.0)
   * - erf.radiation.seb_t_deep_default (real, optional, default 300.0)
   * - erf.radiation.seb_q_deep_default (real, optional, default 0.0)
   * - erf.radiation.seb_diagnostic_enable (bool, optional, default false)
   * - erf.radiation.seb_prognostic_enable (bool, optional, default false)
   * - erf.radiation.seb_surface_heat_capacity (real, optional, default 2.0e4, validated > 0)
   * - erf.radiation.seb_restore_timescale_s (real, optional, default 86400.0, validated > 0)
   * - erf.radiation.seb_moisture_layer_depth_m (real, optional, default 0.1, validated > 0)
   * - erf.radiation.seb_moisture_restore_timescale_s (real, optional, default 86400.0, validated > 0)
   * - erf.radiation.seb_prognostic_t_min_k (real, optional, default 200.0)
   * - erf.radiation.seb_prognostic_t_max_k (real, optional, default 340.0)
   * - erf.radiation.seb_prognostic_q_min (real, optional, default 0.0, clamped to [0, 1])
   * - erf.radiation.seb_prognostic_q_max (real, optional, default 1.0, clamped to [0, 1])
   *
   * Sub-option queries are only performed when the two-stream model is selected.
   */
  /**
   * @brief Bring an out-of-range input back into its valid range and say so.
   * Every erf.radiation.* value with a physical range passes through one of
   * these, so a typo in a deck is reported instead of silently corrected.
   */
  static void clamp_low (const char* key, amrex::Real& value, amrex::Real lo)
  {
    if (value < lo) {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: erf.radiation." << key << " = " << value
                       << " is below its minimum " << lo << "; using " << lo << "\n";
      }
      value = lo;
    }
  }
  static void clamp_high (const char* key, amrex::Real& value, amrex::Real hi)
  {
    if (value > hi) {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: erf.radiation." << key << " = " << value
                       << " is above its maximum " << hi << "; using " << hi << "\n";
      }
      value = hi;
    }
  }
  static void reset_nonpositive (const char* key, amrex::Real& value, amrex::Real fallback)
  {
    if (!std::isfinite(value) || value <= 0.0) {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: erf.radiation." << key << " = " << value
                       << " must be positive; using " << fallback << "\n";
      }
      value = fallback;
    }
  }

  void init_params(bool two_stream_selected, int max_level, std::string pp_prefix)
  {
    amrex::ParmParse pp(pp_prefix);
    enabled = two_stream_selected;

    // A refined hierarchy is supported: every level runs its own column sweep
    // and writes its own heating rates. The one requirement is that each
    // level's grids span that level's domain in z, because the sweep imposes
    // top-of-atmosphere and surface boundary conditions at the ends of a box.
    // That is checked where the grids exist, not here: ERF_RefineBox.cpp for
    // an explicit erf.boxN refinement box, and TwoStreamRadiation::define_level
    // for every path including tagged refinement and regrid.

    // Only query sub-options if the two-stream model is selected
    if (enabled) {
      pp.query("radiation.sw_enabled", sw_enabled);
      pp.query("radiation.lw_enabled", lw_enabled);
      pp.query("radiation.tau_per_layer", tau_per_layer);
      pp.query("radiation.tau_lw_per_layer", tau_lw_per_layer);

      // The sun, the site and the surface temperature come from the keys the
      // RRTMGP interface reads, so the two models share one set of inputs.
      // ParmParse ignores a key nobody reads, so a deck still carrying one of
      // the earlier two-stream-only keys would run with silently different
      // settings: stop and name the replacement instead.
      static const char* retired_keys[][2] = {
          {"radiation.solar_zenith",
           "erf.fixed_solar_zenith_angle, which is the cosine of the angle (60 degrees is 0.5)"},
          {"radiation.solar_zenith_deg",
           "erf.fixed_solar_zenith_angle, which is the cosine of the angle (60 degrees is 0.5)"},
          {"radiation.S0", "erf.fixed_total_solar_irradiance"},
          {"radiation.earth_sun_distance_enable",
           "nothing: leave erf.fixed_total_solar_irradiance unset and the irradiance follows the date"},
          {"radiation.surface_temp_k", "erf.rad_t_sfc"},
          {"radiation.solar_geometry_dynamic_enable",
           "nothing: leave erf.fixed_solar_zenith_angle unset and the sun follows start_datetime"},
          {"radiation.latitude_deg", "erf.rad_cons_lat"},
          {"radiation.longitude_deg", "erf.rad_cons_lon"},
          {"radiation.day_of_year", "start_datetime"},
          {"radiation.time_zone_offset_hours",
           "nothing: start_datetime is UTC and the longitude sets the local solar time"},
      };
      for (const auto& entry : retired_keys) {
          if (pp.contains(entry[0])) {
              amrex::Abort(std::string("erf.") + entry[0] + " is no longer read by the "
                           "two-stream model; use " + entry[1] + ".");
          }
      }
      pp.query("fixed_solar_zenith_angle", fixed_solar_zenith_angle);
      if (fixed_solar_zenith_angle > 1.0) {
          amrex::Abort("erf.fixed_solar_zenith_angle = " + std::to_string(fixed_solar_zenith_angle) +
                       " is the cosine of the solar zenith angle (as for RRTMGP) and cannot "
                       "exceed 1; 60 degrees is 0.5.");
      }
      pp.query("fixed_total_solar_irradiance", fixed_total_solar_irradiance);
      if (!pp.query("rad_t_sfc", rad_t_sfc)) {
          amrex::Abort("erf.radiation_model = TwoStream needs erf.rad_t_sfc, the surface "
                       "temperature [K] used where no land-surface model or surface layer "
                       "supplies one (the same input RRTMGP requires).");
      }
      if (!std::isfinite(rad_t_sfc) || rad_t_sfc <= 0.0) {
          amrex::Abort("erf.rad_t_sfc = " + std::to_string(rad_t_sfc) + " must be a positive temperature [K].");
      }
      pp.query("rad_cons_lat", rad_cons_lat);
      pp.query("rad_cons_lon", rad_cons_lon);
      if (!std::isfinite(rad_cons_lat) || rad_cons_lat < -90.0 || rad_cons_lat > 90.0) {
          amrex::Abort("erf.rad_cons_lat = " + std::to_string(rad_cons_lat) + " must lie in [-90, 90] degrees.");
      }
      if (!std::isfinite(rad_cons_lon) || rad_cons_lon < -180.0 || rad_cons_lon > 180.0) {
          amrex::Abort("erf.rad_cons_lon = " + std::to_string(rad_cons_lon) + " must lie in [-180, 180] degrees.");
      }
      pp.query("rad_orbital_year", rad_orbital_year);
      pp.query("rad_orbital_eccentricity", rad_orbital_eccentricity);
      pp.query("rad_orbital_obliquity", rad_orbital_obliquity);
      pp.query("rad_orbital_mvelp", rad_orbital_mvelp);
      if (pp.contains("rad_freq_in_steps") && amrex::ParallelDescriptor::IOProcessor()) {
          amrex::Print() << "WARNING: erf.rad_freq_in_steps is not read by the two-stream model, "
                            "which runs every step.\n";
      }
      std::string tau_model_str = "per_layer";
      pp.query("radiation.tau_model", tau_model_str);
      if (tau_model_str == "per_layer") {
        tau_model = TauModel::PerLayer;
      } else if (tau_model_str == "mass") {
        tau_model = TauModel::Mass;
      } else {
        // The two models give very different column optical depths, so say so
        // rather than let a misspelling quietly pick the default.
        tau_model = TauModel::PerLayer;
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: unknown erf.radiation.tau_model '"
                         << tau_model_str << "'; valid values are per_layer "
                            "and mass. Using per_layer.\n";
        }
      }
      pp.query("radiation.sw_kabs_dry", sw_kabs_dry);
      pp.query("radiation.sw_kscat_dry", sw_kscat_dry);
      pp.query("radiation.sw_kabs_vapor", sw_kabs_vapor);
      pp.query("radiation.sw_kext_cloud", sw_kext_cloud);
      pp.query("radiation.sw_cloud_omega", sw_cloud_omega);
      pp.query("radiation.sw_cloud_g", sw_cloud_g);
      clamp_low("sw_kabs_dry", sw_kabs_dry, 0.0);
      clamp_low("sw_kscat_dry", sw_kscat_dry, 0.0);
      clamp_low("sw_kabs_vapor", sw_kabs_vapor, 0.0);
      clamp_low("sw_kext_cloud", sw_kext_cloud, 0.0);
      clamp_low("sw_cloud_omega", sw_cloud_omega, 0.0);
      clamp_high("sw_cloud_omega", sw_cloud_omega, 1.0);
      clamp_low("sw_cloud_g", sw_cloud_g, -1.0);
      clamp_high("sw_cloud_g", sw_cloud_g, 1.0);
      pp.query("radiation.lw_mass_absorption_enable", lw_mass_absorption_enable);
      pp.query("radiation.lw_kabs_dry", lw_kabs_dry);
      pp.query("radiation.lw_kabs_vapor", lw_kabs_vapor);
      pp.query("radiation.lw_kabs_cloud", lw_kabs_cloud);
      clamp_low("lw_kabs_dry", lw_kabs_dry, 0.0);
      clamp_low("lw_kabs_vapor", lw_kabs_vapor, 0.0);
      clamp_low("lw_kabs_cloud", lw_kabs_cloud, 0.0);
      pp.query("radiation.v", verbosity);
      pp.query("radiation.diag_file", diag_file);

      // Cloud optical properties (all optional; the defaults leave the
      // column clear-sky)
      std::string tau_profile_type_str = "constant";
      pp.query("radiation.tau_profile_type", tau_profile_type_str);
      if (tau_profile_type_str == "constant") {
        tau_profile_type = TauProfileType::Constant;
      } else if (tau_profile_type_str == "cloud_layer") {
        tau_profile_type = TauProfileType::CloudLayer;
      } else {
        // Falling back to constant drops a requested cloud layer, so say so.
        tau_profile_type = TauProfileType::Constant;
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: unknown erf.radiation.tau_profile_type '"
                         << tau_profile_type_str << "'; valid values are "
                            "constant and cloud_layer. Using constant, so no "
                            "cloud layer is applied.\n";
        }
      }
      pp.query("radiation.cloud_base_height_m", cloud_base_height_m);
      pp.query("radiation.cloud_top_height_m", cloud_top_height_m);
      pp.query("radiation.cloud_tau_per_layer", cloud_tau_per_layer);
      pp.query("radiation.cloud_fraction", cloud_fraction);

      // Parameter validation (R5/R5'): clip cloud_fraction to [0, 1]
      clamp_low("cloud_fraction", cloud_fraction, 0.0);
      clamp_high("cloud_fraction", cloud_fraction, 1.0);
      clamp_low("cloud_tau_per_layer", cloud_tau_per_layer, 0.0);
      if (cloud_top_height_m < cloud_base_height_m) {
        cloud_top_height_m = cloud_base_height_m;
      }

      // Scattering parameters (all optional; the defaults leave the
      // shortwave direct-beam only)
      pp.query("radiation.single_scattering_albedo", single_scattering_albedo);
      pp.query("radiation.asymmetry_factor", asymmetry_factor);
      pp.query("radiation.cloud_single_scattering_albedo", cloud_single_scattering_albedo);
      pp.query("radiation.cloud_asymmetry_factor", cloud_asymmetry_factor);

      // Parameter validation (R5/R5'): clip albedos to [0, 1] and
      // asymmetry factors to [-1, 1]
      clamp_low("single_scattering_albedo", single_scattering_albedo, 0.0);
      clamp_high("single_scattering_albedo", single_scattering_albedo, 1.0);
      clamp_low("asymmetry_factor", asymmetry_factor, -1.0);
      clamp_high("asymmetry_factor", asymmetry_factor, 1.0);
      clamp_low("cloud_single_scattering_albedo", cloud_single_scattering_albedo, 0.0);
      clamp_high("cloud_single_scattering_albedo", cloud_single_scattering_albedo, 1.0);
      clamp_low("cloud_asymmetry_factor", cloud_asymmetry_factor, -1.0);
      clamp_high("cloud_asymmetry_factor", cloud_asymmetry_factor, 1.0);

      // Runtime diagnostics controls (all optional; every channel is on
      // by default)
      pp.query("radiation.diag_enable", diag_enable);
      pp.query("radiation.diag_stdout_enable", diag_stdout_enable);
      pp.query("radiation.diag_tagged_enable", diag_tagged_enable);
      pp.query("radiation.diag_regtest_line_enable", diag_regtest_line_enable);
      pp.query("radiation.diag_csv_enable", diag_csv_enable);
      pp.query("radiation.diag_callsite_mode", diag_callsite_mode);
      pp.query("radiation.diag_dedup_tol", diag_dedup_tol);

      // Validate the call-site mode
      if (diag_callsite_mode != "both" &&
          diag_callsite_mode != "pre_only" &&
          diag_callsite_mode != "post_only") {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: unknown erf.radiation.diag_callsite_mode '"
                         << diag_callsite_mode << "'; valid values are both, "
                            "pre_only and post_only. Using both.\n";
        }
        diag_callsite_mode = "both";
      }
      // Ensure dedup tolerance is non-negative
      clamp_low("diag_dedup_tol", diag_dedup_tol, 0.0);

      // Surface heterogeneity fallback parameters (all optional; only used
      // when hetero
      // LSM/radiation interface fields are absent or invalid)
      pp.query("radiation.surface_albedo_sw", surface_albedo_sw);
      pp.query("radiation.surface_emissivity_lw", surface_emissivity_lw);
      pp.query("radiation.surface_albedo_sw_diffuse", surface_albedo_sw_diffuse);
      clamp_high("surface_albedo_sw_diffuse", surface_albedo_sw_diffuse, 1.0);

      // Clamp albedo/emissivity to [0, 1]
      clamp_low("surface_albedo_sw", surface_albedo_sw, 0.0);
      clamp_high("surface_albedo_sw", surface_albedo_sw, 1.0);
      clamp_low("surface_emissivity_lw", surface_emissivity_lw, 0.0);
      clamp_high("surface_emissivity_lw", surface_emissivity_lw, 1.0);

      // Dynamic optical depth diagnosis from moisture/clouds
      // (all optional; disabled by default)
      pp.query("radiation.tau_sw_dynamic_enable", tau_sw_dynamic_enable);
      pp.query("radiation.tau_lw_dynamic_enable", tau_lw_dynamic_enable);
      pp.query("radiation.tau_sw_coeff_qv", tau_sw_coeff_qv);
      pp.query("radiation.tau_sw_coeff_qc", tau_sw_coeff_qc);
      pp.query("radiation.tau_lw_coeff_qv", tau_lw_coeff_qv);
      pp.query("radiation.tau_lw_coeff_qc", tau_lw_coeff_qc);

      // Clamp the dynamic tau coefficients
      // to nonnegative (negative would give unphysical behavior)
      clamp_low("tau_sw_coeff_qv", tau_sw_coeff_qv, 0.0);
      clamp_low("tau_sw_coeff_qc", tau_sw_coeff_qc, 0.0);
      clamp_low("tau_lw_coeff_qv", tau_lw_coeff_qv, 0.0);
      clamp_low("tau_lw_coeff_qc", tau_lw_coeff_qc, 0.0);

      // Prognostic cloud fraction diagnosis from RH/qc
      // (all optional; disabled by default)
      pp.query("radiation.cloud_fraction_prog_enable", cloud_fraction_prog_enable);
      pp.query("radiation.cloud_fraction_rh_min", cloud_fraction_rh_min);
      pp.query("radiation.cloud_fraction_rh_max", cloud_fraction_rh_max);
      pp.query("radiation.cloud_fraction_qc_scale", cloud_fraction_qc_scale);

      // Clamp the cloud fraction parameters
      // to physically reasonable ranges
      clamp_low("cloud_fraction_rh_min", cloud_fraction_rh_min, 0.0);
      clamp_high("cloud_fraction_rh_min", cloud_fraction_rh_min, 1.0);
      clamp_low("cloud_fraction_rh_max", cloud_fraction_rh_max, 0.0);
      clamp_high("cloud_fraction_rh_max", cloud_fraction_rh_max, 1.0);
      if (cloud_fraction_rh_max < cloud_fraction_rh_min) {
        cloud_fraction_rh_max = cloud_fraction_rh_min;
      }
      clamp_low("cloud_fraction_qc_scale", cloud_fraction_qc_scale, 0.0);

      // Prescribed bulk aerosol/turbidity optical depth
      // (all optional; disabled by default)
      pp.query("radiation.aerosol_enable", aerosol_enable);

      // Parse aerosol profile type as string
      std::string aerosol_profile_type_str = "constant";
      pp.query("radiation.aerosol_profile_type", aerosol_profile_type_str);
      if (aerosol_profile_type_str == "constant") {
        aerosol_profile_type = AerosolProfileType::Constant;
      } else if (aerosol_profile_type_str == "exponential") {
        aerosol_profile_type = AerosolProfileType::Exponential;
      } else if (aerosol_profile_type_str == "table") {
        // The table lookup is not implemented and contributes no optical
        // depth. Keep the selection so the run is not silently given a
        // different profile, and say what it means.
        aerosol_profile_type = AerosolProfileType::Table;
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: erf.radiation.aerosol_profile_type = table is "
                            "not implemented and contributes no aerosol optical "
                            "depth. Use constant or exponential.\n";
        }
      } else {
        // Falling through to constant would hide a misspelling.
        aerosol_profile_type = AerosolProfileType::Constant;
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "WARNING: unknown erf.radiation.aerosol_profile_type '"
                         << aerosol_profile_type_str << "'; valid values are "
                            "constant and exponential. Using constant.\n";
        }
      }

      pp.query("radiation.aerosol_tau_per_layer", aerosol_tau_per_layer);
      pp.query("radiation.aerosol_scale_height_m", aerosol_scale_height_m);
      pp.query("radiation.aerosol_tau_surface", aerosol_tau_surface);

      // Clamp the aerosol parameters
      // to physically reasonable ranges
      clamp_low("aerosol_tau_per_layer", aerosol_tau_per_layer, 0.0);
      clamp_low("aerosol_tau_surface", aerosol_tau_surface, 0.0);
      reset_nonpositive("aerosol_scale_height_m", aerosol_scale_height_m, 2000.0);

      // Simplified SEB infrastructure fallback parameters
      pp.query("radiation.seb_enable", seb_enable);
      pp.query("radiation.seb_use_radiation_fluxes", seb_use_radiation_fluxes);
      pp.query("radiation.seb_sw_flux_default", seb_sw_flux_default);
      pp.query("radiation.seb_lw_flux_default", seb_lw_flux_default);
      pp.query("radiation.seb_hfx_default", seb_hfx_default);
      pp.query("radiation.seb_lh_default", seb_lh_default);
      pp.query("radiation.seb_grdflx_default", seb_grdflx_default);
      pp.query("radiation.seb_q_sfc_default", seb_q_sfc_default);
      pp.query("radiation.seb_t_deep_default", seb_t_deep_default);
      pp.query("radiation.seb_q_deep_default", seb_q_deep_default);

      // SEB diagnostic mode
      pp.query("radiation.seb_diagnostic_enable", seb_diagnostic_enable);

      // Auto-enable seb_enable if seb_diagnostic_enable is requested but seb_enable is false
      if (seb_diagnostic_enable && !seb_enable) {
          seb_enable = true;  // Internal auto-enable for safety/simplicity
      }

      // SEB prognostic mode parameters
      pp.query("radiation.seb_prognostic_enable", seb_prognostic_enable);
      pp.query("radiation.seb_surface_heat_capacity", seb_surface_heat_capacity);
      pp.query("radiation.seb_restore_timescale_s", seb_restore_timescale_s);
      pp.query("radiation.seb_moisture_layer_depth_m", seb_moisture_layer_depth_m);
      pp.query("radiation.seb_moisture_restore_timescale_s", seb_moisture_restore_timescale_s);
      pp.query("radiation.seb_prognostic_t_min_k", seb_prognostic_t_min_k);
      pp.query("radiation.seb_prognostic_t_max_k", seb_prognostic_t_max_k);
      pp.query("radiation.seb_prognostic_q_min", seb_prognostic_q_min);
      pp.query("radiation.seb_prognostic_q_max", seb_prognostic_q_max);

      // The prognostic update needs the SEB fields and the residual, so
      // requesting it turns on the infrastructure it depends on. Report that,
      // since it changes settings the user did not set.
      if (seb_prognostic_enable && (!seb_enable || !seb_diagnostic_enable)) {
          if (amrex::ParallelDescriptor::IOProcessor()) {
              amrex::Print() << "NOTE: erf.radiation.seb_prognostic_enable = true "
                                "also enables seb_enable and "
                                "seb_diagnostic_enable, which it depends on.\n";
          }
          seb_enable = true;
          seb_diagnostic_enable = true;
      }

      if (!std::isfinite(seb_sw_flux_default)) seb_sw_flux_default = 0.0;
      if (!std::isfinite(seb_lw_flux_default)) seb_lw_flux_default = 0.0;
      if (!std::isfinite(seb_hfx_default)) seb_hfx_default = 0.0;
      if (!std::isfinite(seb_lh_default)) seb_lh_default = 0.0;
      if (!std::isfinite(seb_grdflx_default)) seb_grdflx_default = 0.0;
      if (!std::isfinite(seb_q_sfc_default)) seb_q_sfc_default = 0.0;
      reset_nonpositive("seb_t_deep_default", seb_t_deep_default, 300.0);
      if (!std::isfinite(seb_q_deep_default)) seb_q_deep_default = 0.0;

      seb_sw_flux_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_sw_flux_default));
      seb_lw_flux_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_lw_flux_default));
      seb_hfx_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_hfx_default));
      seb_lh_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_lh_default));
      seb_grdflx_default = amrex::max(amrex::Real(-2000.0), amrex::min(amrex::Real(2000.0), seb_grdflx_default));
      seb_q_sfc_default = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), seb_q_sfc_default));
      seb_q_deep_default = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), seb_q_deep_default));


      // Validate the prognostic parameters
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
      seb_prognostic_q_min = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), seb_prognostic_q_min));
      seb_prognostic_q_max = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), seb_prognostic_q_max));
      if (seb_prognostic_q_min > seb_prognostic_q_max) {
          amrex::Real tmp = seb_prognostic_q_min;
          seb_prognostic_q_min = seb_prognostic_q_max;
          seb_prognostic_q_max = tmp;
      }

      // The prognostic surface energy balance is a single-level feature: it evolves the
      // surface temperature that the longwave boundary condition reads, and it runs on
      // level 0 only, so a refined run would give level 0 and its fine levels two
      // different surface boundary conditions for one surface.
      //
      // This belongs here rather than where the fine level is built, because
      // amr.max_level is known now and the combination is wrong regardless of whether a
      // fine level ever materialises: with a refinement criterion that tags nothing the
      // run would otherwise proceed, and with erf.regrid_int > 0 it would abort at the
      // first regrid, hours in.
      if (seb_prognostic_enable && max_level > 0) {
          amrex::Abort("erf.radiation.seb_prognostic_enable = true is supported on a single level "
                       "only: the force-restore surface temperature it evolves is the longwave "
                       "boundary condition, and fine levels have no copy of it. Set "
                       "amr.max_level = 0, or turn off erf.radiation.seb_prognostic_enable to run "
                       "two-stream radiation on a refined hierarchy.");
      }
  }

}

  /**
   * @brief Check if two-stream radiation is active.
   *
   * @return true when erf.radiation_model = TwoStream.
   */
  [[nodiscard]] bool uses_two_stream_radiation() const noexcept
  {
    return enabled;
  }
};

#endif // ERF_RAD_STRUCT_H_
