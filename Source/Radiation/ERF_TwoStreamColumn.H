#ifndef ERF_TWO_STREAM_COLUMN_H_
#define ERF_TWO_STREAM_COLUMN_H_

#include <cmath>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_REAL.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_RadStruct.H>
#include <ERF_TwoStreamSW.H>
#include <ERF_TwoStreamLW.H>
#include <ERF_PrognosticCloudFraction.H>
#include <ERF_AerosolOpticalDepth.H>
#include <ERF_SolarGeometry.H>

/**
 * @file ERF_TwoStreamColumn.H
 * @brief Per-column two-stream radiation kernels.
 *
 * Everything needed to evaluate one (i,j) column of the two-stream model:
 * temperature from the equation of state, per-column surface properties,
 * per-layer optical depth assembly, and vertical_two_stream_sweep(), which
 * performs the SW and LW sweeps and writes per-level heating rates.
 *
 * Vertical orientation follows ERF: k = kmin is the surface layer and
 * k = kmax the top layer. SW sweeps downward from kmax to kmin; LW sweeps
 * downward (TOA -> surface) and then upward (surface -> TOA) on layer
 * interfaces. Layer temperature is obtained from rho*theta through the Exner
 * function.
 *
 * The functions are header-only and GPU-safe so that the column physics can
 * be exercised directly by unit tests (Tests/Unit/Radiation) as well as by
 * the level driver in ERF_AdvanceTwoStreamRadiation.cpp.
 */

/**
 * @brief GPU-safe helper to compute absolute temperature from (rho, rho*theta).
 *
 * ERF stores dry density and dry potential temperature. The pressure follows
 * from the equation of state,
 *   p = p_0 * (R_d * rho * theta_m / p_0)^gamma,   theta_m = theta * (1 + R_v/R_d * qv),
 * and the absolute temperature is recovered through the Exner function,
 *   T = theta * (p / p_0)^(R_d / c_p),
 * which getTgivenRandRTh() evaluates as p / (R_d * rho * (1 + R_v/R_d * qv)).
 *
 * @param[in] rho_theta RhoTheta component [K·kg/m^3]
 * @param[in] rho Density [kg/m^3]
 * @param[in] qv Water-vapor mixing ratio [kg/kg] (0 for dry air)
 * @return Temperature [K], clamped to [100, 400]; 288.15 for unphysical input.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_temperature_from_rhotheta(amrex::Real rho_theta, amrex::Real rho,
                                          amrex::Real qv = 0.0)
{
    if (!(rho > 0.0) || !std::isfinite(rho) ||
        !(rho_theta > 0.0) || !std::isfinite(rho_theta)) {
        return 288.15;  // Defensive: fallback to standard T
    }
    if (!(qv >= 0.0) || !std::isfinite(qv)) {
        qv = 0.0;
    }

    amrex::Real T = getTgivenRandRTh(rho, rho_theta, qv);
    if (!std::isfinite(T) || T <= 0.0) {
        return 288.15;  // Defensive: fallback
    }

    // Defensive clipping to a sensible terrestrial range
    T = std::max(T, 100.0);
    T = std::min(T, 400.0);
    return T;
}

/**
 * @brief Clamp a real value to [min, max], with fallback for invalid (NaN/Inf).
 *
 * If the value is not finite (NaN, Inf), returns fallback. Otherwise returns
 * value clamped to [min, max].
 *
 * @param[in] value Input value to clamp
 * @param[in] min Lower bound
 * @param[in] max Upper bound
 * @param[in] fallback Value to use if input is not finite
 * @return Clamped finite value in [min, max]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real clamp_finite(amrex::Real value, amrex::Real min, amrex::Real max,
                         amrex::Real fallback)
{
    if (!std::isfinite(value)) {
        value = fallback;
    }
    if (value < min) value = min;
    if (value > max) value = max;
    return value;
}

/**
 * @brief Check if a real value is finite and positive.
 *
 * @param[in] value Value to check
 * @return true if finite and > 0, false otherwise
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_finite_positive(amrex::Real value)
{
    return std::isfinite(value) && value > 0.0;
}

/**
 * @brief Resolve per-column shortwave surface albedo from hetero field or fallback.
 *
 * Precedence:
 * 1. If hetero_alb_sw array available and value at (i,j) is finite ∈ [0,1], use it
 * 2. Otherwise, use rad_choice.surface_albedo_sw (already clamped by init_params)
 * 3. Hard default: 0.3
 *
 * @param[in] i,j Column index
 * @param[in] hetero_alb_sw Heterogeneous SW albedo field (may be nullptr)
 * @param[in] rad_choice Radiation parameters with fallback surface_albedo_sw
 * @param[in] has_hetero_alb true if hetero_alb_sw is available
 * @return SW albedo in [0, 1]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real resolve_surface_albedo_sw(
    int i, int j,
    const amrex::Array4<const amrex::Real>* hetero_alb_sw,
    const RadChoice& rad_choice,
    bool has_hetero_alb)
{
    amrex::Real alb = rad_choice.surface_albedo_sw;  // Default fallback (already clamped)
    
    if (has_hetero_alb && hetero_alb_sw != nullptr && hetero_alb_sw->contains(i, j, 0)) {
        amrex::Real hetero_val = (*hetero_alb_sw)(i, j, 0, 0);
        if (std::isfinite(hetero_val)) {
            alb = clamp_finite(hetero_val, 0.0, 1.0, rad_choice.surface_albedo_sw);
        }
    }
    
    return alb;
}

/**
 * @brief Resolve per-column longwave surface emissivity from hetero field or fallback.
 *
 * Precedence:
 * 1. If hetero_emiss_lw array available and value at (i,j) is finite ∈ [0,1], use it
 * 2. Otherwise, use rad_choice.surface_emissivity_lw (already clamped by init_params)
 * 3. Hard default: 0.99
 *
 * @param[in] i,j Column index
 * @param[in] hetero_emiss_lw Heterogeneous LW emissivity field (may be nullptr)
 * @param[in] rad_choice Radiation parameters with fallback surface_emissivity_lw
 * @param[in] has_hetero_emiss true if hetero_emiss_lw is available
 * @return LW emissivity in [0, 1]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real resolve_surface_emissivity_lw(
    int i, int j,
    const amrex::Array4<const amrex::Real>* hetero_emiss_lw,
    const RadChoice& rad_choice,
    bool has_hetero_emiss)
{
    amrex::Real emiss = rad_choice.surface_emissivity_lw;  // Default fallback (already clamped)
    
    if (has_hetero_emiss && hetero_emiss_lw != nullptr && hetero_emiss_lw->contains(i, j, 0)) {
        amrex::Real hetero_val = (*hetero_emiss_lw)(i, j, 0, 0);  // Assume single component
        if (std::isfinite(hetero_val)) {
            emiss = clamp_finite(hetero_val, 0.0, 1.0, rad_choice.surface_emissivity_lw);
        }
    }
    
    return emiss;
}

/**
 * @brief Resolve per-column surface temperature from hetero field or fallback.
 *
 * Precedence:
 * 1. If t_sfc array available and value at (i,j) is finite & positive, use it
 * 2. Otherwise, use rad_choice.surface_temp_k (already validated by init_params)
 * 3. Hard default: 300.0 K
 *
 * @param[in] i,j Column index
 * @param[in] t_sfc Heterogeneous surface temperature field (may be nullptr)
 * @param[in] rad_choice Radiation parameters with fallback surface_temp_k
 * @param[in] has_t_sfc true if t_sfc is available
 * @return Surface temperature [K], strictly positive
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real resolve_surface_temp_k(
    int i, int j,
    const amrex::Array4<const amrex::Real>* t_sfc,
    const RadChoice& rad_choice,
    bool has_t_sfc)
{
    amrex::Real t_surf = rad_choice.surface_temp_k;  // Default fallback (already validated)
    
    if (has_t_sfc && t_sfc != nullptr && t_sfc->contains(i, j, 0)) {
        amrex::Real hetero_val = (*t_sfc)(i, j, 0, 0);  // Assume single component
        if (is_finite_positive(hetero_val)) {
            t_surf = hetero_val;
        }
    }
    
    return t_surf;
}

/**
 * @brief GPU-safe helper to determine whether a layer falls within the
 * cloud band [cloud_base_height_m, cloud_top_height_m].
 *
 * @param[in] z_center Height of the layer center above the surface [m].
 * @param[in] rad_choice Radiation parameters.
 * @return true if this layer is inside the configured cloud band.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_cloud_level(amrex::Real z_center, const RadChoice& rad_choice)
{
    return (z_center >= rad_choice.cloud_base_height_m &&
            z_center <= rad_choice.cloud_top_height_m);
}

/**
 * @brief Water-vapor mixing ratio qv = RhoQv / Rho at (i,j,k).
 *
 * Returns 0 when the state carries no moisture components, or when the
 * stored values are non-finite or negative.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_qv_from_state(int i, int j, int k,
                              const amrex::Array4<const amrex::Real>& state_arr)
{
    if (state_arr.nComp() <= RhoQ1_comp) return 0.0;
    amrex::Real rho    = state_arr(i, j, k, Rho_comp);
    amrex::Real rho_qv = state_arr(i, j, k, RhoQ1_comp);
    if (!(rho > 0.0) || !std::isfinite(rho) || !std::isfinite(rho_qv)) return 0.0;
    amrex::Real qv = rho_qv / rho;
    return (std::isfinite(qv) && qv > 0.0) ? qv : 0.0;
}

/**
 * @brief Cloud-water mixing ratio qc = RhoQc / Rho at (i,j,k).
 *
 * Returns 0 when the state carries no cloud-water component, or when the
 * stored values are non-finite or negative.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_qc_from_state(int i, int j, int k,
                              const amrex::Array4<const amrex::Real>& state_arr)
{
    if (state_arr.nComp() <= RhoQ2_comp) return 0.0;
    amrex::Real rho    = state_arr(i, j, k, Rho_comp);
    amrex::Real rho_qc = state_arr(i, j, k, RhoQ2_comp);
    if (!(rho > 0.0) || !std::isfinite(rho) || !std::isfinite(rho_qc)) return 0.0;
    amrex::Real qc = rho_qc / rho;
    return (std::isfinite(qc) && qc > 0.0) ? qc : 0.0;
}

/**
 * @brief Moisture-dependent (dynamic) per-layer optical depth.
 *
 *   tau = tau_base + coeff_qv * qv + coeff_qc * qc
 *
 * clamped to [0, 100]. With both coefficients zero this returns tau_base
 * unchanged, so the static configuration is reproduced exactly. Used for
 * both the SW and the LW band with the respective coefficient pair.
 *
 * @param[in] tau_base Optical depth before the moisture contribution [unitless].
 * @param[in] qv Water-vapor mixing ratio [kg/kg].
 * @param[in] qc Cloud-water mixing ratio [kg/kg].
 * @param[in] coeff_qv Optical depth per unit qv [unitless per kg/kg].
 * @param[in] coeff_qc Optical depth per unit qc [unitless per kg/kg].
 * @return Optical depth for this layer [unitless].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_dynamic(amrex::Real tau_base, amrex::Real qv, amrex::Real qc,
                                 amrex::Real coeff_qv, amrex::Real coeff_qc)
{
    if (!(qv > 0.0) || !std::isfinite(qv)) qv = 0.0;
    if (!(qc > 0.0) || !std::isfinite(qc)) qc = 0.0;
    amrex::Real tau = tau_base + coeff_qv * qv + coeff_qc * qc;
    if (!std::isfinite(tau)) return tau_base;
    if (tau < 0.0) tau = 0.0;
    if (tau > 100.0) tau = 100.0;
    return tau;
}

/**
 * @brief GPU-safe helper to compute the per-layer optical depth at level k,
 * given the base (clear-sky) optical depth and cloud-layer
 * parameters.
 *
 * When rad_choice.tau_profile_type == Constant, returns tau_base unchanged
 * (byte-identical). When == CloudLayer, adds
 * rad_choice.cloud_tau_per_layer whenever the level height falls within
 * [cloud_base_height_m, cloud_top_height_m].
 *
 * @param[in] z_center Height of the layer center above the surface [m].
 * @param[in] tau_base Clear-sky optical depth per layer.
 * @param[in] rad_choice Radiation parameters.
 * @param[in] apply_cloud If false, always returns tau_base (used for the
 * clear-sky column computation even when cloud_fraction > 0).
 * @return Optical depth for this layer [unitless].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real tau_layer_value(
    amrex::Real z_center, amrex::Real tau_base,
    const RadChoice& rad_choice, bool apply_cloud)
{
    if (!apply_cloud || rad_choice.tau_profile_type != TauProfileType::CloudLayer) {
        return tau_base;
    }
    if (is_cloud_level(z_center, rad_choice)) {
        return tau_base + rad_choice.cloud_tau_per_layer;
    }
    return tau_base;
}

/**
 * @brief GPU-safe helper to select the single-scattering albedo
 * and asymmetry factor to use for level k's diffuse SW calculation.
 *
 * When this column evaluation applies the cloud-layer enhancement
 * (apply_cloud == true, tau_profile_type == CloudLayer, and level k falls
 * within the cloud band), the cloud scattering properties
 * (cloud_single_scattering_albedo, cloud_asymmetry_factor) are used.
 * Otherwise, the clear-sky scattering properties (single_scattering_albedo,
 * asymmetry_factor) are used. Both default to 0.0, so by default this
 * function always yields omega == 0.0, and compute_sw_diffuse_flux()
 * returns exactly 0.0 (direct-beam-only behavior preserved).
 *
 * @param[in] z_center Height of the layer center above the surface [m].
 * @param[in] rad_choice Radiation parameters.
 * @param[in] apply_cloud Same flag passed to tau_layer_value(); true for the
 * cloudy-column evaluation, false for the clear-sky column evaluation.
 * @param[out] omega Selected single-scattering albedo for this level.
 * @param[out] g Selected asymmetry factor for this level.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void select_scattering_props(
    amrex::Real z_center, const RadChoice& rad_choice, bool apply_cloud,
    amrex::Real& omega, amrex::Real& g)
{
    bool use_cloud_props = apply_cloud &&
        rad_choice.tau_profile_type == TauProfileType::CloudLayer &&
        is_cloud_level(z_center, rad_choice);

    if (use_cloud_props) {
        omega = rad_choice.cloud_single_scattering_albedo;
        g = rad_choice.cloud_asymmetry_factor;
    } else {
        omega = rad_choice.single_scattering_albedo;
        g = rad_choice.asymmetry_factor;
    }
}

/**
 * @brief GPU-safe helper to diagnose prognostic cloud fraction
 * from per-level relative humidity and cloud liquid water.
 *
 * Computes cloud fraction from RH and qc using:
 *   cf_rh(k) = linear ramp from 0 at rh_min to 1 at rh_max
 *   cf_qc(k) = qc_scale * qc(k)
 *   cf(k) = min(1, cf_rh + cf_qc)  [saturated blend]
 *
 * Temperature comes from the equation of state (Exner function) and
 * pressure from getPgivenRTh(), both evaluated with the local qv.
 * Returns 0 when cloud_fraction_prog_enable is false.
 *
 * @param[in] i, j, k Grid indices
 * @param[in] state_arr State array proxy (contains Rho, RhoTheta, qv, qc)
 * @param[in] rad_choice Radiation parameters (prognostic cloud fraction settings)
 * @return Diagnosed cloud fraction [0, 1] if enabled; 0 if disabled
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_cloud_fraction_prognostic(
    int i, int j, int k,
    const amrex::Array4<const amrex::Real>& state_arr,
    const RadChoice& rad_choice)
{
    if (!rad_choice.cloud_fraction_prog_enable) {
        return 0.0;
    }

    const amrex::Real qv = get_qv_from_state(i, j, k, state_arr);
    const amrex::Real qc = get_qc_from_state(i, j, k, state_arr);

    amrex::Real rho       = state_arr(i, j, k, Rho_comp);
    amrex::Real rho_theta = state_arr(i, j, k, RhoTheta_comp);
    if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;
    if (rho_theta <= 0.0 || !std::isfinite(rho_theta)) rho_theta = 288.15;

    // Temperature (via the Exner function) and pressure from the equation of state
    amrex::Real T = get_temperature_from_rhotheta(rho_theta, rho, qv);
    amrex::Real P = getPgivenRTh(rho_theta, qv);
    if (P <= 0.0 || !std::isfinite(P)) P = p_0;

    const amrex::Real rh = compute_relative_humidity(qv, T, P);

    return diagnose_cloud_fraction_from_rh_qc(
        rh, qc,
        rad_choice.cloud_fraction_rh_min,
        rad_choice.cloud_fraction_rh_max,
        rad_choice.cloud_fraction_qc_scale);
}

/**
 * @brief Maximum number of vertical levels supported by the
 * fixed-capacity per-column buffers in vertical_two_stream_sweep().
 *
 * The LW heating-rate computation requires the full upward and downward LW
 * flux profile (not just the surface value), so this kernel stores both
 * interface profiles, the per-layer thickness and the per-layer LW optical
 * depth in device-local arrays sized to this capacity. This is a pragmatic,
 * GPU-safe alternative to dynamic per-thread allocation.
 *
 * If bx.length(2) exceeds this value, an AMREX_ALWAYS_ASSERT will fire at
 * runtime. Increase this constant if a taller domain is required; note
 * memory-per-thread scales linearly with it.
 */
constexpr int MAX_RAD_LEVELS = 512;

/**
 * @brief GPU-safe helper to assemble the total per-layer optical depth at
 * level k for either the SW or the LW band.
 *
 * The contributions are applied in the same order for both bands:
 * 1. clear-sky base value, plus the cloud-layer enhancement when the column
 *    is evaluated as "cloudy" and the level lies inside the cloud band;
 * 2. optional moisture-dependent (dynamic) term from qv and qc;
 * 3. optional prognostic cloud fraction, which replaces the cloud-band
 *    enhancement by cf(k) * cloud_tau_per_layer;
 * 4. optional prescribed aerosol term (constant, exponential or table).
 *
 * @param[in] i, j, k Grid indices (k increases upward; kmin is the surface layer).
 * @param[in] dz_layer Thickness of this layer [m].
 * @param[in] z_center Height of the layer center above the surface [m]
 *            (cloud-band detection and the aerosol profile).
 * @param[in] state_arr State array (read-only).
 * @param[in] tau_base Clear-sky optical depth per layer for this band.
 * @param[in] is_sw true for the shortwave band, false for longwave.
 * @param[in] cloudy true for the cloudy-column evaluation.
 * @param[in] rad_choice Radiation parameters.
 * @return Optical depth of layer k [unitless].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_layer_tau(
    int i, int j, int k,
    amrex::Real dz_layer, amrex::Real z_center,
    const amrex::Array4<const amrex::Real>& state_arr,
    amrex::Real tau_base, bool is_sw, bool cloudy,
    const RadChoice& rad_choice)
{
    // Base (clear-sky) value, plus cloud-band enhancement for cloudy columns.
    amrex::Real tau = tau_layer_value(z_center, tau_base, rad_choice, cloudy);

    // Dynamic (moisture-dependent) contribution.
    const bool dynamic_enabled = is_sw ? rad_choice.tau_sw_dynamic_enable
                                       : rad_choice.tau_lw_dynamic_enable;
    if (dynamic_enabled) {
        const amrex::Real qv = get_qv_from_state(i, j, k, state_arr);
        const amrex::Real qc = get_qc_from_state(i, j, k, state_arr);
        tau = is_sw ? diagnose_tau_dynamic(tau, qv, qc, rad_choice.tau_sw_coeff_qv, rad_choice.tau_sw_coeff_qc)
                    : diagnose_tau_dynamic(tau, qv, qc, rad_choice.tau_lw_coeff_qv, rad_choice.tau_lw_coeff_qc);
    }

    // Prognostic cloud fraction: scale the cloud-band enhancement by cf(k).
    if (rad_choice.cloud_fraction_prog_enable && cloudy &&
        rad_choice.tau_profile_type == TauProfileType::CloudLayer &&
        is_cloud_level(z_center, rad_choice)) {
        amrex::Real cf_prog = diagnose_cloud_fraction_prognostic(i, j, k, state_arr, rad_choice);
        amrex::Real tau_clear_k = tau_layer_value(z_center, tau_base, rad_choice, /*apply_cloud=*/false);
        tau = tau_clear_k + cf_prog * rad_choice.cloud_tau_per_layer;
    }

    // Prescribed bulk aerosol contribution (added on top of everything above).
    if (rad_choice.aerosol_enable) {
        amrex::Real tau_aerosol = 0.0;
        if (rad_choice.aerosol_profile_type == AerosolProfileType::Constant) {
            tau_aerosol = diagnose_tau_aerosol_constant(rad_choice.aerosol_tau_per_layer);
        } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Exponential) {
            tau_aerosol = diagnose_tau_aerosol_exponential(z_center, dz_layer,
                                                           rad_choice.aerosol_tau_surface,
                                                           rad_choice.aerosol_scale_height_m);
        } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Table) {
            tau_aerosol = diagnose_tau_aerosol_table(k);
        }
        tau += tau_aerosol;
    }

    return tau;
}

/**
 * @brief GPU-safe per-column vertical integration kernel for two-stream
 * radiation, computing either the clear-sky or cloudy-column fluxes and
 * per-level heating rates, depending on the `cloudy` flag.
 *
 * **Vertical orientation.** ERF's vertical index increases upward: k = kmin
 * is the layer adjacent to the surface and k = kmax is the layer adjacent to
 * the top of the domain (TOA for this model). Layer k spans the interfaces
 * m = k - kmin (bottom) and m = k - kmin + 1 (top), so interface m = 0 is the
 * surface and m = nlev is the TOA.
 *
 * Per (i,j) column:
 * 1. SW: sweep downward from k = kmax to k = kmin, accumulating optical depth
 *    and the direct + diffuse downwelling flux; the per-level SW heating rate
 *    is written to qheating_arr(i,j,k,0).
 * 2. LW: sweep downward from the TOA (F_down = 0) to the surface, then upward
 *    from the surface (F_up = eps * sigma * T_s^4 + (1 - eps) * F_down) to
 *    the TOA, storing both interface profiles. The per-level LW heating rate from the net-flux
 *    divergence is written to qheating_arr(i,j,k,1).
 * 3. Scalar diagnostics (max heating rate, surface fluxes) are returned for
 *    the reduction in the caller.
 *
 * Layer temperature is obtained from (rho, rho*theta, qv) through the
 * equation of state, i.e. including the Exner function, so LW emission uses
 * absolute temperature rather than potential temperature.
 *
 * Integrates per-column heterogeneous surface properties (albedo,
 * emissivity, surface temperature) from optional fields with robust fallback
 * to RadChoice scalar parameters.
 *
 * @param[in] i, j Column indices
 * @param[in] bx Computational box (cell-centered, full vertical extent)
 * @param[in] geom Geometry for this AMR level
 * @param[in] state_arr Array proxy to state data (read-only)
 * @param[in] rad_choice Radiation parameters
 * @param[in] cloudy If true and tau_profile_type == CloudLayer, apply the
 * cloud-layer optical depth enhancement (and cloud scattering properties).
 * @param[out] qheating_arr Component 0 receives the SW heating rate [K/s]
 * and component 1 the LW heating rate [K/s] at every level k in [kmin, kmax].
 * @param[out] max_heating_rate Maximum |Q_sw|+|Q_lw| observed in this column
 * @param[out] sw_surface_flux Absorbed downwelling SW at the surface (direct + diffuse)
 * @param[out] lw_net_surface Net LW (up - down) at the surface
 * @param[in] z_phys_cc Optional cell-centered physical height (nonuniform dz)
 * @param[in] time_utc_seconds UTC seconds within the day (dynamic solar geometry)
 * @param[in] has_hetero_alb_sw true if hetero_alb_sw is available
 * @param[in] hetero_alb_sw Optional per-column SW surface albedo field
 * @param[in] has_hetero_emiss_lw true if hetero_emiss_lw is available
 * @param[in] hetero_emiss_lw Optional per-column LW surface emissivity field
 * @param[in] has_t_sfc true if t_sfc is available
 * @param[in] t_sfc Optional per-column surface temperature field [K]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void vertical_two_stream_sweep(
    int i, int j,
    const amrex::Box& bx,
    const amrex::Geometry& geom,
    const amrex::Array4<const amrex::Real>& state_arr,
    const RadChoice& rad_choice,
    bool cloudy,
    const amrex::Array4<amrex::Real>& qheating_arr,
    amrex::Real& max_heating_rate,
    amrex::Real& sw_surface_flux,
    amrex::Real& lw_net_surface,
    const amrex::Array4<const amrex::Real>& z_phys_cc,
    amrex::Real time_utc_seconds = 0.0,
    bool has_hetero_alb_sw = false,
    const amrex::Array4<const amrex::Real>* hetero_alb_sw = nullptr,
    bool has_hetero_emiss_lw = false,
    const amrex::Array4<const amrex::Real>* hetero_emiss_lw = nullptr,
    bool has_t_sfc = false,
    const amrex::Array4<const amrex::Real>* t_sfc = nullptr)
{
    // Grid bounds: kmin is the surface layer, kmax the top layer.
    const int kmin = bx.smallEnd(2);
    const int kmax = bx.bigEnd(2);
    const int nlev = kmax - kmin + 1;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nlev <= MAX_RAD_LEVELS,
        "vertical_two_stream_sweep: domain vertical extent exceeds "
        "MAX_RAD_LEVELS; increase the constant in ERF_AdvanceTwoStreamRadiation.cpp");

    // Physical constants
    const amrex::Real sigma  = stefan_boltzmann;  // Stefan-Boltzmann [W/(m^2·K^4)]
    const amrex::Real cp_air = Cp_d;              // Dry-air specific heat at constant pressure [J/(kg·K)]

    // ------------------------------------------------------------------
    // Layer thickness dz_level[m] and layer-center height z_level[m] above
    // the surface. Both come from the cell-centered physical heights when
    // z_phys_cc is available (stretched or terrain-following grids), and
    // from the uniform spacing otherwise. z_level locates the cloud band
    // and evaluates the exponential aerosol profile.
    // ------------------------------------------------------------------
    const amrex::Real dz_uniform = geom.CellSize(2);
    amrex::Real dz_level[MAX_RAD_LEVELS];
    amrex::Real z_level[MAX_RAD_LEVELS];

    for (int m = 0; m < nlev; ++m) {
        dz_level[m] = dz_uniform;
    }
    if (z_phys_cc) {
        for (int m = 0; m + 1 < nlev; ++m) {
            amrex::Real dz_computed = z_phys_cc(i, j, kmin + m + 1) - z_phys_cc(i, j, kmin + m);
            if (dz_computed > 0.0 && std::isfinite(dz_computed)) {
                dz_level[m] = dz_computed;
            }
        }
        // Top layer: no level above; reuse the thickness of the layer below.
        if (nlev > 1) dz_level[nlev - 1] = dz_level[nlev - 2];
    }
    if (z_phys_cc) {
        // Surface height = first cell center minus half of the first thickness.
        const amrex::Real z_surface = z_phys_cc(i, j, kmin) - 0.5 * dz_level[0];
        for (int m = 0; m < nlev; ++m) {
            z_level[m] = z_phys_cc(i, j, kmin + m) - z_surface;
        }
    } else {
        for (int m = 0; m < nlev; ++m) {
            z_level[m] = (static_cast<amrex::Real>(m) + 0.5) * dz_uniform;
        }
    }

    // Solar zenith angle (optionally time-varying)
    amrex::Real cos_zenith;
    if (rad_choice.solar_geometry_dynamic_enable) {
        cos_zenith = compute_cos_zenith_angle(
            time_utc_seconds,
            rad_choice.latitude_deg,
            rad_choice.longitude_deg,
            rad_choice.day_of_year,
            rad_choice.time_zone_offset_hours);
    } else {
        amrex::Real zenith_rad = rad_choice.solar_zenith_deg * M_PI / 180.0;
        cos_zenith = std::cos(zenith_rad);
    }

    const amrex::Real S0          = rad_choice.S0;
    const amrex::Real tau_sw_base = rad_choice.tau_per_layer;
    const amrex::Real tau_lw_base = rad_choice.tau_lw_per_layer;

    amrex::Real local_max_heating = 0.0;

    // Zero-initialize this column's heating rate output (covers the
    // sw_enabled=false and lw_enabled=false cases).
    for (int k = kmin; k <= kmax; ++k) {
        qheating_arr(i, j, k, 0) = 0.0;
        qheating_arr(i, j, k, 1) = 0.0;
    }

    // ========================================================================
    // SHORTWAVE: downward sweep from the top layer (k = kmax) to the surface
    // layer (k = kmin). "prev" refers to the interface above the current
    // layer, "curr" to the interface below it.
    // ========================================================================
    amrex::Real tau_sw_cum    = 0.0;  // Cumulative SW optical depth from TOA
    amrex::Real F_sw_dir_prev = 0.0;  // Direct SW flux at the top of the current layer
    amrex::Real F_sw_diff_prev = 0.0; // Diffuse SW flux at the top of the current layer

    if (rad_choice.sw_enabled) {
        if (cos_zenith > 0.0) {
            F_sw_dir_prev = S0 * cos_zenith;  // TOA incident (tau = 0)
        }
        F_sw_diff_prev = 0.0;  // No diffuse flux incident from above TOA

        for (int k = kmax; k >= kmin; --k) {
            const int m = k - kmin;

            amrex::Real rho = state_arr(i, j, k, Rho_comp);
            if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

            amrex::Real tau_sw = diagnose_layer_tau(i, j, k, dz_level[m], z_level[m],
                                                    state_arr, tau_sw_base, /*is_sw=*/true, cloudy,
                                                    rad_choice);
            tau_sw_cum += tau_sw;

            // Direct beam at the bottom of this layer (Beer-Lambert)
            amrex::Real F_sw_dir_curr = compute_sw_direct_flux(tau_sw_cum, S0, cos_zenith);

            // Diffuse flux: transmitted from above plus newly scattered in this layer
            amrex::Real omega = 0.0;
            amrex::Real g = 0.0;
            select_scattering_props(z_level[m], rad_choice, cloudy, omega, g);

            amrex::Real F_sw_diff_layer =
                compute_sw_diffuse_flux(tau_sw, F_sw_dir_prev, cos_zenith, omega, g);
            amrex::Real Tdir_layer = (cos_zenith > 0.0)
                ? std::exp(-tau_sw / cos_zenith)
                : 0.0;
            amrex::Real F_sw_diff_curr = F_sw_diff_prev * Tdir_layer + F_sw_diff_layer;

            // Heating from the divergence of the total downwelling flux
            amrex::Real F_sw_total_top = F_sw_dir_prev + F_sw_diff_prev;
            amrex::Real F_sw_total_bot = F_sw_dir_curr + F_sw_diff_curr;
            amrex::Real Q_sw = compute_sw_heating_rate(F_sw_total_top, F_sw_total_bot,
                                                       dz_level[m], rho, cp_air);
            qheating_arr(i, j, k, 0) = Q_sw;
            local_max_heating = std::max(local_max_heating, std::abs(Q_sw));

            F_sw_dir_prev  = F_sw_dir_curr;
            F_sw_diff_prev = F_sw_diff_curr;
        }
    }
    // After the loop, F_sw_*_prev hold the downwelling fluxes at the surface.

    // ========================================================================
    // LONGWAVE: interface flux profiles, index m = 0 at the surface and
    // m = nlev at the TOA.
    // ========================================================================
    amrex::Real F_lw_up_iface[MAX_RAD_LEVELS + 1];
    amrex::Real F_lw_down_iface[MAX_RAD_LEVELS + 1];
    amrex::Real tau_lw_level[MAX_RAD_LEVELS];

    if (rad_choice.lw_enabled) {
        // Downward sweep: TOA -> surface. Also caches the per-layer LW optical
        // depth for the upward sweep.
        F_lw_down_iface[nlev] = 0.0;  // No incoming LW from space
        for (int k = kmax; k >= kmin; --k) {
            const int m = k - kmin;

            amrex::Real rho       = state_arr(i, j, k, Rho_comp);
            amrex::Real rho_theta = state_arr(i, j, k, RhoTheta_comp);
            if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;
            if (rho_theta <= 0.0 || !std::isfinite(rho_theta)) rho_theta = 288.15;
            amrex::Real qv = get_qv_from_state(i, j, k, state_arr);
            amrex::Real T_layer = get_temperature_from_rhotheta(rho_theta, rho, qv);

            tau_lw_level[m] = diagnose_layer_tau(i, j, k, dz_level[m], z_level[m],
                                                 state_arr, tau_lw_base, /*is_sw=*/false, cloudy,
                                                 rad_choice);

            F_lw_down_iface[m] = compute_lw_flux_down(F_lw_down_iface[m + 1], T_layer, sigma, tau_lw_level[m]);
        }

        // Upward sweep: surface -> TOA. The surface emits eps * sigma * T_s^4
        // and reflects the fraction (1 - eps) of the downwelling flux that
        // reaches it (gray surface, Kirchhoff's law).
        {
            amrex::Real t_surface = resolve_surface_temp_k(i, j, t_sfc, rad_choice, has_t_sfc);
            amrex::Real emiss_lw  = resolve_surface_emissivity_lw(i, j, hetero_emiss_lw, rad_choice, has_hetero_emiss_lw);
            F_lw_up_iface[0] = emiss_lw * compute_thermal_intensity(t_surface, sigma)
                             + (1.0 - emiss_lw) * F_lw_down_iface[0];
        }
        for (int k = kmin; k <= kmax; ++k) {
            const int m = k - kmin;

            amrex::Real rho       = state_arr(i, j, k, Rho_comp);
            amrex::Real rho_theta = state_arr(i, j, k, RhoTheta_comp);
            if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;
            if (rho_theta <= 0.0 || !std::isfinite(rho_theta)) rho_theta = 288.15;
            amrex::Real qv = get_qv_from_state(i, j, k, state_arr);
            amrex::Real T_layer = get_temperature_from_rhotheta(rho_theta, rho, qv);

            F_lw_up_iface[m + 1] = compute_lw_flux_up(F_lw_up_iface[m], T_layer, sigma, tau_lw_level[m]);
        }

        if (rad_choice.isothermal_test) {
            // Isothermal test: every level radiates equally, so the net flux
            // (and hence the heating) vanishes everywhere.
            for (int m = 0; m <= nlev; ++m) {
                F_lw_down_iface[m] = F_lw_up_iface[m];
            }
        }

        // Per-level LW heating rate from the net-flux divergence across
        // each layer (bottom interface m, top interface m+1).
        for (int k = kmin; k <= kmax; ++k) {
            const int m = k - kmin;
            amrex::Real rho = state_arr(i, j, k, Rho_comp);
            if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

            amrex::Real F_net_top = F_lw_up_iface[m + 1] - F_lw_down_iface[m + 1];
            amrex::Real F_net_bot = F_lw_up_iface[m]     - F_lw_down_iface[m];

            amrex::Real Q_lw = compute_lw_heating_rate(F_net_top, F_net_bot, dz_level[m], rho, cp_air);
            qheating_arr(i, j, k, 1) = Q_lw;

            amrex::Real Q_sw_here = qheating_arr(i, j, k, 0);
            local_max_heating = std::max(local_max_heating, std::abs(Q_sw_here) + std::abs(Q_lw));
        }
    }

    // ========================================================================
    // SURFACE AND DIAGNOSTICS
    // ========================================================================
    if (rad_choice.sw_enabled) {
        amrex::Real alb_sw = resolve_surface_albedo_sw(i, j, hetero_alb_sw, rad_choice, has_hetero_alb_sw);

        if (rad_choice.isothermal_test) {
            // Isothermal test: direct beam at the surface, apply albedo
            sw_surface_flux = S0 * std::max(0.0, cos_zenith) * std::exp(-tau_sw_cum) * (1.0 - alb_sw);
        } else {
            // Normal mode: direct plus diffuse, apply albedo
            sw_surface_flux = (F_sw_dir_prev + F_sw_diff_prev) * (1.0 - alb_sw);
        }
    } else {
        sw_surface_flux = 0.0;
    }

    if (rad_choice.lw_enabled) {
        amrex::Real F_lw_up_sfc   = F_lw_up_iface[0];
        amrex::Real F_lw_down_sfc = F_lw_down_iface[0];
        if (rad_choice.isothermal_test) {
            amrex::Real rho_surface       = state_arr(i, j, kmin, Rho_comp);
            amrex::Real rho_theta_surface = state_arr(i, j, kmin, RhoTheta_comp);
            if (rho_surface <= 0.0 || !std::isfinite(rho_surface)) rho_surface = 1.0;
            if (rho_theta_surface <= 0.0 || !std::isfinite(rho_theta_surface)) rho_theta_surface = 288.15;
            amrex::Real qv_surface = get_qv_from_state(i, j, kmin, state_arr);
            amrex::Real T_iso = get_temperature_from_rhotheta(rho_theta_surface, rho_surface, qv_surface);
            amrex::Real I_thermal = compute_thermal_intensity(T_iso, sigma);
            F_lw_up_sfc   = I_thermal;
            F_lw_down_sfc = I_thermal;
            // Isothermal override: net LW flux and heating are exactly zero
            for (int k = kmin; k <= kmax; ++k) {
                qheating_arr(i, j, k, 1) = 0.0;
            }
        }
        lw_net_surface = F_lw_up_sfc - F_lw_down_sfc;
    } else {
        lw_net_surface = 0.0;
    }

    max_heating_rate = local_max_heating;
}

#endif // ERF_TWO_STREAM_COLUMN_H_
