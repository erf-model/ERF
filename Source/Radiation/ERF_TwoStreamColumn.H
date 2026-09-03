#ifndef ERF_TWO_STREAM_COLUMN_H_
#define ERF_TWO_STREAM_COLUMN_H_

#include <cmath>
#include <type_traits>

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
 * @brief Trivially copyable subset of RadChoice used by the column kernels.
 *
 * RadChoice carries std::string members (diagnostics file name, call-site
 * mode) and therefore cannot be captured by value in a GPU device lambda.
 * The level driver builds one TwoStreamParams per call with
 * make_two_stream_params() and passes that to vertical_two_stream_sweep()
 * and its helpers. Field names match RadChoice.
 */
struct TwoStreamParams
{
    bool sw_enabled = true;
    bool lw_enabled = true;
    amrex::Real tau_per_layer = 0.05;
    amrex::Real tau_lw_per_layer = 1.0;
    amrex::Real solar_zenith_deg = 45.0;
    amrex::Real S0 = 1361.0;

    TauProfileType tau_profile_type = TauProfileType::Constant;
    amrex::Real cloud_base_height_m = 500.0;
    amrex::Real cloud_top_height_m = 1000.0;
    amrex::Real cloud_tau_per_layer = 0.5;

    amrex::Real single_scattering_albedo = 0.0;
    amrex::Real asymmetry_factor = 0.0;
    amrex::Real cloud_single_scattering_albedo = 0.0;
    amrex::Real cloud_asymmetry_factor = 0.0;

    amrex::Real surface_albedo_sw = 0.3;
    amrex::Real surface_albedo_sw_diffuse = -1.0;
    amrex::Real surface_emissivity_lw = 0.99;
    amrex::Real surface_temp_k = 300.0;

    bool tau_sw_dynamic_enable = false;
    bool tau_lw_dynamic_enable = false;
    amrex::Real tau_sw_coeff_qv = 0.0;
    amrex::Real tau_sw_coeff_qc = 0.0;
    amrex::Real tau_lw_coeff_qv = 0.0;
    amrex::Real tau_lw_coeff_qc = 0.0;

    bool cloud_fraction_prog_enable = false;
    amrex::Real cloud_fraction_rh_min = 0.0;
    amrex::Real cloud_fraction_rh_max = 1.0;
    amrex::Real cloud_fraction_qc_scale = 1.0e-3;

    bool aerosol_enable = false;
    AerosolProfileType aerosol_profile_type = AerosolProfileType::Constant;
    amrex::Real aerosol_tau_per_layer = 0.0;
    amrex::Real aerosol_scale_height_m = 2000.0;
    amrex::Real aerosol_tau_surface = 0.0;

    TauModel tau_model = TauModel::PerLayer;
    amrex::Real sw_kabs_dry = 4.0e-6;
    amrex::Real sw_kscat_dry = 3.0e-6;
    amrex::Real sw_kabs_vapor = 4.0e-3;
    amrex::Real sw_kext_cloud = 150.0;
    amrex::Real sw_cloud_omega = 0.9999;
    amrex::Real sw_cloud_g = 0.85;

    bool lw_mass_absorption_enable = false;
    amrex::Real lw_kabs_dry = 1.0e-4;
    amrex::Real lw_kabs_vapor = 0.1;
    amrex::Real lw_kabs_cloud = 158.0;

    bool earth_sun_distance_enable = false;
    bool solar_geometry_dynamic_enable = false;
    amrex::Real latitude_deg = 0.0;
    amrex::Real longitude_deg = 0.0;
    amrex::Real day_of_year = 172.0;
    amrex::Real time_zone_offset_hours = 0.0;
};
static_assert(std::is_trivially_copyable<TwoStreamParams>::value,
              "TwoStreamParams must stay trivially copyable so it can be captured by device lambdas");

/**
 * @brief Copy the column-kernel parameters out of a RadChoice.
 */
inline TwoStreamParams make_two_stream_params(const RadChoice& rc)
{
    TwoStreamParams p;
    p.sw_enabled = rc.sw_enabled;
    p.lw_enabled = rc.lw_enabled;
    p.tau_per_layer = rc.tau_per_layer;
    p.tau_lw_per_layer = rc.tau_lw_per_layer;
    p.solar_zenith_deg = rc.solar_zenith_deg;
    p.S0 = rc.S0;
    p.tau_profile_type = rc.tau_profile_type;
    p.cloud_base_height_m = rc.cloud_base_height_m;
    p.cloud_top_height_m = rc.cloud_top_height_m;
    p.cloud_tau_per_layer = rc.cloud_tau_per_layer;
    p.single_scattering_albedo = rc.single_scattering_albedo;
    p.asymmetry_factor = rc.asymmetry_factor;
    p.cloud_single_scattering_albedo = rc.cloud_single_scattering_albedo;
    p.cloud_asymmetry_factor = rc.cloud_asymmetry_factor;
    p.surface_albedo_sw = rc.surface_albedo_sw;
    p.surface_albedo_sw_diffuse = rc.surface_albedo_sw_diffuse;
    p.surface_emissivity_lw = rc.surface_emissivity_lw;
    p.surface_temp_k = rc.surface_temp_k;
    p.tau_sw_dynamic_enable = rc.tau_sw_dynamic_enable;
    p.tau_lw_dynamic_enable = rc.tau_lw_dynamic_enable;
    p.tau_sw_coeff_qv = rc.tau_sw_coeff_qv;
    p.tau_sw_coeff_qc = rc.tau_sw_coeff_qc;
    p.tau_lw_coeff_qv = rc.tau_lw_coeff_qv;
    p.tau_lw_coeff_qc = rc.tau_lw_coeff_qc;
    p.cloud_fraction_prog_enable = rc.cloud_fraction_prog_enable;
    p.cloud_fraction_rh_min = rc.cloud_fraction_rh_min;
    p.cloud_fraction_rh_max = rc.cloud_fraction_rh_max;
    p.cloud_fraction_qc_scale = rc.cloud_fraction_qc_scale;
    p.aerosol_enable = rc.aerosol_enable;
    p.aerosol_profile_type = rc.aerosol_profile_type;
    p.aerosol_tau_per_layer = rc.aerosol_tau_per_layer;
    p.aerosol_scale_height_m = rc.aerosol_scale_height_m;
    p.aerosol_tau_surface = rc.aerosol_tau_surface;
    p.tau_model = rc.tau_model;
    p.sw_kabs_dry = rc.sw_kabs_dry;
    p.sw_kscat_dry = rc.sw_kscat_dry;
    p.sw_kabs_vapor = rc.sw_kabs_vapor;
    p.sw_kext_cloud = rc.sw_kext_cloud;
    p.sw_cloud_omega = rc.sw_cloud_omega;
    p.sw_cloud_g = rc.sw_cloud_g;
    p.lw_mass_absorption_enable = rc.lw_mass_absorption_enable;
    p.lw_kabs_dry = rc.lw_kabs_dry;
    p.lw_kabs_vapor = rc.lw_kabs_vapor;
    p.lw_kabs_cloud = rc.lw_kabs_cloud;
    p.earth_sun_distance_enable = rc.earth_sun_distance_enable;
    p.solar_geometry_dynamic_enable = rc.solar_geometry_dynamic_enable;
    p.latitude_deg = rc.latitude_deg;
    p.longitude_deg = rc.longitude_deg;
    p.day_of_year = rc.day_of_year;
    p.time_zone_offset_hours = rc.time_zone_offset_hours;
    return p;
}

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
    const TwoStreamParams& rad_choice,
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
    const TwoStreamParams& rad_choice,
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
    const TwoStreamParams& rad_choice,
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
bool is_cloud_level(amrex::Real z_center, const TwoStreamParams& rad_choice)
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
 * @brief Inverse Exner function 1/pi at (i,j,k) from the state, with
 * pi = (p / p_0)^(R_d/c_p) evaluated by getExnergivenRTh(). Used to convert
 * a temperature tendency dT/dt into the potential-temperature tendency
 * dtheta/dt = (dT/dt) / pi that the RhoTheta source term expects.
 * Returns 1 for unphysical input.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_inverse_exner(int i, int j, int k,
                              const amrex::Array4<const amrex::Real>& state_arr)
{
    amrex::Real rho_theta = state_arr(i, j, k, RhoTheta_comp);
    if (!(rho_theta > 0.0) || !std::isfinite(rho_theta)) return 1.0;
    const amrex::Real qv = get_qv_from_state(i, j, k, state_arr);
    const amrex::Real pi = getExnergivenRTh(rho_theta, RdoCp, qv);
    return (std::isfinite(pi) && pi > 0.0) ? 1.0 / pi : 1.0;
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
    const TwoStreamParams& rad_choice, bool apply_cloud)
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
 * function always yields omega == 0.0 and the layers neither scatter the
 * direct beam nor reflect diffuse light; the diffuse field then consists of
 * the surface-reflected beam only.
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
    amrex::Real z_center, const TwoStreamParams& rad_choice, bool apply_cloud,
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
    const TwoStreamParams& rad_choice)
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
 * The SW adding method and the LW sweeps need full interface profiles
 * (direct beam, adding-method albedo and source, LW up/down fluxes) and the
 * per-layer thickness, height and optical depths, so this kernel stores them
 * in device-local arrays sized to this capacity (about nine arrays of
 * MAX_RAD_LEVELS reals per thread). This is a pragmatic, GPU-safe
 * alternative to dynamic per-thread allocation.
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
 * @param[in] tau_base Clear-sky optical depth per layer for this band
 *            (for LW, replaced by the mass-path value when
 *            lw_mass_absorption_enable or tau_model = mass is set).
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
    const TwoStreamParams& rad_choice)
{
    // Longwave option: gray optical depth from the layer mass path,
    //   tau = rho dz (k_dry + k_vapor qv + k_cloud qc),
    // which makes the column optical depth independent of the vertical
    // resolution and gives water vapor and cloud water a real greenhouse
    // effect. Replaces the fixed tau_lw_per_layer as the clear-sky base.
    amrex::Real base = tau_base;
    if (!is_sw && (rad_choice.lw_mass_absorption_enable || rad_choice.tau_model == TauModel::Mass)) {
        amrex::Real rho = state_arr(i, j, k, Rho_comp);
        if (!(rho > 0.0) || !std::isfinite(rho)) rho = 0.0;
        const amrex::Real qv = get_qv_from_state(i, j, k, state_arr);
        const amrex::Real qc = get_qc_from_state(i, j, k, state_arr);
        base = rho * dz_layer * (rad_choice.lw_kabs_dry
                                 + rad_choice.lw_kabs_vapor * qv
                                 + rad_choice.lw_kabs_cloud * qc);
        if (!std::isfinite(base) || base < 0.0) base = 0.0;
    }

    // Base (clear-sky) value, plus cloud-band enhancement for cloudy columns.
    amrex::Real tau = tau_layer_value(z_center, base, rad_choice, cloudy);

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
        amrex::Real tau_clear_k = tau_layer_value(z_center, base, rad_choice, /*apply_cloud=*/false);
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
 * @brief Optical depth, single-scattering albedo and asymmetry factor of
 * layer k for one band, for either optical-depth model.
 *
 * Per-layer model (default): the optical depth comes from
 * diagnose_layer_tau() and the scattering properties from the clear-sky or
 * cloud-band inputs through select_scattering_props(), exactly as before.
 *
 * Mass model (tau_model = mass), shortwave: each constituent contributes an
 * extinction optical depth from the layer mass path,
 *   dry absorption   rho dz sw_kabs_dry            (omega = 0)
 *   Rayleigh         rho dz sw_kscat_dry           (omega = 1, g = 0)
 *   water vapor      rho dz sw_kabs_vapor qv       (omega = 0)
 *   cloud water      rho dz sw_kext_cloud qc       (sw_cloud_omega, sw_cloud_g)
 * plus the prescribed cloud-band enhancement (with the cloud-band scattering
 * inputs) and the absorbing moisture-coefficient and aerosol additions. The
 * layer properties are the extinction-weighted mixtures
 *   tau   = sum tau_i,
 *   omega = sum omega_i tau_i / tau,
 *   g     = sum g_i omega_i tau_i / sum omega_i tau_i.
 * Longwave uses diagnose_layer_tau(), whose clear-sky base is the mass path
 * rho dz (lw_kabs_dry + lw_kabs_vapor qv + lw_kabs_cloud qc) in this model.
 *
 * @param[in] i, j, k Grid indices (kmin is the surface layer).
 * @param[in] dz_layer Thickness of this layer [m].
 * @param[in] z_center Height of the layer center above the surface [m].
 * @param[in] state_arr State array (read-only).
 * @param[in] tau_base Per-layer-model clear-sky optical depth for this band.
 * @param[in] is_sw true for the shortwave band, false for longwave.
 * @param[in] cloudy true for the cloudy-column evaluation.
 * @param[in] rad_choice Column-kernel parameters.
 * @param[out] tau Optical depth of the layer [unitless].
 * @param[out] omega Single-scattering albedo of the layer (0 for LW).
 * @param[out] g Asymmetry factor of the layer (0 for LW).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void diagnose_layer_optics(
    int i, int j, int k,
    amrex::Real dz_layer, amrex::Real z_center,
    const amrex::Array4<const amrex::Real>& state_arr,
    amrex::Real tau_base, bool is_sw, bool cloudy,
    const TwoStreamParams& rad_choice,
    amrex::Real& tau, amrex::Real& omega, amrex::Real& g)
{
    omega = 0.0;
    g = 0.0;
    if (!is_sw || rad_choice.tau_model != TauModel::Mass) {
        tau = diagnose_layer_tau(i, j, k, dz_layer, z_center, state_arr, tau_base, is_sw, cloudy, rad_choice);
        if (is_sw) {
            select_scattering_props(z_center, rad_choice, cloudy, omega, g);
        }
        return;
    }

    // Shortwave mass model: constituent extinction optical depths.
    amrex::Real rho = state_arr(i, j, k, Rho_comp);
    if (!(rho > 0.0) || !std::isfinite(rho)) rho = 0.0;
    const amrex::Real qv = get_qv_from_state(i, j, k, state_arr);
    const amrex::Real qc = get_qc_from_state(i, j, k, state_arr);
    const amrex::Real path = rho * dz_layer;

    const amrex::Real tau_dry_abs = path * rad_choice.sw_kabs_dry;
    const amrex::Real tau_rayleigh = path * rad_choice.sw_kscat_dry;
    const amrex::Real tau_vapor = path * rad_choice.sw_kabs_vapor * qv;
    const amrex::Real tau_cloud = path * rad_choice.sw_kext_cloud * qc;

    amrex::Real ext = tau_dry_abs + tau_rayleigh + tau_vapor + tau_cloud;
    amrex::Real sca = tau_rayleigh + rad_choice.sw_cloud_omega * tau_cloud;
    amrex::Real gsca = rad_choice.sw_cloud_g * rad_choice.sw_cloud_omega * tau_cloud;

    // Prescribed cloud band (optionally scaled by the prognostic cloud
    // fraction) with the cloud-band scattering inputs.
    amrex::Real band = tau_layer_value(z_center, 0.0, rad_choice, cloudy);
    if (rad_choice.cloud_fraction_prog_enable && cloudy &&
        rad_choice.tau_profile_type == TauProfileType::CloudLayer &&
        is_cloud_level(z_center, rad_choice)) {
        band = diagnose_cloud_fraction_prognostic(i, j, k, state_arr, rad_choice) * rad_choice.cloud_tau_per_layer;
    }
    if (band > 0.0) {
        ext += band;
        sca += rad_choice.cloud_single_scattering_albedo * band;
        gsca += rad_choice.cloud_asymmetry_factor * rad_choice.cloud_single_scattering_albedo * band;
    }

    // Absorbing additions (moisture coefficients, aerosol): everything that
    // diagnose_layer_tau() adds on top of a zero base, minus the band.
    const amrex::Real additions = diagnose_layer_tau(i, j, k, dz_layer, z_center, state_arr,
                                                     0.0, /*is_sw=*/true, cloudy, rad_choice) - band;
    if (additions > 0.0) ext += additions;

    if (!std::isfinite(ext) || ext <= 0.0) {
        tau = 0.0;
        return;
    }
    tau = ext;
    omega = sca / ext;
    if (omega > 1.0) omega = 1.0;
    if (omega < 0.0) omega = 0.0;
    g = (sca > 0.0) ? gsca / sca : 0.0;
    if (g > 1.0) g = 1.0;
    if (g < -1.0) g = -1.0;
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
 * 1. SW: Beer-Lambert direct beam on the interfaces, then the two-stream
 *    diffuse field (upward and downward streams) by the adding method with
 *    the surface albedo as the lower boundary condition; the per-level SW
 *    heating rate from the net-flux divergence is written to
 *    qheating_arr(i,j,k,0).
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
 * @param[in] rad_choice Column-kernel parameters (TwoStreamParams, built from RadChoice)
 * @param[in] cloudy If true and tau_profile_type == CloudLayer, apply the
 * cloud-layer optical depth enhancement (and cloud scattering properties).
 * @param[out] qheating_arr Component 0 receives the SW and component 1 the
 * LW radiative tendency of potential temperature, dtheta/dt = (dT/dt) / pi
 * [K/s], at every level k in [kmin, kmax]. This is the convention the
 * RhoTheta source term (ERF_MakeSources.cpp) and the RRTMGP path use.
 * @param[out] max_heating_rate Maximum |Q_sw|+|Q_lw| of the stored theta tendencies in this column
 * @param[out] sw_surface_flux SW absorbed by the surface, (1 - alb_dir) F_dir + (1 - alb_dif) F_dif [W/m^2]
 * @param[out] sw_up_toa Upwelling (reflected) SW at the top of the atmosphere [W/m^2]
 * @param[out] lw_net_surface Net LW (up - down) at the surface [W/m^2]
 * @param[out] lw_up_toa Upwelling LW at the top of the atmosphere (outgoing longwave) [W/m^2]
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
    const TwoStreamParams& rad_choice,
    bool cloudy,
    const amrex::Array4<amrex::Real>& qheating_arr,
    amrex::Real& max_heating_rate,
    amrex::Real& sw_surface_flux,
    amrex::Real& sw_up_toa,
    amrex::Real& lw_net_surface,
    amrex::Real& lw_up_toa,
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

    // Solar constant, optionally scaled by the Earth-Sun distance factor.
    const amrex::Real S0 = rad_choice.S0 *
        (rad_choice.earth_sun_distance_enable
             ? compute_earth_sun_distance_factor(rad_choice.day_of_year) : 1.0);
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
    // SHORTWAVE. Interfaces m = 0 (surface) .. nlev (TOA); layer m lies
    // between interfaces m and m+1.
    //  1. Direct beam F_dir(m) by Beer-Lambert from the TOA downward.
    //  2. Adding method, upward pass: A[m] is the albedo of everything below
    //     interface m for diffuse light, S[m] the upward diffuse flux at m
    //     produced by the direct beam illuminating everything below m. The
    //     surface starts the recursion with A[0] = diffuse albedo and
    //     S[0] = direct albedo * F_dir(0).
    //  3. Downward pass from the TOA (no incident diffuse flux): diffuse
    //     downward d(m), diffuse upward u(m) = A[m] d(m) + S[m], and the net
    //     flux F_dir + d - u whose divergence gives the heating rate.
    // ========================================================================
    amrex::Real tau_sw_cum = 0.0;          // Column SW optical depth accumulated from the top
    amrex::Real sw_surface_absorbed = 0.0; // absorbed direct + diffuse at the surface
    amrex::Real sw_up_at_toa = 0.0;        // u(nlev)
    // Direct-beam albedo (per column, LSM field or fallback) and diffuse albedo
    // (surface_albedo_sw_diffuse, or the direct value when negative).
    const amrex::Real alb_dir = resolve_surface_albedo_sw(i, j, hetero_alb_sw, rad_choice, has_hetero_alb_sw);
    const amrex::Real alb_dif = (rad_choice.surface_albedo_sw_diffuse >= 0.0)
                                    ? rad_choice.surface_albedo_sw_diffuse : alb_dir;

    if (rad_choice.sw_enabled && cos_zenith > 0.0) {
        amrex::Real tau_sw_level[MAX_RAD_LEVELS];
        amrex::Real F_dir[MAX_RAD_LEVELS + 1];
        amrex::Real A[MAX_RAD_LEVELS + 1];
        amrex::Real S[MAX_RAD_LEVELS + 1];

        // 1. Per-layer optical depth and the direct beam on the interfaces.
        F_dir[nlev] = S0 * cos_zenith;
        for (int k = kmax; k >= kmin; --k) {
            const int m = k - kmin;
            amrex::Real omega_m = 0.0, g_m = 0.0;
            diagnose_layer_optics(i, j, k, dz_level[m], z_level[m], state_arr, tau_sw_base,
                                  /*is_sw=*/true, cloudy, rad_choice, tau_sw_level[m], omega_m, g_m);
            tau_sw_cum += tau_sw_level[m];
            F_dir[m] = compute_sw_direct_flux(tau_sw_cum, S0, cos_zenith);
        }

        // Two-stream properties of layer m (evaluated in both passes).
        auto layer_props = [&](int m) {
            amrex::Real tau_m = 0.0, omega = 0.0, g = 0.0;
            diagnose_layer_optics(i, j, kmin + m, dz_level[m], z_level[m], state_arr, tau_sw_base,
                                  /*is_sw=*/true, cloudy, rad_choice, tau_m, omega, g);
            return compute_sw_layer_two_stream(tau_sw_level[m], omega, g, cos_zenith);
        };

        // 2. Upward pass of the adding method.
        A[0] = alb_dif;
        S[0] = alb_dir * F_dir[0];
        for (int m = 0; m < nlev; ++m) {
            const TwoStreamLayerSW L = layer_props(m);
            amrex::Real denom = 1.0 - L.R_dif * A[m];
            if (denom < 1.0e-12) denom = 1.0e-12;
            A[m + 1] = L.R_dif + L.T_dif * L.T_dif * A[m] / denom;
            S[m + 1] = L.R_dir * F_dir[m + 1]
                     + L.T_dif * (S[m] + A[m] * L.T_dir * F_dir[m + 1]) / denom;
        }

        // 3. Downward pass: fluxes on the interfaces and layer heating rates.
        amrex::Real d_above = 0.0;          // No diffuse flux incident at the TOA
        amrex::Real u_above = S[nlev];
        sw_up_at_toa = u_above;
        amrex::Real F_net_above = F_dir[nlev] + d_above - u_above;
        for (int k = kmax; k >= kmin; --k) {
            const int m = k - kmin;
            const TwoStreamLayerSW L = layer_props(m);
            amrex::Real denom = 1.0 - L.R_dif * A[m];
            if (denom < 1.0e-12) denom = 1.0e-12;

            const amrex::Real d_m = (L.T_dif * d_above + L.T_dir * F_dir[m + 1] + L.R_dif * S[m]) / denom;
            const amrex::Real u_m = A[m] * d_m + S[m];
            const amrex::Real F_net_m = F_dir[m] + d_m - u_m;

            amrex::Real rho = state_arr(i, j, k, Rho_comp);
            if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

            // dT/dt from the flux divergence, stored as dtheta/dt = (dT/dt)/pi
            // (the RhoTheta source term and RRTMGP use the same convention).
            const amrex::Real Q_sw = compute_sw_heating_rate(F_net_above, F_net_m,
                                                             dz_level[m], rho, cp_air)
                                     * get_inverse_exner(i, j, k, state_arr);
            qheating_arr(i, j, k, 0) = Q_sw;
            local_max_heating = std::max(local_max_heating, std::abs(Q_sw));

            d_above = d_m;
            F_net_above = F_net_m;
            if (m == 0) {
                sw_surface_absorbed = (1.0 - alb_dir) * F_dir[0] + (1.0 - alb_dif) * d_m;
            }
        }
    }

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

        // Per-level LW heating rate from the net-flux divergence across
        // each layer (bottom interface m, top interface m+1).
        for (int k = kmin; k <= kmax; ++k) {
            const int m = k - kmin;
            amrex::Real rho = state_arr(i, j, k, Rho_comp);
            if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

            amrex::Real F_net_top = F_lw_up_iface[m + 1] - F_lw_down_iface[m + 1];
            amrex::Real F_net_bot = F_lw_up_iface[m]     - F_lw_down_iface[m];

            // dT/dt from the net-flux divergence, stored as dtheta/dt = (dT/dt)/pi.
            amrex::Real Q_lw = compute_lw_heating_rate(F_net_top, F_net_bot, dz_level[m], rho, cp_air)
                               * get_inverse_exner(i, j, k, state_arr);
            qheating_arr(i, j, k, 1) = Q_lw;

            amrex::Real Q_sw_here = qheating_arr(i, j, k, 0);
            local_max_heating = std::max(local_max_heating, std::abs(Q_sw_here) + std::abs(Q_lw));
        }
    }

    // ========================================================================
    // SURFACE AND DIAGNOSTICS
    // ========================================================================
    if (rad_choice.sw_enabled) {
        // Absorbed direct plus diffuse (see the SW block above)
        sw_surface_flux = sw_surface_absorbed;
        sw_up_toa = sw_up_at_toa;
    } else {
        sw_surface_flux = 0.0;
        sw_up_toa = 0.0;
    }

    if (rad_choice.lw_enabled) {
        amrex::Real F_lw_up_sfc   = F_lw_up_iface[0];
        amrex::Real F_lw_down_sfc = F_lw_down_iface[0];
        lw_up_toa = F_lw_up_iface[nlev];
        lw_net_surface = F_lw_up_sfc - F_lw_down_sfc;
    } else {
        lw_net_surface = 0.0;
        lw_up_toa = 0.0;
    }

    max_heating_rate = local_max_heating;
}

#endif // ERF_TWO_STREAM_COLUMN_H_
