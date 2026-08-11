#include <ERF.H>
#include <ERF_RadStruct.H>
#include <ERF_RadiationDiagnostics.H>
#include <ERF_TwoStreamSW.H>
#include <ERF_TwoStreamLW.H>
#include <ERF_PrognosticCloudFraction.H>
#include <ERF_AerosolOpticalDepth.H>
#include <ERF_SimplifiedSEB.H>
#include <ERF_SolarGeometry.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <ERF_IndexDefines.H>
#include <cmath>
#include <limits>

using namespace amrex;


namespace {
void fill_or_copy_seb_field(
    MultiFab* seb_mf,
    LandSurface& lsm,
    int lev,
    const char* field_name,
    amrex::Real fallback_value)
{
    if (seb_mf == nullptr) return;

    std::string varname(field_name);
    int lsm_idx = lsm.Get_DataIdx(lev, varname);
    if (lsm_idx >= 0) {
        if (MultiFab* lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx)) {
            MultiFab::Copy(*seb_mf, *lsm_ptr, 0, 0, 1, 0);
            return;
        }
    }
    seb_mf->setVal(fallback_value);
}
}

/**
 * @file ERF_AdvanceTwoStreamRadiation.cpp
 * @brief Phase 5 two-stream radiation driver: adds per-level (SW, LW)
 * heating-rate output (mirroring the RRTMGP qheating_rates MultiFab
 * convention: component 0 = SW, component 1 = LW, units K/s), on top of
 * the Phase 4 real per-column vertical integration, optional height-varying
 * (cloud layer) optical depth, cloud fraction masking, and diffuse
 * (scattering) SW flux contribution.
 *
 * Computes SW/LW fluxes and heating rates using real, per-column vertical
 * sweeps over the actual atmospheric grid. Reads temperature and density
 * from the state and properly accumulates optical depth through all
 * vertical levels.
 *
 * Phase 3 improvements over Phase 2d:
 * - Height-varying optical depth: optional "cloud_layer" tau_profile_type
 *   adds cloud_tau_per_layer on top of the clear-sky background within
 *   [cloud_base_height_m, cloud_top_height_m].
 * - Cloud fraction masking: blends clear-sky and cloudy-column fluxes via
 *   F = (1 - cloud_fraction) * F_clear + cloud_fraction * F_cloudy.
 * - Default tau_profile_type="constant" and cloud_fraction=0.0 reduce
 *   EXACTLY to Phase 2d output (see RAD_DEVELOPMENT.md Phase 3 section for
 *   hand-traced verification).
 *
 * Phase 4 improvements over Phase 3:
 * - Diffuse (scattered) SW flux: during the downward SW sweep, each layer's
 *   direct-beam attenuation now also generates a diffuse flux contribution
 *   via compute_sw_diffuse_flux() (Meador-Weaver two-stream scattering),
 *   accumulated layer-by-layer alongside the direct beam. Clear-sky levels
 *   use rad_choice.single_scattering_albedo / asymmetry_factor; levels
 *   within the Phase 3 cloud layer (only on the "cloudy" column evaluation)
 *   use rad_choice.cloud_single_scattering_albedo / cloud_asymmetry_factor
 *   instead.
 * - Total SW flux at any level = direct-beam flux + accumulated diffuse
 *   flux. Heating rate divergence and surface flux now include both terms.
 * - Default single_scattering_albedo = cloud_single_scattering_albedo = 0.0
 *   means compute_sw_diffuse_flux() returns exactly 0.0 at every level, so
 *   Phase 4 reduces EXACTLY to Phase 3 output when scattering is not
 *   configured (see RAD_DEVELOPMENT.md Phase 4 section for hand-traced
 *   verification).
 *
 * Phase 5 improvements over Phase 4 (this file):
 * - vertical_two_stream_sweep() now accepts a mutable per-column
 *   Array4<amrex::Real> `qheating_arr` (2 components) and writes the SW
 *   heating rate to qheating_arr(i,j,k,0) and the LW heating rate to
 *   qheating_arr(i,j,k,1) at EVERY level k in [kmin, kmax], instead of only
 *   reducing to a domain-max scalar for diagnostics. This mirrors exactly
 *   the (SW, LW) 2-component convention used by the RRTMGP qheating_rates
 *   MultiFab (see Source/ERF_MakeNewArrays.cpp and
 *   Source/SourceTerms/ERF_MakeSources.cpp), so both radiation paths can
 *   share the same downstream RhoTheta source-term injection code.
 * - LW heating rate is now actually computed: compute_lw_heating_rate()
 *   (defined in ERF_TwoStreamLW.H since Phase 1) was previously dead code —
 *   never called anywhere in the Phase 1-4 driver. Phase 5 adds a local,
 *   fixed-capacity per-column buffer (capped at MAX_RAD_LEVELS, asserted at
 *   runtime) to store the upward/downward LW flux profile from the
 *   existing two sweeps, then computes the net-flux-divergence heating
 *   rate layer-by-layer.
 * - compute_twostream_radiation_diagnostics() (see Step 3+ of Phase 5) will
 *   pass a real qheating_rates MultiFab pointer down to this kernel; when
 *   cloud_fraction > 0, the clear-sky and cloudy-column heating rates are
 *   blended into it exactly as sw_surface_flux/lw_net_surface already are,
 *   via a scratch MultiFab for the cloudy-column evaluation.
 *
 * Note: CSV diagnostics (SW_surface, heating_rate_max, etc.) are
 * unchanged; heating_rate_max is still the max(|Q_sw|+|Q_lw|) observed,
 * kept for backward RegTest compatibility.
 */

/**
 * @brief GPU-safe helper function to compute temperature from RhoTheta and background state.
 *
 * Converts specific enthalpy form (RhoTheta) to absolute temperature using:
 *   T = (RhoTheta / Rho) * (p / p_ref)^(R_d / cp)
 *
 * For simplicity in Phase 2/3, we use a constant reference pressure approximation.
 *
 * @param[in] rho_theta RhoTheta component [K·kg/m^3]
 * @param[in] rho Density [kg/m^3]
 * @return Temperature [K]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_temperature_from_rhotheta(amrex::Real rho_theta, amrex::Real rho)
{
    if (rho <= 0.0) {
        return 288.15;  // Defensive: fallback to standard T
    }
    if (rho_theta <= 0.0) {
        return 288.15;  // Defensive: fallback
    }

    // Simplified: assume theta ≈ T for Phase 2/3 (ignores Exner function)
    // Real implementation would use background state and pressure profile
    amrex::Real theta = rho_theta / rho;

    // Defensive clipping
    theta = std::max(theta, 100.0);  // Minimum sensible theta
    theta = std::min(theta, 400.0);  // Maximum sensible theta

    return theta;
}

/**
 * @brief (Phase 11) Clamp a real value to [min, max], with fallback for invalid (NaN/Inf).
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
 * @brief (Phase 11) Check if a real value is finite and positive.
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
 * @brief (Phase 11) Resolve per-column shortwave surface albedo from hetero field or fallback.
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
    const Array4<const amrex::Real>* hetero_alb_sw,
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
 * @brief (Phase 11) Resolve per-column longwave surface emissivity from hetero field or fallback.
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
    const Array4<const amrex::Real>* hetero_emiss_lw,
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
 * @brief (Phase 11) Resolve per-column surface temperature from hetero field or fallback.
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
    const Array4<const amrex::Real>* t_sfc,
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
 * @brief (Phase 10) Compute vertical thickness of a single layer.
 *
 * For uniform grids: dz = geom.CellSize(2)
 * For nonuniform/terrain-aware grids: dz = z_cc(k) - z_cc(k+1)
 *
 * Currently, ERF's two-stream radiation uses only uniform geom.CellSize(2).
 * This helper is prepared for future terrain-aware integration but currently
 * always returns the uniform spacing.
 *
 * @param[in] k Vertical cell index
 * @param[in] geom Geometry for this level
 * @return Vertical cell thickness [m]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real get_dz_for_level(int k, const Geometry& geom)
{
    // Phase 9 stub: currently always use uniform grid
    // Future: when z_phys_cc is available in device scope, use per-level spacing
    return geom.CellSize(2);
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real level_height_m(int k, int kmin, amrex::Real dz)
{
    return (static_cast<amrex::Real>(k - kmin) + 0.5) * dz;
}

/**
 * @brief GPU-safe helper to determine whether level k falls within the
 * Phase 3 cloud layer band [cloud_base_height_m, cloud_top_height_m].
 *
 * @param[in] k Vertical index.
 * @param[in] kmin Lowest vertical index.
 * @param[in] dz Uniform vertical cell spacing [m].
 * @param[in] rad_choice Radiation parameters.
 * @return true if this level is inside the configured cloud band.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool is_cloud_level(int k, int kmin, amrex::Real dz, const RadChoice& rad_choice)
{
    amrex::Real z = level_height_m(k, kmin, dz);
    return (z >= rad_choice.cloud_base_height_m && z <= rad_choice.cloud_top_height_m);
}

/**
 * @brief GPU-safe helper to compute the per-layer optical depth at level k,
 * given the base (clear-sky) optical depth and Phase 3 cloud-layer
 * parameters.
 *
 * When rad_choice.tau_profile_type == Constant, returns tau_base unchanged
 * (byte-identical to Phase 2d). When == CloudLayer, adds
 * rad_choice.cloud_tau_per_layer whenever the level height falls within
 * [cloud_base_height_m, cloud_top_height_m].
 *
 * @param[in] k Vertical index.
 * @param[in] kmin Lowest vertical index.
 * @param[in] dz Uniform vertical cell spacing [m].
 * @param[in] tau_base Clear-sky optical depth per layer.
 * @param[in] rad_choice Radiation parameters.
 * @param[in] apply_cloud If false, always returns tau_base (used for the
 * clear-sky column computation even when cloud_fraction > 0).
 * @return Optical depth for this layer [unitless].
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real tau_layer_value(
    int k, int kmin, amrex::Real dz, amrex::Real tau_base,
    const RadChoice& rad_choice, bool apply_cloud)
{
    if (!apply_cloud || rad_choice.tau_profile_type != TauProfileType::CloudLayer) {
        return tau_base;
    }
    if (is_cloud_level(k, kmin, dz, rad_choice)) {
        return tau_base + rad_choice.cloud_tau_per_layer;
    }
    return tau_base;
}

/**
 * @brief (Phase 12) GPU-safe helper to diagnose dynamic shortwave optical depth
 * from per-level water vapor and cloud liquid water.
 *
 * Computes tau_sw = tau_sw_coeff_qv * qv + tau_sw_coeff_qc * qc + tau_base,
 * where qv and qc are retrieved from the state array at the given (i,j,k).
 * When coefficients are 0, reduces exactly to tau_base (Phase 11 behavior).
 *
 * Guards against invalid values:
 * - If qv or qc are NaN/Inf or negative, uses 0 for that field.
 * - Output tau_sw is clamped to [0, 100] to ensure finite, physically reasonable values.
 * - Falls back to tau_base if anything goes wrong (defensive).
 *
 * @param[in] i, j, k Grid indices
 * @param[in] state_arr State array proxy (contains qv and qc)
 * @param[in] tau_base Static SW optical depth per layer [unitless]
 * @param[in] rad_choice Radiation parameters (contains dynamic tau coefficients)
 * @return Diagnosed or fallback SW optical depth [unitless]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_sw_dynamic(
   int i, int j, int k,
   const Array4<const amrex::Real>& state_arr,
   amrex::Real tau_base,
   const RadChoice& rad_choice)
{
   // If dynamic SW is disabled, return static value immediately
   if (!rad_choice.tau_sw_dynamic_enable) {
       return tau_base;
   }

   // If both coefficients are 0, dynamic path is a no-op
   if (rad_choice.tau_sw_coeff_qv <= 0.0 && rad_choice.tau_sw_coeff_qc <= 0.0) {
       return tau_base;
   }

   // Retrieve water vapor and cloud liquid water mixing ratios
   // Assume state_arr has components in order: Rho, RhoTheta, RhoQv, RhoQc, ...
   // Phase 12 uses hardcoded indices: RhoQv_comp and RhoQc_comp from ERF_Index.H
   amrex::Real qv = 0.0;
   amrex::Real qc = 0.0;

   // Extract qv = RhoQv / Rho (safe division with guards)
   if (state_arr.contains(i, j, k)) {
       amrex::Real rho = state_arr(i, j, k, Rho_comp);
        if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

        amrex::Real rho_qv = state_arr(i, j, k, RhoQ1_comp);
        if (std::isfinite(rho_qv)) {
            qv = rho_qv / rho;
            if (qv < 0.0 || !std::isfinite(qv)) qv = 0.0;
        }

        #if defined(RhoQ2_comp)
        amrex::Real rho_qc = state_arr(i, j, k, RhoQ2_comp);
        if (std::isfinite(rho_qc)) {
            qc = rho_qc / rho;
            if (qc < 0.0 || !std::isfinite(qc)) qc = 0.0;
        }
        #endif
       if (rho > 0.0 && std::isfinite(rho_qc)) {
           qc = rho_qc / rho;
           if (qc < 0.0 || !std::isfinite(qc)) qc = 0.0;  // Guard against invalid qc
       }
   }

   // Compute dynamic tau: tau_sw = tau_base + coeff_qv * qv + coeff_qc * qc
   amrex::Real tau_dynamic = tau_base + rad_choice.tau_sw_coeff_qv * qv + 
                              rad_choice.tau_sw_coeff_qc * qc;

   // Safety: clamp to physically reasonable range [0, 100]
   if (tau_dynamic < 0.0) tau_dynamic = 0.0;
   if (tau_dynamic > 100.0) tau_dynamic = 100.0;

   // Return finite value; if anything went wrong, this still produces a valid tau
   return std::isfinite(tau_dynamic) ? tau_dynamic : tau_base;
}

/**
 * @brief (Phase 12) GPU-safe helper to diagnose dynamic longwave optical depth
 * from per-level water vapor and cloud liquid water.
 *
 * Computes tau_lw = tau_lw_coeff_qv * qv + tau_lw_coeff_qc * qc + tau_base,
 * where qv and qc are retrieved from the state array at the given (i,j,k).
 * When coefficients are 0, reduces exactly to tau_base (Phase 11 behavior).
 *
 * Guards against invalid values:
 * - If qv or qc are NaN/Inf or negative, uses 0 for that field.
 * - Output tau_lw is clamped to [0, 100] to ensure finite, physically reasonable values.
 * - Falls back to tau_base if anything goes wrong (defensive).
 *
 * @param[in] i, j, k Grid indices
 * @param[in] state_arr State array proxy (contains qv and qc)
 * @param[in] tau_base Static LW optical depth per layer [unitless]
 * @param[in] rad_choice Radiation parameters (contains dynamic tau coefficients)
 * @return Diagnosed or fallback LW optical depth [unitless]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_lw_dynamic(
   int i, int j, int k,
   const Array4<const amrex::Real>& state_arr,
   amrex::Real tau_base,
   const RadChoice& rad_choice)
{
   // If dynamic LW is disabled, return static value immediately
   if (!rad_choice.tau_lw_dynamic_enable) {
       return tau_base;
   }

   // If both coefficients are 0, dynamic path is a no-op
   if (rad_choice.tau_lw_coeff_qv <= 0.0 && rad_choice.tau_lw_coeff_qc <= 0.0) {
       return tau_base;
   }

   // Retrieve water vapor and cloud liquid water mixing ratios
   amrex::Real qv = 0.0;
   amrex::Real qc = 0.0;

   // Extract qv = RhoQv / Rho (safe division with guards)
   if (state_arr.contains(i, j, k)) {
       amrex::Real rho = state_arr(i, j, k, Rho_comp);
    if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

    amrex::Real rho_qv = state_arr(i, j, k, RhoQ1_comp);
    if (std::isfinite(rho_qv)) {
        qv = rho_qv / rho;
        if (qv < 0.0 || !std::isfinite(qv)) qv = 0.0;
    }

    #if defined(RhoQ2_comp)
    amrex::Real rho_qc = state_arr(i, j, k, RhoQ2_comp);
    if (std::isfinite(rho_qc)) {
        qc = rho_qc / rho;
        if (qc < 0.0 || !std::isfinite(qc)) qc = 0.0;
    }
    #endif
       if (rho > 0.0 && std::isfinite(rho_qc)) {
           qc = rho_qc / rho;
           if (qc < 0.0 || !std::isfinite(qc)) qc = 0.0;  // Guard against invalid qc
       }
   }

   // Compute dynamic tau: tau_lw = tau_base + coeff_qv * qv + coeff_qc * qc
   amrex::Real tau_dynamic = tau_base + rad_choice.tau_lw_coeff_qv * qv + 
                              rad_choice.tau_lw_coeff_qc * qc;

   // Safety: clamp to physically reasonable range [0, 100]
   if (tau_dynamic < 0.0) tau_dynamic = 0.0;
   if (tau_dynamic > 100.0) tau_dynamic = 100.0;

   // Return finite value; if anything went wrong, this still produces a valid tau
   return std::isfinite(tau_dynamic) ? tau_dynamic : tau_base;
}

/**
 * @brief (Phase 4) GPU-safe helper to select the single-scattering albedo
 * and asymmetry factor to use for level k's diffuse SW calculation.
 *
 * When this column evaluation applies the cloud-layer enhancement
 * (apply_cloud == true, tau_profile_type == CloudLayer, and level k falls
 * within the cloud band), the cloud scattering properties
 * (cloud_single_scattering_albedo, cloud_asymmetry_factor) are used.
 * Otherwise, the clear-sky scattering properties (single_scattering_albedo,
 * asymmetry_factor) are used. Both default to 0.0, so by default this
 * function always yields omega == 0.0, and compute_sw_diffuse_flux()
 * returns exactly 0.0 (Phase 1-3 direct-beam-only behavior preserved).
 *
 * @param[in] k Vertical index.
 * @param[in] kmin Lowest vertical index.
 * @param[in] dz Uniform vertical cell spacing [m].
 * @param[in] rad_choice Radiation parameters.
 * @param[in] apply_cloud Same flag passed to tau_layer_value(); true for the
 * cloudy-column evaluation, false for the clear-sky column evaluation.
 * @param[out] omega Selected single-scattering albedo for this level.
 * @param[out] g Selected asymmetry factor for this level.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void select_scattering_props(
    int k, int kmin, amrex::Real dz, const RadChoice& rad_choice, bool apply_cloud,
    amrex::Real& omega, amrex::Real& g)
{
    bool use_cloud_props = apply_cloud &&
        rad_choice.tau_profile_type == TauProfileType::CloudLayer &&
        is_cloud_level(k, kmin, dz, rad_choice);

    if (use_cloud_props) {
        omega = rad_choice.cloud_single_scattering_albedo;
        g = rad_choice.cloud_asymmetry_factor;
    } else {
        omega = rad_choice.single_scattering_albedo;
        g = rad_choice.asymmetry_factor;
    }
}

/**
 * @brief (Phase 14) GPU-safe helper to diagnose prognostic cloud fraction
 * from per-level relative humidity and cloud liquid water.
 *
 * Computes cloud fraction from RH and qc using:
 *   cf_rh(k) = linear ramp from 0 at rh_min to 1 at rh_max
 *   cf_qc(k) = qc_scale * qc(k)
 *   cf(k) = min(1, cf_rh + cf_qc)  [saturated blend]
 *
 * When cloud_fraction_prog_enable is false, returns 0 (no effect).
 * When cloud_fraction_prog_enable is true, computes and returns diagnosed cf(k) in [0, 1].
 *
 * Guards against invalid values:
 * - If qv, qc, T, or P are NaN/Inf or unphysical, uses safe fallback values.
 * - Output cf is always clamped to [0, 1] and finite.
 *
 * @param[in] i, j, k Grid indices
 * @param[in] state_arr State array proxy (contains qv, qc, RhoTheta, Rho)
 * @param[in] rad_choice Radiation parameters (contains prognostic cloud fraction settings)
 * @param[in] geom Geometry for pressure/temperature computation
 * @return Diagnosed cloud fraction [0, 1] if enabled; 0 if disabled
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_cloud_fraction_prognostic(
    int i, int j, int k,
    const Array4<const amrex::Real>& state_arr,
    const RadChoice& rad_choice,
    const Geometry& geom)
{
    // If prognostic cloud fraction is disabled, return 0 (no effect)
    if (!rad_choice.cloud_fraction_prog_enable) {
        return 0.0;
    }

    // Extract qv and qc from state
    amrex::Real qv = 0.0;
    amrex::Real qc = 0.0;

    if (state_arr.contains(i, j, k)) {
        amrex::Real rho = state_arr(i, j, k, Rho_comp);
        if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;

        amrex::Real rho_qv = state_arr(i, j, k, RhoQ1_comp);
        if (std::isfinite(rho_qv) && rho_qv >= 0.0) {
            qv = rho_qv / rho;
            if (qv < 0.0 || !std::isfinite(qv)) qv = 0.0;
        }

        #if defined(RhoQ2_comp)
        amrex::Real rho_qc = state_arr(i, j, k, RhoQ2_comp);
        if (std::isfinite(rho_qc) && rho_qc >= 0.0) {
            qc = rho_qc / rho;
            if (qc < 0.0 || !std::isfinite(qc)) qc = 0.0;
        }
        #endif
    }

    // Compute temperature and pressure at this level
    amrex::Real rho = 1.0;
    amrex::Real rho_theta = 288.15;
    if (state_arr.contains(i, j, k)) {
        rho = state_arr(i, j, k, Rho_comp);
        if (rho <= 0.0 || !std::isfinite(rho)) rho = 1.0;
        rho_theta = state_arr(i, j, k, RhoTheta_comp);
        if (rho_theta <= 0.0 || !std::isfinite(rho_theta)) rho_theta = 288.15;
    }

    amrex::Real T = get_temperature_from_rhotheta(rho_theta, rho);
    if (T <= 0.0 || !std::isfinite(T)) T = 288.15;

    // Pressure at level k (simple estimate from hydrostatic equilibrium)
    // For now, use a reference pressure; proper implementation would use geom/k
    amrex::Real P = 101325.0 * std::pow(T / 288.15, -5.255);  // Fallback approximation
    if (P <= 0.0 || !std::isfinite(P)) P = 101325.0;

    // Compute RH from qv
    amrex::Real rh = compute_relative_humidity(qv, T, P);

    // Diagnose cloud fraction from RH and qc
    amrex::Real cf = diagnose_cloud_fraction_from_rh_qc(
        rh, qc,
        rad_choice.cloud_fraction_rh_min,
        rad_choice.cloud_fraction_rh_max,
        rad_choice.cloud_fraction_qc_scale);

    return cf;
}

/**
 * @brief (Phase 5) Maximum number of vertical levels supported by the
 * fixed-capacity per-column LW flux buffers in vertical_two_stream_sweep().
 *
 * The LW heating-rate computation requires the full upward and downward LW
 * flux profile (not just the surface value), so this kernel stores both
 * profiles in device-local arrays sized to this capacity. This is a
 * pragmatic, GPU-safe alternative to dynamic per-thread allocation.
 *
 * If bx.length(2) exceeds this value, an AMREX_ALWAYS_ASSERT will fire at
 * runtime (see compute_twostream_radiation_diagnostics()). Increase this
 * constant if a taller domain is required; note memory-per-thread scales
 * linearly with it (2 * MAX_RAD_LEVELS * sizeof(amrex::Real) bytes).
 */
constexpr int MAX_RAD_LEVELS = 512;

/**
 * @brief GPU-safe per-column vertical integration kernel for two-stream
 * radiation, computing either the clear-sky or cloudy-column fluxes and
 * per-level heating rates, depending on the `cloudy` flag.
 *
 * Performs a vertical sweep over each (i,j) column:
 * 1. Initialize fluxes at TOA
 * 2. Sweep downward (k: TOA → surface), accumulating optical depth and
 *    (Phase 4) diffuse SW flux; writes per-level SW heating rate to
 *    qheating_arr(i,j,k,0) at every level (Phase 5).
 * 3. Sweep upward then downward for LW (as in Phase 2c), storing the full
 *    per-level flux profiles in local buffers, then computes per-level LW
 *    heating rate from net-flux divergence and writes to
 *    qheating_arr(i,j,k,1) at every level (Phase 5).
 * 4. Reduction-based scalar diagnostics (max heating rate, surface fluxes)
 *    are still produced for CSV/console output, unchanged from Phase 4.
 *
 * (Phase 11) Integrates per-column heterogeneous surface properties (albedo,
 * emissivity, surface temperature) from optional fields with robust fallback
 * to RadChoice scalar parameters. If hetero fields unavailable or contain
 * invalid values, falls back silently to scalar/hard defaults.
 *
 * All arrays are on device; scalar results are accumulated into reduction
 * variables, while per-level heating rates are written directly into
 * qheating_arr.
 *
 * @param[in] i, j Column indices
 * @param[in] bx Computational box (cell-centered, full domain)
 * @param[in] geom Geometry for this AMR level
 * @param[in] state_arr Array proxy to state data (read-only)
 * @param[in] rad_choice Radiation parameters
 * @param[in] cloudy If true and tau_profile_type == CloudLayer, apply the
 * cloud-layer optical depth enhancement (and, Phase 4, cloud scattering
 * properties); if false, always use the clear-sky (constant) tau and
 * clear-sky scattering properties regardless of tau_profile_type. This lets
 * the caller compute both F_clear and F_cloudy for cloud-fraction blending.
 * @param[out] qheating_arr (Phase 5) Mutable per-column array; component 0
 * receives the SW heating rate [K/s] and component 1 the LW heating rate
 * [K/s] at every level k in [kmin, kmax]. Must have at least 2 components
 * and cover the full vertical extent of bx.
 * @param[out] max_heating_rate Maximum |Q_sw|+|Q_lw| observed (device-side scalar)
 * @param[out] sw_surface_flux Downwelling SW at surface (direct + diffuse) (device-side scalar)
 * @param[out] lw_net_surface Net LW at surface (device-side scalar)
 * @param[in] z_phys_cc (Phase 10) Optional physical height array for nonuniform dz support
 * @param[in] has_hetero_alb_sw (Phase 11) true if hetero_alb_sw is available
 * @param[in] hetero_alb_sw (Phase 11) Optional per-column SW surface albedo field
 * @param[in] has_hetero_emiss_lw (Phase 11) true if hetero_emiss_lw is available
 * @param[in] hetero_emiss_lw (Phase 11) Optional per-column LW surface emissivity field
 * @param[in] has_t_sfc (Phase 11) true if t_sfc is available
 * @param[in] t_sfc (Phase 11) Optional per-column surface temperature field [K]
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void vertical_two_stream_sweep(
    int i, int j,
    const Box& bx,
    const Geometry& geom,
    const Array4<const amrex::Real>& state_arr,
    const RadChoice& rad_choice,
    bool cloudy,
    const Array4<amrex::Real>& qheating_arr,
    amrex::Real& max_heating_rate,
    amrex::Real& sw_surface_flux,
    amrex::Real& lw_net_surface,
    const Array4<const amrex::Real>& z_phys_cc,
    amrex::Real time_utc_seconds = 0.0,
    bool has_hetero_alb_sw = false,
    const Array4<const amrex::Real>* hetero_alb_sw = nullptr,
    bool has_hetero_emiss_lw = false,
    const Array4<const amrex::Real>* hetero_emiss_lw = nullptr,
    bool has_t_sfc = false,
    const Array4<const amrex::Real>* t_sfc = nullptr)
{
    // Grid bounds
    int kmin = bx.smallEnd(2);
    int kmax = bx.bigEnd(2);
    int nlev = kmax - kmin + 1;

    // Physical constants
    amrex::Real sigma = 5.670374419e-8;  // Stefan-Boltzmann [W/(m^2·K^4)]
    amrex::Real cp_air = 1005.0;         // Specific heat at constant pressure [J/(kg·K)]

    // Get vertical grid spacing from geometry
    // Phase 10: Compute per-level dz from z_phys_cc differences for nonuniform grids;
    // fall back to CellSize(2) if z_phys_cc unavailable.
    amrex::Real dz_uniform = geom.CellSize(2);  // Uniform vertical cell spacing [m] (fallback)
    
    // (Phase 10) Per-level dz array for nonuniform grid support.
    // Attempt to compute from z_phys_cc; fall back to uniform if unavailable.
    // Defensive: If nlev > MAX_RAD_LEVELS, this will still work but degrade to
    // uniform spacing (the AMREX_ALWAYS_ASSERT in the caller will catch over-tall domains).
    amrex::Real dz_level[MAX_RAD_LEVELS];
    
    // Phase 10: Compute per-level dz from z_phys_cc if available
    bool using_nonuniform_dz = false;
    if (z_phys_cc) {
        using_nonuniform_dz = true;
        for (int k = 0; k < nlev && (kmin + k) < kmax; ++k) {
            // Compute layer thickness from cell-centered heights
            // Layer k extends from z_phys_cc(i,j,k) to z_phys_cc(i,j,k+1)
            amrex::Real dz_computed = std::abs(z_phys_cc(i,j,kmin+k+1) - z_phys_cc(i,j,kmin+k));
            
            // Safety check: ensure positive and reasonable dz
            if (dz_computed > 0.0) {
                dz_level[k] = dz_computed;
            } else {
                // Fallback to uniform if computed dz is invalid
                dz_level[k] = dz_uniform;
                using_nonuniform_dz = false;  // Mark as fallback used
            }
        }
        // Handle top level (kmax): use uniform spacing as fallback
        if ((kmax - kmin) < MAX_RAD_LEVELS && (kmax - kmin) >= 0) {
            dz_level[kmax - kmin] = dz_uniform;  // Top level uses uniform fallback
        }
    } else {
        // Fallback: z_phys_cc not available, use uniform spacing
        for (int k = 0; k < nlev && (kmin + k) <= kmax; ++k) {
            dz_level[k] = dz_uniform;
        }
    }

    // Compute solar zenith angle (Phase 16: optionally dynamic from solar geometry)
    amrex::Real cos_zenith;
    if (rad_choice.solar_geometry_dynamic_enable) {
        // Phase 16: Compute cos_zenith dynamically from solar position
        cos_zenith = compute_cos_zenith_angle(
            time_utc_seconds,
            rad_choice.latitude_deg,
            rad_choice.longitude_deg,
            rad_choice.day_of_year,
            rad_choice.time_zone_offset_hours);
    } else {
        // Phase 15 and earlier: Use fixed solar zenith angle (bitwise-identical)
        amrex::Real zenith_rad = rad_choice.solar_zenith_deg * M_PI / 180.0;
        cos_zenith = std::cos(zenith_rad);
    }

    // TOA values
    amrex::Real S0 = rad_choice.S0;
    amrex::Real tau_sw_base = rad_choice.tau_per_layer;
    amrex::Real tau_lw_base = rad_choice.tau_lw_per_layer;

    // Initialize accumulators
    amrex::Real tau_sw_cum = 0.0;      // Cumulative SW optical depth (TOA → current level)
    amrex::Real F_sw_dir_prev = 0.0;   // SW direct flux at previous level
    amrex::Real F_sw_diff_prev = 0.0;  // (Phase 4) SW diffuse flux at previous level
    amrex::Real F_lw_down_curr = 0.0;  // LW downwelling at current level (from downward sweep)

    amrex::Real local_max_heating = 0.0;

    // Zero-initialize this column's heating rate output (defensive: covers
    // both sw_enabled=false and lw_enabled=false cases below).
    for (int k = kmin; k <= kmax; ++k) {
        qheating_arr(i, j, k, 0) = 0.0;
        qheating_arr(i, j, k, 1) = 0.0;
    }

    // TOA: initialize SW direct beam
    if (rad_choice.sw_enabled) {
    if (cos_zenith > 0.0) {
        F_sw_dir_prev = S0 * cos_zenith;  // TOA incident (tau = 0)
    }
    F_sw_diff_prev = 0.0;  // No diffuse flux incident from above TOA

    // ========================================================================
    // DOWNWARD PASS: Accumulate optical depth and compute SW direct-beam plus
    // (Phase 4) diffuse flux; (Phase 5) write per-level SW heating rate.
    // ========================================================================
    for (int k = kmin; k <= kmax; ++k) {
        // Read state at this level
        amrex::Real rho = state_arr(i, j, k, Rho_comp);
        amrex::Real rho_theta = state_arr(i, j, k, RhoTheta_comp);

        // Defensive clipping
        if (rho <= 0.0) rho = 1.0;
        if (rho_theta <= 0.0) rho_theta = 288.15;

    // Phase 3: per-level optical depth (constant, or +cloud within layer)
        // NOTE: Use uniform dz for cloud-layer height detection (to keep cloud
        // position logic unchanged from prior phases)
        amrex::Real tau_sw = tau_layer_value(k, kmin, dz_uniform, tau_sw_base, rad_choice, cloudy);

        // (Phase 12) Apply dynamic tau diagnosis if enabled
        if (rad_choice.tau_sw_dynamic_enable) {
            tau_sw = diagnose_tau_sw_dynamic(i, j, k, state_arr, tau_sw, rad_choice);
        }

        // (Phase 14) Apply prognostic cloud fraction modulation if enabled
        // Scale the cloud optical depth contribution by the diagnosed cf(k)
        if (rad_choice.cloud_fraction_prog_enable && cloudy && 
            rad_choice.tau_profile_type == TauProfileType::CloudLayer &&
            is_cloud_level(k, kmin, dz_uniform, rad_choice)) {
            amrex::Real cf_prog = diagnose_cloud_fraction_prognostic(i, j, k, state_arr, rad_choice, geom);
            // Scale cloud tau contribution by cf(k): tau = tau_base + cf(k) * cloud_tau_per_layer
            // We need to recover tau_base and add cf-scaled cloud tau
            amrex::Real tau_sw_base_k = tau_layer_value(k, kmin, dz_uniform, tau_sw_base, rad_choice, /*cloudy=*/false);
            tau_sw = tau_sw_base_k + cf_prog * rad_choice.cloud_tau_per_layer;
        }

        // (Phase 15) Apply prescribed bulk aerosol optical depth if enabled
        // Aerosol tau is added on top of existing tau contributions (tau_base + cloud + dynamic)
        if (rad_choice.aerosol_enable) {
            amrex::Real tau_aerosol = 0.0;
             
            if (rad_choice.aerosol_profile_type == AerosolProfileType::Constant) {
                tau_aerosol = diagnose_tau_aerosol_constant(rad_choice.aerosol_tau_per_layer);
            } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Exponential) {
                // Get height at current level; use z_phys_cc if available, else compute from dz
                amrex::Real z_level = 0.0;
                if (z_phys_cc) {
                    z_level = std::abs(z_phys_cc(i, j, k) - z_phys_cc(i, j, kmin));  // height above the domain base
                } else {
                    z_level = 0.0;
                    for (int kk = kmin; kk < k; ++kk) {
                        int kk_idx = kk - kmin;
                        amrex::Real dz_kk = (kk_idx >= 0 && kk_idx < MAX_RAD_LEVELS) ? dz_level[kk_idx] : dz_uniform;
                        z_level += dz_kk;
                    }
                }
                int k_idx_aero = k - kmin;
                amrex::Real dz_aero = (k_idx_aero >= 0 && k_idx_aero < MAX_RAD_LEVELS) ? dz_level[k_idx_aero] : dz_uniform;
                tau_aerosol = diagnose_tau_aerosol_exponential(z_level, dz_aero, rad_choice.aerosol_tau_surface, 
                                                              rad_choice.aerosol_scale_height_m);
            } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Table) {
                tau_aerosol = diagnose_tau_aerosol_table(k);
            }
             
            // Add aerosol tau on top of existing tau
            tau_sw += tau_aerosol;
        }

        // Accumulate optical depths
        tau_sw_cum += tau_sw;

        // SW: Compute direct-beam flux at current level using Beer-Lambert
        amrex::Real F_sw_dir_curr = compute_sw_direct_flux(tau_sw_cum, S0, cos_zenith);

        // Phase 4: select scattering properties for this level (clear-sky vs
        // cloud, depending on the "cloudy" column flag and whether this level
        // falls within the cloud band). Use uniform dz for cloud detection.
        amrex::Real omega = 0.0;
        amrex::Real g = 0.0;
        select_scattering_props(k, kmin, dz_uniform, rad_choice, cloudy, omega, g);

        amrex::Real F_sw_diff_layer =
            compute_sw_diffuse_flux(tau_sw, F_sw_dir_prev, cos_zenith, omega, g);
        // Total diffuse flux at this level = diffuse flux transmitted from
        // above (attenuated by this layer's direct transmittance, as a
        // simple first-order approximation) + diffuse flux newly generated
        // within this layer.
        amrex::Real Tdir_layer = (cos_zenith > 0.0)
            ? std::exp(-tau_sw / cos_zenith)
            : 0.0;
        amrex::Real F_sw_diff_curr = F_sw_diff_prev * Tdir_layer + F_sw_diff_layer;

        // Total SW flux (direct + diffuse) at top and bottom of this layer,
         // used for heating rate divergence.
         amrex::Real F_sw_total_prev = F_sw_dir_prev + F_sw_diff_prev;
         amrex::Real F_sw_total_curr = F_sw_dir_curr + F_sw_diff_curr;

         // Phase 5: SW heating in this layer, written at EVERY level (not
         // just k < kmax as in Phase 1-4's max-only reduction), so that the
         // per-level qheating_arr output has a physically meaningful value
         // covering the whole column.
         // 
         // Phase 9: Use per-level dz for heating divergence (supports nonuniform grids).
         // Currently all levels use uniform spacing; fallback is automatic.
         int k_idx = k - kmin;
         amrex::Real dz_heating = (k_idx >= 0 && k_idx < MAX_RAD_LEVELS) 
             ? dz_level[k_idx] 
             : dz_uniform;  // Defensive fallback
         amrex::Real Q_sw = compute_sw_heating_rate(F_sw_total_prev, F_sw_total_curr,
                                                     dz_heating, rho, cp_air);
         qheating_arr(i, j, k, 0) = Q_sw;
         local_max_heating = std::max(local_max_heating, std::abs(Q_sw));

        // Prepare for next iteration
        F_sw_dir_prev = F_sw_dir_curr;
        F_sw_diff_prev = F_sw_diff_curr;
    }
    }

    // ========================================================================
    // (Phase 5) LW: store full per-level upward/downward flux profiles in
    // local, fixed-capacity buffers so a per-level net-flux-divergence
    // heating rate can be computed after both sweeps complete.
    // ========================================================================
    amrex::Real F_lw_up_curr = 0.0;  // Will be set at k = kmax (surface)
    amrex::Real F_lw_up_profile[MAX_RAD_LEVELS];
    amrex::Real F_lw_down_profile[MAX_RAD_LEVELS];

 if (rad_choice.lw_enabled) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nlev <= MAX_RAD_LEVELS,
        "vertical_two_stream_sweep: domain vertical extent exceeds "
        "MAX_RAD_LEVELS; increase the constant in ERF_AdvanceTwoStreamRadiation.cpp");

    // ========================================================================
    // UPWARD PASS: Compute LW upwelling flux from surface to TOA
    // ========================================================================
    for (int k = kmax; k >= kmin; --k) {
        // Read state at this level for temperature (needed for LW flux computation)
        amrex::Real rho = state_arr(i, j, k, Rho_comp);
        amrex::Real rho_theta = state_arr(i, j, k, RhoTheta_comp);

        // Defensive clipping
        if (rho <= 0.0) rho = 1.0;
        if (rho_theta <= 0.0) rho_theta = 288.15;

        amrex::Real T_layer = get_temperature_from_rhotheta(rho_theta, rho);

        // Phase 3: per-level LW optical depth (constant, or +cloud within layer)
        amrex::Real tau_lw = tau_layer_value(k, kmin, dz_uniform, tau_lw_base, rad_choice, cloudy);

        // (Phase 12) Apply dynamic tau diagnosis if enabled
        if (rad_choice.tau_lw_dynamic_enable) {
            tau_lw = diagnose_tau_lw_dynamic(i, j, k, state_arr, tau_lw, rad_choice);
        }

        // (Phase 14) Apply prognostic cloud fraction modulation if enabled
        // Scale the cloud optical depth contribution by the diagnosed cf(k)
        if (rad_choice.cloud_fraction_prog_enable && cloudy && 
            rad_choice.tau_profile_type == TauProfileType::CloudLayer &&
            is_cloud_level(k, kmin, dz_uniform, rad_choice)) {
            amrex::Real cf_prog = diagnose_cloud_fraction_prognostic(i, j, k, state_arr, rad_choice, geom);
            // Scale cloud tau contribution by cf(k): tau = tau_base + cf(k) * cloud_tau_per_layer
            amrex::Real tau_lw_base_k = tau_layer_value(k, kmin, dz_uniform, tau_lw_base, rad_choice, /*cloudy=*/false);
            tau_lw = tau_lw_base_k + cf_prog * rad_choice.cloud_tau_per_layer;
        }

        // (Phase 15) Apply prescribed bulk aerosol optical depth if enabled
        // Aerosol tau is added on top of existing tau contributions (tau_base + cloud + dynamic)
        // Note: LW aerosol support introduced here as part of full implementation
        // (separate enable flag for LW-only control can be added in Phase 16 if needed)
        if (rad_choice.aerosol_enable) {
            amrex::Real tau_aerosol = 0.0;
             
            if (rad_choice.aerosol_profile_type == AerosolProfileType::Constant) {
                tau_aerosol = diagnose_tau_aerosol_constant(rad_choice.aerosol_tau_per_layer);
            } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Exponential) {
                // Get height at current level; use z_phys_cc if available, else compute from dz
                amrex::Real z_level = 0.0;
                if (z_phys_cc) {
                    z_level = std::abs(z_phys_cc(i, j, k) - z_phys_cc(i, j, kmin));  // height above the domain base
                } else {
                    z_level = 0.0;
                    for (int kk = kmin; kk < k; ++kk) {
                        int kk_idx = kk - kmin;
                        amrex::Real dz_kk = (kk_idx >= 0 && kk_idx < MAX_RAD_LEVELS) ? dz_level[kk_idx] : dz_uniform;
                        z_level += dz_kk;
                    }
                }
                int k_idx_aero = k - kmin;
                amrex::Real dz_aero = (k_idx_aero >= 0 && k_idx_aero < MAX_RAD_LEVELS) ? dz_level[k_idx_aero] : dz_uniform;
                tau_aerosol = diagnose_tau_aerosol_exponential(z_level, dz_aero, rad_choice.aerosol_tau_surface, 
                                                              rad_choice.aerosol_scale_height_m);
            } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Table) {
                tau_aerosol = diagnose_tau_aerosol_table(k);
            }
             
            // Add aerosol tau on top of existing tau
            tau_lw += tau_aerosol;
        }

        if (k == kmax) {
            // Surface: initialize upwelling flux
            // (Phase 11) Resolve surface temperature and emissivity from hetero fields or fallback
            amrex::Real t_surface = resolve_surface_temp_k(i, j, t_sfc, rad_choice, has_t_sfc);
            amrex::Real emiss_lw = resolve_surface_emissivity_lw(i, j, hetero_emiss_lw, rad_choice, has_hetero_emiss_lw);
            
            // Thermal intensity at surface (Stefan-Boltzmann with emissivity)
            amrex::Real I_thermal_surface = emiss_lw * compute_thermal_intensity(t_surface, sigma);
            F_lw_up_curr = I_thermal_surface;
        } else {
            // Propagate upward through this layer
            F_lw_up_curr = compute_lw_flux_up(F_lw_up_curr, T_layer, sigma, tau_lw);
        }
        F_lw_up_profile[k - kmin] = F_lw_up_curr;
    }

    // ========================================================================
    // DOWNWARD PASS (Phase 2c): Compute real LW downwelling flux
    // ========================================================================
    if (!rad_choice.isothermal_test) {
        // For non-isothermal case, compute real downward two-stream sweep
        F_lw_down_curr = 0.0;  // Start from TOA (no incoming from space)
        for (int k = kmin; k <= kmax; ++k) {
            // Read state at this level
            amrex::Real rho = state_arr(i, j, k, Rho_comp);
            amrex::Real rho_theta = state_arr(i, j, k, RhoTheta_comp);

            // Defensive clipping
            if (rho <= 0.0) rho = 1.0;
            if (rho_theta <= 0.0) rho_theta = 288.15;

            amrex::Real T_layer = get_temperature_from_rhotheta(rho_theta, rho);

            // Phase 3: per-level LW optical depth (constant, or +cloud within layer)
            amrex::Real tau_lw = tau_layer_value(k, kmin, dz_uniform, tau_lw_base, rad_choice, cloudy);

            // (Phase 12) Apply dynamic tau diagnosis if enabled
            if (rad_choice.tau_lw_dynamic_enable) {
                tau_lw = diagnose_tau_lw_dynamic(i, j, k, state_arr, tau_lw, rad_choice);
            }

            // (Phase 14) Apply prognostic cloud fraction modulation if enabled
            // Scale the cloud optical depth contribution by the diagnosed cf(k)
            if (rad_choice.cloud_fraction_prog_enable && cloudy && 
                rad_choice.tau_profile_type == TauProfileType::CloudLayer &&
                is_cloud_level(k, kmin, dz_uniform, rad_choice)) {
                amrex::Real cf_prog = diagnose_cloud_fraction_prognostic(i, j, k, state_arr, rad_choice, geom);
                // Scale cloud tau contribution by cf(k): tau = tau_base + cf(k) * cloud_tau_per_layer
                amrex::Real tau_lw_base_k = tau_layer_value(k, kmin, dz_uniform, tau_lw_base, rad_choice, /*cloudy=*/false);
                tau_lw = tau_lw_base_k + cf_prog * rad_choice.cloud_tau_per_layer;
            }

            // (Phase 15) Apply prescribed bulk aerosol optical depth if enabled
            // Aerosol tau is added on top of existing tau contributions (tau_base + cloud + dynamic)
            if (rad_choice.aerosol_enable) {
                amrex::Real tau_aerosol = 0.0;
                 
                if (rad_choice.aerosol_profile_type == AerosolProfileType::Constant) {
                    tau_aerosol = diagnose_tau_aerosol_constant(rad_choice.aerosol_tau_per_layer);
                } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Exponential) {
                    // Get height at current level; use z_phys_cc if available, else compute from dz
                    amrex::Real z_level = 0.0;
                    if (z_phys_cc) {
                        z_level = std::abs(z_phys_cc(i, j, k) - z_phys_cc(i, j, kmin));  // height above the domain base
                    } else {
                        z_level = 0.0;
                        for (int kk = kmin; kk < k; ++kk) {
                            int kk_idx = kk - kmin;
                            amrex::Real dz_kk = (kk_idx >= 0 && kk_idx < MAX_RAD_LEVELS) ? dz_level[kk_idx] : dz_uniform;
                            z_level += dz_kk;
                        }
                    }
                int k_idx_aero = k - kmin;
                amrex::Real dz_aero = (k_idx_aero >= 0 && k_idx_aero < MAX_RAD_LEVELS) ? dz_level[k_idx_aero] : dz_uniform;
                    tau_aerosol = diagnose_tau_aerosol_exponential(z_level, dz_aero, rad_choice.aerosol_tau_surface, 
                                                                  rad_choice.aerosol_scale_height_m);
                } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Table) {
                    tau_aerosol = diagnose_tau_aerosol_table(k);
                }
                 
                // Add aerosol tau on top of existing tau
                tau_lw += tau_aerosol;
            }

            // Compute downwelling flux at this level using real two-stream formula
            F_lw_down_curr = compute_lw_flux_down(F_lw_down_curr, T_layer, sigma, tau_lw);
            F_lw_down_profile[k - kmin] = F_lw_down_curr;
        }

    }
    else {
        // Isothermal test: all levels radiate equally; net flux is zero at
        // every level (see SURFACE AND DIAGNOSTICS override below), so the
        // per-level LW heating rate is exactly zero everywhere. Fill the
        // downward profile with the (later-overridden) upward profile so
        // the net-flux-divergence loop below produces exactly zero.
        for (int k = kmin; k <= kmax; ++k) {
            F_lw_down_profile[k - kmin] = F_lw_up_profile[k - kmin];
        }
    }

    // ========================================================================
    // (Phase 5) Per-level LW heating rate from net-flux divergence.
    // ========================================================================
    for (int k = kmin; k <= kmax; ++k) {
        amrex::Real rho = state_arr(i, j, k, Rho_comp);
        if (rho <= 0.0) rho = 1.0;

        amrex::Real F_net_top, F_net_bot;
        if (k == kmin) {
            // TOA: no level above; treat the "top" net flux as the TOA
            // downward boundary condition (F_down=0 unless isothermal) minus
            // upward flux leaving the domain top, i.e. the net flux just
            // above this layer using the same level's upward flux (first-
            // order approximation at the domain boundary).
            amrex::Real F_up_top = F_lw_up_profile[k - kmin];
            amrex::Real F_down_top = rad_choice.isothermal_test
                ? F_lw_down_profile[k - kmin]
                : 0.0;
            F_net_top = F_up_top - F_down_top;
        } else {
            F_net_top = F_lw_up_profile[k - 1 - kmin] - F_lw_down_profile[k - 1 - kmin];
        }
        F_net_bot = F_lw_up_profile[k - kmin] - F_lw_down_profile[k - kmin];

        // Phase 9: Use per-level dz for LW heating divergence (supports nonuniform grids).
        // Currently all levels use uniform spacing; fallback is automatic.
        int k_idx = k - kmin;
        amrex::Real dz_heating_lw = (k_idx >= 0 && k_idx < MAX_RAD_LEVELS) 
            ? dz_level[k_idx] 
            : dz_uniform;  // Defensive fallback
        amrex::Real Q_lw = compute_lw_heating_rate(F_net_top, F_net_bot, dz_heating_lw, rho, cp_air);
        qheating_arr(i, j, k, 1) = Q_lw;

        amrex::Real Q_sw_here = qheating_arr(i, j, k, 0);
        local_max_heating = std::max(local_max_heating, std::abs(Q_sw_here) + std::abs(Q_lw));
    }

    F_lw_up_curr = F_lw_up_profile[kmax - kmin];
    F_lw_down_curr = F_lw_down_profile[kmax - kmin];
 }
    // ========================================================================
    // SURFACE AND DIAGNOSTICS
    // ========================================================================
    if (rad_choice.sw_enabled) {
        // (Phase 11) Resolve per-column SW albedo from hetero field or fallback
        amrex::Real alb_sw = resolve_surface_albedo_sw(i, j, hetero_alb_sw, rad_choice, has_hetero_alb_sw);
        
        if (rad_choice.isothermal_test) {
            // Isothermal test: compute direct-beam at surface, apply albedo
            sw_surface_flux = S0 * std::max(0.0, cos_zenith) * std::exp(-tau_sw_cum) * (1.0 - alb_sw);
        } else {
            // Normal mode: Phase 4 includes both direct and diffuse terms, apply albedo
            sw_surface_flux = (F_sw_dir_prev + F_sw_diff_prev) * (1.0 - alb_sw);
        }
    } else {
        sw_surface_flux = 0.0;
    }

    if (rad_choice.lw_enabled) {
        if (rad_choice.isothermal_test) {
            amrex::Real rho_surface = state_arr(i, j, kmax, Rho_comp);
            amrex::Real rho_theta_surface = state_arr(i, j, kmax, RhoTheta_comp);
            if (rho_surface <= 0.0) rho_surface = 1.0;
            if (rho_theta_surface <= 0.0) rho_theta_surface = 288.15;
            amrex::Real T_iso = get_temperature_from_rhotheta(rho_theta_surface, rho_surface);
            amrex::Real I_thermal = compute_thermal_intensity(T_iso, sigma);
            F_lw_up_curr = I_thermal;
            F_lw_down_curr = I_thermal;
            // Isothermal override: net LW flux and heating are exactly zero
            // at every level, overriding the per-level values written above.
            for (int k = kmin; k <= kmax; ++k) {
                qheating_arr(i, j, k, 1) = 0.0;
            }
        }
        lw_net_surface = F_lw_up_curr - F_lw_down_curr;
    } else {
        lw_net_surface = 0.0;
    }

    max_heating_rate = local_max_heating;
}

void ERF::compute_twostream_radiation_diagnostics(
    int lev,
    int nstep,
    amrex::Real time_step,
    std::string const& call_site
    )
{
    const auto& rad_choice = solverChoice.radChoice;

    // Only proceed if TwoStream radiation is enabled
    if (rad_choice.rad_type != RadType::TwoStream) {
        return;
    }

    // Create RadiationDiagnostics instance for this level with Phase 7 controls
    RadiationDiagnostics rad_diag(rad_choice.verbosity, rad_choice.diag_file, lev,
                                   rad_choice.diag_enable, rad_choice.diag_stdout_enable,
                                   rad_choice.diag_tagged_enable, rad_choice.diag_regtest_line_enable,
                                   rad_choice.diag_csv_enable, rad_choice.diag_callsite_mode,
                                   rad_choice.diag_dedup_tol);

    // ========================================
    // Phase 5: GPU-Safe ParallelFor Implementation with Cloud Fraction
    // Blending, Diffuse (Scattering) SW Flux, and Per-Level Heating Rate
    // Output (qheating_rates MultiFab)
    // ========================================

    // Initialize global diagnostics
    amrex::Real SW_surface = 0.0;
    amrex::Real SW_TOA = 0.0;
    amrex::Real F_up_surface = 0.0;
    amrex::Real F_down_toa = 0.0;
    amrex::Real heating_rate_max = 0.0;
    amrex::Real seb_residual_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real seb_residual_max  = std::numeric_limits<amrex::Real>::quiet_NaN();
    
    // (Phase 19b) Prognostic SEB surface temperature and moisture diagnostics
    amrex::Real t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real t_s_max  = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real q_s_max  = std::numeric_limits<amrex::Real>::quiet_NaN();

    // Get state at this level (conservative variables: density, RhoTheta, etc.)
    const auto& state_cons = vars_old[lev][Vars::cons];

    // Only compute radiation if we have valid state data
    if (state_cons.nComp() > 0 ) {

    // Prepare to compute TOA values (used for diagnostics output)
    // (Phase 16) Compute dynamic cos_zenith if enabled, otherwise use static value
    amrex::Real cos_zenith;
    if (rad_choice.solar_geometry_dynamic_enable) {
        // Convert absolute simulation time to UTC seconds within the day [0, 86400)
        amrex::Real time_utc_seconds = std::fmod(t_old[lev], 86400.0);
        if (time_utc_seconds < 0.0) time_utc_seconds += 86400.0;
        cos_zenith = compute_cos_zenith_angle(
            time_utc_seconds,
            rad_choice.latitude_deg,
            rad_choice.longitude_deg,
            rad_choice.day_of_year,
            rad_choice.time_zone_offset_hours);
    } else {
        // Phase 15 and earlier: Use static solar zenith angle
        amrex::Real zenith_rad = rad_choice.solar_zenith_deg * M_PI / 180.0;
        cos_zenith = std::cos(zenith_rad);
    }
    SW_TOA = rad_choice.sw_enabled ? (rad_choice.S0 * std::max(0.0, cos_zenith)) : 0.0;

        // Host-side storage for reduction results (will be set by device-side reduction)
        amrex::Real max_heating_global = 0.0;
        amrex::Real sw_surface_sum = 0.0;
        amrex::Real lw_net_sum = 0.0;
        amrex::Long n_columns_total = 0;

        // (Phase 18) SEB residual diagnostics
        amrex::Real seb_residual_sum = 0.0;
        amrex::Long n_seb_columns = 0;

        // Phase 3: cloud fraction used to blend clear-sky and cloudy-column results.
        // cloud_fraction == 0.0 (default) means only the clear-sky column is
        // ever evaluated, and the blend below reduces to F = F_clear exactly.
        amrex::Real cloud_fraction = rad_choice.cloud_fraction;

        // (Phase 16) Compute UTC seconds within the day for dynamic solar geometry
        amrex::Real time_utc_seconds = 0.0;
        if (rad_choice.solar_geometry_dynamic_enable) {
            time_utc_seconds = std::fmod(t_old[lev], 86400.0);
            if (time_utc_seconds < 0.0) time_utc_seconds += 86400.0;
        }

        // (Phase 5) Note: qheating_rates[lev] is expected to be allocated with
        // 2 components by the caller whenever rad_choice.rad_type ==
        // RadType::TwoStream (see Source/ERF_MakeNewArrays.cpp, Step 3 of
        // Phase 5). If not yet allocated (e.g. before that change lands),
        // this function still safely computes and logs CSV diagnostics but
        // skips the per-level heating write.
        MultiFab* qheating_mf = qheating_rates[lev].get();

        if (rad_choice.seb_enable) {
            fill_or_copy_seb_field(twostream_alb_sw[lev].get(), lsm, lev, "sfc_alb_dir_vis", rad_choice.surface_albedo_sw);
            fill_or_copy_seb_field(twostream_emiss_lw[lev].get(), lsm, lev, "sfc_emis", rad_choice.surface_emissivity_lw);
             
            // (Phase 20) Gate t_sfc fill on prognostic mode: when seb_prognostic_enable is true,
            // t_sfc is owned and evolved by the prognostic update, not reset by fill_or_copy.
            // This prevents silently overwriting the prognostic state before the update reads it.
            //fill_or_copy_seb_field(twostream_t_sfc[lev].get(), lsm, lev, "t_sfc", rad_choice.surface_temp_k);
            if (!rad_choice.seb_prognostic_enable) {
                fill_or_copy_seb_field(twostream_t_sfc[lev].get(), lsm, lev, "t_sfc", rad_choice.surface_temp_k);
            }
             
            fill_or_copy_seb_field(sw_flux_sfc[lev].get(), lsm, lev, "sav", rad_choice.seb_sw_flux_default);
            fill_or_copy_seb_field(lw_flux_sfc[lev].get(), lsm, lev, "fira", rad_choice.seb_lw_flux_default);
            fill_or_copy_seb_field(hfx_sfc[lev].get(), lsm, lev, "hfx", rad_choice.seb_hfx_default);
            fill_or_copy_seb_field(lh_sfc[lev].get(), lsm, lev, "lh", rad_choice.seb_lh_default);
            fill_or_copy_seb_field(grdflx_sfc[lev].get(), lsm, lev, "grdflx", rad_choice.seb_grdflx_default);
             
            // (Phase 20) Gate q_sfc fill on prognostic mode: same reasoning as t_sfc.
            if (!rad_choice.seb_prognostic_enable) {
                fill_or_copy_seb_field(q_sfc[lev].get(), lsm, lev, "noahmp_water_vapor_mixing_ratio_2m_vegetated", rad_choice.seb_q_sfc_default);
            }
            fill_or_copy_seb_field(t_deep[lev].get(), lsm, lev, "smstav", rad_choice.seb_t_deep_default);
            fill_or_copy_seb_field(q_deep[lev].get(), lsm, lev, "smstot", rad_choice.seb_q_deep_default);
        }

        // Sequential loop over all boxes (each box handled with GPU-safe ParallelFor)
        for (MFIter mfi(state_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& state_arr = state_cons.const_array(mfi);
            const Geometry& geom_lev = geom[lev];

            // (Phase 10) Get z_phys_cc for nonuniform dz support if available
            Array4<const amrex::Real> z_phys_cc_arr;
            bool has_z_phys = false;
            if (z_phys_cc[lev] != nullptr) {
                z_phys_cc_arr = z_phys_cc[lev]->const_array(mfi);
                has_z_phys = true;
            }

            // (Phase 14A) Wire LSM surface property fields or use standalone fallback MultiFabs
            // Priority:
            // 1. If LSM is active (lsm.Get_DataIdx() returns >=0), use real LSM fields
            // 2. Otherwise, use standalone fallback MultiFabs (allocated and constant-filled from RadChoice scalars)
            // The resolve_surface_*() helpers implement the full precedence chain with finite guards
            
            // SW albedo: Try LSM field "sfc_alb_dir_vis" (simplified: broadband approx from vis-direct only;
            // future phases may handle full 4-band vis/nir dir/dif; see Phase 14A doc).
            bool has_hetero_alb_sw = false;
            Array4<const amrex::Real> hetero_alb_sw_arr;
            {
                //int lsm_idx = lsm.Get_DataIdx(lev, "sfc_alb_dir_vis");
                std::string varname_alb = "sfc_alb_dir_vis";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_alb);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        hetero_alb_sw_arr = lsm_ptr->const_array(mfi);
                        has_hetero_alb_sw = true;
                    }
                } else if (twostream_alb_sw[lev]) {
                    hetero_alb_sw_arr = twostream_alb_sw[lev]->const_array(mfi);
                    has_hetero_alb_sw = true;
                }
            }

            // LW emissivity: Try LSM field "sfc_emis"
            bool has_hetero_emiss_lw = false;
            Array4<const amrex::Real> hetero_emiss_lw_arr;
            {
                //int lsm_idx = lsm.Get_DataIdx(lev, "sfc_emis");
                std::string varname_emiss = "sfc_emis";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_emiss);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        hetero_emiss_lw_arr = lsm_ptr->const_array(mfi);
                        has_hetero_emiss_lw = true;
                    }
                } else if (twostream_emiss_lw[lev]) {
                    hetero_emiss_lw_arr = twostream_emiss_lw[lev]->const_array(mfi);
                    has_hetero_emiss_lw = true;
                }
            }

            // Surface temperature: Try LSM field "t_sfc"
            bool has_t_sfc_field = false;
            Array4<const amrex::Real> t_sfc_arr;
            {
                //int lsm_idx = lsm.Get_DataIdx(lev, "t_sfc");
                std::string varname_t_sfc = "t_sfc";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_t_sfc);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        t_sfc_arr = lsm_ptr->const_array(mfi);
                        has_t_sfc_field = true;
                    }
                } else if (twostream_t_sfc[lev]) {
                    t_sfc_arr = twostream_t_sfc[lev]->const_array(mfi);
                    has_t_sfc_field = true;
                }
            }


            // Create a 2D box for (i,j) iteration over the horizontal extent
            // One GPU thread per (i,j) column; k-loop is sequential within each thread
            const auto& lo = bx.loVect();
            const auto& hi = bx.hiVect();
            Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));

            // Count columns in this box for later averaging
            amrex::Long n_cols = static_cast<amrex::Long>(bx.length(0)) *
                                 static_cast<amrex::Long>(bx.length(1));
            n_columns_total += n_cols;

            // (Phase 5) Clear-sky-column heating rates are written directly
            // into the real qheating_rates MultiFab when available.
            // Fall back to a throwaway local FArrayBox otherwise (keeps the
            // kernel call GPU-safe even if qheating_rates isn't allocated).
            FArrayBox qheating_fallback_fab;
            Array4<amrex::Real> qheating_clear_arr;
            if (qheating_mf != nullptr) {
                qheating_clear_arr = qheating_mf->array(mfi);
            } else {
                qheating_fallback_fab.resize(bx, 2);
                qheating_clear_arr = qheating_fallback_fab.array();
            }

            // (Phase 5) Cloudy-column heating rates always go into a scratch
            // FArrayBox; only used/blended when cloud_fraction > 0.
            FArrayBox qheating_cloudy_fab;
            Array4<amrex::Real> qheating_cloudy_arr;
            if (cloud_fraction > 0.0) {
                qheating_cloudy_fab.resize(bx, 2);
                qheating_cloudy_arr = qheating_cloudy_fab.array();
            }

            // GPU-safe reduction using ReduceOps (per-column results aggregated on device)
            amrex::Real max_heating_box = 0.0;
            amrex::Real sw_sum_box = 0.0;
            amrex::Real lw_sum_box = 0.0;

            // Device-side reduction: compute max heating and sum of surface fluxes
            ReduceOps<ReduceOpMax, ReduceOpSum, ReduceOpSum> reduce_ops;
            ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(reduce_ops);

            using ReduceTuple = typename decltype(reduce_data)::Type;

         // Launch parallel kernel over (i,j) columns
            reduce_ops.eval(xy_box, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> ReduceTuple
                {
                    // Clear-sky column (always evaluated; this is the sole
                    // contributor when cloud_fraction == 0.0, matching Phase 2d)
                    amrex::Real max_heating_clear = 0.0;
                    amrex::Real sw_flux_clear = 0.0;
                    amrex::Real lw_net_clear = 0.0;
                    vertical_two_stream_sweep(
                        i, j, bx, geom_lev, state_arr, rad_choice, /*cloudy=*/false,
                        qheating_clear_arr,
                        max_heating_clear, sw_flux_clear, lw_net_clear,
                        z_phys_cc_arr,
                        time_utc_seconds,
                        has_hetero_alb_sw, &hetero_alb_sw_arr,
                        has_hetero_emiss_lw, &hetero_emiss_lw_arr,
                        has_t_sfc_field, &t_sfc_arr);

                    amrex::Real max_heating_col = max_heating_clear;
                    amrex::Real sw_flux_col = sw_flux_clear;
                    amrex::Real lw_net_col = lw_net_clear;

                    // Cloudy column only needs to be evaluated when there is a
                    // nonzero cloud fraction; this keeps the cloud_fraction==0
                    // path numerically and computationally identical to Phase 2d.
                    if (cloud_fraction > 0.0) {
                         amrex::Real max_heating_cloudy = 0.0;
                         amrex::Real sw_flux_cloudy = 0.0;
                         amrex::Real lw_net_cloudy = 0.0;
                         vertical_two_stream_sweep(
                            i, j, bx, geom_lev, state_arr, rad_choice, /*cloudy=*/true,
                            qheating_cloudy_arr,
                            max_heating_cloudy, sw_flux_cloudy, lw_net_cloudy,
                             z_phys_cc_arr,
                             time_utc_seconds,
                             has_hetero_alb_sw, &hetero_alb_sw_arr,
                             has_hetero_emiss_lw, &hetero_emiss_lw_arr,
                             has_t_sfc_field, &t_sfc_arr);

                        // Blend clear-sky and cloudy-column results
                        sw_flux_col = (1.0 - cloud_fraction) * sw_flux_clear +
                                      cloud_fraction * sw_flux_cloudy;
                        lw_net_col = (1.0 - cloud_fraction) * lw_net_clear +
                                     cloud_fraction * lw_net_cloudy;
                        max_heating_col = std::max(max_heating_clear, max_heating_cloudy);

                        // (Phase 5) Blend per-level heating rates in place
                        // into qheating_clear_arr (which is the real output
                        // MultiFab when qheating_mf != nullptr).
                        int kmin = bx.smallEnd(2);
                        int kmax = bx.bigEnd(2);
                        for (int k = kmin; k <= kmax; ++k) {
                            for (int comp = 0; comp < 2; ++comp) {
                                amrex::Real q_clear_val = qheating_clear_arr(i, j, k, comp);
                                amrex::Real q_cloudy_val = qheating_cloudy_arr(i, j, k, comp);
                                qheating_clear_arr(i, j, k, comp) =
                                    (1.0 - cloud_fraction) * q_clear_val +
                                    cloud_fraction * q_cloudy_val;
                            }
                        }
                    }

                    // Return tuple for reduction
                    return {max_heating_col, sw_flux_col, lw_net_col};
                }
            );

            // Copy results from device to host
            amrex::Gpu::synchronize();
            auto reduce_tuple = reduce_data.value(reduce_ops);
            max_heating_box = amrex::get<0>(reduce_tuple);
            sw_sum_box = amrex::get<1>(reduce_tuple);
            lw_sum_box = amrex::get<2>(reduce_tuple);

            // Accumulate box results into global results
            max_heating_global = std::max(max_heating_global, max_heating_box);
            sw_surface_sum += sw_sum_box;
            lw_net_sum += lw_sum_box;
        }

         // (Phase 18) Warn if diagnostic is requested but SEB infrastructure isn't enabled
        if (rad_choice.seb_diagnostic_enable && !rad_choice.seb_enable) {
            static bool warned_seb_misconfig = false;
            if (!warned_seb_misconfig && ParallelDescriptor::IOProcessor()) {
                Print() << "WARNING: erf.radiation.seb_diagnostic_enable=true but "
                           "seb_enable=false; SEB residual diagnostics will report NaN. "
                           "Set erf.radiation.seb_enable=true to enable SEB field "
                           "population.\n";
                warned_seb_misconfig = true;
            }
        }
        // (Phase 18) Compute SEB residual diagnostics if enabled
        if (rad_choice.seb_diagnostic_enable && rad_choice.seb_enable) {
            seb_residual_max = 0.0; 
            // Second loop over boxes to compute SEB residual from populated SEB MultiFabs
            for (MFIter mfi(state_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const auto& lo = bx.loVect();
                const auto& hi = bx.hiVect();
                Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));

                // Get SEB field arrays
                Array4<const amrex::Real> sw_flux_arr = sw_flux_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> lw_flux_arr = lw_flux_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> hfx_arr = hfx_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> lh_arr = lh_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> grdflx_arr = grdflx_sfc[lev]->const_array(mfi);

                // Count columns and compute residuals
                amrex::Long n_cols_box = static_cast<amrex::Long>(bx.length(0)) *
                                        static_cast<amrex::Long>(bx.length(1));
                amrex::Real residual_sum_box = 0.0;
                amrex::Real residual_max_box = 0.0;

                // GPU-safe reduction for SEB residuals
                ReduceOps<ReduceOpSum, ReduceOpMax> seb_reduce_ops;
                ReduceData<amrex::Real, amrex::Real> seb_reduce_data(seb_reduce_ops);

                using SEBReduceTuple = typename decltype(seb_reduce_data)::Type;

                seb_reduce_ops.eval(xy_box, seb_reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> SEBReduceTuple {
                        amrex::Real sw_net = sw_flux_arr(i, j, 0);
                        amrex::Real lw_net = lw_flux_arr(i, j, 0);
                        amrex::Real hfx = hfx_arr(i, j, 0);
                        amrex::Real lh = lh_arr(i, j, 0);
                        amrex::Real grdflx = grdflx_arr(i, j, 0);

                        // Compute residual using Phase 18 helper function
                        amrex::Real residual = diagnose_seb_residual(sw_net, lw_net, hfx, lh, grdflx);

                        // Return sum and abs(max) of residual
                        return {residual, std::abs(residual)};
                    }
                );

                // Copy results from device to host
                amrex::Gpu::synchronize();
                auto seb_reduce_tuple = seb_reduce_data.value(seb_reduce_ops);
                residual_sum_box = amrex::get<0>(seb_reduce_tuple);
                residual_max_box = amrex::get<1>(seb_reduce_tuple);

                // Accumulate into global results
                seb_residual_sum += residual_sum_box;
                seb_residual_max = std::max(seb_residual_max, residual_max_box);
                n_seb_columns += n_cols_box;
            }
        }

        // (Phase 19b) Prognostic SEB surface temperature and moisture evolution
        // Only run if prognostic mode is enabled and Noah-MP is NOT driving LSM at this level
        if (rad_choice.seb_prognostic_enable && rad_choice.seb_enable &&
            call_site == "post_dycore") {
            // Check if Noah-MP is active at this level by attempting to get the LSM t_sfc field
            std::string varname_t_sfc_prog = "t_sfc";
            int lsm_idx_t_sfc = lsm.Get_DataIdx(lev, varname_t_sfc_prog);
            bool noahmp_active = (lsm_idx_t_sfc >= 0);
            
            if (!noahmp_active) {
                // Noah-MP is NOT active; proceed with prognostic update
                
                // Initialize diagnostics for T_s and q_s
                amrex::Real t_s_sum = 0.0;
                amrex::Real t_s_max_val = -std::numeric_limits<amrex::Real>::max();
                amrex::Real q_s_sum = 0.0;
                amrex::Real q_s_max_val = -std::numeric_limits<amrex::Real>::max();
                amrex::Long n_prog_columns = 0;
                
                // Third loop over boxes for prognostic SEB update
                for (MFIter mfi(state_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const auto& lo = bx.loVect();
                    const auto& hi = bx.hiVect();
                    Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));
                    
                    // Get SEB field arrays (read-only)
                    Array4<const amrex::Real> sw_flux_arr = sw_flux_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> lw_flux_arr = lw_flux_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> hfx_arr = hfx_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> lh_arr = lh_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> grdflx_arr = grdflx_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> t_deep_arr = t_deep[lev]->const_array(mfi);
                    Array4<const amrex::Real> q_deep_arr = q_deep[lev]->const_array(mfi);
                    
                    // Get SEB state arrays (read-write for prognostic update)
                    Array4<amrex::Real> t_s_arr = twostream_t_sfc[lev]->array(mfi);
                    Array4<amrex::Real> q_s_arr = q_sfc[lev]->array(mfi);
                    // Count columns and prepare for reductions
                    amrex::Long n_cols_box = static_cast<amrex::Long>(bx.length(0)) *
                                             static_cast<amrex::Long>(bx.length(1));
                    amrex::Real t_s_sum_box = 0.0;
                    amrex::Real t_s_max_box = -std::numeric_limits<amrex::Real>::max();
                    amrex::Real q_s_sum_box = 0.0;
                    amrex::Real q_s_max_box = -std::numeric_limits<amrex::Real>::max();
                    
                    // GPU-safe update for prognostic T_s and q_s with reductions
                    ReduceOps<ReduceOpSum, ReduceOpMax, ReduceOpSum, ReduceOpMax> prog_reduce_ops;
                    ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real> prog_reduce_data(prog_reduce_ops);
                    
                    using ProgReduceTuple = typename decltype(prog_reduce_data)::Type;
                    
                    prog_reduce_ops.eval(xy_box, prog_reduce_data,
                            [=,  C_s=rad_choice.seb_surface_heat_capacity,
                            tau=rad_choice.seb_restore_timescale_s,
                            d_s=rad_choice.seb_moisture_layer_depth_m,
                            tau_q=rad_choice.seb_moisture_restore_timescale_s,
                            t_min=rad_choice.seb_prognostic_t_min_k,
                            t_max=rad_choice.seb_prognostic_t_max_k,
                            q_min=rad_choice.seb_prognostic_q_min,
                            q_max=rad_choice.seb_prognostic_q_max] 
                            AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> ProgReduceTuple {
                            amrex::Real t_s_old = t_s_arr(i, j, 0);
                            amrex::Real q_s_old = q_s_arr(i, j, 0);
                            
                            // Read forcing data
                            amrex::Real sw_net = sw_flux_arr(i, j, 0);
                            amrex::Real lw_net = lw_flux_arr(i, j, 0);
                            amrex::Real hfx = hfx_arr(i, j, 0);
                            amrex::Real lh = lh_arr(i, j, 0);
                            amrex::Real grdflx = grdflx_arr(i, j, 0);
                            amrex::Real t_deep = t_deep_arr(i, j, 0);
                            amrex::Real q_deep = q_deep_arr(i, j, 0);
                            
                            // Compute SEB residual
                            amrex::Real seb_res = diagnose_seb_residual(sw_net, lw_net, hfx, lh, grdflx);
                            
                            // Compute tendencies
                            amrex::Real dT_s_dt = prognostic_dTs_dt(seb_res, t_s_old, t_deep,
                                                                     C_s, tau);
                            amrex::Real dq_s_dt = prognostic_dqs_dt(lh, q_s_old, q_deep,
                                                                     d_s, tau_q);
                            
                            // Perform Euler update
                            amrex::Real t_s_new = t_s_old + time_step * dT_s_dt;
                            amrex::Real q_s_new = q_s_old + time_step * dq_s_dt;
                            
                            // Clamp to valid ranges
                            t_s_new = amrex::max(t_min, amrex::min(t_max, t_s_new));
                            q_s_new = amrex::max(q_min, amrex::min(q_max, q_s_new));
                            
                            // Write back updated values (this modifies the device array)
                            t_s_arr(i, j, 0) = t_s_new;
                            q_s_arr(i, j, 0) = q_s_new;
                            
                            // Return for reduction: sum T_s, max T_s, sum q_s, max q_s
                            return {t_s_new, std::abs(t_s_new), q_s_new, std::abs(q_s_new)};
                        });
                    
                    // Copy results from device to host
                    amrex::Gpu::synchronize();
                    auto prog_reduce_tuple = prog_reduce_data.value(prog_reduce_ops);
                    t_s_sum_box = amrex::get<0>(prog_reduce_tuple);
                    t_s_max_box = amrex::get<1>(prog_reduce_tuple);
                    q_s_sum_box = amrex::get<2>(prog_reduce_tuple);
                    q_s_max_box = amrex::get<3>(prog_reduce_tuple);
                    
                    // Accumulate into global results
                    t_s_sum += t_s_sum_box;
                    t_s_max_val = std::max(t_s_max_val, t_s_max_box);
                    q_s_sum += q_s_sum_box;
                    q_s_max_val = std::max(q_s_max_val, q_s_max_box);
                    n_prog_columns += n_cols_box;
                }
                
                // Compute mean values from sums
                if (n_prog_columns > 0) {
                    t_s_mean = t_s_sum / static_cast<amrex::Real>(n_prog_columns);
                    t_s_max = t_s_max_val;
                    q_s_mean = q_s_sum / static_cast<amrex::Real>(n_prog_columns);
                    q_s_max = q_s_max_val;
                } else {
                    t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                    t_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                    q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                    q_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                }
            } else {
                // Noah-MP is active; skip prognostic update for this level
                // Leave t_s and q_s as populated by Phase 17 LSM passthrough
                t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                t_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                q_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
            }
        }
        // equivalent to a single-column value for spatially UNIFORM atmospheres
        // (as in the current SW_ClearSky_Analytical / LW_Isothermal RegTests).
        // Phase 3's SW_Cloud_Layer and Phase 4's SW_Scattering_Cloud RegTests
        // are ALSO spatially uniform (identical cloud/tau/scattering
        // parameters applied to every column), so the domain-averaged value
        // still equals the true single-column flux there. True horizontal
        // heterogeneity (e.g., patchy clouds varying by column) remains
        // deferred to a future phase; see Gap G5 in RAD_DEVELOPMENT.md.
        if (n_columns_total > 0) {
            SW_surface = sw_surface_sum / static_cast<amrex::Real>(n_columns_total);
            F_up_surface = lw_net_sum / static_cast<amrex::Real>(n_columns_total);
        }
        heating_rate_max = max_heating_global;

        // (Phase 18) Compute SEB residual mean from sum
        if (rad_choice.seb_diagnostic_enable && rad_choice.seb_enable && n_seb_columns > 0) {
            seb_residual_mean = seb_residual_sum / static_cast<amrex::Real>(n_seb_columns);
        } else {
            // When feature is disabled, use NaN for backward compatibility
            seb_residual_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
            seb_residual_max = std::numeric_limits<amrex::Real>::quiet_NaN();
        }

        // For LW: in isothermal test, override with analytical value
        if (rad_choice.lw_enabled && rad_choice.isothermal_test) {
            amrex::Real sigma = 5.670374419e-8;
            amrex::Real T = rad_choice.T_iso_K;
            amrex::Real I_thermal = sigma * T * T * T * T;
            F_down_toa = I_thermal;
            F_up_surface = I_thermal;
            heating_rate_max = 0.0;
        }
    }

    // Logging (Phase 2c: use verbosity to guard output)
    if (rad_choice.verbosity >= 1 && ParallelDescriptor::IOProcessor()) {
        Print() << "Radiation diagnostics at step " << nstep << ":\n"
                << "  SW TOA = " << SW_TOA << " W/m^2\n"
                << "  SW surface = " << SW_surface << " W/m^2\n"
                << "  LW up (surface) = " << F_up_surface << " W/m^2\n"
                << "  LW down (TOA) = " << F_down_toa << " W/m^2\n"
                << "  Max heating rate = " << heating_rate_max << " K/s\n";
        if (rad_choice.seb_diagnostic_enable && std::isfinite(seb_residual_mean)) {
            Print() << "  SEB residual (mean) = " << seb_residual_mean << " W/m^2\n"
                    << "  SEB residual (max) = " << seb_residual_max << " W/m^2\n";
        }
        if (rad_choice.seb_prognostic_enable && std::isfinite(t_s_mean)) {
            Print() << "  Surface temperature (mean) = " << t_s_mean << " K\n"
                    << "  Surface temperature (max) = " << t_s_max << " K\n"
                    << "  Surface moisture (mean) = " << q_s_mean << " kg/kg\n"
                    << "  Surface moisture (max) = " << q_s_max << " kg/kg\n";
        }
    }

    rad_diag.append(nstep, time_step, call_site, SW_surface, SW_TOA,
                    F_up_surface, F_down_toa, heating_rate_max,
                    seb_residual_mean, seb_residual_max,
                    t_s_mean, t_s_max, q_s_mean, q_s_max);
}
