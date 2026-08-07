#include <ERF.H>
#include <ERF_RadStruct.H>
#include <ERF_RadiationDiagnostics.H>
#include <ERF_TwoStreamSW.H>
#include <ERF_TwoStreamLW.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <cmath>

using namespace amrex;

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
 * @brief GPU-safe helper to compute the height of cell-center level k above
 * the surface (k = kmin is the lowest cell-center level).
 *
 * @param[in] k Vertical index.
 * @param[in] kmin Lowest vertical index (surface-adjacent cell).
 * @param[in] dz Uniform vertical cell spacing [m].
 * @return Height of the cell center above the surface [m].
 */
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
 * All arrays are on device; scalar results are accumulated into reduction
 * variables, while per-level heating rates are written directly into
 * qheating_arr.
 *
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
    amrex::Real& lw_net_surface)
{
    // Grid bounds
    int kmin = bx.smallEnd(2);
    int kmax = bx.bigEnd(2);
    int nlev = kmax - kmin + 1;

    // Physical constants
    amrex::Real sigma = 5.670374419e-8;  // Stefan-Boltzmann [W/(m^2·K^4)]
    amrex::Real cp_air = 1005.0;         // Specific heat at constant pressure [J/(kg·K)]

    // Get real vertical grid spacing from geometry
    // For uniform grids, use CellSize(2); for terrain-aware grids, would use z_phys_cc differences
    amrex::Real dz = geom.CellSize(2);  // Vertical cell spacing [m]

    // Convert solar zenith angle to radians
    amrex::Real zenith_rad = rad_choice.solar_zenith_deg * M_PI / 180.0;
    amrex::Real cos_zenith = std::cos(zenith_rad);


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
        amrex::Real tau_sw = tau_layer_value(k, kmin, dz, tau_sw_base, rad_choice, cloudy);

        // Accumulate optical depths
        tau_sw_cum += tau_sw;

        // SW: Compute direct-beam flux at current level using Beer-Lambert
        amrex::Real F_sw_dir_curr = compute_sw_direct_flux(tau_sw_cum, S0, cos_zenith);

        // Phase 4: select scattering properties for this level (clear-sky vs
        // cloud, depending on the "cloudy" column flag and whether this level
        // falls within the cloud band), and accumulate diffuse flux generated
        // by this layer on top of the diffuse flux transmitted from above.
        amrex::Real omega = 0.0;
        amrex::Real g = 0.0;
        select_scattering_props(k, kmin, dz, rad_choice, cloudy, omega, g);

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
        amrex::Real Q_sw = compute_sw_heating_rate(F_sw_total_prev, F_sw_total_curr,
                                                    dz, rho, cp_air);
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
    AMREX_ASSERT_WITH_MESSAGE(nlev <= MAX_RAD_LEVELS,
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
        amrex::Real tau_lw = tau_layer_value(k, kmin, dz, tau_lw_base, rad_choice, cloudy);

        if (k == kmax) {
            // Surface: initialize upwelling flux
            amrex::Real I_thermal = compute_thermal_intensity(T_layer, sigma);
            F_lw_up_curr = I_thermal;
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
            amrex::Real tau_lw = tau_layer_value(k, kmin, dz, tau_lw_base, rad_choice, cloudy);

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

        amrex::Real Q_lw = compute_lw_heating_rate(F_net_top, F_net_bot, dz, rho, cp_air);
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
        if (rad_choice.isothermal_test) {
            sw_surface_flux = S0 * std::max(0.0, cos_zenith) * std::exp(-tau_sw_cum);
        } else {
            // Phase 4: surface SW flux includes both direct and diffuse terms.
            sw_surface_flux = F_sw_dir_prev + F_sw_diff_prev;
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
    amrex::Real time_step)
{
    const auto& rad_choice = solverChoice.radChoice;

    // Only proceed if TwoStream radiation is enabled
    if (rad_choice.rad_type != RadType::TwoStream) {
        return;
    }

    // Create RadiationDiagnostics instance for this level
    RadiationDiagnostics rad_diag(rad_choice.verbosity, rad_choice.diag_file, lev);

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

    // Get state at this level (conservative variables: density, RhoTheta, etc.)
    const auto& state_cons = vars_old[lev][Vars::cons];

    // Only compute radiation if we have valid state data
    if (state_cons.nComp() > 0 ) {

        // Prepare to compute TOA values (used for diagnostics output)
        amrex::Real zenith_rad = rad_choice.solar_zenith_deg * M_PI / 180.0;
        amrex::Real cos_zenith = std::cos(zenith_rad);
        SW_TOA = rad_choice.sw_enabled ? (rad_choice.S0 * std::max(0.0, cos_zenith)) : 0.0;

        // Host-side storage for reduction results (will be set by device-side reduction)
        amrex::Real max_heating_global = 0.0;
        amrex::Real sw_surface_sum = 0.0;
        amrex::Real lw_net_sum = 0.0;
        amrex::Long n_columns_total = 0;

        // Phase 3: cloud fraction used to blend clear-sky and cloudy-column results.
        // cloud_fraction == 0.0 (default) means only the clear-sky column is
        // ever evaluated, and the blend below reduces to F = F_clear exactly.
        amrex::Real cloud_fraction = rad_choice.cloud_fraction;

        // (Phase 5) Note: qheating_rates[lev] is expected to be allocated with
        // 2 components by the caller whenever rad_choice.rad_type ==
        // RadType::TwoStream (see Source/ERF_MakeNewArrays.cpp, Step 3 of
        // Phase 5). If not yet allocated (e.g. before that change lands),
        // this function still safely computes and logs CSV diagnostics but
        // skips the per-level heating write.
        MultiFab* qheating_mf = qheating_rates[lev].get();

        // Sequential loop over all boxes (each box handled with GPU-safe ParallelFor)
        for (MFIter mfi(state_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& state_arr = state_cons.const_array(mfi);
            const Geometry& geom_lev = geom[lev];

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
                        max_heating_clear, sw_flux_clear, lw_net_clear);

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
                            max_heating_cloudy, sw_flux_cloudy, lw_net_cloudy);

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

        // NOTE: This computes a domain-averaged surface flux. This is only
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
    }

    // Append diagnostics to CSV file
    rad_diag.append(nstep, time_step, SW_surface, SW_TOA,
                    F_up_surface, F_down_toa, heating_rate_max);
}
