#ifndef ERF_CROWN_FIRE_H
#define ERF_CROWN_FIRE_H

#include <AMReX_MultiFab.H>
#include <cmath>

/**
 * @file ERF_CrownFire.H
 * @brief Crown fire initiation and active crown fire rate of spread
 *
 * Implements multiple crown fire models from published literature:
 *   1. Van Wagner (1977) criterion for surface-to-crown fire transition
 *   2. Cruz et al. (2005) model for active crown fire rate of spread
 *   3. Rothermel (1991) crown ROS as multiple of surface ROS
 *   4. Van Wagner proxy model for crown ROS
 *   5. Van Wagner passive-to-active blending
 *   6. Crown Fraction Burned (CFB) diagnostic
 *
 * Assumptions:
 *   - Canopy parameters are uniform across the domain; spatially varying canopy
 *     parameters are deferred to Phase 10
 *   - Cruz et al. (2005) was calibrated for North American conifer stands and is
 *     not appropriate for grass or shrub fuel models FM1–FM3, FM5–FM6
 *   - Crown fire cells do not return to surface fire state once activated
 *
 * References:
 *   Van Wagner, C.E. (1977). Conditions for the start and spread of crown fire.
 *     Canadian Journal of Forest Research, 7(1), 23–34.
 *   Cruz, M.G., Alexander, M.E. & Wakimoto, R.H. (2005). Development and testing
 *     of models for predicting crown fire rate of spread in conifer forest stands.
 *     Canadian Journal of Forest Research, 35(7), 1626–1639.
 *   Rothermel, R.C. (1991). Predicting behavior and size of crown fires in the
 *     northern Rocky Mountain conifers. USDA Forest Service INT-RP-438.
 *   Scott, J.H. & Reinhardt, E.D. (2001). Assessing crown fire potential by linking
 *     models of surface and crown fire behavior. USDA Forest Service RMRS-RP-29.
 */

/**
 * @brief Van Wagner (1977) critical fireline intensity for crown fire initiation
 *
 * Formula used here:
 *   I_B_crit = 0.010 * C_BH * (460 + 25.9 * FMC_excess_pct)
 * where FMC_excess_pct = 100 * (foliar_moisture - M_c).
 *
 * If foliar_moisture <= M_c the transition is treated as unattainable and a
 * large sentinel value is returned.
 *
 * @param[in] C_BH              Canopy base height [m]
 * @param[in] foliar_moisture   Foliar moisture [fraction, 0-1]
 * @param[in] M_c               Critical foliar moisture [fraction]
 * @return Critical intensity [kW/m], or 1.0e12 kW/m when foliar_moisture <= M_c
 *
 * Reference: Van Wagner, C.E. (1977). Conditions for the start and spread of
 *   crown fire. Canadian Journal of Forest Research, 7(1), 23–34.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real van_wagner_critical_intensity(amrex::Real C_BH, amrex::Real foliar_moisture,
                                           amrex::Real M_c) noexcept
{
    constexpr amrex::Real LARGE_INTENSITY = 1.0e12;

    if (foliar_moisture <= M_c) {
        return LARGE_INTENSITY;
    }

    C_BH = amrex::max(C_BH, 0.1);
    amrex::Real FMC_excess_pct = amrex::max((foliar_moisture - M_c) * 100.0, 0.0);
    return 0.010 * C_BH * (460.0 + 25.9 * FMC_excess_pct);
}

/**
 * @brief Cruz et al. (2005) active crown fire rate of spread
 *
 * Formula: R [m/min] = 11.02 * U_10^0.9 * C_BD^0.19 * exp(-0.17 * MC10_pct)
 * Returns value in m/s (input formula is m/min, so divide by 60).
 *
 * Note: moisture_10hr is converted from fraction to percent for the formula.
 *
 * @param[in] U_10_ms       10-meter wind speed [m/s] (converted to km/h internally)
 * @param[in] C_BD          Canopy bulk density [kg/m³]
 * @param[in] moisture_10hr 10-hr dead fuel moisture [fraction, 0-1]
 * @return Crown fire rate of spread [m/s], or 0 when U_10_ms <= 0
 *
 * Reference: Cruz, M.G., Alexander, M.E. & Wakimoto, R.H. (2005). Development
 *   and testing of models for predicting crown fire rate of spread in conifer
 *   forest stands. Canadian Journal of Forest Research, 35(7), 1626–1639.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cruz_crown_ros(amrex::Real U_10_ms, amrex::Real C_BD, amrex::Real moisture_10hr) noexcept
{
    // Zero or negative wind produces zero crown ROS
    if (U_10_ms <= 0.0) {
        return 0.0;
    }

    // Convert wind from m/s to km/h
    amrex::Real U_kmh = U_10_ms * 3.6;
    U_kmh = amrex::max(U_kmh, 1.0);  // Minimum 1 km/h (model not calibrated at calm)

    // Convert moisture from fraction to percent
    amrex::Real MC10_pct = moisture_10hr * 100.0;

    // Clamp inputs to valid ranges
    C_BD = amrex::max(C_BD, 0.01);
    MC10_pct = amrex::max(MC10_pct, 0.0);

    // Cruz formula: R [m/min] = 11.02 * U^0.90 * CBD^0.19 * exp(-0.17 * MC10)
    constexpr amrex::Real A_COEFF  = 11.02;   // leading coefficient
    constexpr amrex::Real W_EXP    = 0.90;    // wind exponent
    constexpr amrex::Real CBD_EXP  = 0.19;    // CBD exponent
    constexpr amrex::Real MC_COEFF = 0.17;    // moisture coefficient

    amrex::Real R_m_min = A_COEFF
                        * std::pow(U_kmh, W_EXP)
                        * std::pow(C_BD, CBD_EXP)
                        * std::exp(-MC_COEFF * MC10_pct);

    // Convert m/min to m/s
    amrex::Real R_m_s = R_m_min / 60.0;
    return amrex::max(0.0, R_m_s);
}

/**
 * @brief Rothermel (1991) active crown fire rate of spread
 *
 * Simple empirical formula: R_crown = 3.34 * R_surface
 *
 * @param[in] R_surface_ms  Surface fire rate of spread [m/s]
 * @return Crown fire rate of spread [m/s]
 *
 * Reference: Rothermel, R.C. (1991). Predicting behavior and size of crown fires
 *   in the northern Rocky Mountain conifers. USDA Forest Service INT-RP-438.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_rothermel_1991_crown_ros(amrex::Real R_surface_ms) noexcept
{
    constexpr amrex::Real ROTHERMEL_MULT = 3.34;
    return amrex::max(ROTHERMEL_MULT * R_surface_ms, 0.0);
}

/**
 * @brief Van Wagner (1977) proxy active crown fire rate of spread
 *
 * Formula: R_crown = (3.0 / CBD) * moisture_factor / 60.0
 * where moisture_factor = max(0.3, min(1.0, 1.0 - (FMC - 100) / 200))
 *
 * @param[in] C_BD              Canopy bulk density [kg/m³]
 * @param[in] foliar_moisture   Foliar moisture [fraction, 0-1]
 * @return Crown fire rate of spread [m/s]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_van_wagner_proxy_ros(amrex::Real C_BD, amrex::Real foliar_moisture) noexcept
{
    // Convert moisture to percent
    amrex::Real FMC_pct = foliar_moisture * 100.0;

    // Clamp CBD
    C_BD = amrex::max(C_BD, 0.01);

    // Base crown ROS [m/min] = 3.0 / CBD
    amrex::Real R_m_min = 3.0 / C_BD;

    // Moisture factor (reduction for higher moisture)
    amrex::Real mf = 1.0 - (FMC_pct - 100.0) / 200.0;
    mf = amrex::max(0.3, amrex::min(1.0, mf));

    R_m_min *= mf;

    // Convert m/min to m/s
    return amrex::max(0.0, R_m_min / 60.0);
}

/**
 * @brief Van Wagner (1977) passive-to-active crown fire blending
 *
 * A continuous crowning fraction CF blends between surface and active crown ROS:
 *   CF = min((I_B / I_o)^(2/3), 1.0)
 *   R_total = R_surface * (1 - CF) + R_crown_active * CF
 *
 * @param[in] R_surface_ms       Surface fire ROS [m/s]
 * @param[in] R_crown_active_ms  Active crown fire ROS [m/s]
 * @param[in] I_B_kwm            Byram fireline intensity [kW/m]
 * @param[in] I_o_kwm            Van Wagner critical surface intensity [kW/m]
 * @return Blended ROS [m/s]
 *
 * Reference: Van Wagner, C.E. (1977). Conditions for the start and spread of
 *   crown fire. Canadian Journal of Forest Research, 7(1), 23–34.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_van_wagner_passive_blend(
    amrex::Real R_surface_ms,
    amrex::Real R_crown_active_ms,
    amrex::Real I_B_kwm,
    amrex::Real I_o_kwm) noexcept
{
    if (I_o_kwm < 1.0e-6 || I_B_kwm <= 0.0) return R_surface_ms;

    amrex::Real ratio = amrex::max(I_B_kwm / I_o_kwm, 0.0);
    amrex::Real CF = amrex::min(std::pow(ratio, 2.0 / 3.0), 1.0);
    return R_surface_ms * (1.0 - CF) + amrex::max(R_crown_active_ms, R_surface_ms) * CF;
}

/**
 * @brief Crown Fraction Burned (CFB) diagnostic
 *
 * Computes the crown fraction burned, which distinguishes between passive
 * (torching) and active crown fire. Used for emissions and fire effects modeling.
 *
 * Formula: CFB = (R_total - R_surface) / (R_active_crown - R_surface)
 * Clipped to range [0, 1].
 *
 * @param[in] R_total            Total (observed or modeled) ROS [m/s]
 * @param[in] R_surface          Surface fire ROS [m/s]
 * @param[in] R_active_crown     Active crown fire ROS [m/s]
 * @return CFB value [0-1]
 *
 * Reference: Scott, J.H. & Reinhardt, E.D. (2001). Assessing crown fire potential
 *   by linking models of surface and crown fire behavior. USDA Forest Service
 *   Research Paper RMRS-RP-29.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_crown_fraction_burned(
    amrex::Real R_total,
    amrex::Real R_surface,
    amrex::Real R_active_crown) noexcept
{
    // Ensure positive rates
    R_total       = amrex::max(R_total,       0.0);
    R_surface     = amrex::max(R_surface,     0.0);
    R_active_crown = amrex::max(R_active_crown, 0.0);

    // If active crown ROS not greater than surface, no crown fire
    if (R_active_crown <= R_surface + 1.0e-8) {
        return 0.0;
    }

    // Scott & Reinhardt (2001) CFB formula
    amrex::Real CFB = (R_total - R_surface) / (R_active_crown - R_surface);

    // Clamp to valid range [0, 1]
    CFB = amrex::max(0.0, amrex::min(1.0, CFB));

    return CFB;
}

/**
 * @brief Update crown fire status and rate of spread
 *
 * GPU kernel over fire grid cells. For each burned cell (phi < 0):
 *   - If fireline intensity >= I_B_crit, activate crown fire (set flag to 1.0)
 *   - If already crown-active or just activated, replace ROS with crown value
 *   - Once activated, the flag is never reset to 0.0
 *
 * @param[out] fire_crown_active  Crown fire status: 1.0=active, 0.0=surface
 * @param[out] fire_ros           Rate of spread field [m/s]
 * @param[in]  fire_phi           Level-set field (<0 = burned)
 * @param[in]  fireline_intensity Byram intensity [kW/m]
 * @param[in]  I_B_crit           Critical intensity for initiation [kW/m]
 * @param[in]  R_crown_ms         Crown fire ROS [m/s]
 * @param[in]  fire_debug         Print debug diagnostics if true
 */
inline void update_crown_fire_status(amrex::MultiFab& fire_crown_active,
                                      amrex::MultiFab& fire_ros,
                                      const amrex::MultiFab& fire_phi,
                                      const amrex::MultiFab& fireline_intensity,
                                      amrex::Real I_B_crit,
                                      amrex::Real R_crown_ms,
                                      bool fire_debug = false)
{
    using namespace amrex;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fire_crown_active); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto ca  = fire_crown_active[mfi].array();
        auto ros = fire_ros[mfi].array();
        auto phi = fire_phi[mfi].const_array();
        auto I_B = fireline_intensity[mfi].const_array();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Only update burned cells
            if (phi(i, j, k) < 0.0) {
                // Check if crown fire should be initiated
                if (ca(i, j, k) < 0.5 && I_B(i, j, k) >= I_B_crit) {
                    ca(i, j, k) = 1.0;
                }
                // If crown fire is active, use crown ROS
                if (ca(i, j, k) >= 0.5) {
                    ros(i, j, k) = R_crown_ms;
                }
            }
        });
    }

    if (fire_debug) {
        amrex::Real n_crown = fire_crown_active.sum(0);
        amrex::Real ros_max = fire_ros.max(0);
        amrex::Print() << "[FIRE CROWN] active_crown_cells=" << static_cast<long>(n_crown)
                       << "  I_B_crit=" << I_B_crit << " kW/m"
                       << "  R_crown=" << R_crown_ms << " m/s"
                       << "  ros_max_after=" << ros_max << " m/s\n";
    }
}

/**
 * @brief Compute and accumulate crown fire heat flux
 *
 * GPU kernel. For each cell where crown fire is active (crown_active = 1.0) and
 * cell is burned (phi < 0):
 *   - On first activation, initialize crown load: w = C_BD * canopy_depth [kg/m²]
 *   - Compute crown fuel burnout residence time: tau = canopy_depth / max(R_crown, 1e-10)
 *   - Heat flux from crown: Q = crown_load * h_crown / tau
 *   - Accumulate crown flux into existing surface flux
 *   - Deplete crown load exponentially: w_new = w * exp(-dt / tau)
 *
 * @param[out] fire_heat_flux       Heat flux field [W/m²] (accumulated)
 * @param[out] fire_crown_load      Crown fuel load [kg/m²] (initialized on first call, then depleted)
 * @param[in]  fire_crown_active    Crown fire status (1.0 = active)
 * @param[in]  fire_phi             Level-set field (<0 = burned)
 * @param[in]  canopy_depth         Canopy depth [m]
 * @param[in]  C_BD                 Canopy bulk density [kg/m³]; used to initialize crown load
 * @param[in]  h_crown_J_per_kg     Crown heat content [J/kg]
 * @param[in]  R_crown_ms           Crown fire ROS [m/s]
 * @param[in]  dt_s                 Timestep [s]
 */
inline void fill_crown_heat_flux(amrex::MultiFab& fire_heat_flux,
                                  amrex::MultiFab& fire_crown_load,
                                  const amrex::MultiFab& fire_crown_active,
                                  const amrex::MultiFab& fire_phi,
                                  amrex::Real canopy_depth,
                                  amrex::Real C_BD,
                                  amrex::Real h_crown_J_per_kg,
                                  amrex::Real R_crown_ms,
                                  amrex::Real dt_s)
{
    using namespace amrex;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fire_heat_flux); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto Q       = fire_heat_flux[mfi].array();
        auto w_crown = fire_crown_load[mfi].array();
        auto ca      = fire_crown_active[mfi].const_array();
        auto phi     = fire_phi[mfi].const_array();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Only process cells where crown fire is active and burned
            if (ca(i, j, k) >= 0.5 && phi(i, j, k) < 0.0) {

                // Initialize crown load on first activation:
                // load [kg/m²] = bulk density [kg/m³] × canopy depth [m]
                if (w_crown(i, j, k) <= 0.0) {
                    w_crown(i, j, k) = C_BD * canopy_depth;
                }

                // Crown burnout residence time [s]
                amrex::Real R_safe    = amrex::max(R_crown_ms, 1.0e-10);
                amrex::Real tau_crown = amrex::max(canopy_depth / R_safe, 1.0);

                // Heat flux from crown fuel [W/m²]
                amrex::Real Q_crown = w_crown(i, j, k) * h_crown_J_per_kg / tau_crown;
                Q(i, j, k) += Q_crown;

                // Exponential depletion of crown load
                amrex::Real decay_factor = std::exp(-dt_s / tau_crown);
                w_crown(i, j, k) *= decay_factor;

                // Floor at zero
                w_crown(i, j, k) = amrex::max(0.0, w_crown(i, j, k));
            }
        });
    }
}

#endif  // ERF_CROWN_FIRE_H