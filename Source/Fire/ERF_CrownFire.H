#ifndef ERF_CROWN_FIRE_H
#define ERF_CROWN_FIRE_H

#include <AMReX_MultiFab.H>
#include <cmath>

/**
 * @file ERF_CrownFire.H
 * @brief Crown fire initiation and active crown fire rate of spread
 *
 * Implements two models:
 *   1. Van Wagner (1977) criterion for surface-to-crown fire transition
 *   2. Cruz et al. (2005) model for active crown fire rate of spread
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
 */

/**
 * @brief Van Wagner (1977) critical fireline intensity for crown fire initiation
 *
 * Formula: I_B_crit = 0.010 * (C_BH * (M_fold - M_c))^1.5
 *
 * @param[in] C_BH      Canopy base height [m]
 * @param[in] M_fold    Foliar moisture [fraction]
 * @param[in] M_c       Critical foliar moisture [fraction]
 * @return Critical intensity [kW/m], or 1.0e12 kW/m when M_fold <= M_c
 *
 * Reference: Van Wagner, C.E. (1977). Conditions for the start and spread of
 *   crown fire. Canadian Journal of Forest Research, 7(1), 23–34.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real van_wagner_critical_intensity(amrex::Real C_BH, amrex::Real M_fold,
                                           amrex::Real M_c) noexcept
{
    // Crown fire thermally impossible when foliar moisture is at or below critical
    if (M_fold <= M_c) {
        return 1.0e12_rt;
    }
    amrex::Real delta_m = M_fold - M_c;
    amrex::Real I_B_crit = 0.010_rt * std::pow(C_BH * delta_m, 1.5_rt);
    return I_B_crit;
}

/**
 * @brief Cruz et al. (2005) active crown fire rate of spread
 *
 * Formula: R [m/min] = 11.02 * U_10^0.9 * C_BD^0.19 * exp(-17.08 * M_fold)
 * Returns value in m/s (input formula is m/min, so divide by 60).
 *
 * @param[in] U_10_kmh  10-meter wind speed [km/h]
 * @param[in] C_BD      Canopy bulk density [kg/m³]
 * @param[in] M_fold    Foliar moisture [fraction]
 * @return Crown fire rate of spread [m/s], or 0 when U_10_kmh <= 0
 *
 * Reference: Cruz, M.G., Alexander, M.E. & Wakimoto, R.H. (2005). Development
 *   and testing of models for predicting crown fire rate of spread in conifer
 *   forest stands. Canadian Journal of Forest Research, 35(7), 1626–1639.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cruz_crown_ros(amrex::Real U_10_kmh, amrex::Real C_BD, amrex::Real M_fold) noexcept
{
    // Zero wind produces zero crown ROS
    if (U_10_kmh <= 0.0_rt) {
        return 0.0_rt;
    }
    amrex::Real R_m_min = 11.02_rt * std::pow(U_10_kmh, 0.9_rt)
                        * std::pow(C_BD, 0.19_rt)
                        * std::exp(-17.08_rt * M_fold);
    // Convert m/min to m/s
    amrex::Real R_m_s = R_m_min / 60.0_rt;
    return amrex::max(0.0_rt, R_m_s);
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
 * @param[out] fire_ros            Rate of spread field [m/s]
 * @param[in]  fire_phi            Level-set field (<0 = burned)
 * @param[in]  fireline_intensity  Byram intensity [kW/m]
 * @param[in]  I_B_crit            Critical intensity for initiation [kW/m]
 * @param[in]  R_crown_ms          Crown fire ROS [m/s]
 */
inline void update_crown_fire_status(amrex::MultiFab& fire_crown_active,
                                      amrex::MultiFab& fire_ros,
                                      const amrex::MultiFab& fire_phi,
                                      const amrex::MultiFab& fireline_intensity,
                                      amrex::Real I_B_crit,
                                      amrex::Real R_crown_ms)
{
    using namespace amrex;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fire_crown_active); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto& ca = fire_crown_active[mfi];
        auto& ros = fire_ros[mfi];
        const auto& phi = fire_phi[mfi];
        const auto& I_B = fireline_intensity[mfi];

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Only update burned cells
            if (phi(i, j, k) < 0.0_rt) {
                // Check if crown fire should be initiated
                if (ca(i, j, k) < 0.5_rt && I_B(i, j, k) >= I_B_crit) {
                    ca(i, j, k) = 1.0_rt;
                }
                // If crown fire is active, use crown ROS
                if (ca(i, j, k) >= 0.5_rt) {
                    ros(i, j, k) = R_crown_ms;
                }
            }
        });
    }
}

/**
 * @brief Compute and accumulate crown fire heat flux
 *
 * GPU kernel. For each cell where crown fire is active (crown_active = 1.0) and
 * cell is burned (phi < 0):
 *   - Compute crown fuel burnout residence time: tau = canopy_depth / max(R_crown, 1e-10)
 *   - Heat flux from crown: Q = crown_load * h_crown / tau
 *   - Accumulate crown flux into existing surface flux
 *   - Deplete crown load exponentially: w_new = w * exp(-dt / tau)
 *
 * @param[out] fire_heat_flux       Heat flux field [W/m²] (accumulated)
 * @param[out] fire_crown_load      Crown fuel load [kg/m²] (depleted)
 * @param[in]  fire_crown_active    Crown fire status (1.0 = active)
 * @param[in]  fire_phi             Level-set field (<0 = burned)
 * @param[in]  canopy_depth         Canopy depth [m]
 * @param[in]  C_BD                 Canopy bulk density [kg/m³]
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
        auto& Q = fire_heat_flux[mfi];
        auto& w_crown = fire_crown_load[mfi];
        const auto& ca = fire_crown_active[mfi];
        const auto& phi = fire_phi[mfi];

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Only process cells where crown fire is active and burned
            if (ca(i, j, k) >= 0.5_rt && phi(i, j, k) < 0.0_rt) {
                // Crown burnout residence time [s]
                amrex::Real R_safe = amrex::max(R_crown_ms, 1.0e-10_rt);
                amrex::Real tau_crown = amrex::max(canopy_depth / R_safe, 1.0_rt);

                // Heat flux from crown fuel [W/m²]
                amrex::Real Q_crown = w_crown(i, j, k) * h_crown_J_per_kg / tau_crown;
                Q(i, j, k) += Q_crown;

                // Exponential depletion of crown load
                amrex::Real decay_factor = std::exp(-dt_s / tau_crown);
                w_crown(i, j, k) *= decay_factor;

                // Floor at zero
                w_crown(i, j, k) = amrex::max(0.0_rt, w_crown(i, j, k));
            }
        });
    }
}

#endif
