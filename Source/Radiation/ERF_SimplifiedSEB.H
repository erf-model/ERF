#ifndef ERF_SIMPLIFIED_SEB_H_
#define ERF_SIMPLIFIED_SEB_H_

#include <AMReX_GpuControl.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

/**
 * @file ERF_SimplifiedSEB.H
 * @brief Phase 18: Simplified SEB — Diagnostic mode residual computation.
 *
 * Implements diagnostic-only surface energy balance residual diagnosis.
 * No prognostic update to surface temperature or fluxes; residual is computed
 * and reported for validation and diagnostic purposes only.
 *
 * Surface Energy Balance Equation:
 *   R_net(i,j) = SW_net(i,j) + LW_net(i,j)
 *              = [sw_flux_sfc(i,j)] + [lw_flux_sfc(i,j)]
 *
 *   SEB_residual(i,j) = R_net(i,j) - hfx_sfc(i,j) - lh_sfc(i,j) - grdflx_sfc(i,j)
 *
 * Where:
 * - sw_flux_sfc: net shortwave flux at surface [W/m^2]
 * - lw_flux_sfc: net longwave flux at surface [W/m^2]
 * - hfx_sfc: sensible heat flux (H) [W/m^2]
 * - lh_sfc: latent heat flux (LE) [W/m^2]
 * - grdflx_sfc: ground heat flux (G) [W/m^2]
 *
 * A perfectly closed budget gives SEB_residual = 0; non-zero residual
 * indicates energy not accounted for by the four main SEB terms.
 *
 * References:
 * -----------
 * - Phase 18 Radiation Development notes (Source/Radiation/RAD_DEVELOPMENT.md)
 * - Oke, T. R., 1987: Boundary Layer Climates (2nd ed.), Routledge.
 */

/**
 * @brief Diagnose the surface energy balance residual.
 *
 * Computes the SEB residual from net radiative flux and turbulent/ground
 * heat fluxes. All inputs are guarded with std::isfinite() to safely handle
 * NaN/Inf; if any input is non-finite, returns 0.0 (safe no-op).
 *
 * **GPU Safety**:
 * - AMREX_GPU_DEVICE AMREX_FORCE_INLINE
 * - No host-side I/O
 * - No dynamic allocation
 * - Safe for use in device reduction kernels
 *
 * @param[in] sw_net Shortwave net flux at surface [W/m^2]
 * @param[in] lw_net Longwave net flux at surface [W/m^2]
 * @param[in] hfx Sensible heat flux [W/m^2]
 * @param[in] lh Latent heat flux [W/m^2]
 * @param[in] grdflx Ground heat flux [W/m^2]
 *
 * @return SEB residual [W/m^2] = (sw_net + lw_net) - hfx - lh - grdflx
 *         Returns 0.0 if any input is non-finite (NaN or Inf)
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_seb_residual(amrex::Real sw_net,
                                   amrex::Real lw_net,
                                   amrex::Real hfx,
                                   amrex::Real lh,
                                   amrex::Real grdflx)
{
    // Guard all inputs against non-finite values (NaN, Inf)
    if (!std::isfinite(sw_net) || !std::isfinite(lw_net) ||
        !std::isfinite(hfx) || !std::isfinite(lh) || !std::isfinite(grdflx)) {
        return 0.0;  // Safe no-op: return zero residual if any input is invalid
    }

    // Compute net radiation
    amrex::Real r_net = sw_net + lw_net;

    // Compute SEB residual
    amrex::Real seb_residual = r_net - hfx - lh - grdflx;

    return seb_residual;
}

/**
 * @brief Compute prognostic surface temperature tendency using force-restore formulation.
 *
 * Computes dT_s/dt from SEB residual and restoring term toward deep soil temperature:
 *   dT_s/dt = SEB_residual / C_s - (2*pi/tau) * (T_s - T_deep)
 *
 * All inputs are guarded with std::isfinite() to safely handle NaN/Inf; if any input
 * is non-finite, any parameter is non-positive (C_s <= 0, tau <= 0), returns 0.0 (safe no-op).
 *
 * **GPU Safety**:
 * - AMREX_GPU_DEVICE AMREX_FORCE_INLINE
 * - No host-side I/O
 * - No dynamic allocation
 * - Safe for use in device kernels
 *
 * @param[in] seb_residual SEB residual [W/m^2] = R_net - H - LE - G
 * @param[in] t_s Current surface temperature [K]
 * @param[in] t_deep Deep soil temperature [K]
 * @param[in] c_s Effective surface heat capacity [J/(m^2*K)]
 * @param[in] tau Force-restore timescale [s]
 *
 * @return Temperature tendency dT_s/dt [K/s]
 *         Returns 0.0 if any input is non-finite or any parameter is non-positive
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real prognostic_dTs_dt(amrex::Real seb_residual,
                               amrex::Real t_s,
                               amrex::Real t_deep,
                               amrex::Real c_s,
                               amrex::Real tau)
{
    // Guard all inputs against non-finite values
    if (!std::isfinite(seb_residual) || !std::isfinite(t_s) ||
        !std::isfinite(t_deep) || !std::isfinite(c_s) || !std::isfinite(tau)) {
        return 0.0;  // Safe no-op
    }

    // Validate parameter ranges
    if (c_s <= 0.0 || tau <= 0.0) {
        return 0.0;  // Safe no-op for invalid parameters
    }

    // Compute tendency: SEB_residual / C_s - (2*pi/tau) * (T_s - T_deep)
    amrex::Real seb_term = seb_residual / c_s;
    amrex::Real restore_term = (2.0 * M_PI / tau) * (t_s - t_deep);
    amrex::Real dt_s_dt = seb_term - restore_term;

    return dt_s_dt;
}

/**
 * @brief Compute prognostic surface moisture tendency using force-restore formulation.
 *
 * Computes dq_s/dt from latent heat flux and restoring term toward deep soil moisture:
 *   dq_s/dt = -(LE / (L_v * rho_w * d_s)) - (1/tau_q) * (q_s - q_deep)
 *
 * Constants:
 * - L_v = 2.5e6 J/kg (latent heat of vaporization)
 * - rho_w = 1000.0 kg/m^3 (water density)
 *
 * All inputs are guarded with std::isfinite() to safely handle NaN/Inf; if any input
 * is non-finite, any parameter is non-positive (d_s <= 0, tau_q <= 0), returns 0.0 (safe no-op).
 *
 * **GPU Safety**:
 * - AMREX_GPU_DEVICE AMREX_FORCE_INLINE
 * - No host-side I/O
 * - No dynamic allocation
 * - Safe for use in device kernels
 *
 * @param[in] le Latent heat flux [W/m^2]
 * @param[in] q_s Current surface moisture [kg/kg]
 * @param[in] q_deep Deep soil moisture [kg/kg]
 * @param[in] d_s Effective surface moisture layer depth [m]
 * @param[in] tau_q Moisture force-restore timescale [s]
 *
 * @return Moisture tendency dq_s/dt [kg/kg/s]
 *         Returns 0.0 if any input is non-finite or any parameter is non-positive
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real prognostic_dqs_dt(amrex::Real le,
                               amrex::Real q_s,
                               amrex::Real q_deep,
                               amrex::Real d_s,
                               amrex::Real tau_q)
{
    // Guard all inputs against non-finite values
    if (!std::isfinite(le) || !std::isfinite(q_s) ||
        !std::isfinite(q_deep) || !std::isfinite(d_s) || !std::isfinite(tau_q)) {
        return 0.0;  // Safe no-op
    }

    // Validate parameter ranges
    if (d_s <= 0.0 || tau_q <= 0.0) {
        return 0.0;  // Safe no-op for invalid parameters
    }

    // Physical constants
    const amrex::Real L_v = 2.5e6;    // Latent heat of vaporization [J/kg]
    const amrex::Real rho_w = 1000.0; // Water density [kg/m^3]

    // Compute tendency: -(LE / (L_v * rho_w * d_s)) - (1/tau_q) * (q_s - q_deep)
    amrex::Real le_term = -(le / (L_v * rho_w * d_s));
    amrex::Real restore_term = (1.0 / tau_q) * (q_s - q_deep);
    amrex::Real dq_s_dt = le_term - restore_term;

    return dq_s_dt;
}

#endif // ERF_SIMPLIFIED_SEB_H_
