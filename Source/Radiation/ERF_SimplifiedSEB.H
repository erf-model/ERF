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

#endif // ERF_SIMPLIFIED_SEB_H_
