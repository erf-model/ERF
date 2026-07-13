/**
 * @file ERF_DustThreshold.H
 * @brief Threshold friction velocity computation for the ERF-Dust emission layer.
 *
 * Provides inline GPU-device functions and a MultiFab kernel for computing
 * the threshold friction velocity u*_t at each dust grid cell.
 *
 * Emission begins when u* (from MRF PBL, Phase 9) exceeds u*_t.
 *
 * Base threshold (Bagnold 1941):
 *   u*_t = A * sqrt(rho_p * g * d / rho_a)
 *
 * Modifier chain:
 *   f_chem  = (1 - alpha_c * crust_index) * (1 - alpha_e * efflorescence)
 *   f_moist = 1 + beta_m  * surface_moisture   (raises u*_t; inhibits emission)
 *   f_supp  = 1 + gamma_s * suppression         (raises u*_t; inhibits emission)
 *
 *   u*_t_final = u*_t_base * f_chem / (f_moist * f_supp)
 *   clamped to [USTAR_T_MIN, USTAR_T_MAX].
 *
 * References:
 *   Bagnold, R.A. (1941). The Physics of Blown Sand and Desert Dunes.
 *     Methuen, London.
 *   Marticorena, B., & Bergametti, G. (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Shao, Y., & Lu, H. (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 *   Fecan, F., Marticorena, B., & Bergametti, G. (1999). Ann. Geophys., 17, 149.
 *     https://doi.org/10.1007/s00585-999-0149-7
 */

#ifndef ERF_DUST_THRESHOLD_H
#define ERF_DUST_THRESHOLD_H

#ifdef ERF_USE_DUST

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelFor.H>
#include <AMReX_MFIter.H>
#include <cmath>

namespace DustThresholdConst {
    constexpr amrex::Real USTAR_T_MIN  = 0.05;   ///< Floor for u*_t [m/s].
    constexpr amrex::Real USTAR_T_MAX  = 5.0;    ///< Ceiling for u*_t [m/s].
    constexpr amrex::Real RHO_AIR_STD  = 1.225;  ///< Standard air density [kg/m^3].
    constexpr amrex::Real G_ACCEL      = 9.81;   ///< Gravitational acceleration [m/s^2].
    constexpr amrex::Real BETA_MOIST   = 4.0;    ///< Moisture inhibition coefficient [-].
                                                  ///<   Fecan et al. (1999).
    constexpr amrex::Real GAMMA_SUPP   = 6.0;    ///< Suppression elevation coefficient [-].
}

/**
 * @brief Compute Bagnold base threshold friction velocity.
 *
 * u*_t = A * sqrt(rho_p * g * d / rho_a)
 *
 * @param[in] A        Bagnold coefficient (typically 0.0123) [-].
 * @param[in] rho_p    Particle density [kg/m^3].
 * @param[in] d_m      Particle diameter [m].
 * @param[in] rho_a    Air density [kg/m^3]. If <= 0, RHO_AIR_STD is used.
 *
 * @return Bagnold u*_t [m/s].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_ustar_t_bagnold(amrex::Real A, amrex::Real rho_p,
                                     amrex::Real d_m, amrex::Real rho_a) noexcept
{
    amrex::Real rho = (rho_a > 0.0) ? rho_a : DustThresholdConst::RHO_AIR_STD;
    return A * std::sqrt(rho_p * DustThresholdConst::G_ACCEL * d_m / rho);
}

/**
 * @brief Compute moisture inhibition factor f_moist.
 *
 * f_moist = 1 + beta_m * w; w clamped to [0,1].
 * Reference: Fecan et al. (1999), https://doi.org/10.1007/s00585-999-0149-7
 *
 * @param[in] w Surface moisture [0,1]. Values outside this range are clamped.
 *
 * @return f_moist, always >= 1.0.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_moisture_inhibition(amrex::Real w) noexcept
{
    w = amrex::min(amrex::max(w, 0.0), 1.0);
    return 1.0 + DustThresholdConst::BETA_MOIST * w;
}

/**
 * @brief Compute suppression elevation factor f_supp.
 *
 * f_supp = 1 + gamma_s * s; s clamped to [0,1].
 *
 * @param[in] s Suppression agent coverage [0,1]. Values outside this range are clamped.
 *
 * @return f_supp, always >= 1.0.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_suppression_factor(amrex::Real s) noexcept
{
    s = amrex::min(amrex::max(s, 0.0), 1.0);
    return 1.0 + DustThresholdConst::GAMMA_SUPP * s;
}

/**
 * @brief Compute final threshold friction velocity with all modifiers.
 *
 * Combines Bagnold base with chemistry (crust/efflorescence), moisture, and suppression:
 *
 *   f_chem  = (1 - alpha_crust * crust_index) * (1 - alpha_efflor * efflorescence)
 *   f_moist = 1 + beta_m * surface_moisture
 *   f_supp  = 1 + gamma_s * suppression
 *
 *   u*_t = u*_t_base * f_chem / (f_moist * f_supp)
 *   clamped to [USTAR_T_MIN, USTAR_T_MAX].
 *
 * @param[in] ustar_t_base       Bagnold base threshold [m/s].
 * @param[in] crust_index        Crust strength [0,1].
 * @param[in] efflorescence      Efflorescence coverage [0,1].
 * @param[in] alpha_crust        Crust reduction coefficient [-].
 * @param[in] alpha_efflor       Efflorescence reduction coefficient [-].
 * @param[in] surface_moisture   Surface moisture [0,1].
 * @param[in] suppression        Suppression coverage [0,1].
 *
 * @return Final u*_t clamped to [USTAR_T_MIN, USTAR_T_MAX] [m/s].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_ustar_t_full(amrex::Real ustar_t_base,
                                  amrex::Real crust_index,
                                  amrex::Real efflorescence,
                                  amrex::Real alpha_crust,
                                  amrex::Real alpha_efflor,
                                  amrex::Real surface_moisture,
                                  amrex::Real suppression) noexcept
{
    // Chemistry reduction from PHREEQC (Phase 4).
    // Marticorena & Bergametti (1995), https://doi.org/10.1029/95JD00690
    amrex::Real f_chem = (1.0 - alpha_crust * amrex::min(amrex::max(crust_index,  0.0), 1.0))
                       * (1.0 - alpha_efflor * amrex::min(amrex::max(efflorescence,0.0), 1.0));
    f_chem = amrex::max(f_chem, 0.01);

    amrex::Real f_moist = compute_moisture_inhibition(surface_moisture);
    amrex::Real f_supp  = compute_suppression_factor(suppression);

    amrex::Real ut = ustar_t_base * f_chem / (f_moist * f_supp);
    return amrex::min(amrex::max(ut, DustThresholdConst::USTAR_T_MIN),
                      DustThresholdConst::USTAR_T_MAX);
}

/**
 * @brief Recompute dust_ustar_t at all dust grid cells.
 *
 * GPU ParallelFor over dust grid tiles. Applies the full modifier chain:
 * crust/efflorescence (PHREEQC, Phase 4), moisture inhibition, suppression.
 * Called from DustLayer::initialize and DustLayer::advance.
 *
 * @param[in,out] dust_ustar_t      Threshold friction velocity [m/s].
 * @param[in]     dust_ustar_base   Base Bagnold threshold (read-only) [m/s].
 * @param[in]     dust_crust_index  Crust strength index [0,1].
 * @param[in]     dust_efflor       Efflorescence fraction [0,1].
 * @param[in]     dust_surf_moist   Surface moisture [0,1].
 * @param[in]     dust_suppression  Suppression coverage [0,1].
 * @param[in]     alpha_crust       Crust reduction coefficient [-].
 * @param[in]     alpha_efflor      Efflorescence reduction coefficient [-].
 */
inline void recompute_dust_ustar_t(
    amrex::MultiFab&       dust_ustar_t,
    const amrex::MultiFab& dust_ustar_base,
    const amrex::MultiFab& dust_crust_index,
    const amrex::MultiFab& dust_efflor,
    const amrex::MultiFab& dust_surf_moist,
    const amrex::MultiFab& dust_suppression,
    amrex::Real            alpha_crust,
    amrex::Real            alpha_efflor)
{
    for (amrex::MFIter mfi(dust_ustar_t, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        auto ut_arr     = dust_ustar_t.array(mfi);
        auto base_arr   = dust_ustar_base.const_array(mfi);
        auto crust_arr  = dust_crust_index.const_array(mfi);
        auto efflor_arr = dust_efflor.const_array(mfi);
        auto moist_arr  = dust_surf_moist.const_array(mfi);
        auto supp_arr   = dust_suppression.const_array(mfi);
        const amrex::Real ac = alpha_crust;
        const amrex::Real ae = alpha_efflor;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ut_arr(i,j,k) = compute_ustar_t_full(
                base_arr(i,j,k), crust_arr(i,j,k), efflor_arr(i,j,k),
                ac, ae, moist_arr(i,j,k), supp_arr(i,j,k));
        });
    }
}

#endif // ERF_USE_DUST

#endif // ERF_DUST_THRESHOLD_H
