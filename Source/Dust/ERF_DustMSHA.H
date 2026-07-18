/**
 * @file ERF_DustMSHA.H
 * @brief MSHA worker dust exposure: 8h TWA dose accumulator and exceedance flag.
 *
 * TWA model: dose += PM_resp * dt/3600, TWA_8h = dose/8
 * PM_resp = dust_pm10[ug/m3] * 1e-3  (converts to mg/m3)
 * MSHA PEL: 5 mg/m3, 30 CFR 56.5001.
 * Shift reset: dose zeroed when floor(t/T_shift) increases.
 *
 * References:
 *   30 CFR Part 56: https://www.ecfr.gov/current/title-30/chapter-I/subchapter-K/part-56
 *   30 CFR Part 57: https://www.ecfr.gov/current/title-30/chapter-I/subchapter-K/part-57
 *   MSHA: https://www.msha.gov/rules-regs/dust
 *   Stayner et al. (1998). Am J Respir Crit Care Med 157:69.
 *     https://doi.org/10.1164/ajrccm.157.1.9703022
 *   Attfield & Morring (1992). Am J Ind Med 22:417.
 *     https://doi.org/10.1002/ajim.4700220311
 */

#ifdef ERF_USE_DUST

#ifndef ERF_DUST_MSHA_H
#define ERF_DUST_MSHA_H

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>

namespace DustMSHAConst {
    constexpr amrex::Real UG_TO_MG  = 1.0e-3;     ///< Convert µg/m³ to mg/m³
    constexpr amrex::Real S_TO_H    = 1.0/3600.0; ///< Convert seconds to hours
    constexpr amrex::Real TWA_HOURS = 8.0;        ///< 8-hour time-weighted average
    constexpr amrex::Real PEL_DEF   = 5.0;        ///< Default PEL [mg/m³], 30 CFR 56.5001
}

/**
 * @brief Update MSHA dose accumulator and compute 8h TWA.
 *
 * Computes:
 *   dose(i,j) += max(pm10(i,j), 0) [ug/m3] * 1e-3 * dt/3600 [mg/m3/h]
 *   twa(i,j)   = dose(i,j) / 8.0 [mg/m3]
 *
 * Called each timestep from compute_msha_exposure.
 *
 * @param[in,out] dose  Cumulative dose MultiFab [mg/m³/h]
 * @param[in,out] twa   8h TWA MultiFab [mg/m³]
 * @param[in]     pm10  Instantaneous PM10 [µg/m³]
 * @param[in]     dt    Timestep [s]
 */
AMREX_FORCE_INLINE
void update_msha_dose(amrex::MultiFab& dose, amrex::MultiFab& twa,
                      const amrex::MultiFab& pm10, amrex::Real dt)
{
    using namespace amrex;
    const Real dth = dt * DustMSHAConst::S_TO_H;
    for (MFIter mfi(dose, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto d = dose.array(mfi);
        auto tw = twa.array(mfi);
        auto pm = pm10.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            d(i,j,k)  += amrex::max(pm(i,j,k), 0.0) * DustMSHAConst::UG_TO_MG * dth;
            tw(i,j,k)  = d(i,j,k) / DustMSHAConst::TWA_HOURS;
        });
    }
}

/**
 * @brief Compute MSHA PEL exceedance flag.
 *
 * Sets exceedance(i,j) = 1 where twa(i,j) > pel, else 0.
 *
 * @param[in,out] exceed  Exceedance flag MultiFab [0/1]
 * @param[in]     twa     8h TWA [mg/m³]
 * @param[in]     pel     PEL threshold [mg/m³]
 */
AMREX_FORCE_INLINE
void compute_msha_exceed(amrex::MultiFab& exceed, const amrex::MultiFab& twa,
                         amrex::Real pel)
{
    using namespace amrex;
    for (MFIter mfi(exceed, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto e = exceed.array(mfi);
        auto tw = twa.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            e(i,j,k) = (tw(i,j,k) > pel) ? 1.0 : 0.0;
        });
    }
}

#endif // ERF_DUST_MSHA_H
#endif // ERF_USE_DUST
