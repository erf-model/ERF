/**
 * @file ERF_DustPM.H
 * @brief Particulate matter size fractionation for EPA NAAQS compliance.
 *
 * Classifies dust size bins into PM2.5 and PM10 fractions based on
 * aerodynamic diameter. The classification follows US EPA convention:
 *
 *   PM2.5: aerodynamic diameter <= 2.5 µm
 *   PM10:  aerodynamic diameter <= 10 µm (includes PM2.5)
 *
 * A bin with diameter d [m] contributes its full mass to the PM class
 * if d <= threshold. No partial fractional attribution is applied
 * within a bin at this phase.
 *
 * References:
 *   US EPA (2024). NAAQS for Particulate Matter.
 *     https://www.epa.gov/pm-pollution/
 *     national-ambient-air-quality-standards-naaqs-particulate-matter-pm
 *   40 CFR Part 50.
 *     https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-50
 *   Watson, J.G. et al. (1994). Atmos. Environ., 28, 2493.
 *     https://doi.org/10.1016/1352-2310(94)90400-6
 */

#ifndef ERF_DUST_PM_H
#define ERF_DUST_PM_H

#include <vector>
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>

namespace DustPMConst {
    constexpr amrex::Real PM25_THRESHOLD_M  = 2.5e-6; ///< 2.5 µm in meters
    constexpr amrex::Real PM10_THRESHOLD_M  = 10.0e-6; ///< 10 µm in meters
    constexpr amrex::Real UG_PER_KG         = 1.0e9;  ///< kg/m^3 -> µg/m^3
    constexpr amrex::Real PM25_24H_NAAQS    = 35.0;   ///< µg/m^3, 40 CFR 50.18
    constexpr amrex::Real PM10_24H_NAAQS    = 150.0;  ///< µg/m^3, 40 CFR 50.6
}

/**
 * @brief Returns true when bin diameter contributes to PM2.5.
 * @param d_m  Bin aerodynamic diameter [m].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_pm25(amrex::Real d_m) noexcept {
    return d_m <= DustPMConst::PM25_THRESHOLD_M;
}

/**
 * @brief Returns true when bin diameter contributes to PM10.
 * @param d_m  Bin aerodynamic diameter [m].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_pm10(amrex::Real d_m) noexcept {
    return d_m <= DustPMConst::PM10_THRESHOLD_M;
}

/**
 * @brief Accumulate PM2.5 and PM10 concentrations from per-bin dust_conc_sfc.
 *
 * Loops over active bins and sums mass into pm25 and pm10 MultiFabs [µg/m^3].
 * Both MultiFabs are zeroed on entry and filled by this function.
 *
 * @param[out] pm25         PM2.5 concentration [µg/m^3] on dust grid.
 * @param[out] pm10         PM10 concentration [µg/m^3] on dust grid.
 * @param[in]  dust_conc    Near-surface dust concentration [kg/m^3] on dust grid.
 *                          When transport_bins_separately=false this is a single
 *                          component representing the total dust (all bins).
 *                          When transport_bins_separately=true, comp b holds bin b.
 * @param[in]  bin_diameters  Vector of bin diameters [m].
 * @param[in]  n_active     Number of active bins (1 or n_size_bins).
 */
inline void compute_pm_concentrations(
    amrex::MultiFab&       pm25,
    amrex::MultiFab&       pm10,
    const amrex::MultiFab& dust_conc,
    const std::vector<amrex::Real>& bin_diameters,
    int n_active)
{
    pm25.setVal(0.0);
    pm10.setVal(0.0);

    for (int b = 0; b < n_active; ++b) {
        int d_idx = (b < (int)bin_diameters.size())
                  ? b : (int)bin_diameters.size() - 1;
        amrex::Real d_m = bin_diameters[d_idx];

        bool contrib_pm25 = is_pm25(d_m);
        bool contrib_pm10 = is_pm10(d_m);
        if (!contrib_pm25 && !contrib_pm10) continue;

        const amrex::Real scale = DustPMConst::UG_PER_KG;

        for (amrex::MFIter mfi(dust_conc, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto conc = dust_conc.const_array(mfi);
            auto p25  = pm25.array(mfi);
            auto p10  = pm10.array(mfi);
            const int src_comp = b;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real c = amrex::max(conc(i,j,k,src_comp), 0.0) * scale;
                if (contrib_pm25) p25(i,j,k) += c;
                if (contrib_pm10) p10(i,j,k) += c;
            });
        }
    }
}

/**
 * @brief Update running time average using exponential moving average.
 *
 * C_avg(t) = C_avg(t-dt) * (T - dt) / T  +  C_now * dt / T
 *
 * At t=0 or when C_avg has never been set, initialise to C_now.
 *
 * @param[in,out] C_avg   Running average MultiFab [µg/m^3].
 * @param[in]     C_now   Current concentration [µg/m^3].
 * @param[in]     dt      Timestep [s].
 * @param[in]     T_avg   Averaging period [s] (3600 or 86400).
 * @param[in]     nstep   Current step number (used for initialisation).
 */
inline void update_running_average(
    amrex::MultiFab& C_avg,
    const amrex::MultiFab& C_now,
    amrex::Real dt, amrex::Real T_avg, int nstep)
{
    amrex::Real alpha = dt / T_avg;
    amrex::Real beta  = 1.0 - alpha;
    if (beta < 0.0) beta = 0.0;

    if (nstep == 0) {
        amrex::MultiFab::Copy(C_avg, C_now, 0, 0, 1, 0);
        return;
    }
    for (amrex::MFIter mfi(C_avg, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto avg = C_avg.array(mfi);
        auto now = C_now.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            avg(i,j,k) = beta * avg(i,j,k) + alpha * now(i,j,k);
        });
    }
}

/**
 * @brief Compute exceedance flag: 1 where C_avg > threshold, 0 otherwise.
 *
 * @param[out] flag      Exceedance MultiFab [0/1].
 * @param[in]  C_avg     Running average [µg/m^3].
 * @param[in]  threshold NAAQS threshold [µg/m^3].
 */
inline void compute_exceedance_flag(
    amrex::MultiFab& flag,
    const amrex::MultiFab& C_avg,
    amrex::Real threshold)
{
    for (amrex::MFIter mfi(flag, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto fl  = flag.array(mfi);
        auto avg = C_avg.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            fl(i,j,k) = (avg(i,j,k) > threshold) ? 1.0 : 0.0;
        });
    }
}

#endif  // ERF_DUST_PM_H