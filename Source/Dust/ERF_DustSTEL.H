/**
 * @file ERF_DustSTEL.H
 * @brief 15-minute short-term exposure limit (STEL) for total dust.
 *
 * STEL_15min [mg/m³] = running average of PM10 [µg/m³]*1e-3 over 900 s.
 * Uses update_running_average() from ERF_DustPM.H with T_avg=900 s.
 * OSHA STEL nuisance dust: 10 mg/m³.
 *
 * References:
 *   OSHA 29 CFR 1910.1000 Table Z-1.
 *   ACGIH TLV-STEL guidelines.
 */

#ifdef ERF_USE_DUST

#ifndef ERF_DUST_STEL_H
#define ERF_DUST_STEL_H

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <ERF_DustPM.H>
#include <fstream>
#include <string>

/**
 * @brief Update 15-minute STEL running average from PM10 concentration.
 *
 * Uses update_running_average() from ERF_DustPM.H with T_avg=900 s.
 * Converts PM10 from µg/m³ to mg/m³ (factor of 1e-3).
 *
 * @param[in,out] stel_avg_mf   Running average STEL MultiFab [mg/m³]
 * @param[in]     pm10_mf       Instantaneous PM10 [µg/m³]
 * @param[in]     dt            Timestep [s]
 * @param[in]     T_avg         Averaging period [s], typically 900 s
 */
inline void update_stel_average(
    amrex::MultiFab&       stel_avg_mf,
    const amrex::MultiFab& pm10_mf,
    amrex::Real            dt,
    amrex::Real            T_avg)
{
    using namespace amrex;

    // Create temporary MultiFab for PM10 in mg/m³
    MultiFab pm10_mg_m3(pm10_mf.boxArray(), pm10_mf.DistributionMap(), 1, 0);

    // Convert PM10 from µg/m³ to mg/m³: multiply by 1e-3
    for (MFIter mfi(pm10_mg_m3, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto pm10_mg = pm10_mg_m3.array(mfi);
        auto pm10 = pm10_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pm10_mg(i,j,k) = pm10(i,j,k) * 1.0e-3;
        });
    }

    // Update STEL running average using exponential moving average
    for (MFIter mfi(stel_avg_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto stel = stel_avg_mf.array(mfi);
        auto pm10_mg = pm10_mg_m3.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Running average: stel = stel * (1 - dt/T_avg) + pm10_mg * (dt/T_avg)
            Real alpha = dt / T_avg;
            stel(i,j,k) = stel(i,j,k) * (1.0 - alpha) + pm10_mg(i,j,k) * alpha;
        });
    }
}

/**
 * @brief Append STEL diagnostics to CSV file.
 *
 * Header format:
 *   step,time_s,stel_max_mg_m3,stel_mean_mg_m3,stel_exceeded_cells
 *
 * Only IOProcessor writes; all ranks compute reductions.
 *
 * @param[in] step           Current step number
 * @param[in] time_s         Current time [s]
 * @param[in] stel_avg_mf    STEL running average MultiFab [mg/m³]
 * @param[in] stel_threshold STEL threshold [mg/m³]
 * @param[in] filename       Output CSV filename
 */
inline void write_stel_stats(int step, amrex::Real time_s,
                             const amrex::MultiFab& stel_avg_mf,
                             amrex::Real stel_threshold,
                             const std::string& filename)
{
    using namespace amrex;

    // Compute domain reductions on all ranks
    Real stel_max = stel_avg_mf.max(0);
    Real stel_sum = stel_avg_mf.sum(0);
    Real stel_mean = stel_sum / stel_avg_mf.boxArray().numPts();

    // Count cells exceeding STEL threshold using MultiFab mask
    MultiFab exceed_mask(stel_avg_mf.boxArray(), stel_avg_mf.DistributionMap(), 1, 0);
    exceed_mask.setVal(0.0);
    for (MFIter mfi(exceed_mask, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto mask = exceed_mask.array(mfi);
        auto stel = stel_avg_mf.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            mask(i,j,k) = (stel(i,j,k) > stel_threshold) ? 1.0 : 0.0;
        });
    }
    Real n_exceed = exceed_mask.sum(0);

    if (!ParallelDescriptor::IOProcessor()) return;

    bool exists = std::ifstream(filename).good();
    std::ofstream out(filename, std::ios::app);
    if (!out.is_open()) {
        Print() << "[DUST STEL] WARNING: Cannot open file: " << filename << "\n";
        return;
    }

    if (!exists) {
        out << "# ERF-Dust 15-minute STEL diagnostics\n"
            << "step,time_s,stel_max_mg_m3,stel_mean_mg_m3,stel_exceeded_cells\n";
    }

    out << step << "," << time_s << "," << stel_max << "," << stel_mean << ","
        << (long)n_exceed << "\n";
}

#endif // ERF_DUST_STEL_H
#endif // ERF_USE_DUST
