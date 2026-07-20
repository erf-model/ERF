/**
 * @file ERF_DustSilica.H
 * @brief Respirable crystalline silica (RCS) concentration from PM10.
 *
 * RCS [mg/m³] = PM10 [µg/m³] * 1e-3 * silica_fraction
 * OSHA PEL (quartz): 0.05 mg/m³, 29 CFR 1910.1053.
 * NIOSH REL: 0.05 mg/m³ TWA.
 *
 * References:
 *   29 CFR 1910.1053: https://www.osha.gov/silica-crystalline
 *   Checkoway et al. (1997). Am J Ind Med 31:432.
 */

#ifdef ERF_USE_DUST

#ifndef ERF_DUST_SILICA_H
#define ERF_DUST_SILICA_H

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <string>

/**
 * @brief Compute respirable crystalline silica (RCS) concentration from PM10.
 *
 * RCS [mg/m³] = PM10 [µg/m³] * 1e-3 * silica_fraction
 *
 * @param[out] rcs_mf              RCS concentration MultiFab [mg/m³]
 * @param[in]  pm10_mf             PM10 concentration [µg/m³] on dust grid
 * @param[in]  silica_fraction     Silica mass fraction of dust [-]
 */
inline void compute_silica_concentration(
    amrex::MultiFab&       rcs_mf,
    const amrex::MultiFab& pm10_mf,
    amrex::Real            silica_fraction)
{
    using namespace amrex;
    for (MFIter mfi(rcs_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto rcs = rcs_mf.array(mfi);
        auto pm10 = pm10_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // PM10 is in µg/m³, convert to mg/m³: multiply by 1e-3
            // RCS = PM10 [µg/m³] * 1e-3 [mg/µg] * silica_fraction [-]
            rcs(i,j,k) = amrex::max(pm10(i,j,k), 0.0) * 1.0e-3 * silica_fraction;
        });
    }
}

/**
 * @brief Append silica RCS diagnostics to CSV file.
 *
 * Header format:
 *   step,time_s,rcs_max_mg_m3,rcs_mean_mg_m3,osha_pel_exceeded_cells
 *
 * Only IOProcessor writes; all ranks compute reductions.
 *
 * @param[in] step       Current step number
 * @param[in] time_s     Current time [s]
 * @param[in] rcs_mf     RCS concentration MultiFab [mg/m³]
 * @param[in] osha_pel   OSHA PEL threshold [mg/m³]
 * @param[in] filename   Output CSV filename
 */
inline void write_silica_stats(int step, amrex::Real time_s,
                               const amrex::MultiFab& rcs_mf,
                               amrex::Real osha_pel,
                               const std::string& filename)
{
    using namespace amrex;

    // Compute domain reductions on all ranks
    Real rcs_max = rcs_mf.max(0);
    Real rcs_sum = rcs_mf.sum(0);
    Real rcs_mean = rcs_sum / rcs_mf.boxArray().numPts();

    // Count cells exceeding OSHA PEL using MultiFab mask
    MultiFab exceed_mask(rcs_mf.boxArray(), rcs_mf.DistributionMap(), 1, 0);
    exceed_mask.setVal(0.0);
    for (MFIter mfi(exceed_mask, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto mask = exceed_mask.array(mfi);
        auto rcs = rcs_mf.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            mask(i,j,k) = (rcs(i,j,k) > osha_pel) ? 1.0 : 0.0;
        });
    }
    Real n_exceed = exceed_mask.sum(0);

    if (!ParallelDescriptor::IOProcessor()) return;

    bool exists = std::ifstream(filename).good();
    std::ofstream out(filename, std::ios::app);
    if (!out.is_open()) {
        Print() << "[DUST SILICA] WARNING: Cannot open file: " << filename << "\n";
        return;
    }

    if (!exists) {
        out << "# ERF-Dust silica RCS diagnostics\n"
            << "step,time_s,rcs_max_mg_m3,rcs_mean_mg_m3,osha_pel_exceeded_cells\n";
    }

    out << step << "," << time_s << "," << rcs_max << "," << rcs_mean << ","
        << (long)n_exceed << "\n";
}

#endif // ERF_DUST_SILICA_H
#endif // ERF_USE_DUST
