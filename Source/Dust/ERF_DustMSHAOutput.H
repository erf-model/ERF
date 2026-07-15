/**
 * @file ERF_DustMSHAOutput.H
 * @brief CSV writers for MSHA exposure, shift summary, and receptor samples.
 *
 * Three output files:
 *   msha_exposure.csv     - per-step: step,time_s,TWA_max,n_exceed,dose_max
 *   msha_shift_summary.csv- per-shift-end: shift,time_s,TWA_peak,n_exceed
 *   msha_receptor_<N>.csv - per-step: step,time_s,PM10_ug_m3 at receptor N
 *
 * Pattern: ERF_FireStatsOutput.H in Source/Fire/:
 *   https://github.com/hgopalan/ERF/blob/ERF-Fire/Source/Fire/ERF_FireStatsOutput.H
 *
 * References:
 *   30 CFR Part 56: https://www.ecfr.gov/current/title-30/chapter-I/subchapter-K/part-56
 *   MSHA dust: https://www.msha.gov/rules-regs/dust
 */

#ifdef ERF_USE_DUST

#ifndef ERF_DUST_MSHA_OUTPUT_H
#define ERF_DUST_MSHA_OUTPUT_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <string>

/**
 * @brief Append one row to per-step MSHA exposure CSV.
 *
 * Header format:
 *   # MSHA 8h TWA. PEL=5 mg/m3, 30 CFR 56.5001
 *   step,time_s,TWA_max_mg_m3,n_cells_exceed,dose_max_mg_m3_h
 *
 * Only IOProcessor writes; all ranks compute reductions.
 *
 * @param[in] step      Current step number
 * @param[in] t         Current time [s]
 * @param[in] f         CSV filename
 * @param[in] twa       8h TWA MultiFab [mg/m³]
 * @param[in] ex        Exceedance flag MultiFab [0/1]
 * @param[in] dose      Cumulative dose MultiFab [mg/m³/h]
 */
inline void append_msha_stats(int step, amrex::Real t,
    const std::string& f, const amrex::MultiFab& twa,
    const amrex::MultiFab& ex, const amrex::MultiFab& dose)
{
    amrex::Real tmax = twa.max(0);
    amrex::Real nex  = ex.sum(0);
    amrex::Real dmax = dose.max(0);

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    bool exists = std::ifstream(f).good();
    std::ofstream out(f, std::ios::app);

    if (!exists) {
        out << "# MSHA 8h TWA. PEL=5 mg/m3, 30 CFR 56.5001\n"
            << "step,time_s,TWA_max_mg_m3,n_cells_exceed,dose_max_mg_m3_h\n";
    }

    out << step << "," << t << "," << tmax << "," << (long)nex << "," << dmax << "\n";
}

/**
 * @brief Write end-of-shift summary to CSV.
 *
 * Header format:
 *   # MSHA end-of-shift summary. 30 CFR Part 56
 *   shift_number,shift_end_time_s,TWA_peak_mg_m3,n_cells_exceed
 *
 * Only IOProcessor writes; all ranks compute reductions.
 *
 * @param[in] shift     Shift number (1-indexed)
 * @param[in] t         Shift end time [s]
 * @param[in] f         CSV filename
 * @param[in] twa       8h TWA MultiFab at shift end [mg/m³]
 * @param[in] ex        Exceedance flag MultiFab at shift end [0/1]
 */
inline void write_msha_shift_summary(int shift, amrex::Real t,
    const std::string& f, const amrex::MultiFab& twa, const amrex::MultiFab& ex)
{
    amrex::Real tp = twa.max(0);
    amrex::Real nex = ex.sum(0);

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    bool exists = std::ifstream(f).good();
    std::ofstream out(f, std::ios::app);

    if (!exists) {
        out << "# MSHA end-of-shift summary. 30 CFR Part 56\n"
            << "shift_number,shift_end_time_s,TWA_peak_mg_m3,n_cells_exceed\n";
    }

    out << shift << "," << t << "," << tp << "," << (long)nex << "\n";
}

/**
 * @brief Append receptor point sample to CSV.
 *
 * Finds the nearest dust grid cell to (rx, ry) and writes PM10 value.
 * No interpolation: uses nearest-cell-only.
 *
 * Header format (created on first call):
 *   # Receptor: <name> x=<rx> y=<ry>
 *   # 30 CFR Part 56. MSHA: https://www.msha.gov/rules-regs/dust
 *   step,time_s,PM10_ug_m3
 *
 * Only IOProcessor writes; all ranks participate in nearest-cell search and allreduce.
 *
 * @param[in] step      Current step number
 * @param[in] t         Current time [s]
 * @param[in] f         CSV filename (e.g., "msha_receptor_<name>.csv")
 * @param[in] name      Receptor name (for header comment)
 * @param[in] rx        Receptor x-coordinate [m]
 * @param[in] ry        Receptor y-coordinate [m]
 * @param[in] pm10      PM10 MultiFab [µg/m³]
 * @param[in] geom      Dust grid geometry
 */
inline void append_receptor_sample(int step, amrex::Real t,
    const std::string& f, const std::string& name,
    amrex::Real rx, amrex::Real ry,
    const amrex::MultiFab& pm10, const amrex::Geometry& geom)
{
    using namespace amrex;

    // Compute cell indices from coordinates
    const auto& dx = geom.CellSize();
    const auto& lo = geom.ProbLo();
    const Box& dom = geom.Domain();

    int ic = amrex::max(dom.smallEnd(0),
             amrex::min(dom.bigEnd(0), (int)((rx - lo[0]) / dx[0])));
    int jc = amrex::max(dom.smallEnd(1),
             amrex::min(dom.bigEnd(1), (int)((ry - lo[1]) / dx[1])));

    // Find value in the nearest cell
    Real val = 0.0;
    IntVect iv(ic, jc, 0);

    for (MFIter mfi(pm10); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(iv)) {
            val = pm10[mfi](iv, 0);
        }
    }

    // Allreduce to get value (only one rank will have it)
    ParallelDescriptor::ReduceRealSum(val);

    if (!ParallelDescriptor::IOProcessor()) return;

    bool exists = std::ifstream(f).good();
    std::ofstream out(f, std::ios::app);

    if (!exists) {
        out << "# Receptor: " << name << " x=" << rx << " y=" << ry << "\n"
            << "# 30 CFR Part 56. MSHA: https://www.msha.gov/rules-regs/dust\n"
            << "step,time_s,PM10_ug_m3\n";
    }

    out << step << "," << t << "," << val << "\n";
}

#endif // ERF_DUST_MSHA_OUTPUT_H
#endif // ERF_USE_DUST
