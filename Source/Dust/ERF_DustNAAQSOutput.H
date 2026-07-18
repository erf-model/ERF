/**
 * @file ERF_DustNAAQSOutput.H
 * @brief CSV output for EPA NAAQS PM2.5/PM10 compliance diagnostics.
 *
 * Pattern: ERF_DustStatsOutput.H (Phase 16) and ERF_FireStatsOutput.H:
 *   https://github.com/hgopalan/ERF/blob/ERF-Fire/Source/Fire/ERF_FireStatsOutput.H
 *
 * Writes per-step domain statistics to dust_naaqs.csv:
 *   step, time_s, PM25_max_ug_m3, PM25_24h_max_ug_m3, PM10_max_ug_m3,
 *   PM10_24h_max_ug_m3, n_cells_PM25_exceed, n_cells_PM10_exceed
 *
 * References:
 *   US EPA NAAQS PM2.5 (9 µg/m^3 annual, 35 µg/m^3 24-hour):
 *     https://www.epa.gov/pm-pollution/
 *     national-ambient-air-quality-standards-naaqs-particulate-matter-pm
 *   40 CFR Part 50: https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-50
 */

#ifndef ERF_DUST_NAAQS_OUTPUT_H
#define ERF_DUST_NAAQS_OUTPUT_H

#ifdef ERF_USE_DUST
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <string>

inline void append_naaqs_stats(
    int step, amrex::Real time_s,
    const std::string& filename,
    const amrex::MultiFab& pm25,
    const amrex::MultiFab& pm25_24h,
    const amrex::MultiFab& pm10,
    const amrex::MultiFab& pm10_24h,
    const amrex::MultiFab& flag_pm25,
    const amrex::MultiFab& flag_pm10)
{
    amrex::Real pm25_max    = pm25.max(0);
    amrex::Real pm25_24h_mx = pm25_24h.max(0);
    amrex::Real pm10_max    = pm10.max(0);
    amrex::Real pm10_24h_mx = pm10_24h.max(0);
    amrex::Real n_pm25_ex   = flag_pm25.sum(0);
    amrex::Real n_pm10_ex   = flag_pm10.sum(0);

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    bool file_exists = std::ifstream(filename).good();
    std::ofstream out(filename, std::ios::app);
    if (!out.is_open()) {
        amrex::Print() << "[DUST] WARNING: Cannot open NAAQS file: "
                       << filename << "\n";
        return;
    }
    if (!file_exists) {
        out << "# EPA NAAQS PM2.5/PM10 compliance output\n"
            << "# PM2.5 24h NAAQS: 35 ug/m^3 (40 CFR 50.18)\n"
            << "# PM10  24h NAAQS: 150 ug/m^3 (40 CFR 50.6)\n"
            << "step,time_s,PM25_max_ug_m3,PM25_24h_max_ug_m3,"
            << "PM10_max_ug_m3,PM10_24h_max_ug_m3,"
            << "n_cells_PM25_exceed,n_cells_PM10_exceed\n";
    }
    out << step         << "," << time_s      << ","
        << pm25_max     << "," << pm25_24h_mx << ","
        << pm10_max     << "," << pm10_24h_mx << ","
        << (long)n_pm25_ex << "," << (long)n_pm10_ex << "\n";
}
#endif

#endif  // ERF_DUST_NAAQS_OUTPUT_H
