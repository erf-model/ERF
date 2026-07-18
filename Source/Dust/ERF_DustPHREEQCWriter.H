/**
 * @file ERF_DustPHREEQCWriter.H
 * @brief Writes dust deposition accumulator to PHREEQC-compatible input files.
 *
 * Writes the dust_deposition_rate MultiFab [kg/m^2] to a whitespace-delimited
 * text file at user-specified intervals. Format:
 *   x_centre[m]  y_centre[m]  deposition[kg/m^2]
 * one row per dust grid cell, sorted by (j, i) index order.
 *
 * Also writes a per-site summary CSV appending one block per write:
 *   time_s, site_index, site_name, total_deposition_kg
 *
 * The deposition file is overwritten each interval (not appended).
 * The summary CSV is appended.
 *
 * Fire reference (periodic write pattern):
 *   ERF_FirePrerequisites.H/cpp at commit f264ce7dd0f5e727c71c939659a6a5caf001d21b
 *   https://github.com/hgopalan/ERF/blob/f264ce7dd0f5e727c71c939659a6a5caf001d21b/
 *   Source/Fire/ERF_FirePrerequisites.H
 *
 * Fire stats CSV pattern:
 *   ERF_FireStatsOutput.H at same commit.
 *   https://github.com/hgopalan/ERF/blob/f264ce7dd0f5e727c71c939659a6a5caf001d21b/
 *   Source/Fire/ERF_FireStatsOutput.H
 *
 * References:
 *   Parkhurst & Appelo (2013). PHREEQC Version 3.
 *     https://pubs.usgs.gov/tm/06/a43/
 *   Marticorena & Bergametti (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Shao & Lu (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 */
#ifdef ERF_USE_DUST
#ifndef ERF_DUST_PHREEQC_WRITER_H_
#define ERF_DUST_PHREEQC_WRITER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>
#include <fstream>
#include <string>
#include <iomanip>

/**
 * @brief Write dust_deposition_rate to PHREEQC input file.
 *
 * File format (whitespace-delimited, one row per cell):
 *   x_centre[m]  y_centre[m]  deposition[kg/m^2]
 * File is overwritten (not appended). Only IOProcessor writes.
 * Rows ordered by i (x-fast), j (y-slow).
 *
 * @param filename       Output file path.
 * @param dep_rate       dust_deposition_rate MultiFab [kg/m^2].
 * @param geom_dust      Dust grid geometry.
 * @param cur_time       Simulation time [s] (written as header comment).
 * @param nstep          Current step (written as header comment).
 */
inline void write_phreeqc_deposition_file(
    const std::string&      filename,
    const amrex::MultiFab&  dep_rate,
    const amrex::Geometry&  geom_dust,
    amrex::Real             cur_time,
    int                     nstep)
{
    // Collect cell data to IOProcessor via ReduceReal.
    // For each cell, gather (x_c, y_c, dep) on host.
    // Pattern: ERF_FirePrerequisites.cpp LoopOnCpu + IOProcessor write.
    const auto& dx  = geom_dust.CellSize();
    const auto& lo  = geom_dust.ProbLo();
    const amrex::Box& domain = geom_dust.Domain();
    const int ni = domain.length(0);
    const int nj = domain.length(1);

    // Collect deposition values into a host vector (index = j*ni + i).
    amrex::Vector<amrex::Real> dep_host(ni * nj, 0.0);

    for (amrex::MFIter mfi(dep_rate); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto dep_arr = dep_rate.const_array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
            int idx = j * ni + i;
            dep_host[idx] = dep_arr(i,j,k);
        });
    }
    // Reduce across MPI ranks (sum = correct since each cell owned by one rank).
    amrex::ParallelDescriptor::ReduceRealSum(dep_host.data(), ni * nj);

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        amrex::Print() << "[DUST] WARNING: Cannot open PHREEQC file: "
                       << filename << "\n";
        return;
    }
    // Header comments for PHREEQC batch script identification.
    out << "# ERF-Dust deposition feedback file\n"
        << "# time_s=" << std::fixed << std::setprecision(3) << cur_time
        << "  step=" << nstep << "\n"
        << "# x_centre_m  y_centre_m  deposition_kg_m2\n";
    out << std::scientific << std::setprecision(6);
    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            amrex::Real xc = lo[0] + (i + 0.5) * dx[0];
            amrex::Real yc = lo[1] + (j + 0.5) * dx[1];
            amrex::Real dep = dep_host[j * ni + i];
            out << xc << "  " << yc << "  " << dep << "\n";
        }
    }
}

/**
 * @brief Append per-site deposition summary to CSV.
 *
 * One row per site per call:
 *   time_s, site_index, site_name, total_deposition_kg
 * Site 0 = unassigned cells; sites 1..N = Phase 20 named sites.
 * total_deposition_kg = sum(dep_rate * cell_area) over site cells.
 *
 * Pattern: ERF_FireStatsOutput.H append_fire_stats.
 *
 * @param filename       CSV file path.
 * @param dep_rate       dust_deposition_rate [kg/m^2] on dust grid.
 * @param site_id        dust_site_id MultiFab from Phase 20 (may be null).
 * @param geom_dust      Dust grid geometry (for cell area).
 * @param cur_time       Simulation time [s].
 * @param site_names     Phase 20 site name vector (may be empty).
 */
inline void append_phreeqc_site_summary(
    const std::string&                  filename,
    const amrex::MultiFab&              dep_rate,
    const amrex::MultiFab*              site_id,
    const amrex::Geometry&              geom_dust,
    amrex::Real                         cur_time,
    const amrex::Vector<std::string>&   site_names)
{
    const auto& dx = geom_dust.CellSize();
    const amrex::Real cell_area = dx[0] * dx[1];
    int n_sites = (int)site_names.size();

    // Sum deposition per site (site 0 = unassigned).
    amrex::Vector<amrex::Real> site_dep(n_sites + 1, 0.0);

    for (amrex::MFIter mfi(dep_rate); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto dep_arr = dep_rate.const_array(mfi);
        if (site_id) {
            auto sid_arr = site_id->const_array(mfi);
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
                int s = amrex::max(0, amrex::min(n_sites, (int)sid_arr(i,j,k)));
                site_dep[s] += dep_arr(i,j,k) * cell_area;
            });
        } else {
            // No site registry: accumulate all into site 0.
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
                site_dep[0] += dep_arr(i,j,k) * cell_area;
            });
        }
    }
    for (int s = 0; s <= n_sites; ++s)
        amrex::ParallelDescriptor::ReduceRealSum(site_dep[s]);

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    bool file_exists = std::ifstream(filename).good();
    std::ofstream out(filename, std::ios::app);
    if (!out.is_open()) {
        amrex::Print() << "[DUST] WARNING: Cannot open site summary file: "
                       << filename << "\n";
        return;
    }
    if (!file_exists) {
        out << "# PHREEQC per-site deposition summary\n"
            << "# References: Parkhurst & Appelo (2013)."
            << " https://pubs.usgs.gov/tm/06/a43/\n"
            << "time_s,site_index,site_name,total_deposition_kg\n";
    }
    for (int s = 0; s <= n_sites; ++s) {
        std::string sname = (s == 0)
            ? "unassigned"
            : (s <= n_sites ? site_names[s-1] : "unknown");
        out << std::fixed << std::setprecision(3) << cur_time << ","
            << s << "," << sname << ","
            << std::scientific << std::setprecision(6) << site_dep[s] << "\n";
    }
}
#endif
#endif
