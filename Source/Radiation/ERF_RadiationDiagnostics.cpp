#include "ERF_RadiationDiagnostics.H"
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <iomanip>

/**
 * @file ERF_RadiationDiagnostics.cpp
 * @brief Implementation of RadiationDiagnostics CSV logger and debug output.
 */

RadiationDiagnostics::RadiationDiagnostics(int verbosity,
                                           const std::string& diag_file,
                                           int amr_level)
  : m_verbosity(verbosity), m_diag_file(diag_file), m_amr_level(amr_level)
{
  // Constructor: nothing special needed
  // write_header_if_needed() is called on first append()
}

RadiationDiagnostics::~RadiationDiagnostics()
{
  // Destructor: nothing special needed
  // File handles are closed by ofstream RAII
}

void RadiationDiagnostics::write_header_if_needed()
{
  if (m_header_written) {
    return;
  }

  if (!amrex::ParallelDescriptor::IOProcessor()) {
    m_header_written = true;
    return;
  }

  // Check if file already exists (append vs. new)
  std::ifstream infile(m_diag_file);
  bool file_exists = infile.good();
  infile.close();

  std::ofstream outfile(m_diag_file, std::ios::app);
  if (!outfile.good()) {
    amrex::Warning("RadiationDiagnostics: Could not open file " + m_diag_file);
    m_header_written = true;
    return;
  }

  // Only write header if file is new
  if (!file_exists) {
    outfile << "step,time,SW_surface,SW_TOA,F_up_surface,F_down_toa,heating_rate_max"
            << "\n";
  }
  outfile.close();
  m_header_written = true;
}

void RadiationDiagnostics::append(int step, amrex::Real time,
                                  amrex::Real SW_surface, amrex::Real SW_TOA,
                                  amrex::Real F_up_surface, amrex::Real F_down_toa,
                                  amrex::Real heating_rate_max)
{
  // Guard against duplicate writes for the same step
  if (step == m_last_write_step) {
    return;
  }
  m_last_write_step = step;

  write_header_if_needed();

  // Print debug output (if verbosity >= 1, and IOProcessor only)
  //
  // NOTE: The bracketed tag below is intentionally NOT hardcoded to a
  // specific phase (e.g., "[Phase1]"). This diagnostics module is shared
  // by every phase of the radiation development (Phase 1 through the
  // current Phase 4 scattering work and beyond); using a static phase
  // label here was a latent bug — it silently kept printing "[Phase1]"
  // even after Phase 2/3/4 functionality was added, misleading anyone
  // grepping logs by phase. Use the module-generic "[RAD]" tag plus the
  // call-site tag, which remains accurate and grepable across all phases.
  if (m_verbosity >= 1 && amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "[RAD][RadiationDiagnostics::append] step=" << step
                   << " time=" << time << " SW_surface=" << SW_surface
                   << " SW_TOA=" << SW_TOA << " F_up_surface=" << F_up_surface
                   << " F_down_toa=" << F_down_toa
                   << " heating_rate_max=" << heating_rate_max << "\n";
  }

  // Print RADIATION_DIAG line (always, for RegTest scripts)
  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "RADIATION_DIAG: step=" << step << " time=" << time
                   << " SW_surface=" << std::scientific << std::setprecision(6)
                   << SW_surface << " SW_TOA=" << SW_TOA
                   << " F_up_surface=" << F_up_surface << " F_down_toa=" << F_down_toa
                   << " heating_rate_max=" << heating_rate_max << "\n";
  }

  // Append CSV row
  if (!amrex::ParallelDescriptor::IOProcessor()) {
    return;
  }

  std::ofstream outfile(m_diag_file, std::ios::app);
  if (!outfile.good()) {
    amrex::Warning("RadiationDiagnostics: Could not open file " + m_diag_file +
                   " for append");
    return;
  }

  // Write CSV row with scientific notation for fluxes
  outfile << step << "," << std::scientific << std::setprecision(6) << time << ","
          << SW_surface << "," << SW_TOA << "," << F_up_surface << "," << F_down_toa
          << "," << heating_rate_max << "\n";
  outfile.close();
}
