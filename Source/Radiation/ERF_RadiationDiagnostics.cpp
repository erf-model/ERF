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
                                           int amr_level,
                                           bool diag_enable,
                                           bool diag_stdout_enable,
                                           bool diag_tagged_enable,
                                           bool diag_regtest_line_enable,
                                           bool diag_csv_enable,
                                           const std::string& diag_callsite_mode,
                                           amrex::Real diag_dedup_tol)
  : m_verbosity(verbosity), m_diag_file(diag_file), m_amr_level(amr_level),
    m_diag_enable(diag_enable), m_diag_stdout_enable(diag_stdout_enable),
    m_diag_tagged_enable(diag_tagged_enable),
    m_diag_regtest_line_enable(diag_regtest_line_enable),
    m_diag_csv_enable(diag_csv_enable), m_diag_callsite_mode(diag_callsite_mode),
    m_diag_dedup_tol(diag_dedup_tol)
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
  if (m_header_written || !m_diag_csv_enable || !m_diag_enable) {
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
    // Base columns (Phase 1-7): always present
    outfile << "step,time,call_site,SW_surface,SW_TOA,F_up_surface,F_down_toa,heating_rate_max";
    // Phase 18: SEB diagnostic columns (added at end for backward compatibility)
    outfile << ",SEB_residual_mean,SEB_residual_max";
    outfile << "\n";
  }
  outfile.close();
  m_header_written = true;
}

void RadiationDiagnostics::append(int step, amrex::Real time, const std::string& call_site,
                                  amrex::Real SW_surface, amrex::Real SW_TOA,
                                  amrex::Real F_up_surface, amrex::Real F_down_toa,
                                  amrex::Real heating_rate_max,
                                  amrex::Real seb_residual_mean,
                                  amrex::Real seb_residual_max)
{
  // Phase 7: Master enable gate
  if (!m_diag_enable) {
    return;
  }

  // Phase 7: Call-site mode filtering
  bool should_emit_pre = (m_diag_callsite_mode == "both" || m_diag_callsite_mode == "pre_only");
  bool should_emit_post = (m_diag_callsite_mode == "both" || m_diag_callsite_mode == "post_only");
  
  bool is_pre_site = (call_site.find("pre") != std::string::npos);
  bool is_post_site = (call_site.find("post") != std::string::npos);

  // Check if this call site should be filtered based on mode
  if ((is_pre_site && !should_emit_pre) || (is_post_site && !should_emit_post)) {
    return;
  }

  // Phase 9: Guard against duplicate writes using 3-tuple identity:
  //   (step, call_site, time)
  // This ensures:
  // 1. Accidental repeated calls at same (step, call_site, time) are suppressed.
  // 2. Legitimate pre_dycore + post_dycore entries at same step both retained
  //    (different call_site values, so identity tuple differs).
  // 3. Mode filtering above already runs first, so pre_only/post_only modes
  //    naturally prevent unwanted entries from reaching this dedup logic.
  //
  // Time tolerance (m_diag_dedup_tol) accounts for floating-point rounding when
  // multiple diagnostics functions are called at effectively the same time.
  if (step == m_last_write_step &&
      call_site == m_last_write_call_site &&
      std::abs(time - m_last_write_time) < m_diag_dedup_tol) {
    return;  // Duplicate detected; skip write
  }
  m_last_write_step = step;
  m_last_write_call_site = call_site;
  m_last_write_time = time;

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
  if (m_verbosity >= 1 && amrex::ParallelDescriptor::IOProcessor() &&
      m_diag_tagged_enable && m_diag_stdout_enable) {
    amrex::Print() << "[RAD][RadiationDiagnostics::append] step=" << step
                   << " time=" << time 
                   << " call_site=" << call_site
                   << " SW_surface=" << SW_surface
                   << " SW_TOA=" << SW_TOA << " F_up_surface=" << F_up_surface
                   << " F_down_toa=" << F_down_toa
                   << " heating_rate_max=" << heating_rate_max << "\n";
  }

  // Print RADIATION_DIAG line (for RegTest scripts)
  if (amrex::ParallelDescriptor::IOProcessor() && m_diag_regtest_line_enable &&
      m_diag_stdout_enable) {
    amrex::Print() << "RADIATION_DIAG: step=" << step << " time=" << time
                   << " call_site=" << call_site
                   << " SW_surface=" << std::scientific << std::setprecision(6)
                   << SW_surface << " SW_TOA=" << SW_TOA
                   << " F_up_surface=" << F_up_surface << " F_down_toa=" << F_down_toa
                   << " heating_rate_max=" << heating_rate_max << "\n";
  }

  // Append CSV row
  if (!amrex::ParallelDescriptor::IOProcessor() || !m_diag_csv_enable) {
    return;
  }

  write_header_if_needed();

  std::ofstream outfile(m_diag_file, std::ios::app);
  if (!outfile.good()) {
    amrex::Warning("RadiationDiagnostics: Could not open file " + m_diag_file +
                   " for append");
    return;
  }

  // Write CSV row with scientific notation for fluxes
  outfile << step << "," << std::scientific << std::setprecision(6) << time << ","
          << call_site << ","
          << SW_surface << "," << SW_TOA << "," << F_up_surface << "," << F_down_toa
          << "," << heating_rate_max;

  // Phase 18: Append SEB residual columns if finite (backward compatible: write NaN if not available)
  outfile << "," << seb_residual_mean << "," << seb_residual_max;
  outfile << "\n";
  outfile.close();
}
