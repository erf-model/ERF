#ifndef ERF_RADIATION_DIAGNOSTICS_H_
#define ERF_RADIATION_DIAGNOSTICS_H_

#include <string>
#include <fstream>
#include <AMReX.H>
#include <AMReX_REAL.H>

/**
 * @file ERF_RadiationDiagnostics.H
 * @brief Radiation module diagnostic output and CSV logging.
 *
 * Follows the debug-print and diagnostic-output conventions established in
 * the ERF-SLUCM branch (see Source/UrbanCanopy/ERF_UCMDiagnostics.H for
 * reference pattern).
 *
 * TAGGED DEBUG FORMAT:
 * -------------------
 * All debug prints use the bracketed-tag prefix format:
 *   [RAD][Phase1][FunctionName] step={nstep} time={time} ...
 *
 * Tags:
 * - [RAD]       — module tag (parallel to [UCM])
 * - [Phase1]    — phase tag (will extend to [Phase2], [Phase3], etc.)
 * - [ClassName::function] — call-site tag for grepability
 *
 * VERBOSITY CONTROL:
 * ------------------
 * Debug prints only emit when erf.radiation.v >= 1 (default 0).
 * However, RADIATION_DIAG: lines always print (required for RegTest scripts).
 *
 * MPI SAFETY:
 * -----------
 * All host-side debug output is guarded by amrex::ParallelDescriptor::IOProcessor(),
 * matching the UCMDiagnostics pattern.
 *
 * CSV OUTPUT:
 * -----------
 * The RadiationDiagnostics class appends rows to a CSV file containing
 * simulation-level radiation quantities (e.g., surface SW/LW, TOA fluxes).
 * RegTest check scripts parse the RADIATION_DIAG: prefix to extract values.
 */

/**
 * @class RadiationDiagnostics
 * @brief CSV logger and debug-print manager for radiation module.
 *
 * Manages output of Phase 1 radiation diagnostics including:
 * - Tagged debug prints (when verbosity >= 1)
 * - CSV rows for flux/heating diagnostics
 * - RADIATION_DIAG: lines for RegTest parsing
 *
 * **Phase 9 Dedup Contract:**
 *
 * Duplicate detection uses a 3-tuple identity (step, call_site, time):
 * - step: timestep number (int)
 * - call_site: diagnostic call-site label (e.g., "pre_dycore", "post_dycore")
 * - time: simulation time with tolerance m_diag_dedup_tol
 *
 * This ensures:
 * 1. Accidental repeated calls at the same (step, call_site, time) are suppressed.
 * 2. Legitimate pre_dycore + post_dycore entries at the same step are both retained
 *    (they have different call_site values, so dedup identity differs).
 * 3. Mode filtering (pre_only, post_only, both) preserves this contract by checking
 *    whether to emit BEFORE the dedup guard runs.
 *
 * Follows the duplicate-write-guard pattern from UCMDiagnostics.
 */
class RadiationDiagnostics
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] verbosity Debug verbosity level (0=off, 1+=print).
   * @param[in] diag_file Output CSV file path (relative to working directory).
   * @param[in] amr_level AMR level index.
   * @param[in] diag_enable (Phase 7) Master switch for all diagnostics.
   * @param[in] diag_stdout_enable (Phase 7) Enable human-readable stdout block.
   * @param[in] diag_tagged_enable (Phase 7) Enable tagged debug lines.
   * @param[in] diag_regtest_line_enable (Phase 7) Enable RADIATION_DIAG: lines.
   * @param[in] diag_csv_enable (Phase 7) Enable CSV file writes.
   * @param[in] diag_callsite_mode (Phase 7) Call-site filter: "both", "pre_only", "post_only".
   * @param[in] diag_dedup_tol (Phase 7) Tolerance for time equality in duplicate guard.
   */
  RadiationDiagnostics(int verbosity, const std::string& diag_file, int amr_level,
                      bool diag_enable = true, bool diag_stdout_enable = true,
                      bool diag_tagged_enable = true, bool diag_regtest_line_enable = true,
                      bool diag_csv_enable = true, const std::string& diag_callsite_mode = "both",
                      amrex::Real diag_dedup_tol = 1.0e-12);

  /**
   * @brief Destructor.
   */
  ~RadiationDiagnostics();

  /**
   * @brief Append a diagnostic row to the CSV file.
   *
   * Records a row of Phase 1 radiation diagnostics for the current timestep:
   *
   * CSV columns (order):
   *   step, time, SW_surface, SW_TOA, F_up_surface, F_down_toa, heating_rate_max
   *
   * Also prints tagged debug output and RADIATION_DIAG: line for RegTest parsing.
   *
   * Guards against duplicate writes (uses m_last_write_step).
   *
   * @param[in] step           Current timestep number.
   * @param[in] time           Current simulation time [s].
   * @param[in] SW_surface     Downwelling SW flux at surface [W/m^2].
   * @param[in] SW_TOA         Downwelling SW flux at TOA (should match input) [W/m^2].
   * @param[in] F_up_surface   Upwelling LW flux at surface [W/m^2].
   * @param[in] F_down_toa     Downwelling LW flux at TOA [W/m^2].
   * @param[in] heating_rate_max Maximum heating rate magnitude [K/s].
   */
  /*void append(int step, amrex::Real time, amrex::Real SW_surface,
              amrex::Real SW_TOA, amrex::Real F_up_surface,
              amrex::Real F_down_toa, amrex::Real heating_rate_max);*/
  void append(int step, amrex::Real time, std::string const& call_site,
            amrex::Real sw_surface, amrex::Real sw_toa,
            amrex::Real f_up_surface, amrex::Real f_down_toa,
            amrex::Real heating_rate_max);

private:
  /**
   * @brief Write CSV header if file is new or doesn't exist.
   */
  void write_header_if_needed();

  int m_verbosity;              ///< Debug verbosity level
  std::string m_diag_file;      ///< Output CSV file path
  int m_amr_level;              ///< AMR level index
  int m_last_write_step = -1;   ///< Guard against duplicate writes
  bool m_header_written = false; ///< Track whether CSV header was written
  std::string m_last_write_call_site = "";
  amrex::Real m_last_write_time = -1.0;

  // Phase 7 diagnostics controls
  bool m_diag_enable = true;
  bool m_diag_stdout_enable = true;
  bool m_diag_tagged_enable = true;
  bool m_diag_regtest_line_enable = true;
  bool m_diag_csv_enable = true;
  std::string m_diag_callsite_mode = "both";
  amrex::Real m_diag_dedup_tol = 1.0e-12;
};

#endif // ERF_RADIATION_DIAGNOSTICS_H_
