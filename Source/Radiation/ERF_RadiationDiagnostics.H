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
 * Follows the duplicate-write-guard pattern from UCMDiagnostics (m_last_write_step).
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
   */
  RadiationDiagnostics(int verbosity, const std::string& diag_file, int amr_level);

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
};

#endif // ERF_RADIATION_DIAGNOSTICS_H_
