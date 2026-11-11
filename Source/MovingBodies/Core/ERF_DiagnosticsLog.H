#ifndef ERF_DIAGNOSTICS_LOG_H
#define ERF_DIAGNOSTICS_LOG_H

/** Opening of the moving-bodies diagnostics files for a fresh or a restarted run. */

#include <fstream>
#include <string>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

namespace erf_actuator {

/**
 * Open a diagnostics text file for writing (the caller writes CSV or space-separated rows),
 * creating its directory. With truncate (a fresh run) the file starts empty; without it (a
 * restarted run) the rows are appended to what the run it continues wrote. Returns whether the
 * caller must write the header: always after truncating, and on append only when the file is
 * new or empty (a restart in a clean directory).
 */
inline bool open_log (std::ofstream& out, const std::string& fname, bool truncate)
{
    const auto slash = fname.rfind('/');
    if (slash != std::string::npos) { amrex::UtilCreateDirectory(fname.substr(0, slash), 0755); }
    bool need_header = truncate;
    if (!truncate) {
        std::ifstream probe(fname, std::ios::ate);
        need_header = !probe || probe.tellg() <= 0;
    }
    out.open(fname, truncate ? std::ios::trunc : std::ios::app);
    if (!out) { amrex::Abort("cannot open moving-bodies diagnostics file '" + fname + "'"); }
    return need_header;
}

/**
 * Whether fname exists, as the I/O rank sees it, on every rank: a checkpoint read that every rank
 * enters (ParallelDescriptor::ReadAndBcastFile) must be decided alike everywhere, whatever a
 * parallel file system shows each node.
 */
inline bool file_exists_everywhere (const std::string& fname)
{
    int exists = amrex::ParallelDescriptor::IOProcessor() ? (amrex::FileExists(fname) ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&exists, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return exists != 0;
}

} // namespace erf_actuator

#endif
