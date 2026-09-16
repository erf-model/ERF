#ifndef ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_
#define ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_

#include <istream>
#include <ostream>
#include <string>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>

namespace erf_checkpoint_surface_temperature {

constexpr int contract_version = 1;

inline void
write_contract_version (std::ostream& output)
{
    output << contract_version << "\n";
}

/**
 * Return the first AMR level containing a legacy SST or TSK surface array.
 *
 * The arrays are identified by their AMReX header files.  A return value of
 * -1 means that neither legacy array is present on any level.
 */
inline int
first_legacy_surface_temperature_level (const std::string& checkpointname,
                                         const int finest_level)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        const bool has_sst = amrex::FileExists(
            amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "SST_0_H"));
        const bool has_tsk = amrex::FileExists(
            amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "TSK_0_H"));
        if (has_sst || has_tsk) {
            return lev;
        }
    }
    return -1;
}

enum class ContractReadStatus {
    Valid,
    Malformed,
    UnknownVersion
};

/**
 * Parse the integer surface-temperature contract marker from a checkpoint.
 */
inline ContractReadStatus
read_contract_version (std::istream& input, int& version)
{
    if (!(input >> version)) {
        return ContractReadStatus::Malformed;
    }

    std::string extra;
    if (input >> extra) {
        return ContractReadStatus::Malformed;
    }

    return (version == contract_version) ? ContractReadStatus::Valid
                                         : ContractReadStatus::UnknownVersion;
}

/**
 * Legacy Metgrid SST/TSK arrays are absolute temperature, whereas current
 * Metgrid initialization stores the SurfaceLayer representation (theta).
 * Without the marker there is no pressure in a checkpoint with which to
 * reinterpret those legacy arrays safely.
 */
inline bool
legacy_metgrid_surface_temperature_is_unsafe (const bool marker_present,
                                              const bool is_metgrid,
                                              const bool has_sst,
                                              const bool has_tsk)
{
    return !marker_present && is_metgrid && (has_sst || has_tsk);
}

} // namespace erf_checkpoint_surface_temperature

#endif
