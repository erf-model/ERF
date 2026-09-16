#ifndef ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_
#define ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_

#include <istream>
#include <ostream>
#include <string>

namespace erf_checkpoint_surface_temperature {

constexpr int contract_version = 1;

inline void
write_contract_version (std::ostream& output)
{
    output << contract_version << "\n";
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
