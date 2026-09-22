#ifndef ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_
#define ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_

#include <cctype>
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

enum class LegacyInitType {
    WRFInput,
    Metgrid,
    Unknown
};

enum class LegacySurfaceTemperatureCompatibility {
    Compatible,
    UnsafeMetgrid,
    UnknownProvenance
};

inline LegacyInitType
parse_legacy_init_type_value (std::string value)
{
    for (char& character : value) {
        character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    }
    if (value == "wrfinput") {
        return LegacyInitType::WRFInput;
    }
    if (value == "metgrid") {
        return LegacyInitType::Metgrid;
    }
    return LegacyInitType::Unknown;
}

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
inline LegacyInitType
parse_legacy_init_type_from_job_info (std::istream& input)
{
    constexpr const char* key = "erf.init_type";
    constexpr std::size_t key_length = 13;
    LegacyInitType result = LegacyInitType::Unknown;
    std::string line;

    while (std::getline(input, line)) {
        std::size_t pos = 0;
        while (pos < line.size() &&
               std::isspace(static_cast<unsigned char>(line[pos]))) {
            ++pos;
        }
        if (line.compare(pos, key_length, key) != 0) {
            continue;
        }
        pos += key_length;
        if (pos < line.size() && !std::isspace(static_cast<unsigned char>(line[pos])) &&
            line[pos] != '=') {
            continue;
        }
        while (pos < line.size() &&
               std::isspace(static_cast<unsigned char>(line[pos]))) {
            ++pos;
        }
        if (pos == line.size() || line[pos] != '=') {
            result = LegacyInitType::Unknown;
            continue;
        }
        ++pos;
        while (pos < line.size() &&
               std::isspace(static_cast<unsigned char>(line[pos]))) {
            ++pos;
        }

        std::string value;
        bool malformed = (pos == line.size());
        if (!malformed && line[pos] == '"') {
            ++pos;
            const std::size_t end_quote = line.find('"', pos);
            if (end_quote == std::string::npos) {
                malformed = true;
            } else {
                value = line.substr(pos, end_quote - pos);
                pos = end_quote + 1;
            }
        } else if (!malformed) {
            const std::size_t value_begin = pos;
            while (pos < line.size() &&
                   !std::isspace(static_cast<unsigned char>(line[pos]))) {
                ++pos;
            }
            value = line.substr(value_begin, pos - value_begin);
        }

        while (pos < line.size() &&
               std::isspace(static_cast<unsigned char>(line[pos]))) {
            ++pos;
        }
        if (malformed || pos != line.size()) {
            result = LegacyInitType::Unknown;
            continue;
        }

        result = parse_legacy_init_type_value(value);
    }

    return result;
}

inline LegacySurfaceTemperatureCompatibility
classify_legacy_surface_temperature_checkpoint (
    const bool marker_present,
    const bool has_surface_arrays,
    const LegacyInitType legacy_init_type)
{
    if (marker_present || !has_surface_arrays) {
        return LegacySurfaceTemperatureCompatibility::Compatible;
    }
    if (legacy_init_type == LegacyInitType::WRFInput) {
        return LegacySurfaceTemperatureCompatibility::Compatible;
    }
    if (legacy_init_type == LegacyInitType::Metgrid) {
        return LegacySurfaceTemperatureCompatibility::UnsafeMetgrid;
    }
    return LegacySurfaceTemperatureCompatibility::UnknownProvenance;
}

} // namespace erf_checkpoint_surface_temperature

#endif
