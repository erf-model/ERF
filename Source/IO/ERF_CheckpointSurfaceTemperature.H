#ifndef ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_
#define ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_

#include <cctype>
#include <istream>
#include <ostream>
#include <string>

#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>

#include "ERF_EOS.H"

namespace erf_checkpoint_surface_temperature {

/**
 * Version of the surface_temperature_contract marker written with a checkpoint, which says
 * what its SST_0 and TSK_0 arrays hold:
 *
 *   1: potential temperature, referenced to the surface pressure of the input they were
 *      read from
 *   2: absolute temperature [K], which the surface layer converts at the current surface
 *      pressure
 */
constexpr int contract_version = 2;

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

/**
 * What the SST_0 and TSK_0 arrays of a checkpoint hold
 */
enum class StoredSurfaceTemperature {
    Absolute,          ///< absolute temperature [K], restored as it is
    Potential,         ///< potential temperature, converted to absolute temperature on restart
    UnknownProvenance  ///< no marker, and the init type that wrote them is not known
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
 * Parse the integer surface-temperature contract marker from a checkpoint.  Every version
 * up to contract_version is understood.
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

    return (version >= 1 && version <= contract_version) ? ContractReadStatus::Valid
                                                          : ContractReadStatus::UnknownVersion;
}

/**
 * Read the init type a checkpoint without a contract marker was written by, from the last
 * erf.init_type assignment in its job_info.  Before the marker, WRFInput stored SST and TSK
 * as potential temperature and Metgrid as absolute temperature, so this tells which of the
 * two the checkpoint holds.
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

/**
 * Classify the SST_0 and TSK_0 arrays of a checkpoint.  marker_version is the version of its
 * contract marker, or 0 if it has none; a checkpoint without one is classified by the init
 * type that wrote it (see parse_legacy_init_type_from_job_info).
 */
inline StoredSurfaceTemperature
classify_stored_surface_temperature (const int marker_version,
                                     const bool has_surface_arrays,
                                     const LegacyInitType legacy_init_type)
{
    if (!has_surface_arrays || marker_version >= 2) {
        return StoredSurfaceTemperature::Absolute;
    }
    if (marker_version == 1 || legacy_init_type == LegacyInitType::WRFInput) {
        return StoredSurfaceTemperature::Potential;
    }
    if (legacy_init_type == LegacyInitType::Metgrid) {
        return StoredSurfaceTemperature::Absolute;
    }
    return StoredSurfaceTemperature::UnknownProvenance;
}

/**
 * Convert a surface field (SST or TSK) from potential to absolute temperature with the
 * surface pressure it was converted with, T = theta (p/p0)^(R_d/c_p).  A ghost cell takes the
 * pressure of the nearest valid cell of its box; the ghost cells are then filled from the
 * neighboring boxes where there are any.
 *
 * @param[in,out] field    potential temperature in, absolute temperature [K] out
 * @param[in]     pressure surface pressure [Pa] on the valid cells of the same 2D layout
 * @param[in]     rdOcp    R_d/c_p
 * @param[in]     period   periodicity of the level
 */
inline void
potential_to_absolute_temperature (amrex::MultiFab& field,
                                   const amrex::MultiFab& pressure,
                                   const amrex::Real rdOcp,
                                   const amrex::Periodicity& period)
{
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const amrex::Box& vbx = mfi.validbox();
        const amrex::IntVect lo = vbx.smallEnd();
        const amrex::IntVect hi = vbx.bigEnd();
        const amrex::Array4<amrex::Real> f = field.array(mfi);
        const amrex::Array4<const amrex::Real> p = pressure.const_array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int li = amrex::min(amrex::max(i, lo[0]), hi[0]);
            const int lj = amrex::min(amrex::max(j, lo[1]), hi[1]);
            f(i,j,k) = getTgivenPandTh(p(li,lj,lo[2]), f(i,j,k), rdOcp);
        });
    }
    field.FillBoundary(period);
}

} // namespace erf_checkpoint_surface_temperature

#endif
