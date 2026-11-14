#ifndef ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_
#define ERF_CHECKPOINT_SURFACE_TEMPERATURE_H_

#include <cctype>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Periodicity.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include "ERF_EOS.H"
#include "ERF_SurfaceTemperature.H"

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

/** Name of the marker file in the checkpoint directory */
constexpr const char* contract_file_name = "surface_temperature_contract";

/** Parameter (in the erf namespace) asserting the init type that wrote a checkpoint */
constexpr const char* legacy_init_type_key = "legacy_surface_temperature_init_type";

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

inline const char*
legacy_init_type_name (const LegacyInitType type)
{
    switch (type) {
    case LegacyInitType::WRFInput:
        return "WRFInput";
    case LegacyInitType::Metgrid:
        return "Metgrid";
    default:
        return "unknown";
    }
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
 * Read the init type that wrote a checkpoint, from the last erf.init_type assignment in its
 * job_info.  Before the marker, WRFInput stored SST and TSK as potential temperature and
 * Metgrid as absolute temperature, so for a checkpoint without one this tells which of the
 * two it holds.  For one with a version 1 marker, which holds potential temperature, it
 * tells which surface pressure they were converted with.
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
 * Where the init type that wrote the SST_0 and TSK_0 of a checkpoint was learned from
 */
enum class ProvenanceOrigin {
    JobInfo,    ///< the erf.init_type recorded in the job_info of the checkpoint
    Assertion,  ///< erf.legacy_surface_temperature_init_type, in the inputs of the restart
    None        ///< not known
};

/**
 * What the SST_0 and TSK_0 arrays of a checkpoint hold, and which init type wrote them
 */
struct SurfaceTemperatureProvenance {
    int marker_version = 0;   ///< version of the contract marker, 0 if there is none
    int first_level    = -1;  ///< first AMR level holding SST_0 or TSK_0, -1 if none does
    LegacyInitType source   = LegacyInitType::Unknown;  ///< init type that wrote them
    ProvenanceOrigin origin = ProvenanceOrigin::None;   ///< where source was learned from
    LegacyInitType asserted = LegacyInitType::Unknown;  ///< init type the restart inputs assert
    StoredSurfaceTemperature stored = StoredSurfaceTemperature::Absolute;
};

/**
 * Establish what the SST_0 and TSK_0 of a checkpoint hold, and which init type wrote them.
 *
 * The erf.init_type recorded in the job_info of the checkpoint is its provenance.  The
 * erf.legacy_surface_temperature_init_type assertion stands in for it only where the
 * checkpoint records none.  The init type of the restart plays no part: a run can be
 * restarted with inputs other than those that wrote the checkpoint.
 *
 * @param marker_version  version of the contract marker, 0 if there is none
 * @param first_level     first AMR level holding SST_0 or TSK_0, -1 if none does
 * @param job_info_source init type recorded in the job_info of the checkpoint
 * @param asserted_source init type asserted by erf.legacy_surface_temperature_init_type
 */
inline SurfaceTemperatureProvenance
resolve_surface_temperature_provenance (const int marker_version,
                                        const int first_level,
                                        const LegacyInitType job_info_source,
                                        const LegacyInitType asserted_source)
{
    SurfaceTemperatureProvenance result;
    result.marker_version = marker_version;
    result.first_level = first_level;
    result.asserted = asserted_source;
    if (job_info_source != LegacyInitType::Unknown) {
        result.source = job_info_source;
        result.origin = ProvenanceOrigin::JobInfo;
    } else if (asserted_source != LegacyInitType::Unknown) {
        result.source = asserted_source;
        result.origin = ProvenanceOrigin::Assertion;
    }
    result.stored = classify_stored_surface_temperature(marker_version, first_level >= 0,
                                                        result.source);
    return result;
}

/**
 * How the SST_0 and TSK_0 of a checkpoint are brought to the absolute temperature they are
 * now kept as
 */
enum class SurfaceTemperatureMigration {
    None,             ///< they hold absolute temperature: restored as they are
    InputPressure,    ///< exact: converted back with the wrfinput PSFC they were converted with
    RestartPressure,  ///< approximate: converted with the surface pressure of the restarted state
    Reject            ///< they cannot be migrated: the restart aborts, giving the reason
};

struct SurfaceTemperatureMigrationPlan {
    SurfaceTemperatureMigration method = SurfaceTemperatureMigration::None;
    SurfaceTemperatureProvenance provenance;  ///< what the decision was made from
    std::string reason;                       ///< for Reject: what is missing, and how to supply it
};

/**
 * Decide how to bring the SST_0 and TSK_0 of a checkpoint to absolute temperature.
 *
 * Arrays that hold potential temperature were converted with the surface pressure of the
 * input they were read from, which only the init type that wrote them identifies.  A WRFInput
 * run converted them with the PSFC of its wrfinput files: a WRFInput restart reads that PSFC
 * again and recovers them exactly, provided its inputs name the wrfinput files of the original
 * run.  A Metgrid run converted them with the surface pressure of its met_em files, which a
 * restart does not read: they are migrated with the surface pressure of the restarted state
 * instead.  That keeps the surface potential temperature the run had at the restart, not the
 * SST of the input: the absolute temperature it gives carries the factor
 * (p_restart/p_input)^(R_d/c_p) between the two pressures.
 *
 * @param provenance             what the checkpoint holds and which init type wrote it
 * @param restart_reads_wrfinput whether the restart is a WRFInput run, in a build with NetCDF
 * @param checkpointname         name of the checkpoint, for the reason of a rejection
 */
inline SurfaceTemperatureMigrationPlan
plan_surface_temperature_migration (const SurfaceTemperatureProvenance& provenance,
                                    const bool restart_reads_wrfinput,
                                    const std::string& checkpointname)
{
    SurfaceTemperatureMigrationPlan plan;
    plan.provenance = provenance;
    if (provenance.stored == StoredSurfaceTemperature::Absolute) {
        return plan;
    }

    const std::string holds = "Restart checkpoint '" + checkpointname + "' holds SST_0/TSK_0 " +
                              "(from AMR level " + std::to_string(provenance.first_level) + ")";
    const std::string set_assertion = "Set erf." + std::string(legacy_init_type_key) +
                                      " to the erf.init_type of the run that wrote it: ";
    plan.method = SurfaceTemperatureMigration::Reject;
    if (provenance.stored == StoredSurfaceTemperature::UnknownProvenance) {
        plan.reason = holds + " without a surface-temperature contract marker, and its job_info "
                      "does not record the erf.init_type that wrote them. " + set_assertion +
                      "WRFInput, which stored them as potential temperature, or Metgrid, which "
                      "stored them as absolute temperature.";
        return plan;
    }

    switch (provenance.source) {
    case LegacyInitType::WRFInput:
        if (restart_reads_wrfinput) {
            plan.method = SurfaceTemperatureMigration::InputPressure;
        } else {
            plan.reason = holds + " as potential temperature, which a WRFInput run converted with "
                          "the surface pressure (PSFC) of its wrfinput files. Recovering them "
                          "needs that PSFC: restart with erf.init_type = WRFInput, in a build "
                          "with NetCDF, and the erf.nc_init_file_<lev> of the original run.";
        }
        break;
    case LegacyInitType::Metgrid:
        plan.method = SurfaceTemperatureMigration::RestartPressure;
        break;
    case LegacyInitType::Unknown:
        plan.reason = holds + " as potential temperature (surface-temperature contract version " +
                      std::to_string(provenance.marker_version) + "), but its job_info does not "
                      "record the erf.init_type that converted them, so the surface pressure "
                      "they were converted with is not known. " + set_assertion +
                      "WRFInput, to recover them exactly with the PSFC of its wrfinput files, or "
                      "Metgrid, to migrate them approximately with the surface pressure of the "
                      "restarted state.";
        break;
    }
    return plan;
}

/**
 * Where the wrfinput surface pressure that a level's SST and TSK were converted with is
 * recovered from
 */
enum class LegacySurfacePressureSource {
    OwnInitFile,   ///< the PSFC of the level's own (first) init file
    CoarserLevel,  ///< the coarser level's, interpolated: the level took its SST and TSK from it
    Unavailable    ///< it cannot be recovered
};

/**
 * Where the wrfinput surface pressure of a level is recovered from, on restart from a WRFInput
 * checkpoint that holds the SST and TSK as potential temperature.  It is read from where the
 * level read its SST and TSK, and converted them with:
 *
 *  - level 0 read them from its first init file, the only one init_from_wrfinput reads there,
 *    however many erf.nc_init_file_0 lists;
 *  - a finer level read them from its own init file, unless it has none, or level 0 reads a
 *    wrflowinp file and the level has none of its own (see init_fine_low_slices): it then
 *    took them from the coarser level, interpolated, with the coarser level's PSFC;
 *  - a finer level with several init files read them from all of them, each over its own
 *    footprint, which is not recovered.
 *
 * @param lev                     AMR level
 * @param num_own_init_files      number of init files of the level (erf.nc_init_file_<lev>)
 * @param level0_reads_low_file   whether level 0 reads a wrflowinp file (erf.nc_low_file)
 * @param has_own_low_file        whether the level has a wrflowinp file of its own
 * @param coarser_level_recovered whether the surface pressure of the coarser level was recovered
 */
inline LegacySurfacePressureSource
legacy_wrfinput_pressure_source (const int lev,
                                 const int num_own_init_files,
                                 const bool level0_reads_low_file,
                                 const bool has_own_low_file,
                                 const bool coarser_level_recovered)
{
    if (lev == 0) {
        return (num_own_init_files > 0) ? LegacySurfacePressureSource::OwnInitFile
                                        : LegacySurfacePressureSource::Unavailable;
    }
    const bool from_coarser_level = (num_own_init_files == 0) ||
                                    (level0_reads_low_file && !has_own_low_file);
    if (from_coarser_level) {
        return coarser_level_recovered ? LegacySurfacePressureSource::CoarserLevel
                                       : LegacySurfacePressureSource::Unavailable;
    }
    return (num_own_init_files == 1) ? LegacySurfacePressureSource::OwnInitFile
                                     : LegacySurfacePressureSource::Unavailable;
}

/**
 * What decides, at each cell of a level, whether its surface layer takes the SST, the TSK
 * or the prescribed land temperature there (erf_surface_temperature::surface_temperature_source)
 */
struct SurfaceTemperatureUse {
    const amrex::iMultiFab* lmask = nullptr;  ///< land mask; without one, every cell is land
    bool have_tsk   = false;                  ///< whether the level has TSK
    bool ignore_sst = false;                  ///< whether its surface layer ignores the SST
};

/**
 * What decides it at a level on restart: the land mask restored from the checkpoint, and
 * whether the checkpoint holds TSK at the level.  The surface layer ignores the SST of a
 * level with TSK whatever most.ignore_sst says, and that parameter selects nothing at a
 * level without TSK, so it does not matter here.
 *
 * @param lmask    the restored land mask of the level, or nullptr
 * @param have_tsk whether the level has TSK
 */
inline SurfaceTemperatureUse
restart_surface_temperature_use (const amrex::iMultiFab* lmask, const bool have_tsk)
{
    return SurfaceTemperatureUse{lmask, have_tsk,
                                 erf_surface_temperature::sst_ignored(false, have_tsk)};
}

/** The name of a surface field: "SST" or "TSK" */
inline const char*
surface_field_name (const erf_surface_temperature::SurfaceTemperatureSource source)
{
    return (source == erf_surface_temperature::SurfaceTemperatureSource::tsk) ? "TSK" : "SST";
}

/**
 * The valid cells of a surface field that were not converted on this rank: those that could
 * not be, how many and the first one recorded, and the placeholders left as they are
 */
struct LegacyConversionFailure {
    int count     = 0;  ///< valid cells that could not be converted
    int preserved = 0;  ///< valid cells left as they are: placeholders where the field is not used
    int claimed   = 0;  ///< whether the first cell that could not be converted is described below
    int i = 0;
    int j = 0;
    int land = 0;       ///< its land mask
    int used = 0;       ///< whether the surface layer takes the field there
    amrex::Real theta    = amrex::Real(0.0);
    amrex::Real pressure = amrex::Real(0.0);
};

/**
 * Convert a surface field (SST or TSK) from potential to absolute temperature with the
 * surface pressure it was converted with, T = theta (p/p0)^(R_d/c_p).  A ghost cell takes the
 * pressure of the nearest valid cell of its box; the ghost cells are then filled from the
 * neighboring boxes where there are any.
 *
 * Every cell that holds a finite, positive potential temperature is converted, whether or not
 * the surface layer takes the field there, so that no value the run may come to use (after a
 * regrid or an interpolation, say) is left as a potential temperature; that takes a finite,
 * positive surface pressure.  A cell that holds none is left as it is.  Of the valid cells
 * left unconverted:
 *
 *  - one where the surface layer takes the field (see SurfaceTemperatureUse), or that holds a
 *    potential temperature it has no pressure to convert with, is counted, and the first one
 *    described, for abort_on_legacy_conversion_failure;
 *  - one where the surface layer does not take the field, and that holds no temperature, is
 *    a placeholder, such as a zero SST over land: it is only counted as preserved.  It is not
 *    made to look like a temperature, so should the field ever be taken there, the surface
 *    layer rejects it then.
 *
 * A ghost cell left unconverted is not counted.
 *
 * @param[in,out] field    potential temperature in, absolute temperature [K] out
 * @param[in]     source   which field it is: SurfaceTemperatureSource::sst or ::tsk
 * @param[in]     use      where the surface layer takes it; the land mask, if any, on the
 *                         layout of field
 * @param[in]     pressure surface pressure [Pa] on the valid cells of the layout of field
 * @param[in]     rdOcp    R_d/c_p
 * @param[in]     period   periodicity of the level
 * @return the valid cells of this rank that were not converted
 */
inline LegacyConversionFailure
potential_to_absolute_temperature (amrex::MultiFab& field,
                                   const erf_surface_temperature::SurfaceTemperatureSource source,
                                   const SurfaceTemperatureUse& use,
                                   const amrex::MultiFab& pressure,
                                   const amrex::Real rdOcp,
                                   const amrex::Periodicity& period)
{
    AMREX_ALWAYS_ASSERT(source != erf_surface_temperature::SurfaceTemperatureSource::prescribed);
    AMREX_ALWAYS_ASSERT(!use.lmask || (use.lmask->boxArray() == field.boxArray() &&
                                       use.lmask->DistributionMap() == field.DistributionMap()));
    amrex::Gpu::DeviceScalar<LegacyConversionFailure> d_fail(LegacyConversionFailure{});
    LegacyConversionFailure* fail = d_fail.dataPtr();
    const bool have_tsk = use.have_tsk;
    const bool ignore_sst = use.ignore_sst;
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const amrex::Box& vbx = mfi.validbox();
        const amrex::IntVect lo = vbx.smallEnd();
        const amrex::IntVect hi = vbx.bigEnd();
        const amrex::Array4<amrex::Real> f = field.array(mfi);
        const amrex::Array4<const amrex::Real> p = pressure.const_array(mfi);
        const amrex::Array4<const int> mask = (use.lmask) ? use.lmask->const_array(mfi)
                                                          : amrex::Array4<const int>{};
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int li = amrex::min(amrex::max(i, lo[0]), hi[0]);
            const int lj = amrex::min(amrex::max(j, lo[1]), hi[1]);
            const amrex::Real theta = f(i,j,k);
            const amrex::Real p_s = p(li,lj,lo[2]);
            amrex::Real temperature = theta;
            if (erf_surface_temperature::theta_to_temperature(theta, p_s, rdOcp, temperature)) {
                f(i,j,k) = temperature;
                return;
            }
            if (li != i || lj != j) {
                return;
            }
            const int land = (mask) ? mask(i,j,k) : 1;
            const bool used = erf_surface_temperature::surface_temperature_source(
                                  land != 0, have_tsk, ignore_sst) == source;
            if (!used && !erf_surface_temperature::valid_temperature(theta)) {
                amrex::Gpu::Atomic::Add(&fail->preserved, 1);
                return;
            }
            amrex::Gpu::Atomic::Add(&fail->count, 1);
            if (amrex::Gpu::Atomic::CAS(&fail->claimed, 0, 1) == 0) {
                fail->i = i;
                fail->j = j;
                fail->land = land;
                fail->used = used ? 1 : 0;
                fail->theta = theta;
                fail->pressure = p_s;
            }
        });
    }
    amrex::Gpu::streamSynchronize();
    field.FillBoundary(period);
    return d_fail.dataValue();
}

/**
 * Say how many placeholders of a converted surface field were left as they are, and abort the
 * restart if any rank could not convert a valid cell, after each rank that could not
 * describes its first such cell.  Collective: every rank must call this.
 *
 * @param[in] fail   this rank's record, from potential_to_absolute_temperature
 * @param[in] source which field was converted: SurfaceTemperatureSource::sst or ::tsk
 * @param[in] what   what was being converted, to open the messages with
 */
inline void
abort_on_legacy_conversion_failure (const LegacyConversionFailure& fail,
                                    const erf_surface_temperature::SurfaceTemperatureSource source,
                                    const std::string& what)
{
    const std::string name = surface_field_name(source);
    int counts[2] = {fail.count, fail.preserved};
    amrex::ParallelDescriptor::ReduceIntSum(counts, 2);
    const int count = counts[0];
    const int preserved = counts[1];
    if (preserved > 0) {
        amrex::Print() << what << ": kept " << preserved << " cells where the surface layer "
                       << "does not take the " << name << " as they are, placeholders without "
                       << "a finite, positive temperature\n";
    }
    if (count == 0) {
        return;
    }
    if (fail.count > 0) {
        // Assemble the whole report before writing it, so that the ranks do not interleave
        std::ostringstream msg;
        msg << std::setprecision(10) << what << " [rank " << amrex::ParallelDescriptor::MyProc()
            << "]: " << fail.count << " cells cannot be converted; the first, (" << fail.i
            << "," << fail.j << "), a " << ((fail.land != 0) ? "land" : "water") << " cell "
            << "where the surface layer " << ((fail.used != 0) ? "takes" : "does not take")
            << " the " << name << ", has theta = " << fail.theta << ", surface pressure = "
            << fail.pressure << " Pa\n";
        amrex::AllPrint() << msg.str() << std::flush;
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::Abort(what + ": " + std::to_string(count) + " cells could not be converted to "
                 "absolute temperature (see above). A cell where the surface layer takes the " +
                 name + " needs a finite, positive potential temperature, and any cell that "
                 "holds one a finite, positive surface pressure to convert it with");
}

/**
 * Give each cell of a planar field that lies outside the domain, in a direction that is not
 * periodic, the value of the nearest cell of the domain, so that a finer level can interpolate
 * the field up to the domain boundary.  The cells of the domain must be filled.
 *
 * @param[in,out] mf   planar field, cell-centered
 * @param[in]     geom geometry of its level
 */
inline void
fill_outside_domain_from_nearest (amrex::MultiFab& mf, const amrex::Geometry& geom)
{
    const amrex::Box& domain = geom.Domain();
    const amrex::IntVect dlo = domain.smallEnd();
    const amrex::IntVect dhi = domain.bigEnd();
    const bool periodic_x = geom.isPeriodic(0);
    const bool periodic_y = geom.isPeriodic(1);
    const int ncomp = mf.nComp();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Array4<amrex::Real> a = mf.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ni = periodic_x ? i : amrex::min(amrex::max(i, dlo[0]), dhi[0]);
            const int nj = periodic_y ? j : amrex::min(amrex::max(j, dlo[1]), dhi[1]);
            if (ni != i || nj != j) {
                for (int n = 0; n < ncomp; ++n) {
                    a(i,j,k,n) = a(ni,nj,k,n);
                }
            }
        });
    }
}

/**
 * Read what the SST_0 and TSK_0 arrays of a checkpoint hold, and which init type wrote them.
 *
 * The contract marker says what they hold.  Where it does not settle that they hold absolute
 * temperature (no marker, or version 1), the init type that wrote them is read from the
 * erf.init_type in the job_info of the checkpoint or, where that records none, from the
 * erf.legacy_surface_temperature_init_type assertion.  Aborts on a malformed or unknown marker
 * or an invalid assertion.
 *
 * @param checkpointname name of the checkpoint directory
 * @param finest_level   finest AMR level of the checkpoint
 */
inline SurfaceTemperatureProvenance
read_surface_temperature_provenance (const std::string& checkpointname,
                                     const int finest_level)
{
    int marker_version = 0;
    const std::string marker_name = checkpointname + "/" + contract_file_name;
    if (amrex::FileExists(marker_name)) {
        amrex::Vector<char> marker_chars;
        amrex::ParallelDescriptor::ReadAndBcastFile(marker_name, marker_chars);
        std::istringstream marker_stream(std::string(marker_chars.dataPtr()),
                                         std::istringstream::in);
        const auto status = read_contract_version(marker_stream, marker_version);
        if (status == ContractReadStatus::Malformed) {
            amrex::Abort("Malformed surface-temperature contract marker in '" + marker_name + "'");
        }
        if (status == ContractReadStatus::UnknownVersion) {
            amrex::Abort("Unsupported surface-temperature contract version " +
                         std::to_string(marker_version) + " in '" + marker_name + "'");
        }
    }
    const int first_level = first_legacy_surface_temperature_level(checkpointname, finest_level);

    // Which init type wrote the arrays matters only where they may hold potential temperature
    LegacyInitType job_info_source = LegacyInitType::Unknown;
    LegacyInitType asserted_source = LegacyInitType::Unknown;
    if (classify_stored_surface_temperature(marker_version, first_level >= 0,
                                            LegacyInitType::Unknown) !=
        StoredSurfaceTemperature::Absolute)
    {
        const std::string job_info_name = checkpointname + "/job_info";
        if (amrex::FileExists(job_info_name)) {
            amrex::Vector<char> job_info_chars;
            amrex::ParallelDescriptor::ReadAndBcastFile(job_info_name, job_info_chars);
            std::istringstream job_info_stream(std::string(job_info_chars.dataPtr()),
                                               std::istringstream::in);
            job_info_source = parse_legacy_init_type_from_job_info(job_info_stream);
        }

        amrex::ParmParse pp_erf("erf");
        if (pp_erf.countval(legacy_init_type_key) > 0) {
            std::string value;
            pp_erf.get(legacy_init_type_key, value);
            asserted_source = parse_legacy_init_type_value(value);
            if (asserted_source == LegacyInitType::Unknown) {
                amrex::Abort("Invalid erf." + std::string(legacy_init_type_key) + " value '" +
                             value + "'; accepted values are WRFInput and Metgrid "
                             "(case-insensitive).");
            }
        }
    }

    const auto provenance = resolve_surface_temperature_provenance(
        marker_version, first_level, job_info_source, asserted_source);
    if (provenance.origin == ProvenanceOrigin::Assertion) {
        amrex::Print() << "WARNING: checkpoint '" << checkpointname << "' does not record the "
                       << "erf.init_type that wrote its SST_0/TSK_0; trusting the explicit erf."
                       << legacy_init_type_key << " = " << legacy_init_type_name(asserted_source)
                       << " assertion.\n";
    } else if (asserted_source != LegacyInitType::Unknown && asserted_source != job_info_source) {
        amrex::Print() << "WARNING: ignoring erf." << legacy_init_type_key << " = "
                       << legacy_init_type_name(asserted_source) << ": checkpoint '"
                       << checkpointname << "' records erf.init_type = "
                       << legacy_init_type_name(job_info_source) << " in its job_info.\n";
    }
    return provenance;
}

} // namespace erf_checkpoint_surface_temperature

#endif
