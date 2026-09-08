/**
 * \file ERF_Provenance.cpp
 */
#include "ERF_Provenance.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace
{

using erf_provenance::ArtifactType;
using erf_provenance::LineageStatus;
using erf_provenance::ProvenanceReadStatus;

bool is_hex (char c)
{
    return (c >= '0' && c <= '9') ||
           (c >= 'a' && c <= 'f') ||
           (c >= 'A' && c <= 'F');
}

bool is_canonical_uuid_layout (std::string_view uuid)
{
    if (uuid.size() != 36) {
        return false;
    }
    for (std::size_t i = 0; i < uuid.size(); ++i) {
        const bool hyphen_position = i == 8 || i == 13 || i == 18 || i == 23;
        if (hyphen_position ? uuid[i] != '-' : !is_hex(uuid[i])) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> split_lineage (std::string_view value)
{
    std::vector<std::string> result;
    std::size_t start = 0;
    while (start <= value.size()) {
        const std::size_t comma = value.find(',', start);
        const std::size_t end = comma == std::string_view::npos ? value.size() : comma;
        if (end == start) {
            return {};
        }
        result.emplace_back(value.substr(start, end - start));
        if (comma == std::string_view::npos) break;
        start = comma + 1;
    }
    return result;
}

bool parse_int (std::string_view text, int& value)
{
    try {
        std::size_t used = 0;
        const std::string input(text);
        const long parsed = std::stol(input, &used, 10);
        if (used != input.size() || parsed < std::numeric_limits<int>::min() ||
            parsed > std::numeric_limits<int>::max()) {
            return false;
        }
        value = static_cast<int>(parsed);
        return true;
    } catch (...) {
        return false;
    }
}

bool parse_double (std::string_view text, double& value)
{
    try {
        std::size_t used = 0;
        const std::string input(text);
        value = std::stod(input, &used);
        return used == input.size();
    } catch (...) {
        return false;
    }
}

bool parse_bool (std::string_view text, bool& value)
{
    if (text == "true") {
        value = true;
        return true;
    }
    if (text == "false") {
        value = false;
        return true;
    }
    return false;
}

bool parse_artifact_type (std::string_view text, ArtifactType& type)
{
    if (text == "checkpoint") {
        type = ArtifactType::Checkpoint;
    } else if (text == "plotfile_2d") {
        type = ArtifactType::Plotfile2D;
    } else if (text == "plotfile_3d") {
        type = ArtifactType::Plotfile3D;
    } else {
        return false;
    }
    return true;
}

bool parse_lineage_status (std::string_view text, LineageStatus& status)
{
    if (text == "complete") status = LineageStatus::Complete;
    else if (text == "incomplete_ancestor") status = LineageStatus::IncompleteAncestor;
    else if (text == "missing_job_info") status = LineageStatus::MissingJobInfo;
    else if (text == "missing_provenance_block") status = LineageStatus::MissingProvenanceBlock;
    else if (text == "malformed_provenance") status = LineageStatus::MalformedProvenance;
    else if (text == "unsupported_schema") status = LineageStatus::UnsupportedSchema;
    else if (text == "artifact_type_mismatch") status = LineageStatus::ArtifactTypeMismatch;
    else return false;
    return true;
}

LineageStatus failure_lineage_status (ProvenanceReadStatus status)
{
    switch (status) {
    case ProvenanceReadStatus::MissingJobInfo: return LineageStatus::MissingJobInfo;
    case ProvenanceReadStatus::MissingProvenanceBlock: return LineageStatus::MissingProvenanceBlock;
    case ProvenanceReadStatus::MalformedProvenance: return LineageStatus::MalformedProvenance;
    case ProvenanceReadStatus::UnsupportedSchema: return LineageStatus::UnsupportedSchema;
    case ProvenanceReadStatus::ArtifactTypeMismatch: return LineageStatus::ArtifactTypeMismatch;
    case ProvenanceReadStatus::Valid: break;
    }
    return LineageStatus::MalformedProvenance;
}

void broadcast_string (std::string& value)
{
    const int io_rank = amrex::ParallelDescriptor::IOProcessorNumber();
    int length = amrex::ParallelDescriptor::IOProcessor()
               ? static_cast<int>(value.size()) : 0;
    amrex::ParallelDescriptor::Bcast(&length, 1, io_rank);
    value.resize(static_cast<std::size_t>(length));
    if (length > 0) {
        amrex::ParallelDescriptor::Bcast(value.data(), length, io_rank);
    }
}

std::string strip_cr (std::string line)
{
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return line;
}

bool validate_record (const erf_provenance::ProvenanceRecord& record, std::string& diagnostic)
{
    const auto& execution = record.execution;
    const auto& artifact = record.artifact;
    if (execution.schema_version != erf_provenance::schema_version) {
        diagnostic = "unsupported provenance schema";
        return false;
    }
    if (!erf_provenance::is_valid_uuid_v4(execution.simulation_uuid) ||
        !erf_provenance::is_valid_uuid_v4(execution.execution_uuid) ||
        (!execution.parent_execution_uuid.empty() &&
         !erf_provenance::is_valid_uuid_v4(execution.parent_execution_uuid)) ||
        !erf_provenance::is_valid_uuid_v4(artifact.artifact_uuid)) {
        diagnostic = "provenance contains an invalid UUID";
        return false;
    }
    if (execution.execution_lineage.empty()) {
        diagnostic = "execution lineage is empty";
        return false;
    }
    for (const auto& uuid : execution.execution_lineage) {
        if (!erf_provenance::is_valid_uuid_v4(uuid)) {
            diagnostic = "execution lineage contains an invalid UUID";
            return false;
        }
    }
    if (execution.execution_lineage.back() != execution.execution_uuid) {
        diagnostic = "lineage does not end at execution UUID";
        return false;
    }
    if (execution.lineage_complete) {
        if (execution.execution_lineage.front() != execution.simulation_uuid ||
            execution.restart_generation != static_cast<int>(execution.execution_lineage.size()) - 1 ||
            (execution.execution_lineage.size() == 1 && !execution.parent_execution_uuid.empty()) ||
            (execution.execution_lineage.size() > 1 &&
             execution.parent_execution_uuid != execution.execution_lineage[execution.execution_lineage.size()-2]) ||
            execution.lineage_status != LineageStatus::Complete) {
            diagnostic = "complete lineage invariants failed";
            return false;
        }
    } else {
        if (execution.restart_generation != -1 ||
            execution.lineage_status == LineageStatus::Complete ||
            (execution.execution_lineage.size() == 1 && !execution.parent_execution_uuid.empty()) ||
            (execution.execution_lineage.size() > 1 &&
             execution.parent_execution_uuid != execution.execution_lineage[execution.execution_lineage.size()-2])) {
            diagnostic = "incomplete lineage invariants failed";
            return false;
        }
    }
    if (!execution.source_checkpoint_uuid.empty() &&
        !erf_provenance::is_valid_uuid_v4(execution.source_checkpoint_uuid)) {
        diagnostic = "source checkpoint UUID is invalid";
        return false;
    }
    if (execution.execution_lineage.size() > 1 &&
        execution.source_checkpoint_uuid.empty()) {
        diagnostic = "restart lineage has no source checkpoint UUID";
        return false;
    }
    if (!std::isfinite(artifact.artifact_time_seconds)) {
        diagnostic = "artifact time is not finite";
        return false;
    }
    if (artifact.artifact_step < 0 || execution.execution_start_utc.empty() ||
        artifact.artifact_created_utc.empty()) {
        diagnostic = "required provenance value is empty or invalid";
        return false;
    }
    return true;
}

erf_provenance::ProvenanceParseResult failure (ProvenanceReadStatus status,
                                               std::string diagnostic)
{
    erf_provenance::ProvenanceParseResult result;
    result.status = status;
    result.diagnostic = std::move(diagnostic);
    return result;
}

} // namespace

namespace erf_provenance
{

std::string uuid_v4_from_bytes (const std::array<std::uint8_t, 16>& bytes)
{
    static constexpr char digits[] = "0123456789abcdef";
    std::string result;
    result.reserve(36);
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        if (i == 4 || i == 6 || i == 8 || i == 10) result.push_back('-');
        result.push_back(digits[(bytes[i] >> 4) & 0x0f]);
        result.push_back(digits[bytes[i] & 0x0f]);
    }
    return result;
}

bool is_valid_uuid_v4 (std::string_view uuid)
{
    if (!is_canonical_uuid_layout(uuid) || uuid[14] != '4') return false;
    const char variant = uuid[19];
    return variant == '8' || variant == '9' || variant == 'a' || variant == 'b' ||
           variant == 'A' || variant == 'B';
}

std::string generate_uuid_v4 ()
{
    std::array<std::uint8_t, 16> bytes{};
    bool generated = false;
    try {
        std::random_device random;
        for (auto& byte : bytes) byte = static_cast<std::uint8_t>(random());
        generated = true;
    } catch (...) {
        generated = false;
    }
    const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const auto steady = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto address = reinterpret_cast<std::uintptr_t>(&bytes);
    if (!generated) {
        for (std::size_t i = 0; i < bytes.size(); ++i) {
            bytes[i] = static_cast<std::uint8_t>((now >> ((i % 8) * 8)) ^
                                                 (steady >> (((i + 3) % 8) * 8)) ^
                                                 (address >> ((i % sizeof(std::uintptr_t)) * 8)));
        }
    } else {
        for (std::size_t i = 0; i < bytes.size(); ++i) {
            bytes[i] ^= static_cast<std::uint8_t>((now >> ((i % 8) * 8)) ^
                                                   (steady >> (((i + 3) % 8) * 8)) ^
                                                   (address >> ((i % sizeof(std::uintptr_t)) * 8)));
        }
    }
    bytes[6] = static_cast<std::uint8_t>((bytes[6] & 0x0f) | 0x40);
    bytes[8] = static_cast<std::uint8_t>((bytes[8] & 0x3f) | 0x80);
    return uuid_v4_from_bytes(bytes);
}

std::string format_utc (std::time_t time_value)
{
    std::tm utc{};
#if defined(_WIN32)
    gmtime_s(&utc, &time_value);
#else
    gmtime_r(&time_value, &utc);
#endif
    std::ostringstream out;
    out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

std::string current_utc ()
{
    return format_utc(std::time(nullptr));
}

const char* lineage_status_token (LineageStatus status) noexcept
{
    switch (status) {
    case LineageStatus::Complete: return "complete";
    case LineageStatus::IncompleteAncestor: return "incomplete_ancestor";
    case LineageStatus::MissingJobInfo: return "missing_job_info";
    case LineageStatus::MissingProvenanceBlock: return "missing_provenance_block";
    case LineageStatus::MalformedProvenance: return "malformed_provenance";
    case LineageStatus::UnsupportedSchema: return "unsupported_schema";
    case LineageStatus::ArtifactTypeMismatch: return "artifact_type_mismatch";
    }
    return "malformed_provenance";
}

const char* artifact_type_token (ArtifactType type) noexcept
{
    switch (type) {
    case ArtifactType::Checkpoint: return "checkpoint";
    case ArtifactType::Plotfile2D: return "plotfile_2d";
    case ArtifactType::Plotfile3D: return "plotfile_3d";
    }
    return "checkpoint";
}

std::string serialize_provenance_block (const ProvenanceRecord& record)
{
    std::ostringstream out;
    out << provenance_begin << '\n';
    out << "schema_version=" << record.execution.schema_version << '\n';
    out << "simulation_uuid=" << record.execution.simulation_uuid << '\n';
    out << "execution_uuid=" << record.execution.execution_uuid << '\n';
    out << "parent_execution_uuid=" << record.execution.parent_execution_uuid << '\n';
    out << "execution_lineage=";
    for (std::size_t i = 0; i < record.execution.execution_lineage.size(); ++i) {
        if (i != 0) out << ',';
        out << record.execution.execution_lineage[i];
    }
    out << '\n';
    out << "restart_generation=" << record.execution.restart_generation << '\n';
    out << "lineage_complete=" << (record.execution.lineage_complete ? "true" : "false") << '\n';
    out << "lineage_status=" << lineage_status_token(record.execution.lineage_status) << '\n';
    out << "source_checkpoint_uuid=" << record.execution.source_checkpoint_uuid << '\n';
    out << "source_checkpoint_path=" << record.execution.source_checkpoint_path << '\n';
    out << "execution_start_utc=" << record.execution.execution_start_utc << '\n';
    out << "artifact_uuid=" << record.artifact.artifact_uuid << '\n';
    out << "artifact_type=" << artifact_type_token(record.artifact.artifact_type) << '\n';
    out << "artifact_step=" << record.artifact.artifact_step << '\n';
    out << "artifact_time_seconds=" << std::setprecision(17)
        << record.artifact.artifact_time_seconds << '\n';
    out << "artifact_created_utc=" << record.artifact.artifact_created_utc << '\n';
    out << provenance_end << '\n';
    return out.str();
}

ProvenanceParseResult parse_provenance_block (std::string_view text)
{
    enum class BlockState { BeforeBlock, InsideBlock, AfterBlock };
    BlockState state = BlockState::BeforeBlock;
    std::vector<std::string> block_lines;
    std::istringstream input{std::string(text)};
    std::string line;
    while (std::getline(input, line)) {
        line = strip_cr(std::move(line));
        const bool is_begin = line == provenance_begin;
        const bool is_end = line == provenance_end;
        if (state == BlockState::BeforeBlock) {
            if (is_end) {
                return failure(ProvenanceReadStatus::MalformedProvenance,
                               "provenance end marker precedes begin marker");
            }
            if (is_begin) state = BlockState::InsideBlock;
        } else if (state == BlockState::InsideBlock) {
            if (is_begin) {
                return failure(ProvenanceReadStatus::MalformedProvenance,
                               "provenance block has a nested begin marker");
            }
            if (is_end) {
                state = BlockState::AfterBlock;
            } else {
                block_lines.push_back(line);
            }
        } else if (is_begin || is_end) {
            return failure(ProvenanceReadStatus::MalformedProvenance,
                           "job_info has multiple provenance blocks");
        }
    }
    if (state == BlockState::BeforeBlock) {
        return failure(ProvenanceReadStatus::MissingProvenanceBlock,
                       "job_info has no exact ERF provenance block");
    }
    if (state == BlockState::InsideBlock) {
        return failure(ProvenanceReadStatus::MalformedProvenance,
                       "provenance block has no matching end marker");
    }

    std::unordered_map<std::string, std::string> fields;
    static const std::unordered_set<std::string> required = {
        "schema_version", "simulation_uuid", "execution_uuid", "parent_execution_uuid",
        "execution_lineage", "restart_generation", "lineage_complete", "lineage_status",
        "source_checkpoint_uuid", "source_checkpoint_path", "execution_start_utc",
        "artifact_uuid", "artifact_type", "artifact_step", "artifact_time_seconds",
        "artifact_created_utc"
    };
    for (const auto& block_line : block_lines) {
        if (block_line.empty()) continue;
        const std::size_t equals = block_line.find('=');
        if (equals == std::string::npos || equals == 0) {
            return failure(ProvenanceReadStatus::MalformedProvenance,
                           "provenance line is not a key=value pair");
        }
        const std::string key = block_line.substr(0, equals);
        if (required.count(key) != 0 && fields.count(key) != 0) {
            return failure(ProvenanceReadStatus::MalformedProvenance,
                           "duplicate required provenance key: " + key);
        }
        if (required.count(key) != 0) fields.emplace(key, block_line.substr(equals + 1));
    }
    ProvenanceParseResult result;
    int schema = 0;
    if (fields.count("schema_version") == 0 ||
        !parse_int(fields.at("schema_version"), schema)) {
        return failure(ProvenanceReadStatus::MalformedProvenance, "invalid provenance schema_version");
    }
    if (schema != erf_provenance::schema_version) {
        return failure(ProvenanceReadStatus::UnsupportedSchema, "unsupported provenance schema_version");
    }
    for (const auto& key : required) {
        if (fields.count(key) == 0) {
            return failure(ProvenanceReadStatus::MalformedProvenance,
                           "missing required provenance key: " + key);
        }
    }
    result.record.execution.schema_version = schema;
    result.record.execution.simulation_uuid = fields.at("simulation_uuid");
    result.record.execution.execution_uuid = fields.at("execution_uuid");
    result.record.execution.parent_execution_uuid = fields.at("parent_execution_uuid");
    result.record.execution.execution_lineage = split_lineage(fields.at("execution_lineage"));
    if (result.record.execution.execution_lineage.empty()) {
        return failure(ProvenanceReadStatus::MalformedProvenance, "invalid execution_lineage");
    }
    if (!parse_int(fields.at("restart_generation"), result.record.execution.restart_generation) ||
        !parse_bool(fields.at("lineage_complete"), result.record.execution.lineage_complete) ||
        !parse_lineage_status(fields.at("lineage_status"), result.record.execution.lineage_status) ||
        !parse_artifact_type(fields.at("artifact_type"), result.record.artifact.artifact_type) ||
        !parse_int(fields.at("artifact_step"), result.record.artifact.artifact_step) ||
        !parse_double(fields.at("artifact_time_seconds"), result.record.artifact.artifact_time_seconds)) {
        return failure(ProvenanceReadStatus::MalformedProvenance, "invalid typed provenance field");
    }
    result.record.execution.source_checkpoint_uuid = fields.at("source_checkpoint_uuid");
    result.record.execution.source_checkpoint_path = fields.at("source_checkpoint_path");
    result.record.execution.execution_start_utc = fields.at("execution_start_utc");
    result.record.artifact.artifact_uuid = fields.at("artifact_uuid");
    result.record.artifact.artifact_created_utc = fields.at("artifact_created_utc");

    std::string diagnostic;
    if (!validate_record(result.record, diagnostic)) {
        return failure(ProvenanceReadStatus::MalformedProvenance, std::move(diagnostic));
    }
    result.status = ProvenanceReadStatus::Valid;
    return result;
}

ProvenanceReadResult read_job_info_file (const std::string& path)
{
    int read_status = static_cast<int>(ProvenanceReadStatus::Valid);
    std::string contents;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream input(path, std::ios::in | std::ios::binary);
        if (!input.good()) {
            read_status = static_cast<int>(ProvenanceReadStatus::MissingJobInfo);
        } else {
            contents.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
            if (!input.good() && !input.eof()) {
                read_status = static_cast<int>(ProvenanceReadStatus::MalformedProvenance);
            }
        }
    }
    const int io_rank = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&read_status, 1, io_rank);
    broadcast_string(contents);
    if (read_status != static_cast<int>(ProvenanceReadStatus::Valid)) {
        return failure(static_cast<ProvenanceReadStatus>(read_status),
                       read_status == static_cast<int>(ProvenanceReadStatus::MissingJobInfo)
                           ? "job_info does not exist" : "job_info could not be read");
    }
    return parse_provenance_block(contents);
}

ExecutionProvenance make_cold_start_provenance (const std::string& execution_uuid,
                                                const std::string& execution_start_utc)
{
    ExecutionProvenance result;
    result.simulation_uuid = execution_uuid;
    result.execution_uuid = execution_uuid;
    result.execution_lineage = {execution_uuid};
    result.execution_start_utc = execution_start_utc;
    return result;
}

ExecutionProvenance make_restart_provenance (const ExecutionProvenance& current_invocation,
                                             const ProvenanceRecord& parent_checkpoint,
                                             const std::string& checkpoint_path)
{
    ExecutionProvenance result = current_invocation;
    const auto& parent = parent_checkpoint.execution;
    result.simulation_uuid = parent.simulation_uuid;
    result.parent_execution_uuid = parent.execution_uuid;
    result.execution_lineage = parent.execution_lineage;
    result.execution_lineage.push_back(current_invocation.execution_uuid);
    result.source_checkpoint_uuid = parent_checkpoint.artifact.artifact_uuid;
    result.source_checkpoint_path = checkpoint_path;
    if (parent.lineage_complete) {
        result.restart_generation = parent.restart_generation + 1;
        result.lineage_complete = true;
        result.lineage_status = LineageStatus::Complete;
    } else {
        result.restart_generation = -1;
        result.lineage_complete = false;
        result.lineage_status = LineageStatus::IncompleteAncestor;
    }
    return result;
}

ExecutionProvenance make_incomplete_restart_provenance (const ExecutionProvenance& current_invocation,
                                                         ProvenanceReadStatus failure_status,
                                                         const std::string& checkpoint_path)
{
    ExecutionProvenance result = current_invocation;
    result.simulation_uuid = current_invocation.execution_uuid;
    result.parent_execution_uuid.clear();
    result.execution_lineage = {current_invocation.execution_uuid};
    result.restart_generation = -1;
    result.lineage_complete = false;
    result.lineage_status = failure_lineage_status(failure_status);
    result.source_checkpoint_uuid.clear();
    result.source_checkpoint_path = checkpoint_path;
    return result;
}

ExecutionProvenance initialize_execution_provenance ()
{
    static bool initialized = false;
    static ExecutionProvenance invocation;
    if (!initialized) {
        std::string execution_uuid;
        std::string execution_start_utc;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            execution_uuid = generate_uuid_v4();
            execution_start_utc = current_utc();
        }
        broadcast_string(execution_uuid);
        broadcast_string(execution_start_utc);
        invocation = make_cold_start_provenance(execution_uuid, execution_start_utc);
        initialized = true;
    }
    return invocation;
}

} // namespace erf_provenance
