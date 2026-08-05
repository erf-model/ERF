/**
 * \file ERF_Provenance.H
 */
#ifndef ERF_PROVENANCE_H_
#define ERF_PROVENANCE_H_

#include <array>
#include <cstdint>
#include <ctime>
#include <string>
#include <string_view>
#include <vector>

namespace erf_provenance
{

constexpr int schema_version = 1;
constexpr const char* provenance_begin = "ERF_PROVENANCE_BEGIN";
constexpr const char* provenance_end = "ERF_PROVENANCE_END";

/**
 * The reason a lineage is incomplete is part of the persisted contract.
 * Missing or invalid metadata never prevents ERF from reading physical state.
 */
enum class LineageStatus
{
    Complete,
    IncompleteAncestor,
    MissingJobInfo,
    MissingProvenanceBlock,
    MalformedProvenance,
    UnsupportedSchema,
    ArtifactTypeMismatch
};

enum class ArtifactType
{
    Checkpoint,
    Plotfile2D,
    Plotfile3D
};

enum class ProvenanceReadStatus
{
    Valid,
    MissingJobInfo,
    MissingProvenanceBlock,
    MalformedProvenance,
    UnsupportedSchema,
    ArtifactTypeMismatch
};

/**
 * Identity and known ancestry for one ERF process invocation.
 *
 * Simulation identity names one rooted restart family, which may branch.
 * Execution identity names one invocation. Artifact identity is kept separate
 * because one invocation can write many checkpoints and plotfiles. The known
 * ancestry path is stored in every artifact so old directories are not needed
 * to reconstruct repeated restarts. The parent UUID gives the immediate edge;
 * the source checkpoint UUID selects the exact checkpoint artifact. These UUIDs
 * identify records; they do not authenticate contents or provide cryptographic
 * proof.
 */
struct ExecutionProvenance
{
    int schema_version = erf_provenance::schema_version;
    std::string simulation_uuid;
    std::string execution_uuid;
    std::string parent_execution_uuid;
    std::vector<std::string> execution_lineage;
    int restart_generation = 0;
    bool lineage_complete = true;
    LineageStatus lineage_status = LineageStatus::Complete;
    std::string source_checkpoint_uuid;
    std::string source_checkpoint_path;
    std::string execution_start_utc;
};

struct ArtifactProvenance
{
    std::string artifact_uuid;
    ArtifactType artifact_type = ArtifactType::Checkpoint;
    int artifact_step = -1;
    double artifact_time_seconds = 0.0;
    std::string artifact_created_utc;
};

struct ProvenanceRecord
{
    ExecutionProvenance execution;
    ArtifactProvenance artifact;
};

struct ProvenanceParseResult
{
    ProvenanceReadStatus status = ProvenanceReadStatus::MissingProvenanceBlock;
    ProvenanceRecord record;
    std::string diagnostic;

    [[nodiscard]] bool valid () const noexcept
    {
        return status == ProvenanceReadStatus::Valid;
    }
};

using ProvenanceReadResult = ProvenanceParseResult;

/** Format supplied bytes as a lowercase canonical RFC UUIDv4 string. */
std::string uuid_v4_from_bytes (const std::array<std::uint8_t, 16>& bytes);

/** Validate a lowercase or uppercase canonical RFC UUIDv4 string. */
bool is_valid_uuid_v4 (std::string_view uuid);

/** Generate a UUID without using AMReX's scientific random-number stream. */
std::string generate_uuid_v4 ();

/** Format a supplied clock value as an ISO 8601 UTC timestamp. */
std::string format_utc (std::time_t time_value);

/** Return the current time in the persisted UTC format. */
std::string current_utc ();

/** Serialize exactly one schema-1 machine-readable block. */
std::string serialize_provenance_block (const ProvenanceRecord& record);

/** Parse only delimited provenance markers from a human-readable job_info. */
ProvenanceParseResult parse_provenance_block (std::string_view job_info_text);

/** Read and broadcast one job_info before parsing the same bytes on every rank. */
ProvenanceReadResult read_job_info_file (const std::string& job_info_path);

/** Create the root record for one cold-start invocation. */
ExecutionProvenance make_cold_start_provenance (const std::string& execution_uuid,
                                                const std::string& execution_start_utc);

/** Create a child record from a valid checkpoint record, including its artifact UUID. */
ExecutionProvenance make_restart_provenance (const ExecutionProvenance& current_invocation,
                                             const ProvenanceRecord& parent_checkpoint,
                                             const std::string& checkpoint_path);

/** Start a new known lineage segment after a nonfatal provenance read failure. */
ExecutionProvenance make_incomplete_restart_provenance (const ExecutionProvenance& current_invocation,
                                                         ProvenanceReadStatus failure,
                                                         const std::string& checkpoint_path);

/** Return the stable persisted token for a lineage status. */
const char* lineage_status_token (LineageStatus status) noexcept;

/** Return the stable persisted token for an artifact type. */
const char* artifact_type_token (ArtifactType type) noexcept;

/** Initialize one process-wide invocation record and share it across multiblock ERF objects. */
ExecutionProvenance initialize_execution_provenance ();

} // namespace erf_provenance

#endif
