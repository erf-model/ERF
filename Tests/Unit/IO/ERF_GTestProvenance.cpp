#include "ERF_Provenance.H"

#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace
{

using namespace erf_provenance;

const std::string uuid_a = "11111111-1111-4111-8111-111111111111";
const std::string uuid_b = "22222222-2222-4222-8222-222222222222";
const std::string uuid_c = "33333333-3333-4333-8333-333333333333";
const std::string uuid_d = "44444444-4444-4444-8444-444444444444";
const std::string artifact_a = "aaaaaaaa-aaaa-4aaa-8aaa-aaaaaaaaaaaa";
const std::string artifact_b = "bbbbbbbb-bbbb-4bbb-8bbb-bbbbbbbbbbbb";
const std::string artifact_c = "cccccccc-cccc-4ccc-8ccc-cccccccccccc";

ProvenanceRecord checkpoint_record (const ExecutionProvenance& execution,
                                    const std::string& artifact_uuid,
                                    int step = 4)
{
    ProvenanceRecord record;
    record.execution = execution;
    record.artifact.artifact_uuid = artifact_uuid;
    record.artifact.artifact_type = ArtifactType::Checkpoint;
    record.artifact.artifact_step = step;
    record.artifact.artifact_time_seconds = 12.5;
    record.artifact.artifact_created_utc = "2026-07-11T18:42:17Z";
    return record;
}

std::string replace_once (std::string text, const std::string& from, const std::string& to)
{
    const std::size_t position = text.find(from);
    EXPECT_NE(position, std::string::npos);
    if (position != std::string::npos) text.replace(position, from.size(), to);
    return text;
}

} // namespace

TEST(ERFProvenance, UUIDFormatting)
{
    // Motivation: Provenance identifiers must use one canonical UUIDv4 form so checkpoints and analysis tools can validate and compare them reliably.
    const std::array<std::uint8_t, 16> bytes{
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x4a, 0x77,
        0x8f, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
    const std::string uuid = uuid_v4_from_bytes(bytes);
    EXPECT_EQ(uuid, "00112233-4455-4a77-8f99-aabbccddeeff");
    EXPECT_EQ(uuid.size(), 36u);
    EXPECT_EQ(uuid[8], '-');
    EXPECT_EQ(uuid[13], '-');
    EXPECT_EQ(uuid[18], '-');
    EXPECT_EQ(uuid[23], '-');
}

TEST(ERFProvenance, UUIDValidation)
{
    // Motivation: Restart metadata must reject malformed identifiers without rejecting the physical checkpoint that contains them.
    EXPECT_TRUE(is_valid_uuid_v4(uuid_a));
    EXPECT_TRUE(is_valid_uuid_v4("AAAAAAAA-AAAA-4AAA-8AAA-AAAAAAAAAAAA"));
    EXPECT_FALSE(is_valid_uuid_v4("11111111-1111-4111-8111-11111111111"));
    EXPECT_FALSE(is_valid_uuid_v4("111111111111-4111-8111-111111111111"));
    EXPECT_FALSE(is_valid_uuid_v4("11111111-1111-5111-8111-111111111111"));
    EXPECT_FALSE(is_valid_uuid_v4("11111111-1111-4111-7111-111111111111"));
    EXPECT_FALSE(is_valid_uuid_v4("11111111-1111-4111-8111-11111111111g"));
}

TEST(ERFProvenance, ColdStart)
{
    // Motivation: A cold start defines the root of a simulation lineage and must not claim a parent or a restart generation.
    const auto execution = make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z");
    EXPECT_EQ(execution.simulation_uuid, uuid_a);
    EXPECT_EQ(execution.execution_uuid, uuid_a);
    EXPECT_TRUE(execution.parent_execution_uuid.empty());
    ASSERT_EQ(execution.execution_lineage.size(), 1u);
    EXPECT_EQ(execution.execution_lineage.front(), uuid_a);
    EXPECT_EQ(execution.restart_generation, 0);
    EXPECT_TRUE(execution.lineage_complete);
    EXPECT_EQ(execution.lineage_status, LineageStatus::Complete);
    EXPECT_TRUE(execution.source_checkpoint_uuid.empty());
    EXPECT_TRUE(execution.source_checkpoint_path.empty());
}

TEST(ERFProvenance, CompleteRestart)
{
    // Motivation: A restart must preserve the simulation identity, identify its parent execution and source checkpoint, and append exactly one execution.
    const auto parent = make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z");
    const auto child = make_restart_provenance(make_cold_start_provenance(uuid_b, "2026-07-11T19:00:00Z"),
                                               checkpoint_record(parent, artifact_a), "/tmp/chk00004");
    EXPECT_EQ(child.simulation_uuid, uuid_a);
    EXPECT_EQ(child.execution_uuid, uuid_b);
    EXPECT_EQ(child.parent_execution_uuid, uuid_a);
    EXPECT_EQ(child.execution_lineage, (std::vector<std::string>{uuid_a, uuid_b}));
    EXPECT_EQ(child.restart_generation, 1);
    EXPECT_TRUE(child.lineage_complete);
    EXPECT_EQ(child.lineage_status, LineageStatus::Complete);
    EXPECT_EQ(child.source_checkpoint_uuid, artifact_a);
    EXPECT_EQ(child.source_checkpoint_path, "/tmp/chk00004");
}

TEST(ERFProvenance, MultipleRestartGenerations)
{
    // Motivation: Full lineage must survive repeated checkpoint restarts without requiring older checkpoint directories to remain available.
    const auto a = make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z");
    const auto b = make_restart_provenance(make_cold_start_provenance(uuid_b, "2026-07-11T19:00:00Z"),
                                           checkpoint_record(a, artifact_a), "/tmp/chk-a");
    const auto c = make_restart_provenance(make_cold_start_provenance(uuid_c, "2026-07-11T20:00:00Z"),
                                           checkpoint_record(b, artifact_b), "/tmp/chk-b");
    EXPECT_EQ(c.execution_lineage, (std::vector<std::string>{uuid_a, uuid_b, uuid_c}));
    EXPECT_EQ(c.restart_generation, 2);
    EXPECT_EQ(c.source_checkpoint_uuid, artifact_b);
}

TEST(ERFProvenance, BranchingRestarts)
{
    // Motivation: Two executions may restart from the same checkpoint, so they must share a parent and simulation UUID while retaining distinct identities.
    const auto a = make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z");
    const auto b = make_restart_provenance(make_cold_start_provenance(uuid_b, "2026-07-11T19:00:00Z"),
                                           checkpoint_record(a, artifact_a), "/tmp/chk-a");
    const auto c = make_restart_provenance(make_cold_start_provenance(uuid_c, "2026-07-11T20:00:00Z"),
                                           checkpoint_record(b, artifact_b), "/tmp/chk-b");
    const auto d = make_restart_provenance(make_cold_start_provenance(uuid_d, "2026-07-11T20:30:00Z"),
                                           checkpoint_record(b, artifact_b), "/tmp/chk-b");
    EXPECT_EQ(c.simulation_uuid, d.simulation_uuid);
    EXPECT_EQ(c.parent_execution_uuid, d.parent_execution_uuid);
    EXPECT_NE(c.execution_uuid, d.execution_uuid);
    EXPECT_EQ(c.execution_lineage.back(), uuid_c);
    EXPECT_EQ(d.execution_lineage.back(), uuid_d);
}

TEST(ERFProvenance, IncompleteLegacyLineage)
{
    // Motivation: Legacy checkpoints must remain restartable while clearly marking that ancestry before the current execution is unknown.
    const auto child = make_incomplete_restart_provenance(
        make_cold_start_provenance(uuid_b, "2026-07-11T19:00:00Z"),
        ProvenanceReadStatus::MissingJobInfo, "/tmp/legacy");
    EXPECT_EQ(child.simulation_uuid, uuid_b);
    EXPECT_EQ(child.execution_lineage, (std::vector<std::string>{uuid_b}));
    EXPECT_EQ(child.restart_generation, -1);
    EXPECT_FALSE(child.lineage_complete);
    EXPECT_EQ(child.lineage_status, LineageStatus::MissingJobInfo);
    EXPECT_TRUE(child.source_checkpoint_uuid.empty());
    EXPECT_EQ(child.source_checkpoint_path, "/tmp/legacy");
}

TEST(ERFProvenance, IncompleteLineagePropagation)
{
    // Motivation: Once ancestry is incomplete, later restarts must preserve the known suffix without falsely claiming that the full history was recovered.
    auto parent_execution = make_incomplete_restart_provenance(
        make_cold_start_provenance(uuid_b, "2026-07-11T19:00:00Z"),
        ProvenanceReadStatus::MissingProvenanceBlock, "/tmp/legacy");
    auto parent = checkpoint_record(parent_execution, artifact_a);
    const auto child = make_restart_provenance(make_cold_start_provenance(uuid_c, "2026-07-11T20:00:00Z"),
                                               parent, "/tmp/chk-b");
    EXPECT_EQ(child.simulation_uuid, uuid_b);
    EXPECT_EQ(child.execution_lineage, (std::vector<std::string>{uuid_b, uuid_c}));
    EXPECT_EQ(child.parent_execution_uuid, uuid_b);
    EXPECT_EQ(child.restart_generation, -1);
    EXPECT_FALSE(child.lineage_complete);
    EXPECT_EQ(child.lineage_status, LineageStatus::IncompleteAncestor);
    EXPECT_EQ(child.source_checkpoint_uuid, artifact_a);
}

TEST(ERFProvenance, SerializationRoundTrip)
{
    // Motivation: The checkpoint reader and job_info writer must agree on every schema-1 field so a written checkpoint can reconstruct its lineage.
    auto execution = make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z");
    execution.source_checkpoint_path = "/tmp/checkpoint=a";
    const auto original = checkpoint_record(execution, artifact_a);
    const auto parsed = parse_provenance_block(serialize_provenance_block(original));
    ASSERT_TRUE(parsed.valid()) << parsed.diagnostic;
    EXPECT_EQ(serialize_provenance_block(parsed.record), serialize_provenance_block(original));
}

TEST(ERFProvenance, HumanTextAroundBlock)
{
    // Motivation: job_info contains build and input text, so the parser must read only the delimited provenance block.
    const auto record = checkpoint_record(make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z"), artifact_a);
    const std::string text = "input=value=with=equals UUID-like=" + uuid_d + "\n" +
                             serialize_provenance_block(record) +
                             "build hash=" + uuid_c + "\n";
    const auto parsed = parse_provenance_block(text);
    ASSERT_TRUE(parsed.valid()) << parsed.diagnostic;
    EXPECT_EQ(parsed.record.artifact.artifact_uuid, artifact_a);
}

TEST(ERFProvenance, MissingBlock)
{
    // Motivation: Checkpoints written before provenance support contain no block and must produce a nonfatal, explicit status.
    const auto parsed = parse_provenance_block("old job_info\ninputs=foo\n");
    EXPECT_EQ(parsed.status, ProvenanceReadStatus::MissingProvenanceBlock);
    EXPECT_FALSE(parsed.valid());
}

TEST(ERFProvenance, MalformedAndDuplicateFields)
{
    // Motivation: Ambiguous provenance must be rejected rather than silently constructing a false lineage.
    const auto serialized = serialize_provenance_block(
        checkpoint_record(make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z"), artifact_a));
    EXPECT_EQ(parse_provenance_block(replace_once(serialized, "artifact_step=4\n", "artifact_step=bad\n")).status,
              ProvenanceReadStatus::MalformedProvenance);
    EXPECT_EQ(parse_provenance_block(replace_once(serialized, "artifact_uuid=" + artifact_a + "\n",
                                                  "artifact_uuid=" + artifact_a + "\nartifact_uuid=" + artifact_b + "\n")).status,
              ProvenanceReadStatus::MalformedProvenance);
    EXPECT_EQ(parse_provenance_block(replace_once(serialized, "lineage_complete=true", "lineage_complete=TRUE")).status,
              ProvenanceReadStatus::MalformedProvenance);
    EXPECT_EQ(parse_provenance_block(replace_once(serialized, "execution_uuid=" + uuid_a,
                                                  "execution_uuid=not-a-uuid")).status,
              ProvenanceReadStatus::MalformedProvenance);
    EXPECT_EQ(parse_provenance_block(replace_once(serialized, "parent_execution_uuid=\n",
                                                  "parent_execution_uuid=" + uuid_b + "\n")).status,
              ProvenanceReadStatus::MalformedProvenance);
    EXPECT_EQ(parse_provenance_block(serialized + serialized).status,
              ProvenanceReadStatus::MalformedProvenance);
}

TEST(ERFProvenance, UnsupportedSchema)
{
    // Motivation: New ERF versions must not make older physical checkpoint data unreadable when provenance schema support differs.
    auto serialized = serialize_provenance_block(
        checkpoint_record(make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z"), artifact_a));
    serialized = replace_once(serialized, "schema_version=1", "schema_version=99");
    const auto parsed = parse_provenance_block(serialized);
    EXPECT_EQ(parsed.status, ProvenanceReadStatus::UnsupportedSchema);
    EXPECT_FALSE(parsed.valid());
}

TEST(ERFProvenance, ArtifactIdentity)
{
    // Motivation: Multiple outputs from one execution must share execution provenance while each output keeps a distinct artifact identity.
    auto execution = make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z");
    const auto checkpoint = checkpoint_record(execution, artifact_a);
    const auto plotfile = checkpoint_record(execution, artifact_b);
    EXPECT_EQ(checkpoint.execution.execution_uuid, plotfile.execution.execution_uuid);
    EXPECT_NE(checkpoint.artifact.artifact_uuid, plotfile.artifact.artifact_uuid);
}

TEST(ERFProvenance, CRLFAndEqualsInPath)
{
    // Motivation: Provenance parsing must tolerate common line endings and paths that contain equals signs without corrupting checkpoint identity.
    auto execution = make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z");
    execution.source_checkpoint_path = "/scratch/run=a/chk00004";
    const auto serialized = serialize_provenance_block(checkpoint_record(execution, artifact_a));
    std::string crlf;
    for (const char character : serialized) {
        if (character == '\n') crlf += '\r';
        crlf += character;
    }
    const auto parsed = parse_provenance_block(crlf);
    ASSERT_TRUE(parsed.valid()) << parsed.diagnostic;
    EXPECT_EQ(parsed.record.execution.source_checkpoint_path, "/scratch/run=a/chk00004");
}

TEST(ERFProvenance, FileReadStatuses)
{
    // Motivation: File-level restart metadata failures must return structured nonfatal statuses without relying on an abort or death test.
    const auto directory = std::filesystem::temp_directory_path() / "erf-provenance-file-test";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    const auto path = directory / "job_info";
    const auto missing = read_job_info_file(path.string());
    EXPECT_EQ(missing.status, ProvenanceReadStatus::MissingJobInfo);
    {
        std::ofstream output(path);
        output << "legacy text\n";
    }
    EXPECT_EQ(read_job_info_file(path.string()).status, ProvenanceReadStatus::MissingProvenanceBlock);
    {
        std::ofstream output(path, std::ios::trunc);
        output << "ERF_PROVENANCE_BEGIN\nartifact_step=bad\nERF_PROVENANCE_END\n";
    }
    EXPECT_EQ(read_job_info_file(path.string()).status, ProvenanceReadStatus::MalformedProvenance);
    {
        std::ofstream output(path, std::ios::trunc);
        auto block = serialize_provenance_block(
            checkpoint_record(make_cold_start_provenance(uuid_a, "2026-07-11T18:42:17Z"), artifact_a));
        output << replace_once(block, "schema_version=1", "schema_version=88");
    }
    EXPECT_EQ(read_job_info_file(path.string()).status, ProvenanceReadStatus::UnsupportedSchema);
    std::filesystem::remove_all(directory);
}
