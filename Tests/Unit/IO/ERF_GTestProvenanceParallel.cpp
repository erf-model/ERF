#include "ERF_Provenance.H"

#include <AMReX_ParallelDescriptor.H>
#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace
{

void broadcast_reference_string (std::string& value)
{
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    int length = amrex::ParallelDescriptor::IOProcessor()
               ? static_cast<int>(value.size()) : 0;
    amrex::ParallelDescriptor::Bcast(&length, 1, root);
    value.resize(static_cast<std::size_t>(length));
    if (length > 0) {
        amrex::ParallelDescriptor::Bcast(value.data(), length, root);
    }
}

} // namespace

TEST(ERFProvenanceParallel, BroadcastExecutionIdentity)
{
    // Motivation: Every MPI rank must attach the same execution identity and
    // lineage to collectively written ERF outputs. Valid but different
    // per-rank UUIDs must fail this test.
    const auto execution = erf_provenance::initialize_execution_provenance();
    const bool is_io = amrex::ParallelDescriptor::IOProcessor();
    std::string reference_execution_uuid = is_io ? execution.execution_uuid : "";
    std::string reference_simulation_uuid = is_io ? execution.simulation_uuid : "";
    std::string reference_start_utc = is_io ? execution.execution_start_utc : "";
    std::string reference_parent_uuid = is_io ? execution.parent_execution_uuid : "";
    std::string reference_lineage_status = is_io
        ? erf_provenance::lineage_status_token(execution.lineage_status) : "";
    broadcast_reference_string(reference_execution_uuid);
    broadcast_reference_string(reference_simulation_uuid);
    broadcast_reference_string(reference_start_utc);
    broadcast_reference_string(reference_parent_uuid);
    broadcast_reference_string(reference_lineage_status);

    int reference_generation = is_io ? execution.restart_generation : 0;
    int reference_lineage_complete = is_io && execution.lineage_complete ? 1 : 0;
    int reference_lineage_size = is_io
        ? static_cast<int>(execution.execution_lineage.size()) : 0;
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&reference_generation, 1, root);
    amrex::ParallelDescriptor::Bcast(&reference_lineage_complete, 1, root);
    amrex::ParallelDescriptor::Bcast(&reference_lineage_size, 1, root);
    std::vector<std::string> reference_lineage(static_cast<std::size_t>(reference_lineage_size));
    for (int i = 0; i < reference_lineage_size; ++i) {
        reference_lineage[static_cast<std::size_t>(i)] = is_io
            ? execution.execution_lineage[static_cast<std::size_t>(i)] : "";
        broadcast_reference_string(reference_lineage[static_cast<std::size_t>(i)]);
    }

    bool lineage_equal = execution.execution_lineage.size() == reference_lineage.size();
    if (lineage_equal) {
        for (std::size_t i = 0; i < reference_lineage.size(); ++i) {
            lineage_equal = lineage_equal &&
                execution.execution_lineage[i] == reference_lineage[i];
        }
    }
    int local_equal = erf_provenance::is_valid_uuid_v4(execution.execution_uuid) &&
                      execution.execution_uuid == reference_execution_uuid &&
                      execution.simulation_uuid == reference_simulation_uuid &&
                      execution.execution_start_utc == reference_start_utc &&
                      execution.parent_execution_uuid == reference_parent_uuid &&
                      lineage_equal &&
                      execution.restart_generation == reference_generation &&
                      (execution.lineage_complete ? 1 : 0) == reference_lineage_complete &&
                      std::string(erf_provenance::lineage_status_token(execution.lineage_status)) ==
                          reference_lineage_status &&
                      execution.simulation_uuid == execution.execution_uuid &&
                      execution.parent_execution_uuid.empty() &&
                      execution.execution_lineage == std::vector<std::string>{execution.execution_uuid} &&
                      execution.restart_generation == 0 &&
                      execution.lineage_complete &&
                      execution.lineage_status == erf_provenance::LineageStatus::Complete ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntSum(local_equal);
    EXPECT_EQ(local_equal, amrex::ParallelDescriptor::NProcs());
}
