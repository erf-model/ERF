#include "ERF_Provenance.H"

#include <AMReX_ParallelDescriptor.H>
#include <gtest/gtest.h>

TEST(ERFProvenanceParallel, BroadcastExecutionIdentity)
{
    // Motivation: Every MPI rank must attach the same execution identity and lineage to collectively written ERF outputs.
    const auto execution = erf_provenance::initialize_execution_provenance();
    int local_equal = erf_provenance::is_valid_uuid_v4(execution.execution_uuid) &&
                      execution.simulation_uuid == execution.execution_uuid &&
                      execution.execution_lineage.size() == 1u ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntSum(local_equal);
    EXPECT_EQ(local_equal, amrex::ParallelDescriptor::NProcs());
}
