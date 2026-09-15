#ifndef ERF_SBM_RESTART_H_
#define ERF_SBM_RESTART_H_

#include "ERF_SBMLayout.H"

#include <AMReX_REAL.H>

#include <string>

namespace erf_sbm {

struct SBMCheckpointSchema {
    std::string schema_version;
    std::string population_ids;
    std::string moment_modes;
    std::string grid_identity;
    std::string property_identity;
    std::string constraint_policy;
    std::string projection_identity;
    std::string transport_identity;
    std::string numerical_policy;
};

[[nodiscard]] SBMCheckpointSchema
make_checkpoint_schema(const SBMLayout& layout,
                       const std::string& constraint_policy,
                       const std::string& transport_identity,
                       const std::string& numerical_policy);

[[nodiscard]] std::string compare_checkpoint_schema(const SBMCheckpointSchema& expected,
                                                    const SBMCheckpointSchema& actual);

void write_checkpoint_schema(const std::string& path, const SBMCheckpointSchema& schema);
[[nodiscard]] SBMCheckpointSchema read_checkpoint_schema(const std::string& path);

/** Scale-aware equality used for checkpointed bulk projections. */
[[nodiscard]] bool compare_projection(amrex::Real checkpointed,
                                      amrex::Real reconstructed,
                                      amrex::Real scale,
                                      int operation_count,
                                      amrex::Real* absolute_error = nullptr) noexcept;

} // namespace erf_sbm

#endif
