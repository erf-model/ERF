#ifndef ERF_SBM_RESTART_H_
#define ERF_SBM_RESTART_H_

#include "ERF_SBMBulkProjection.H"

#include <AMReX_MultiFab.H>

#include <string>

namespace erf_sbm {

/** Exact M1 semantic schema; this deliberately identifies the bounded fixture. */
[[nodiscard]] std::string restart_schema(const SBMLayout& layout);

/** Compare persisted and expected schemas byte-for-byte. */
[[nodiscard]] bool restart_schema_matches(const SBMLayout& layout,
                                          const std::string& persisted);

/** Validate every valid-cell value against the layout's runtime constraints. */
[[nodiscard]] bool authoritative_state_admissible(
    const amrex::MultiFab& spectrum, const SBMLayout& layout, int level,
    std::string* diagnostic = nullptr);

/** Compare the projected spectrum to checkpointed qc/qr without mutating core. */
[[nodiscard]] bool restart_projection_matches(const amrex::MultiFab& spectrum,
                                               const amrex::MultiFab& persisted_core,
                                               const SBMBulkProjection& projection,
                                               int qc_component, int qr_component,
                                               amrex::Real tolerance_scale = amrex::Real(128.0));

} // namespace erf_sbm

#endif
