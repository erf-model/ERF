#ifndef ERF_SBM_RESTART_H_
#define ERF_SBM_RESTART_H_

#include "ERF_SBMBulkProjection.H"

#include <AMReX_MultiFab.H>

#include <cstdint>
#include <memory>
#include <string>

namespace erf_sbm {

/** Exact semantic restart schema; no implicit checkpoint migration is performed. */
[[nodiscard]] std::string restart_schema (const SBMLayout& layout);

/** Compare persisted and expected schemas byte-for-byte. */
[[nodiscard]] bool restart_schema_matches (const SBMLayout& layout,
                                           const std::string& persisted);

/** Validate every valid-cell value against the layout's runtime constraints. */
[[nodiscard]] bool authoritative_state_admissible (
    const amrex::MultiFab& spectrum, const SBMLayout& layout, int level,
    std::string* diagnostic = nullptr);

/** Persistent device-side equivalent of the M3 full-state admission gate.
 * Define one scalar mask per active level layout and reuse all semantic
 * descriptors across candidate checks. The detailed host oracle is invoked
 * only after a device rejection.
 */
class SBMStateAdmission {
public:
    explicit SBMStateAdmission (const SBMLayout& layout);
    ~SBMStateAdmission ();
    SBMStateAdmission (const SBMStateAdmission&) = delete;
    SBMStateAdmission& operator= (const SBMStateAdmission&) = delete;

    void define_level (int level, const amrex::BoxArray& grids,
                       const amrex::DistributionMapping& mapping);
    void destroy_level (int level);
    [[nodiscard]] bool admissible (const amrex::MultiFab& spectrum, int level,
                                   std::string& diagnostic) const;

    struct Counters {
        std::uint64_t device_admissions{0};
        std::uint64_t host_scans{0};
        std::uint64_t host_scan_fab_bytes{0};
        std::uint64_t host_d2h_bytes{0};
        std::uint64_t successful_admissions{0};
        std::uint64_t mpi_reductions{0};
    };
    [[nodiscard]] Counters counters () const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

/** Compare the projected spectrum to checkpointed qc/qr without mutating core. */
[[nodiscard]] bool restart_projection_matches (const amrex::MultiFab& spectrum,
                                               const amrex::MultiFab& persisted_core,
                                               const SBMBulkProjection& projection,
                                               int qc_component, int qr_component,
                                               amrex::Real tolerance_scale = amrex::Real(128.0));

} // namespace erf_sbm

#endif
