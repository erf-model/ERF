#ifndef ERF_SBM_TRANSPORT_H_
#define ERF_SBM_TRANSPORT_H_

#include "ERF_SBMStateManager.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMAdvectionBoundary.H"
#include "AuxiliaryState/ERF_AuxiliaryStage.H"

#include <AMReX_Geometry.H>

#include <memory>
#include <string>
#include <vector>

namespace erf_sbm {

/**
 * Bounded mapped advection for the authoritative SBM spectrum.
 *
 * The state manager remains the owner of prognostic state. This object owns
 * per-level stage scratch, the static mapped measure, the compact projected
 * qc/qr completed-step face ledger, and an optional all-component completed
 * spectral face ledger for AMR-capable instances. Operational AMR execution
 * remains disabled by the ERF runtime guard.
 */
class SBMTransport
{
  public:
    explicit SBMTransport (const SBMLayout& layout,
                           int number_of_levels,
                           int max_groups_per_chunk = 16);
    ~SBMTransport ();
    SBMTransport (const SBMTransport&) = delete;
    SBMTransport& operator= (const SBMTransport&) = delete;

    void define (int level,
                 const amrex::BoxArray& cell_ba,
                 const amrex::DistributionMapping& dm);
    bool rebuild_static_measure (int level,
                                 const amrex::MultiFab& detJ,
                                 const amrex::MultiFab& mx,
                                 const amrex::MultiFab& my,
                                 std::string& diagnostic);

    /**
     * ERF slow-RHS seam. The conserved-state views are ordered as the
     * old-step anchor, the actual host predictor/input, and the target being
     * written. Keep this forwarding point shared by the ERF call site and
     * seam-level tests so density-time pairing cannot drift silently.
     */
    void advance_stage_from_host (int level,
                                  erf_auxiliary::HostIntegrator method,
                                  int stage,
                                  double step_old_time,
                                  double input_time,
                                  double target_time,
                                  double host_stage_interval,
                                  SBMStateManager& state_manager,
                                  const amrex::MultiFab& S_old_cons,
                                  const amrex::MultiFab& S_new_cons,
                                  amrex::MultiFab& S_data_cons,
                                  const amrex::MultiFab& avg_xmom,
                                  const amrex::MultiFab& avg_ymom,
                                  const amrex::MultiFab& avg_zmom,
                                  const amrex::Geometry& geometry,
                                  int qc_component,
                                  int qr_component);

    void advance_stage (int level,
                        erf_auxiliary::HostIntegrator method,
                        int stage,
                        double step_old_time,
                        double input_time,
                        double target_time,
                        double host_stage_interval,
                        SBMStateManager& state_manager,
                        const amrex::MultiFab& conserved_anchor,
                        const amrex::MultiFab& conserved_input,
                        amrex::MultiFab& conserved_target,
                        const amrex::MultiFab& avg_xmom,
                        const amrex::MultiFab& avg_ymom,
                        const amrex::MultiFab& avg_zmom,
                        const amrex::Geometry& geometry,
                        int qc_component,
                        int qr_component);

    void destroy (int level);
    [[nodiscard]] bool is_defined (int level) const;
    [[nodiscard]] bool measure_is_ready (int level) const;
    [[nodiscard]] const amrex::MultiFab& static_measure (int level) const;
    [[nodiscard]] const erf_auxiliary::CompletedStepFluxLedger&
    projected_ledger (int level) const;
    /** Authoritative all-component accepted transfer, defined for AMR-capable instances only. */
    [[nodiscard]] const erf_auxiliary::CompletedStepFluxLedger& completed_spectral_ledger (int level) const;

  private:
    struct LevelStorage;
    struct DeviceMetadata;

    const SBMLayout& m_layout;
    SBMBulkProjection m_projection;
    std::vector<ConstraintGroup> m_groups;
    FlattenedConstraintSet m_flat_constraints;
    std::vector<AttachedPropertySupportDescriptor> m_property_support;
    std::vector<ConstraintClosureChunk> m_chunks;
    std::vector<int> m_max_constraints_per_group;
    std::vector<std::unique_ptr<LevelStorage>> m_levels;
    std::unique_ptr<DeviceMetadata> m_device_metadata;
};

} // namespace erf_sbm

#endif
