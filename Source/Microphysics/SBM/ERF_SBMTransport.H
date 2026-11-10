#ifndef ERF_SBM_TRANSPORT_H_
#define ERF_SBM_TRANSPORT_H_

#include "ERF_SBMStateManager.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMAdvectionBoundary.H"
#include "AuxiliaryState/ERF_AuxiliaryStage.H"

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace erf_sbm {

class SBMAMRFluxRegister;
class SBMStateAdmission;

/** Optional device-side probe used by the small-grid G1 ghost oracle. */
struct SBMAMRStageGhostProbe {
    amrex::IntVect fine_cell{0};
    int global_component{0};
    int requested_stage{-1};
    int requested_role{0}; // input=0, limiter-base=1, anchor=2
    amrex::Real* sampled_intensive_device{nullptr};
    amrex::Real* sampled_endpoint_device{nullptr};
    int* hit_device{nullptr};
    int* endpoint_hit_device{nullptr};
    // Optional A4 test-only observation: three physical (M,N,property)
    // components across the positive-x WENO stencil and the accepted face
    // rates at one selected interface face. The caller owns device storage.
    amrex::IntVect sampled_face{-1};
    int sampled_face_direction{-1};
    int sampled_stage{-1};
    int sampled_mass_component{-1};
    int sampled_number_component{-1};
    int sampled_property_component{-1};
    amrex::Real* sampled_physical_stencil_device{nullptr};
    amrex::Real* sampled_endpoint_stencil_device{nullptr};
    amrex::Real* sampled_mass_rate_device{nullptr};
    amrex::Real* sampled_number_rate_device{nullptr};
    amrex::Real* sampled_property_rate_device{nullptr};
    amrex::Real* sampled_high_mass_rate_device{nullptr};
    amrex::Real* sampled_high_number_rate_device{nullptr};
    amrex::Real* sampled_high_property_rate_device{nullptr};
    amrex::Real* sampled_low_mass_rate_device{nullptr};
    amrex::Real* sampled_low_number_rate_device{nullptr};
    amrex::Real* sampled_low_property_rate_device{nullptr};
    amrex::Real* sampled_lambda_device{nullptr};
    int* sampled_rate_hit_device{nullptr};
    int skipped_global_component{-1};
    bool seed_finite_stale_value{false};
    // Store the skipped component but omit its write-tally increment.
    bool skip_write_tally_only{false};
    amrex::Real stale_value{amrex::Real(0.25)};
};

/** Coarse host density endpoints paired with the state-manager's coarse SBM
 * endpoints during a fine-level stage. The caller supplies the actual ERF
 * owner times and the fine physical-step bound used for endpoint roundoff.
 */
struct SBMAMRStageParent {
    const amrex::MultiFab* conserved_old{nullptr};
    const amrex::MultiFab* conserved_new{nullptr};
    double old_time{0.0};
    double new_time{0.0};
    double fine_step_dt{0.0};
    amrex::IntVect refinement_ratio{1};
    SBMAMRStageGhostProbe* ghost_probe{nullptr};
    bool profile_parent_fill{false};
    // The bounded pack reduced the ParallelCopy call fanout but regressed
    // measured fine-step time. Keep the lower-cost componentwise path as the
    // production default and enable packing only in its equivalence oracle.
    bool componentwise_parent_fill{true};
    // Optional second face observation for the A4 stencil control. This is
    // consumed only by test probes and adds no persistent transport storage.
    SBMAMRStageGhostProbe* secondary_ghost_probe{nullptr};
    // Qualification-only distributed-value check. It compares received
    // ghost data against values read directly from the owning valid FAB.
    bool validate_split_ghost_values{false};
    // Qualification-only negative test: rank 0 skips one coarse-only ghost store.
    bool inject_parent_store_fault{false};
};

/**
 * Bounded mapped advection for the authoritative SBM spectrum.
 *
 * The state manager remains the owner of prognostic state. This object owns
 * per-level stage scratch, the static mapped measure, the compact projected
 * qc/qr completed-step face ledger, and (in reference mode) an optional
 * all-component completed spectral face ledger. Production AMR instances add
 * accepted chunk-local stage rates directly to their pending interface YA
 * register and do not retain that full spectral ledger. Operational AMR
 * execution remains disabled by the ERF runtime guard until lifecycle gates
 * are qualified.
 */
class SBMTransport
{
  public:
    static constexpr int default_max_groups_per_chunk = 16;

    explicit SBMTransport (const SBMLayout& layout,
                           int number_of_levels,
                           int max_groups_per_chunk = default_max_groups_per_chunk,
                           bool keep_spectral_ledger_reference = true);
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

    /** Bind the adjacent interface registers after final level layouts exist. */
    void set_amr_flux_registers (int level,
                                 SBMAMRFluxRegister* as_coarse,
                                 SBMAMRFluxRegister* as_fine);
    bool begin_physical_step (int level, double step_begin, double step_end,
                              std::string& diagnostic);
    bool accept_physical_step (int level, double step_begin, double step_end,
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
                                  int qr_component,
                                  const SBMAMRStageParent* parent = nullptr);

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
                        int qr_component,
                        const SBMAMRStageParent* parent = nullptr);

    void destroy (int level);
    [[nodiscard]] int nlevels () const noexcept
    { return static_cast<int>(m_levels.size()); }
    [[nodiscard]] bool is_defined (int level) const;
    [[nodiscard]] bool measure_is_ready (int level) const;
    [[nodiscard]] SBMStateAdmission& admission () noexcept;
    [[nodiscard]] const amrex::MultiFab& static_measure (int level) const;
    /** Complete constraint-group chunks shared by transport and AMR transfer. */
    [[nodiscard]] int max_groups_per_chunk () const noexcept
    { return m_max_groups_per_chunk; }
    [[nodiscard]] const std::vector<ConstraintClosureChunk>& chunks () const noexcept
    { return m_chunks; }
    [[nodiscard]] int max_chunk_components () const noexcept
    { return m_max_chunk_components; }
    [[nodiscard]] const erf_auxiliary::CompletedStepFluxLedger&
    projected_ledger (int level) const;
    /** Full all-component accepted-integral reference, present only when enabled. */
    [[nodiscard]] const erf_auxiliary::CompletedStepFluxLedger& completed_spectral_ledger (int level) const;
    // Qualification diagnostics; populated only in ERF_SBM_AMR_QUALIFICATION builds.
    void enable_zero_flow_diagnostics (bool enabled) noexcept
    { m_zero_flow_diagnostics_enabled = enabled; }
    [[nodiscard]] int interval_limiter_active_groups () const noexcept
    { return m_interval_limiter_active_groups; }
    [[nodiscard]] bool interval_limiter_changed_rate () const noexcept
    { return m_interval_limiter_changed_rate; }
    [[nodiscard]] amrex::Real interval_carrier_linf () const noexcept
    { return m_interval_carrier_linf; }
    [[nodiscard]] amrex::Real interval_accepted_rate_linf () const noexcept
    { return m_interval_accepted_rate_linf; }

  private:
    struct LevelStorage;
    struct DeviceMetadata;

    const SBMLayout& m_layout;
    SBMBulkProjection m_projection;
    std::vector<ConstraintGroup> m_groups;
    FlattenedConstraintSet m_flat_constraints;
    std::vector<AttachedPropertySupportDescriptor> m_property_support;
    std::vector<ConstraintClosureChunk> m_chunks;
    int m_max_groups_per_chunk{default_max_groups_per_chunk};
    std::vector<int> m_max_constraints_per_group;
    int m_max_chunk_components{0};
    std::vector<std::unique_ptr<LevelStorage>> m_levels;
    std::unique_ptr<DeviceMetadata> m_device_metadata;
    std::unique_ptr<SBMStateAdmission> m_admission;
    bool m_keep_spectral_ledger_reference{true};
    bool m_zero_flow_diagnostics_enabled{false};
    // Split qualification checks parent ghosts once at the coarse start and
    // once strictly inside the coarse interval.
    std::array<bool, 2> m_split_ghost_values_checked{};
    bool m_host_integrator_reported{false};
    int m_interval_limiter_active_groups{0};
    bool m_interval_limiter_changed_rate{false};
    amrex::Real m_interval_carrier_linf{amrex::Real(0.0)};
    amrex::Real m_interval_accepted_rate_linf{amrex::Real(0.0)};
};

} // namespace erf_sbm

#endif
