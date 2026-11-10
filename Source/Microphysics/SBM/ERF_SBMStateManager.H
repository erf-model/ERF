#ifndef ERF_SBM_STATE_MANAGER_H_
#define ERF_SBM_STATE_MANAGER_H_

#include "ERF_SBMBulkProjection.H"

#include <AMReX_MultiFab.H>

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace erf_sbm {

/** Owns the authoritative old/new spectral state and its semantic times. */
class SBMStateManager {
public:
    SBMStateManager (SBMLayout layout, int nlevels);

    [[nodiscard]] const SBMLayout& layout () const noexcept { return m_layout; }
    [[nodiscard]] int nlevels () const noexcept { return static_cast<int>(m_levels.size()); }
    [[nodiscard]] bool is_defined (int lev) const;

    void define (int lev, const amrex::BoxArray& grids,
                 const amrex::DistributionMapping& mapping,
                 double current_time, int ngrow = 0);
    void destroy (int lev);

    /** Accepted state at the beginning of an active physical step. */
    [[nodiscard]] const amrex::MultiFab& old_state (int lev) const;
    [[nodiscard]] double old_time (int lev) const;
    [[nodiscard]] bool old_valid (int lev) const;

    /** Accepted/current state. Throws while the target buffer is invalid. */
    [[nodiscard]] const amrex::MultiFab& new_state (int lev) const;
    [[nodiscard]] double new_time (int lev) const;
    [[nodiscard]] bool new_valid (int lev) const;

    /** Mutable accepted state for initialization and synchronized restart I/O. */
    [[nodiscard]] amrex::MultiFab& new_state_for_initialization (int lev);

    /** Mutable target storage for transport; semantic reads remain invalid until accepted. */
    [[nodiscard]] amrex::MultiFab& new_target_storage (int lev);

    /** Rotate accepted new state into old state at the ERF::Advance seam. */
    /** A fine-level caller may supply its actual dt to admit bounded
     * scheduler-roundoff between consecutive substep endpoint expressions.
     * Without it, owner time identity remains exact.
     */
    bool begin_step (int lev, double step_old_time, std::string& diagnostic,
                     std::optional<double> dt_min = std::nullopt);

    /** Publish a candidate accepted by the host stage sequence. */
    bool accept_stage_target (int lev, double target_time,
                              bool physical_step_complete,
                              std::string& diagnostic);
    /** Publish an admitted same-time AMR synchronization candidate. */
    bool accept_synchronized_state (int lev,
                                   const amrex::MultiFab& candidate,
                                   double semantic_time,
                                   std::string& diagnostic);
    bool can_accept_synchronized_state (int lev,
                                        const amrex::MultiFab& candidate,
                                        double semantic_time,
                                        std::string& diagnostic) const;
    /** Commit a candidate after collective can_accept_synchronized_state succeeds. */
    void publish_synchronized_state (int lev,
                                     const amrex::MultiFab& candidate);
    [[nodiscard]] bool step_active (int lev) const;

    /** One-way projection from the spectrum into ERF's compact qc and qr lanes. */
    void project_to_core (int lev, amrex::MultiFab& core,
                          int qc_component, int qr_component) const;

private:
    struct LevelState {
        std::unique_ptr<amrex::MultiFab> old_state;
        std::unique_ptr<amrex::MultiFab> new_state;
        double old_time{0.0};
        double new_time{0.0};
        bool old_valid{false};
        bool new_valid{false};
        bool step_active{false};
    };

    [[nodiscard]] LevelState& level (int lev);
    [[nodiscard]] const LevelState& level (int lev) const;

    SBMLayout m_layout;
    SBMBulkProjection m_projection;
    std::vector<LevelState> m_levels;
};

} // namespace erf_sbm

#endif
