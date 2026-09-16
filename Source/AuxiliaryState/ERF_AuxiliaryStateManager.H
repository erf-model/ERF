#ifndef ERF_AUXILIARY_STATE_MANAGER_H_
#define ERF_AUXILIARY_STATE_MANAGER_H_

#include "ERF_AuxiliaryStateLayout.H"
#include "ERF_AuxiliaryFaceTransfer.H"

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace erf_auxiliary {

/**
 * Owns auxiliary old/evaluation/output storage outside ERF's compact cons
 * state.  The class is generic and has no dependency on any physics provider.
 */
class AuxiliaryStateManager {
public:
    AuxiliaryStateManager() = default;
    explicit AuxiliaryStateManager(AuxiliaryStateLayout layout)
        : m_layout(std::move(layout)) {}

    void define_level(int level, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm, int ngrow = 1,
                      int scratch_ncomp = -1);
    void remake_level(int level, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm, int ngrow = 1,
                      const amrex::Periodicity& periodicity = amrex::Periodicity::NonPeriodic(),
                      int scratch_ncomp = -1);
    void remake_level_from_coarse(int level, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm, int ngrow,
                                  const amrex::Periodicity& periodicity,
                                  int coarse_level,
                                  const amrex::Geometry& coarse_geometry,
                                  const amrex::Geometry& fine_geometry,
                                  const amrex::IntVect& ref_ratio,
                                  double time, int scratch_ncomp = -1);
    void destroy_level(int level);
    void average_down_to(int coarse_level, int fine_level,
                         const amrex::IntVect& ref_ratio);
    void prolong_from_coarse(int coarse_level, int fine_level,
                             const amrex::Geometry& coarse_geometry,
                             const amrex::Geometry& fine_geometry,
                             const amrex::IntVect& ref_ratio,
                             double time = -1.0);
    void fill_stage_from_coarse(int coarse_level, int fine_level, double time,
                                const amrex::Geometry& coarse_geometry,
                                const amrex::Geometry& fine_geometry,
                                const amrex::IntVect& ref_ratio);
    // `reset_face_transfer=false` starts another AMR fine substep while
    // retaining the accepted ledger accumulated over the parent coarse step.
    void begin_step(int level, double old_time = 0.0, bool reset_face_transfer = true);
    // Restore the temporal meaning of the three views after an external
    // checkpoint read.  A restart may regrid before the next begin_step call,
    // so the coarse-time bracket must be available during that remake.
    void set_time_views(int level, double old_time, double evaluation_time,
                        double output_time);
    void accept_stage(int level, double stage_time = 0.0);
    void record_stage_face_transfer(int level, const StageContext& context,
                                    const AuxiliaryFaceTransfer& stage_flux);

    [[nodiscard]] bool active() const noexcept { return !m_layout.schema_id().empty(); }
    [[nodiscard]] bool has_level(const int level) const noexcept {
        return level >= 0 && static_cast<std::size_t>(level) < m_output.size() &&
               m_output[static_cast<std::size_t>(level)] != nullptr;
    }
    [[nodiscard]] int nlevels() const noexcept { return static_cast<int>(m_output.size()); }
    [[nodiscard]] const AuxiliaryStateLayout& layout() const noexcept { return m_layout; }
    [[nodiscard]] amrex::MultiFab& old(int level) { return *m_old.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] amrex::MultiFab& evaluation(int level) { return *m_evaluation.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] amrex::MultiFab& output(int level) { return *m_output.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] amrex::MultiFab& scratch(int level) {
        if (!m_scratch.at(static_cast<std::size_t>(level))) {
            throw std::logic_error("auxiliary scratch was disabled for this production manager");
        }
        return *m_scratch.at(static_cast<std::size_t>(level));
    }
    [[nodiscard]] const amrex::MultiFab& old(int level) const { return *m_old.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] const amrex::MultiFab& evaluation(int level) const { return *m_evaluation.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] const amrex::MultiFab& output(int level) const { return *m_output.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] double old_time(int level) const { return m_old_time.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] double evaluation_time(int level) const { return m_evaluation_time.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] double output_time(int level) const { return m_output_time.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] std::size_t state_resident_bytes() const noexcept { return m_state_resident_bytes; }
    [[nodiscard]] std::size_t face_transfer_resident_bytes() const noexcept {
        return m_face_transfer_resident_bytes;
    }
    [[nodiscard]] std::size_t resident_bytes() const noexcept { return m_resident_bytes; }
    [[nodiscard]] const AuxiliaryFaceTransferLedger& face_transfer_ledger(int level) const {
        return *m_face_ledgers.at(static_cast<std::size_t>(level));
    }
    [[nodiscard]] AuxiliaryFaceTransferLedger& face_transfer_ledger(int level) {
        return *m_face_ledgers.at(static_cast<std::size_t>(level));
    }

private:
    void recompute_resident_bytes() noexcept;
    AuxiliaryStateLayout m_layout;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_old;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_evaluation;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_output;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_scratch;
    std::vector<std::unique_ptr<AuxiliaryFaceTransferLedger>> m_face_ledgers;
    std::vector<double> m_old_time;
    std::vector<double> m_evaluation_time;
    std::vector<double> m_output_time;
    std::size_t m_state_resident_bytes{0};
    std::size_t m_face_transfer_resident_bytes{0};
    std::size_t m_resident_bytes{0};
};

} // namespace erf_auxiliary

#endif
