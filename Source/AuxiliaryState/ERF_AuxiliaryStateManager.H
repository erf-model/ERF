#ifndef ERF_AUXILIARY_STATE_MANAGER_H_
#define ERF_AUXILIARY_STATE_MANAGER_H_

#include "ERF_AuxiliaryStateLayout.H"
#include "ERF_AuxiliaryFaceTransfer.H"

#include <AMReX_MultiFab.H>

#include <memory>
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
                      const amrex::DistributionMapping& dm, int ngrow = 1);
    void begin_step(int level);
    void accept_stage(int level);
    void record_stage_face_transfer(int level, const StageContext& context,
                                    const AuxiliaryFaceTransfer& stage_flux);

    [[nodiscard]] bool active() const noexcept { return !m_layout.schema_id().empty(); }
    [[nodiscard]] const AuxiliaryStateLayout& layout() const noexcept { return m_layout; }
    [[nodiscard]] amrex::MultiFab& old(int level) { return *m_old.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] amrex::MultiFab& evaluation(int level) { return *m_evaluation.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] amrex::MultiFab& output(int level) { return *m_output.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] amrex::MultiFab& scratch(int level) { return *m_scratch.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] const amrex::MultiFab& old(int level) const { return *m_old.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] const amrex::MultiFab& evaluation(int level) const { return *m_evaluation.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] const amrex::MultiFab& output(int level) const { return *m_output.at(static_cast<std::size_t>(level)); }
    [[nodiscard]] std::size_t resident_bytes() const noexcept { return m_resident_bytes; }
    [[nodiscard]] const AuxiliaryFaceTransferLedger& face_transfer_ledger(int level) const {
        return *m_face_ledgers.at(static_cast<std::size_t>(level));
    }
    [[nodiscard]] AuxiliaryFaceTransferLedger& face_transfer_ledger(int level) {
        return *m_face_ledgers.at(static_cast<std::size_t>(level));
    }

private:
    AuxiliaryStateLayout m_layout;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_old;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_evaluation;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_output;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_scratch;
    std::vector<std::unique_ptr<AuxiliaryFaceTransferLedger>> m_face_ledgers;
    std::size_t m_resident_bytes{0};
};

} // namespace erf_auxiliary

#endif
