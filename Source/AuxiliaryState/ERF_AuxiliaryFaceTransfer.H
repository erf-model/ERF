#ifndef ERF_AUXILIARY_FACE_TRANSFER_H_
#define ERF_AUXILIARY_FACE_TRANSFER_H_

#include "ERF_AuxiliaryStageContext.H"

#include <AMReX_MultiFab.H>

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace erf_auxiliary {

/**
 * @brief Stage-local face flux storage for an auxiliary conserved state.
 *
 * Each direction is a face-centered MultiFab with the same box decomposition
 * as the cell-centered auxiliary state.  The class contains no physics or
 * component-name knowledge, so an AMR reflux consumer can use it later.
 */
class AuxiliaryFaceTransfer {
public:
    AuxiliaryFaceTransfer() = default;

    void define(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                int ncomp, int ngrow = 0);
    void setVal(amrex::Real value);

    [[nodiscard]] bool defined() const noexcept { return m_face[0] != nullptr; }
    [[nodiscard]] int ncomp() const noexcept { return m_ncomp; }
    [[nodiscard]] std::size_t resident_bytes() const noexcept { return m_resident_bytes; }

    [[nodiscard]] amrex::MultiFab& direction(int dir) { return *m_face.at(static_cast<std::size_t>(dir)); }
    [[nodiscard]] const amrex::MultiFab& direction(int dir) const { return *m_face.at(static_cast<std::size_t>(dir)); }
    [[nodiscard]] amrex::MultiFab& x() { return direction(0); }
    [[nodiscard]] amrex::MultiFab& y() { return direction(1); }
    [[nodiscard]] amrex::MultiFab& z() { return direction(2); }
    [[nodiscard]] const amrex::MultiFab& x() const { return direction(0); }
    [[nodiscard]] const amrex::MultiFab& y() const { return direction(1); }
    [[nodiscard]] const amrex::MultiFab& z() const { return direction(2); }

private:
    std::array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM> m_face;
    int m_ncomp{0};
    std::size_t m_resident_bytes{0};
};

/**
 * @brief Retains actual stage fluxes and the accepted full-step transfer.
 *
 * `record_stage` copies the numerical face flux used by a stage into its
 * stage slot and accumulates the ERF temporal ledger.  It intentionally keeps
 * only the current three-stage RK/Heun history and one accepted transfer.
 */
class AuxiliaryFaceTransferLedger {
public:
    void define(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                int ncomp, int nstages = 3, int ngrow = 0);
    void begin_step();
    void record_stage(const StageContext& context, const AuxiliaryFaceTransfer& stage_flux);

    [[nodiscard]] AuxiliaryFaceTransfer& stage(int index) { return *m_stages.at(static_cast<std::size_t>(index)); }
    [[nodiscard]] const AuxiliaryFaceTransfer& stage(int index) const { return *m_stages.at(static_cast<std::size_t>(index)); }
    [[nodiscard]] const AuxiliaryFaceTransfer& accepted() const { return *m_accepted; }
    [[nodiscard]] int nstages() const noexcept { return static_cast<int>(m_stages.size()); }
    [[nodiscard]] std::size_t resident_bytes() const noexcept { return m_resident_bytes; }

private:
    std::vector<std::unique_ptr<AuxiliaryFaceTransfer>> m_stages;
    std::unique_ptr<AuxiliaryFaceTransfer> m_accepted;
    std::vector<bool> m_recorded;
    int m_ncomp{0};
    std::size_t m_resident_bytes{0};
};

} // namespace erf_auxiliary

#endif
