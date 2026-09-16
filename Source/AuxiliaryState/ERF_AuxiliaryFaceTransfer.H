#ifndef ERF_AUXILIARY_FACE_TRANSFER_H_
#define ERF_AUXILIARY_FACE_TRANSFER_H_

#include "ERF_AuxiliaryStageContext.H"

#include <AMReX_MultiFab.H>

#include <array>
#include <cstddef>
#include <memory>

namespace erf_auxiliary {

/**
 * @brief Logical MultiFab payload size, including each FAB's allocated ghosts.
 *
 * This reports the data payload rather than allocator metadata.  Counting the
 * grown BoxArray entries makes the result independent of MPI ownership and
 * matches the FAB boxes allocated by AMReX.
 */
[[nodiscard]] std::size_t allocated_payload_bytes(const amrex::MultiFab& state) noexcept;

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
    [[nodiscard]] bool compatible_with(const amrex::BoxArray& ba,
                                       const amrex::DistributionMapping& dm,
                                       int ncomp) const noexcept;
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
 * @brief Reuses one stage flux and retains the accepted full-step transfer.
 *
 * `stage()` is overwritten by each production stage.  `record_stage` accumulates
 * the exact stage flux using `StageContext`: compressible `I = dt * F(stage 2)`;
 * anelastic `I = dt/2 * (F(stage 0) + F(stage 1))`.  No full stage history is
 * retained, so storage is two face-transfer objects regardless of stage count.
 */
class AuxiliaryFaceTransferLedger {
public:
    void define(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                int ncomp, int nstages = 3, int ngrow = 0);
    void begin_step(bool reset_accepted = true);
    void record_stage(const StageContext& context, const AuxiliaryFaceTransfer& stage_flux);

    [[nodiscard]] AuxiliaryFaceTransfer& stage() { return *m_stage; }
    [[nodiscard]] const AuxiliaryFaceTransfer& stage() const { return *m_stage; }
    [[nodiscard]] const AuxiliaryFaceTransfer& accepted() const { return *m_accepted; }
    [[nodiscard]] int nstages() const noexcept { return m_nstages; }
    [[nodiscard]] std::size_t stage_resident_bytes() const noexcept {
        return m_stage ? m_stage->resident_bytes() : 0;
    }
    [[nodiscard]] std::size_t accepted_resident_bytes() const noexcept {
        return m_accepted ? m_accepted->resident_bytes() : 0;
    }
    [[nodiscard]] std::size_t resident_bytes() const noexcept { return m_resident_bytes; }

private:
    std::unique_ptr<AuxiliaryFaceTransfer> m_stage;
    std::unique_ptr<AuxiliaryFaceTransfer> m_accepted;
    std::array<bool, 3> m_recorded{{false, false, false}};
    int m_nstages{0};
    int m_ncomp{0};
    std::size_t m_resident_bytes{0};
};

} // namespace erf_auxiliary

#endif
