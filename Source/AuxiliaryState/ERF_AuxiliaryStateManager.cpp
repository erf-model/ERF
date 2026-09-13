#include "ERF_AuxiliaryStateManager.H"

#include <AMReX_MultiFabUtil.H>

#include <stdexcept>

namespace erf_auxiliary {

void AuxiliaryStateManager::define_level(const int level,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm,
                                         const int ngrow)
{
    if (level < 0 || m_layout.ncomp() <= 0 || ngrow < 0) {
        throw std::invalid_argument("invalid auxiliary state definition");
    }
    const auto n = static_cast<std::size_t>(level + 1);
    m_old.resize(n);
    m_evaluation.resize(n);
    m_output.resize(n);
    m_scratch.resize(n);
    m_face_ledgers.resize(n);

    m_old[static_cast<std::size_t>(level)] = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    m_evaluation[static_cast<std::size_t>(level)] = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    m_output[static_cast<std::size_t>(level)] = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    m_scratch[static_cast<std::size_t>(level)] = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    m_old[static_cast<std::size_t>(level)]->setVal(0.0);
    m_evaluation[static_cast<std::size_t>(level)]->setVal(0.0);
    m_output[static_cast<std::size_t>(level)]->setVal(0.0);
    m_scratch[static_cast<std::size_t>(level)]->setVal(0.0);
    m_face_ledgers[static_cast<std::size_t>(level)] = std::make_unique<AuxiliaryFaceTransferLedger>();
    m_face_ledgers[static_cast<std::size_t>(level)]->define(ba, dm, m_layout.ncomp(), 3, 0);

    const std::size_t cells = static_cast<std::size_t>(ba.numPts());
    const std::size_t nboxes = static_cast<std::size_t>(ba.size());
    const std::size_t state_bytes = nboxes == 0 ? 0 : cells * static_cast<std::size_t>(m_layout.ncomp()) *
        static_cast<std::size_t>(4) * sizeof(amrex::Real);
    m_resident_bytes = state_bytes + m_face_ledgers[static_cast<std::size_t>(level)]->resident_bytes();
}

void AuxiliaryStateManager::begin_step(const int level)
{
    amrex::MultiFab::Copy(old(level), output(level), 0, 0, m_layout.ncomp(), output(level).nGrowVect());
    amrex::MultiFab::Copy(evaluation(level), output(level), 0, 0, m_layout.ncomp(), output(level).nGrowVect());
    face_transfer_ledger(level).begin_step();
}

void AuxiliaryStateManager::accept_stage(const int level)
{
    // Keep output as the accepted state and publish it as the evaluation
    // state for the next callback.  The old full-step baseline remains
    // untouched throughout all explicit stages.
    amrex::MultiFab::Copy(evaluation(level), output(level), 0, 0, m_layout.ncomp(), output(level).nGrowVect());
}

void AuxiliaryStateManager::record_stage_face_transfer(const int level,
                                                       const StageContext& context,
                                                       const AuxiliaryFaceTransfer& stage_flux)
{
    face_transfer_ledger(level).record_stage(context, stage_flux);
}

} // namespace erf_auxiliary
