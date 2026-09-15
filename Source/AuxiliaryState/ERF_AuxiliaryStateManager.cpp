#include "ERF_AuxiliaryStateManager.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
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

    recompute_resident_bytes();
}

void AuxiliaryStateManager::remake_level(const int level, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm,
                                         const int ngrow, const amrex::Periodicity& periodicity)
{
    if (level < 0 || m_layout.ncomp() <= 0 || ngrow < 0) {
        throw std::invalid_argument("invalid auxiliary state remake");
    }
    const bool had_level = has_level(level);
    const auto n = static_cast<std::size_t>(level + 1);
    m_old.resize(n); m_evaluation.resize(n); m_output.resize(n);
    m_scratch.resize(n); m_face_ledgers.resize(n);
    auto new_old = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    auto new_evaluation = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    auto new_output = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    auto new_scratch = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    new_old->setVal(0.0); new_evaluation->setVal(0.0); new_output->setVal(0.0); new_scratch->setVal(0.0);
    if (had_level) {
        new_old->ParallelCopy(*m_old[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                              new_old->nGrowVect(), m_old[static_cast<std::size_t>(level)]->nGrowVect(), periodicity);
        new_evaluation->ParallelCopy(*m_evaluation[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                                     new_evaluation->nGrowVect(), m_evaluation[static_cast<std::size_t>(level)]->nGrowVect(), periodicity);
        new_output->ParallelCopy(*m_output[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                                 new_output->nGrowVect(), m_output[static_cast<std::size_t>(level)]->nGrowVect(), periodicity);
        new_old->FillBoundary(periodicity); new_evaluation->FillBoundary(periodicity); new_output->FillBoundary(periodicity);
    }
    m_old[static_cast<std::size_t>(level)] = std::move(new_old);
    m_evaluation[static_cast<std::size_t>(level)] = std::move(new_evaluation);
    m_output[static_cast<std::size_t>(level)] = std::move(new_output);
    m_scratch[static_cast<std::size_t>(level)] = std::move(new_scratch);
    m_face_ledgers[static_cast<std::size_t>(level)] = std::make_unique<AuxiliaryFaceTransferLedger>();
    m_face_ledgers[static_cast<std::size_t>(level)]->define(ba, dm, m_layout.ncomp(), 3, 0);
    recompute_resident_bytes();
}

void AuxiliaryStateManager::destroy_level(const int level)
{
    if (level < 0 || !has_level(level)) return;
    m_old[static_cast<std::size_t>(level)].reset();
    m_evaluation[static_cast<std::size_t>(level)].reset();
    m_output[static_cast<std::size_t>(level)].reset();
    m_scratch[static_cast<std::size_t>(level)].reset();
    m_face_ledgers[static_cast<std::size_t>(level)].reset();
    recompute_resident_bytes();
}

void AuxiliaryStateManager::average_down_to(const int coarse_level, const int fine_level,
                                            const amrex::IntVect& ref_ratio)
{
    if (!has_level(coarse_level) || !has_level(fine_level) || coarse_level >= fine_level) {
        throw std::invalid_argument("invalid auxiliary average-down levels");
    }
    for (auto* coarse : {m_old[static_cast<std::size_t>(coarse_level)].get(),
                         m_evaluation[static_cast<std::size_t>(coarse_level)].get(),
                         m_output[static_cast<std::size_t>(coarse_level)].get()}) {
        const amrex::MultiFab* fine = coarse == m_old[static_cast<std::size_t>(coarse_level)].get()
            ? m_old[static_cast<std::size_t>(fine_level)].get()
            : (coarse == m_evaluation[static_cast<std::size_t>(coarse_level)].get()
                ? m_evaluation[static_cast<std::size_t>(fine_level)].get()
                : m_output[static_cast<std::size_t>(fine_level)].get());
        amrex::average_down(*fine, *coarse, 0, m_layout.ncomp(), ref_ratio);
        coarse->FillBoundary(amrex::Periodicity::NonPeriodic());
    }
}

void AuxiliaryStateManager::prolong_from_coarse(const int coarse_level, const int fine_level,
                                                const amrex::Geometry& coarse_geometry,
                                                const amrex::Geometry& fine_geometry,
                                                const amrex::IntVect& ref_ratio)
{
    if (!has_level(coarse_level) || !has_level(fine_level) || coarse_level >= fine_level ||
        ref_ratio.min() <= 0) {
        throw std::invalid_argument("invalid auxiliary prolongation levels or refinement ratio");
    }
    // P2 supports only fully periodic static Cartesian geometry.  The
    // no-physical-BC overload is consequently the precise piecewise-constant
    // injection required here and does not invent a provider boundary rule.
    if (!coarse_geometry.isAllPeriodic() || !fine_geometry.isAllPeriodic()) {
        throw std::invalid_argument("SBM auxiliary prolongation requires periodic geometry");
    }
    amrex::Vector<amrex::BCRec> bcs(static_cast<std::size_t>(m_layout.ncomp()));
    auto prolong = [this, coarse_level, fine_level, &coarse_geometry, &fine_geometry,
                    &ref_ratio, &bcs](amrex::MultiFab& fine, const amrex::MultiFab& coarse) {
        amrex::InterpFromCoarseLevel(fine, fine.nGrowVect(), amrex::IntVect(0), coarse,
                                     0, 0, m_layout.ncomp(), coarse_geometry, fine_geometry,
                                     ref_ratio, &amrex::pc_interp, bcs, 0);
        fine.FillBoundary(fine_geometry.periodicity());
    };
    prolong(output(fine_level), output(coarse_level));
    prolong(old(fine_level), old(coarse_level));
    prolong(evaluation(fine_level), evaluation(coarse_level));
    scratch(fine_level).setVal(0.0);
    scratch(fine_level).FillBoundary(fine_geometry.periodicity());
}

void AuxiliaryStateManager::recompute_resident_bytes() noexcept
{
    m_state_resident_bytes = 0;
    m_face_transfer_resident_bytes = 0;
    for (std::size_t level = 0; level < m_output.size(); ++level) {
        if (m_old[level]) {
            m_state_resident_bytes += allocated_payload_bytes(*m_old[level]) +
                allocated_payload_bytes(*m_evaluation[level]) + allocated_payload_bytes(*m_output[level]) +
                allocated_payload_bytes(*m_scratch[level]);
        }
        if (m_face_ledgers[level]) m_face_transfer_resident_bytes += m_face_ledgers[level]->resident_bytes();
    }
    m_resident_bytes = m_state_resident_bytes + m_face_transfer_resident_bytes;
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
