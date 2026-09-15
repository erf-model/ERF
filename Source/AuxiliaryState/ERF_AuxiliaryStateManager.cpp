#include "ERF_AuxiliaryStateManager.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <cmath>
#include <limits>
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
    m_old_time.resize(n, 0.0);
    m_evaluation_time.resize(n, 0.0);
    m_output_time.resize(n, 0.0);

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
    m_old_time.resize(n, 0.0); m_evaluation_time.resize(n, 0.0); m_output_time.resize(n, 0.0);
    const double saved_old_time = had_level ? m_old_time[static_cast<std::size_t>(level)] : 0.0;
    const double saved_evaluation_time = had_level ? m_evaluation_time[static_cast<std::size_t>(level)] : 0.0;
    const double saved_output_time = had_level ? m_output_time[static_cast<std::size_t>(level)] : 0.0;
    auto new_old = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    auto new_evaluation = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    auto new_output = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    auto new_scratch = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
    new_old->setVal(0.0); new_evaluation->setVal(0.0); new_output->setVal(0.0); new_scratch->setVal(0.0);
    if (had_level) {
        new_old->ParallelCopy(*m_old[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                              amrex::IntVect(0), amrex::IntVect(0),
                              amrex::Periodicity::NonPeriodic());
        new_evaluation->ParallelCopy(*m_evaluation[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                                     amrex::IntVect(0), amrex::IntVect(0),
                                     amrex::Periodicity::NonPeriodic());
        new_output->ParallelCopy(*m_output[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                                 amrex::IntVect(0), amrex::IntVect(0),
                                 amrex::Periodicity::NonPeriodic());
        new_old->FillBoundary(periodicity); new_evaluation->FillBoundary(periodicity); new_output->FillBoundary(periodicity);
    }
    m_old[static_cast<std::size_t>(level)] = std::move(new_old);
    m_evaluation[static_cast<std::size_t>(level)] = std::move(new_evaluation);
    m_output[static_cast<std::size_t>(level)] = std::move(new_output);
    m_scratch[static_cast<std::size_t>(level)] = std::move(new_scratch);
    m_face_ledgers[static_cast<std::size_t>(level)] = std::make_unique<AuxiliaryFaceTransferLedger>();
    m_face_ledgers[static_cast<std::size_t>(level)]->define(ba, dm, m_layout.ncomp(), 3, 0);
    m_old_time[static_cast<std::size_t>(level)] = saved_old_time;
    m_evaluation_time[static_cast<std::size_t>(level)] = saved_evaluation_time;
    m_output_time[static_cast<std::size_t>(level)] = saved_output_time;
    recompute_resident_bytes();
}

void AuxiliaryStateManager::remake_level_from_coarse(
    const int level, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
    const int ngrow, const amrex::Periodicity& periodicity, const int coarse_level,
    const amrex::Geometry& coarse_geometry, const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ref_ratio, const double time)
{
    if (!has_level(coarse_level) || level <= coarse_level || ref_ratio.min() <= 0 ||
        !coarse_geometry.isAllPeriodic() || !fine_geometry.isAllPeriodic()) {
        throw std::invalid_argument("invalid auxiliary coarse-filled remake contract");
    }
    const double coarse_old_time = old_time(coarse_level);
    const double coarse_new_time = output_time(coarse_level);
    const double time_scale = 1.0 + std::max(std::abs(coarse_old_time), std::abs(coarse_new_time));
    if (time < coarse_old_time - 128.0 * std::numeric_limits<double>::epsilon() * time_scale ||
        time > coarse_new_time + 128.0 * std::numeric_limits<double>::epsilon() * time_scale) {
        throw std::invalid_argument("auxiliary remake time is outside the authoritative coarse time bracket");
    }

    const bool had_level = has_level(level);
    const auto n = static_cast<std::size_t>(level + 1);
    m_old.resize(n); m_evaluation.resize(n); m_output.resize(n); m_scratch.resize(n);
    m_face_ledgers.resize(n); m_old_time.resize(n, 0.0); m_evaluation_time.resize(n, 0.0); m_output_time.resize(n, 0.0);
    const double saved_old_time = had_level ? m_old_time[static_cast<std::size_t>(level)] : time;
    const double saved_evaluation_time = had_level ? m_evaluation_time[static_cast<std::size_t>(level)] : time;
    const double saved_output_time = had_level ? m_output_time[static_cast<std::size_t>(level)] : time;

    auto make_state = [&]() {
        auto state = std::make_unique<amrex::MultiFab>(ba, dm, m_layout.ncomp(), ngrow);
        state->setVal(0.0);
        return state;
    };
    auto new_old = make_state();
    auto new_evaluation = make_state();
    auto new_output = make_state();
    auto new_scratch = make_state();
    amrex::Vector<amrex::BCRec> bcs(static_cast<std::size_t>(m_layout.ncomp()));

    auto coarse_at_time = [&](const amrex::MultiFab& coarse_old,
                              const amrex::MultiFab& coarse_output) {
        auto interpolated = std::make_unique<amrex::MultiFab>(
            coarse_old.boxArray(), coarse_old.DistributionMap(), m_layout.ncomp(), coarse_old.nGrowVect());
        const double denominator = coarse_new_time - coarse_old_time;
        if (std::abs(denominator) <= 128.0 * std::numeric_limits<double>::epsilon() * time_scale) {
            amrex::MultiFab::Copy(*interpolated, coarse_output, 0, 0, m_layout.ncomp(), coarse_output.nGrowVect());
        } else {
            const amrex::Real theta = static_cast<amrex::Real>((time - coarse_old_time) / denominator);
            amrex::MultiFab::LinComb(*interpolated, amrex::Real(1.0) - theta, coarse_old,
                                     0, theta, coarse_output, 0, 0, m_layout.ncomp(), coarse_old.nGrowVect());
        }
        interpolated->FillBoundary(coarse_geometry.periodicity());
        return interpolated;
    };
    auto fill_from_coarse = [&](amrex::MultiFab& fine, const amrex::MultiFab& coarse) {
        amrex::InterpFromCoarseLevel(fine, fine.nGrowVect(), amrex::IntVect(0), coarse,
                                     0, 0, m_layout.ncomp(), coarse_geometry, fine_geometry,
                                     ref_ratio, &amrex::pc_interp, bcs, 0);
        fine.FillBoundary(periodicity);
    };

    // New coverage is created at one regrid time, so every newly introduced
    // semantic state is initialized from the same temporally interpolated
    // authoritative coarse spectrum.  Existing fine valid cells are overlaid
    // below, preserving their old/evaluation/output meanings.
    auto coarse_state = coarse_at_time(old(coarse_level), output(coarse_level));
    fill_from_coarse(*new_old, *coarse_state);
    fill_from_coarse(*new_evaluation, *coarse_state);
    fill_from_coarse(*new_output, *coarse_state);
    if (had_level) {
        new_old->ParallelCopy(*m_old[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                              amrex::IntVect(0), amrex::IntVect(0),
                              amrex::Periodicity::NonPeriodic());
        new_evaluation->ParallelCopy(*m_evaluation[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                                     amrex::IntVect(0), amrex::IntVect(0),
                                     amrex::Periodicity::NonPeriodic());
        new_output->ParallelCopy(*m_output[static_cast<std::size_t>(level)], 0, 0, m_layout.ncomp(),
                                 amrex::IntVect(0), amrex::IntVect(0),
                                 amrex::Periodicity::NonPeriodic());
        new_old->FillBoundary(periodicity); new_evaluation->FillBoundary(periodicity); new_output->FillBoundary(periodicity);
    }
    m_old[static_cast<std::size_t>(level)] = std::move(new_old);
    m_evaluation[static_cast<std::size_t>(level)] = std::move(new_evaluation);
    m_output[static_cast<std::size_t>(level)] = std::move(new_output);
    m_scratch[static_cast<std::size_t>(level)] = std::move(new_scratch);
    m_face_ledgers[static_cast<std::size_t>(level)] = std::make_unique<AuxiliaryFaceTransferLedger>();
    m_face_ledgers[static_cast<std::size_t>(level)]->define(ba, dm, m_layout.ncomp(), 3, 0);
    m_old_time[static_cast<std::size_t>(level)] = had_level ? saved_old_time : time;
    m_evaluation_time[static_cast<std::size_t>(level)] = had_level ? saved_evaluation_time : time;
    m_output_time[static_cast<std::size_t>(level)] = had_level ? saved_output_time : time;
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
                                                const amrex::IntVect& ref_ratio,
                                                const double requested_time)
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
    const double coarse_old_time = old_time(coarse_level);
    const double coarse_new_time = output_time(coarse_level);
    const double time = requested_time >= 0.0 ? requested_time : coarse_new_time;
    const double time_scale = 1.0 + std::max(std::abs(coarse_old_time), std::abs(coarse_new_time));
    if (time < coarse_old_time - 128.0 * std::numeric_limits<double>::epsilon() * time_scale ||
        time > coarse_new_time + 128.0 * std::numeric_limits<double>::epsilon() * time_scale) {
        throw std::invalid_argument("auxiliary prolongation time is outside the authoritative coarse time bracket");
    }
    auto coarse_state = std::make_unique<amrex::MultiFab>(
        old(coarse_level).boxArray(), old(coarse_level).DistributionMap(), m_layout.ncomp(), old(coarse_level).nGrowVect());
    const double denominator = coarse_new_time - coarse_old_time;
    if (std::abs(denominator) <= 128.0 * std::numeric_limits<double>::epsilon() * time_scale) {
        amrex::MultiFab::Copy(*coarse_state, output(coarse_level), 0, 0, m_layout.ncomp(), output(coarse_level).nGrowVect());
    } else {
        const amrex::Real theta = static_cast<amrex::Real>((time - coarse_old_time) / denominator);
        amrex::MultiFab::LinComb(*coarse_state, amrex::Real(1.0) - theta, old(coarse_level),
                                 0, theta, output(coarse_level), 0, 0, m_layout.ncomp(), old(coarse_level).nGrowVect());
    }
    coarse_state->FillBoundary(coarse_geometry.periodicity());
    amrex::Vector<amrex::BCRec> bcs(static_cast<std::size_t>(m_layout.ncomp()));
    auto prolong = [this, &coarse_geometry, &fine_geometry, &ref_ratio, &bcs,
                    &coarse_state](amrex::MultiFab& fine) {
        amrex::InterpFromCoarseLevel(fine, fine.nGrowVect(), amrex::IntVect(0), *coarse_state,
                                     0, 0, m_layout.ncomp(), coarse_geometry, fine_geometry,
                                     ref_ratio, &amrex::pc_interp, bcs, 0);
        fine.FillBoundary(fine_geometry.periodicity());
    };
    prolong(output(fine_level));
    prolong(old(fine_level));
    prolong(evaluation(fine_level));
    scratch(fine_level).setVal(0.0);
    scratch(fine_level).FillBoundary(fine_geometry.periodicity());
}

void AuxiliaryStateManager::fill_stage_from_coarse(
    const int coarse_level, const int fine_level, const double time,
    const amrex::Geometry& coarse_geometry, const amrex::Geometry& fine_geometry,
    const amrex::IntVect& ref_ratio)
{
    if (!has_level(coarse_level) || !has_level(fine_level) || coarse_level >= fine_level ||
        ref_ratio.min() <= 0 || !coarse_geometry.isAllPeriodic() || !fine_geometry.isAllPeriodic()) {
        throw std::invalid_argument("invalid auxiliary stage FillPatch levels or geometry");
    }
    const double coarse_old_time = old_time(coarse_level);
    const double coarse_new_time = output_time(coarse_level);
    const double time_scale = 1.0 + std::max(std::abs(coarse_old_time), std::abs(coarse_new_time));
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * time_scale;
    if (time < coarse_old_time - tolerance || time > coarse_new_time + tolerance) {
        throw std::invalid_argument("auxiliary stage time is outside the authoritative coarse time bracket");
    }
    old(coarse_level).FillBoundary(coarse_geometry.periodicity());
    output(coarse_level).FillBoundary(coarse_geometry.periodicity());
    evaluation(fine_level).FillBoundary(fine_geometry.periodicity());
    amrex::Vector<amrex::MultiFab*> coarse_states{&old(coarse_level), &output(coarse_level)};
    amrex::Vector<amrex::Real> coarse_times{static_cast<amrex::Real>(coarse_old_time),
                                            static_cast<amrex::Real>(coarse_new_time)};
    amrex::Vector<amrex::MultiFab*> fine_states{&evaluation(fine_level), &evaluation(fine_level)};
    amrex::Vector<amrex::Real> fine_times{static_cast<amrex::Real>(time), static_cast<amrex::Real>(time)};
    amrex::Vector<amrex::BCRec> bcs(static_cast<std::size_t>(m_layout.ncomp()));
    amrex::FillPatchTwoLevels(evaluation(fine_level), evaluation(fine_level).nGrowVect(),
                              amrex::IntVect(0), static_cast<amrex::Real>(time),
                              coarse_states, coarse_times, fine_states, fine_times,
                              0, 0, m_layout.ncomp(), coarse_geometry, fine_geometry,
                              ref_ratio, &amrex::pc_interp, bcs, 0);
    evaluation(fine_level).FillBoundary(fine_geometry.periodicity());
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

void AuxiliaryStateManager::begin_step(const int level, const double old_time_value)
{
    amrex::MultiFab::Copy(old(level), output(level), 0, 0, m_layout.ncomp(), output(level).nGrowVect());
    amrex::MultiFab::Copy(evaluation(level), output(level), 0, 0, m_layout.ncomp(), output(level).nGrowVect());
    face_transfer_ledger(level).begin_step();
    m_old_time[static_cast<std::size_t>(level)] = old_time_value;
    m_evaluation_time[static_cast<std::size_t>(level)] = old_time_value;
    m_output_time[static_cast<std::size_t>(level)] = old_time_value;
}

void AuxiliaryStateManager::accept_stage(const int level, const double stage_time)
{
    // Keep output as the accepted state and publish it as the evaluation
    // state for the next callback.  The old full-step baseline remains
    // untouched throughout all explicit stages.
    amrex::MultiFab::Copy(evaluation(level), output(level), 0, 0, m_layout.ncomp(), output(level).nGrowVect());
    m_evaluation_time[static_cast<std::size_t>(level)] = stage_time;
    m_output_time[static_cast<std::size_t>(level)] = stage_time;
}

void AuxiliaryStateManager::record_stage_face_transfer(const int level,
                                                       const StageContext& context,
                                                       const AuxiliaryFaceTransfer& stage_flux)
{
    face_transfer_ledger(level).record_stage(context, stage_flux);
}

} // namespace erf_auxiliary
