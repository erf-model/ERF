#include "ERF_SBMAMRFluxRegister.H"
#include "ERF_SBMAMRTime.H"

#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace erf_sbm {
namespace {

} // namespace

void
SBMAMRFluxRegister::reset_side (SideState& side) const
{
    side = {};
    side.components_seen.assign(static_cast<std::size_t>(m_ncomp), 0);
}

void
SBMAMRFluxRegister::define (const amrex::BoxArray& fine_ba,
                            const amrex::BoxArray& coarse_ba,
                            const amrex::DistributionMapping& fine_dm,
                            const amrex::DistributionMapping& coarse_dm,
                            const amrex::Geometry& fine_geom,
                            const amrex::Geometry& coarse_geom,
                            const amrex::IntVect& ratio,
                            const int fine_level,
                            const int ncomp)
{
    clear();
    AMREX_ALWAYS_ASSERT(fine_level > 0 && ncomp > 0);
    m_fine_ba = fine_ba;
    m_coarse_ba = coarse_ba;
    m_fine_dm = fine_dm;
    m_coarse_dm = coarse_dm;
    m_fine_geom = fine_geom;
    m_coarse_geom = coarse_geom;
    m_ratio = ratio;
    m_fine_level = fine_level;
    m_ncomp = ncomp;
    m_reflux_components_seen.assign(static_cast<std::size_t>(m_ncomp), 0);
    m_register = std::make_unique<amrex::YAFluxRegister>(
        m_fine_ba, m_coarse_ba, m_fine_dm, m_coarse_dm, m_fine_geom,
        m_coarse_geom, m_ratio, m_fine_level, m_ncomp);
    amrex::MultiFab fine_coverage(m_fine_ba, m_fine_dm, 1, 0);
    m_covered_parent_mask.define(m_coarse_ba, m_coarse_dm, 1, 1);
    fine_coverage.setVal(amrex::Real(1.0));
    m_covered_parent_mask.setVal(amrex::Real(0.0));
    amrex::average_down(fine_coverage, m_covered_parent_mask, 0, 1, m_ratio);
    m_covered_parent_mask.FillBoundary(m_coarse_geom.periodicity());
    reset_side(m_coarse);
    reset_side(m_fine);
}

void
SBMAMRFluxRegister::clear ()
{
    if (m_register) { m_register->clear(); }
    m_register.reset();
    m_covered_parent_mask.clear();
    m_fine_ba.clear();
    m_coarse_ba.clear();
    m_fine_dm = amrex::DistributionMapping();
    m_coarse_dm = amrex::DistributionMapping();
    m_fine_geom = amrex::Geometry();
    m_coarse_geom = amrex::Geometry();
    m_ratio = amrex::IntVect(1);
    m_fine_level = 0;
    m_ncomp = 0;
    m_interval_active = false;
    m_coarse_step_accepted = false;
    m_reflux_applied = false;
    m_reflux_components_seen.clear();
    m_consumed = true;
    m_interval_begin = 0.0;
    m_interval_end = 0.0;
    m_min_fine_dt = 0.0;
    m_expected_fine_steps = 0;
    m_accepted_fine_steps = 0;
    m_reflux_call_count = 0;
    m_reflux_sync_count = 0;
    m_coarse = {};
    m_fine = {};
}

int
SBMAMRFluxRegister::stage_count (
    const erf_auxiliary::HostIntegrator method) const noexcept
{
    return erf_auxiliary::HostIntegratorStageCount(method);
}

bool
SBMAMRFluxRegister::collective_check (
    const bool local_ok, const std::string& local_diagnostic,
    const char* remote_diagnostic, std::string& diagnostic) const
{
    int vote = local_ok ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntMin(vote);
    if (vote != 0) {
        diagnostic.clear();
        return true;
    }
    diagnostic = local_ok ? remote_diagnostic : local_diagnostic;
    return false;
}

bool
SBMAMRFluxRegister::begin_interval (const double coarse_begin,
                                    const double coarse_end,
                                    const int expected_fine_steps,
                                    std::string& diagnostic)
{
    diagnostic.clear();
    const bool local_ok = m_register != nullptr && !m_interval_active &&
        m_consumed && std::isfinite(coarse_begin) &&
        std::isfinite(coarse_end) && coarse_end > coarse_begin &&
        std::isfinite(coarse_end - coarse_begin) && expected_fine_steps > 0 &&
        std::isfinite((coarse_end - coarse_begin) /
                      static_cast<double>(expected_fine_steps)) &&
        (coarse_end - coarse_begin) / static_cast<double>(expected_fine_steps) > 0.0;
    const std::string local_diagnostic = m_register == nullptr
        ? "SBM interface register is undefined"
        : m_interval_active || !m_consumed
            ? "SBM interface register still owns an unconsumed interval"
            : "SBM interface interval has invalid times or fine-substep count";
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM interface interval is invalid on another MPI rank",
                          diagnostic)) {
        return false;
    }
    m_register->reset();
    m_interval_begin = coarse_begin;
    m_interval_end = coarse_end;
    m_expected_fine_steps = expected_fine_steps;
    m_min_fine_dt = (coarse_end - coarse_begin) /
                    static_cast<double>(expected_fine_steps);
    m_reflux_call_count = 0;
    m_reflux_sync_count = 0;
    m_accepted_fine_steps = 0;
    m_coarse_step_accepted = false;
    m_reflux_applied = false;
    std::fill(m_reflux_components_seen.begin(),
              m_reflux_components_seen.end(), 0);
    m_interval_active = true;
    m_consumed = false;
    reset_side(m_coarse);
    reset_side(m_fine);
    return true;
}

bool
SBMAMRFluxRegister::begin_level_step (const bool coarse_side,
                                      const double step_begin,
                                      const double step_end,
                                      std::string& diagnostic)
{
    diagnostic.clear();
    SideState& side = coarse_side ? m_coarse : m_fine;
    bool local_ok = m_interval_active && !m_consumed && !side.step_active &&
        std::isfinite(step_begin) && std::isfinite(step_end) &&
        step_end > step_begin;
    std::string local_diagnostic = !m_interval_active || m_consumed
        ? "SBM level step arrived outside an active interface interval"
        : side.step_active ? "SBM interface level step is already active"
                           : "SBM interface level step has invalid times";
    if (local_ok && coarse_side) {
        local_ok = !m_coarse_step_accepted &&
                   semantic_time_matches(step_begin, m_interval_begin,
                                         m_min_fine_dt) &&
                   semantic_time_matches(step_end, m_interval_end,
                                         m_min_fine_dt);
        if (!local_ok) {
            local_diagnostic = "coarse SBM step does not match its interface interval";
        }
    } else if (local_ok) {
        const double fine_dt = (m_interval_end - m_interval_begin) /
                               static_cast<double>(m_expected_fine_steps);
        const double expected_begin = m_interval_begin +
            static_cast<double>(m_accepted_fine_steps) * fine_dt;
        const double expected_end = expected_begin + fine_dt;
        local_ok = m_coarse_step_accepted &&
            m_accepted_fine_steps < m_expected_fine_steps &&
            semantic_time_matches(step_begin, expected_begin, m_min_fine_dt) &&
            semantic_time_matches(step_end, expected_end, m_min_fine_dt);
        if (!local_ok) {
            local_diagnostic = "fine SBM substep is missing, duplicated, or outside the coarse interval";
        }
    }
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM interface level step is invalid on another MPI rank",
                          diagnostic)) {
        return false;
    }
    side.step_active = true;
    side.stage_active = false;
    side.step_accepted = false;
    side.step_begin = step_begin;
    side.step_end = step_end;
    side.current_stage_time = step_begin;
    side.next_stage = 0;
    side.components_seen.assign(static_cast<std::size_t>(m_ncomp), 0);
    return true;
}

bool
SBMAMRFluxRegister::begin_stage (
    const bool coarse_side, const erf_auxiliary::HostIntegrator method,
    const int stage, const double step_begin, std::string& diagnostic)
{
    diagnostic.clear();
    SideState& side = coarse_side ? m_coarse : m_fine;
    const int count = stage_count(method);
    const bool local_ok = m_interval_active && side.step_active &&
        !side.stage_active && stage == side.next_stage && count > 0 &&
        stage < count && semantic_time_matches(side.step_begin, step_begin,
                                               m_min_fine_dt) &&
        (stage == 0 || side.method == method);
    std::string local_diagnostic = !m_interval_active || !side.step_active
        ? "SBM interface stage arrived without an active physical step"
        : side.stage_active ? "SBM interface stage is already active"
        : count == 0 ? "SBM interface does not support this host integrator"
        : stage != side.next_stage
            ? "SBM interface stage is duplicate, skipped, or out of order"
                : !semantic_time_matches(side.step_begin, step_begin,
                                         m_min_fine_dt)
                ? "SBM interface stage changed its physical-step start time"
                : "SBM interface host integrator changed during a step";
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM interface stage is invalid on another MPI rank",
                          diagnostic)) {
        return false;
    }
    if (stage == 0) { side.method = method; }
    side.stage_active = true;
    side.components_seen.assign(static_cast<std::size_t>(m_ncomp), 0);
    return true;
}

bool
SBMAMRFluxRegister::add_stage_rates (
    const bool coarse_side,
    const erf_auxiliary::MappedFaceFluxRate& accepted_rate,
    const std::vector<int>& register_components,
    const double completed_ledger_time,
    std::string& diagnostic)
{
    diagnostic.clear();
    SideState& side = coarse_side ? m_coarse : m_fine;
    bool local_ok = m_register != nullptr && side.stage_active &&
        std::isfinite(completed_ledger_time) && completed_ledger_time >= 0.0 &&
        erf_auxiliary::MappedFaceLayoutMatchesCellLayout(
            accepted_rate, coarse_side ? m_coarse_ba : m_fine_ba,
            coarse_side ? m_coarse_dm : m_fine_dm);
    std::string local_diagnostic = m_register == nullptr
        ? "SBM interface register is undefined"
        : !side.stage_active ? "accepted SBM face rates arrived outside a stage"
        : !std::isfinite(completed_ledger_time) || completed_ledger_time < 0.0
            ? "SBM accepted face rates have an invalid completed-step coefficient"
            : "accepted SBM face-rate layout does not match the interface grid";
    if (local_ok && (register_components.empty() || register_components.size() >
            static_cast<std::size_t>(accepted_rate.nComp()))) {
        local_ok = false;
        local_diagnostic = "SBM chunk component map does not match its accepted-rate FAB";
    }
    if (local_ok) {
        std::vector<unsigned char> local_seen(static_cast<std::size_t>(m_ncomp), 0);
        for (const int component : register_components) {
            if (component < 0 || component >= m_ncomp ||
                side.components_seen[static_cast<std::size_t>(component)] != 0 ||
                local_seen[static_cast<std::size_t>(component)] != 0) {
                local_ok = false;
                local_diagnostic = "SBM chunk contains an invalid or duplicate global component";
                break;
            }
            local_seen[static_cast<std::size_t>(component)] = 1;
        }
    }
    if (!collective_check(local_ok, local_diagnostic,
                          "accepted SBM face-rate chunk is invalid on another MPI rank",
                          diagnostic)) {
        return false;
    }

    if (completed_ledger_time != 0.0) {
        const auto& grids = coarse_side ? m_coarse_ba : m_fine_ba;
        const auto& mapping = coarse_side ? m_coarse_dm : m_fine_dm;
        amrex::MultiFab cell_layout(grids, mapping, 1, 0,
                                    amrex::MFInfo().SetAlloc(false));
        const auto dx = (coarse_side ? m_coarse_geom : m_fine_geom).CellSizeArray();
        const amrex::RunOn run_on = amrex::Gpu::inLaunchRegion()
            ? amrex::RunOn::Gpu : amrex::RunOn::Host;
        for (amrex::MFIter mfi(cell_layout, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            std::array<amrex::FArrayBox const*, AMREX_SPACEDIM> fluxes{};
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                fluxes[static_cast<std::size_t>(dir)] =
                    &accepted_rate.dir(dir)[mfi];
            }
            for (std::size_t first = 0; first < register_components.size();) {
                std::size_t end = first + 1U;
                while (end < register_components.size() &&
                       register_components[end] ==
                           register_components[end - 1U] + 1) {
                    ++end;
                }
                const int local = static_cast<int>(first);
                const int global_component = register_components[first];
                const int count = static_cast<int>(end - first);
                if (coarse_side) {
                    m_register->CrseAdd(mfi, fluxes, dx.data(),
                        static_cast<amrex::Real>(completed_ledger_time),
                        local, global_component, count, run_on);
                } else {
                    m_register->FineAdd(mfi, fluxes, dx.data(),
                        static_cast<amrex::Real>(completed_ledger_time),
                        local, global_component, count, run_on);
                }
                first = end;
            }
        }
        // accepted_rate is reused by the next chunk; finish all YA readers
        // before that storage can be overwritten on an asynchronous backend.
        amrex::Gpu::streamSynchronize();
    }
    for (const int component : register_components) {
        side.components_seen[static_cast<std::size_t>(component)] = 1;
    }
    return true;
}

bool
SBMAMRFluxRegister::finish_stage (const bool coarse_side,
                                  std::string& diagnostic)
{
    diagnostic.clear();
    SideState& side = coarse_side ? m_coarse : m_fine;
    const bool all_components_seen = std::all_of(
        side.components_seen.begin(), side.components_seen.end(),
        [](const unsigned char seen) { return seen != 0; });
    const bool local_ok = side.stage_active && all_components_seen;
    const std::string local_diagnostic = !side.stage_active
        ? "SBM interface stage finish requires an open stage"
        : "SBM interface stage omitted one or more global spectral components";
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM interface stage is incomplete on another MPI rank",
                          diagnostic)) {
        return false;
    }
    side.stage_active = false;
    ++side.next_stage;
    return true;
}

bool
SBMAMRFluxRegister::accept_level_step (const bool coarse_side,
                                       const double step_begin,
                                       const double step_end,
                                       std::string& diagnostic)
{
    diagnostic.clear();
    SideState& side = coarse_side ? m_coarse : m_fine;
    const bool stages_complete = side.next_stage == stage_count(side.method);
    bool local_ok = side.step_active && !side.stage_active && stages_complete &&
        semantic_time_matches(step_begin, side.step_begin, m_min_fine_dt) &&
        semantic_time_matches(step_end, side.step_end, m_min_fine_dt);
    std::string local_diagnostic = !side.step_active
        ? "SBM interface level-step acceptance has no active step"
        : side.stage_active ? "SBM interface level-step acceptance found an open stage"
        : !stages_complete ? "SBM interface level step omitted one or more host stages"
        : "SBM interface level-step acceptance changed its semantic times";
    if (local_ok && coarse_side) {
        local_ok = semantic_time_matches(step_begin, m_interval_begin,
                                         m_min_fine_dt) &&
                   semantic_time_matches(step_end, m_interval_end,
                                         m_min_fine_dt) &&
                   !m_coarse_step_accepted;
        if (!local_ok) {
            local_diagnostic = "accepted coarse SBM step does not match its interval";
        }
    } else if (local_ok) {
        const double fine_dt = (m_interval_end - m_interval_begin) /
                               static_cast<double>(m_expected_fine_steps);
        const double expected_begin = m_interval_begin +
            static_cast<double>(m_accepted_fine_steps) * fine_dt;
        local_ok = m_coarse_step_accepted &&
            semantic_time_matches(step_begin, expected_begin, m_min_fine_dt) &&
            semantic_time_matches(step_end, expected_begin + fine_dt,
                                  m_min_fine_dt) &&
            m_accepted_fine_steps < m_expected_fine_steps;
        if (!local_ok) {
            local_diagnostic = "accepted fine SBM step is missing, duplicated, or out of order";
        }
    }
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM interface level step was rejected on another MPI rank",
                          diagnostic)) {
        return false;
    }
    side.step_active = false;
    side.step_accepted = true;
    if (coarse_side) {
        m_coarse_step_accepted = true;
    } else {
        ++m_accepted_fine_steps;
    }
    return true;
}

bool
SBMAMRFluxRegister::ready_to_reflux () const noexcept
{
    return m_interval_active && !m_consumed && m_coarse_step_accepted &&
        m_accepted_fine_steps == m_expected_fine_steps &&
        !m_coarse.step_active && !m_fine.step_active &&
        !m_coarse.stage_active && !m_fine.stage_active;
}

bool
SBMAMRFluxRegister::has_outstanding_amount () const noexcept
{
    return m_interval_active && !m_consumed;
}

bool
SBMAMRFluxRegister::reflux (amrex::MultiFab& mapped_candidate,
                            std::string& diagnostic)
{
    diagnostic.clear();
    const bool candidate_matches = m_register != nullptr &&
        mapped_candidate.nComp() == m_ncomp &&
        mapped_candidate.boxArray() == m_coarse_ba &&
        mapped_candidate.DistributionMap() == m_coarse_dm;
    const bool no_component_was_refluxed = std::none_of(
        m_reflux_components_seen.begin(), m_reflux_components_seen.end(),
        [] (const unsigned char seen) { return seen != 0; });
    const bool local_ok = ready_to_reflux() && !m_reflux_applied &&
                          candidate_matches && no_component_was_refluxed;
    const std::string local_diagnostic = !ready_to_reflux()
        ? "SBM YA reflux requires one accepted coarse step and every accepted fine substep"
        : m_reflux_applied || !no_component_was_refluxed
            ? "SBM YA reflux candidate was already applied for this pending interval"
        : "SBM YA reflux candidate does not match the coarse full-component layout";
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM YA reflux preflight failed on another MPI rank",
                          diagnostic)) {
        return false;
    }
    m_register->Reflux(mapped_candidate, 0, 0, m_ncomp);
    amrex::Gpu::streamSynchronize();
    std::fill(m_reflux_components_seen.begin(),
              m_reflux_components_seen.end(), 1);
    m_reflux_applied = true;
    return true;
}

bool SBMAMRFluxRegister::reflux_components (
    amrex::MultiFab& mapped_correction,
    const std::vector<int>& global_components,
    std::string& diagnostic)
{
    diagnostic.clear();
    bool mapping_ok = m_register != nullptr && !global_components.empty() &&
        mapped_correction.nComp() >= static_cast<int>(global_components.size()) &&
        mapped_correction.boxArray() == m_coarse_ba &&
        mapped_correction.DistributionMap() == m_coarse_dm &&
        m_reflux_components_seen.size() == static_cast<std::size_t>(m_ncomp);
    std::vector<int> seen_now(static_cast<std::size_t>(m_ncomp), 0);
    for (const int component : global_components) {
        if (component < 0 || component >= m_ncomp ||
            (component >= 0 && component < m_ncomp &&
             (m_reflux_components_seen[static_cast<std::size_t>(component)] != 0 ||
              seen_now[static_cast<std::size_t>(component)] != 0))) {
            mapping_ok = false;
        } else {
            seen_now[static_cast<std::size_t>(component)] = 1;
        }
    }
    const bool local_ok = ready_to_reflux() && !m_reflux_applied && mapping_ok;
    const std::string local_diagnostic = !ready_to_reflux()
        ? "SBM YA reflux requires one accepted coarse step and every accepted fine substep"
        : m_reflux_applied
            ? "SBM YA reflux candidate was already applied for this pending interval"
        : "SBM YA correction chunk has an invalid layout or repeated component mapping";
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM YA reflux preflight failed on another MPI rank",
                          diagnostic)) {
        return false;
    }

    std::size_t first = 0;
    while (first < global_components.size()) {
        std::size_t end = first + 1U;
        while (end < global_components.size() &&
               global_components[end] == global_components[end - 1U] + 1) {
            ++end;
        }
        const int source_component = global_components[first];
        const int destination_component = static_cast<int>(first);
        const int count = static_cast<int>(end - first);
        m_register->Reflux(mapped_correction, source_component,
                           destination_component, count);
        ++m_reflux_call_count;
        amrex::Gpu::streamSynchronize();
        ++m_reflux_sync_count;
        for (std::size_t index = first; index < end; ++index) {
            m_reflux_components_seen[static_cast<std::size_t>(
                global_components[index])] = 1;
        }
        first = end;
    }
    m_reflux_applied = std::all_of(m_reflux_components_seen.begin(),
                                   m_reflux_components_seen.end(),
                                   [] (const unsigned char seen) { return seen != 0; });
    return true;
}

bool
SBMAMRFluxRegister::discard_interval (std::string& diagnostic)
{
    diagnostic.clear();
    const bool local_ok = m_register != nullptr && m_interval_active &&
                          !m_consumed;
    const std::string local_diagnostic = m_register == nullptr
        ? "SBM interface register is undefined"
        : "SBM interface interval has no pending amount to discard";
    if (!collective_check(local_ok, local_diagnostic,
                          "SBM rejected interface interval is invalid on another MPI rank",
                          diagnostic)) {
        return false;
    }
    m_register->reset();
    m_interval_active = false;
    m_consumed = true;
    m_reflux_applied = false;
    m_coarse_step_accepted = false;
    m_min_fine_dt = 0.0;
    m_accepted_fine_steps = 0;
    std::fill(m_reflux_components_seen.begin(),
              m_reflux_components_seen.end(), 0);
    reset_side(m_coarse);
    reset_side(m_fine);
    return true;
}

bool
SBMAMRFluxRegister::can_accept_reflux (std::string& diagnostic) const
{
    diagnostic.clear();
    const bool all_components_refluxed = std::all_of(
        m_reflux_components_seen.begin(), m_reflux_components_seen.end(),
        [] (const unsigned char seen) { return seen != 0; });
    const bool local_ok = ready_to_reflux() && m_reflux_applied &&
                          all_components_refluxed;
    const std::string local_diagnostic = !ready_to_reflux()
        ? "SBM YA amount cannot be consumed before all level steps are accepted"
        : "SBM YA amount cannot be consumed before a reflux candidate was formed";
    return collective_check(local_ok, local_diagnostic,
                            "SBM YA amount consumption was rejected on another MPI rank",
                            diagnostic);
}

bool
SBMAMRFluxRegister::accept_reflux (std::string& diagnostic)
{
    if (!can_accept_reflux(diagnostic)) { return false; }
    m_register->reset();
    m_consumed = true;
    m_interval_active = false;
    m_reflux_applied = false;
    m_coarse_step_accepted = false;
    m_min_fine_dt = 0.0;
    m_accepted_fine_steps = 0;
    std::fill(m_reflux_components_seen.begin(),
              m_reflux_components_seen.end(), 0);
    reset_side(m_coarse);
    reset_side(m_fine);
    return true;
}

} // namespace erf_sbm
