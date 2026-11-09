#include "ERF_SBMStateManager.H"
#include "ERF_SBMAMRTime.H"
#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

namespace erf_sbm {

SBMStateManager::SBMStateManager (SBMLayout layout, const int nlevels)
    : m_layout(std::move(layout)), m_projection(m_layout)
{
    if (nlevels <= 0) throw std::invalid_argument("SBM state manager needs at least one level");
    m_levels.resize(static_cast<std::size_t>(nlevels));
}

SBMStateManager::LevelState& SBMStateManager::level (const int lev)
{
    if (lev < 0 || lev >= nlevels()) throw std::out_of_range("SBM level is outside the manager");
    return m_levels[static_cast<std::size_t>(lev)];
}

const SBMStateManager::LevelState& SBMStateManager::level (const int lev) const
{
    if (lev < 0 || lev >= nlevels()) throw std::out_of_range("SBM level is outside the manager");
    return m_levels[static_cast<std::size_t>(lev)];
}

bool SBMStateManager::is_defined (const int lev) const
{
    const auto& state = level(lev);
    return state.old_state != nullptr && state.new_state != nullptr;
}

void SBMStateManager::define (const int lev, const amrex::BoxArray& grids,
                              const amrex::DistributionMapping& mapping,
                              const double current_time, const int ngrow)
{
    if (ngrow < 0) throw std::invalid_argument("SBM ghost width cannot be negative");
    if (!std::isfinite(current_time)) {
        throw std::invalid_argument("SBM current semantic time must be finite");
    }
    auto& state = level(lev);
    if (state.old_state || state.new_state) {
        throw std::logic_error("SBM level state is already defined");
    }

    state.old_state = std::make_unique<amrex::MultiFab>(
        grids, mapping, m_layout.ncomp(), ngrow);
    state.new_state = std::make_unique<amrex::MultiFab>(
        grids, mapping, m_layout.ncomp(), ngrow);
    // Both buffers are deterministic. Zero is a valid physical spectrum;
    // validity is carried only by the semantic flags below.
    state.old_state->setVal(amrex::Real(0.0));
    state.new_state->setVal(amrex::Real(0.0));
    state.old_valid = false;
    state.new_valid = true;
    state.old_time = 0.0;
    state.new_time = current_time;
    state.step_active = false;
}

void SBMStateManager::destroy (const int lev)
{
    level(lev) = LevelState{};
}

const amrex::MultiFab& SBMStateManager::old_state (const int lev) const
{
    const auto& state = level(lev);
    if (!state.old_state || !state.old_valid) {
        throw std::logic_error("SBM old spectral view is invalid at this lifecycle point");
    }
    return *state.old_state;
}

double SBMStateManager::old_time (const int lev) const
{
    static_cast<void>(old_state(lev));
    return level(lev).old_time;
}

bool SBMStateManager::old_valid (const int lev) const
{
    return level(lev).old_state != nullptr && level(lev).old_valid;
}

const amrex::MultiFab& SBMStateManager::new_state (const int lev) const
{
    const auto& state = level(lev);
    if (!state.new_state || !state.new_valid) {
        throw std::logic_error("SBM new spectral view is invalid until a stage is accepted");
    }
    return *state.new_state;
}

double SBMStateManager::new_time (const int lev) const
{
    static_cast<void>(new_state(lev));
    return level(lev).new_time;
}

bool SBMStateManager::new_valid (const int lev) const
{
    return level(lev).new_state != nullptr && level(lev).new_valid;
}

amrex::MultiFab& SBMStateManager::new_state_for_initialization (const int lev)
{
    auto& state = level(lev);
    if (!state.new_state || !state.new_valid || state.old_valid || state.step_active) {
        throw std::logic_error(
            "SBM accepted new state can be edited only before the first physical step");
    }
    return *state.new_state;
}

amrex::MultiFab& SBMStateManager::new_target_storage (const int lev)
{
    auto& state = level(lev);
    if (!state.new_state || !state.step_active) {
        throw std::logic_error("SBM new target storage requires an active physical step");
    }
    return *state.new_state;
}

bool SBMStateManager::begin_step (const int lev, const double step_old_time,
                                  std::string& diagnostic,
                                  const std::optional<double> dt_min)
{
    diagnostic.clear();
    auto& state = level(lev);
    if (!state.old_state || !state.new_state) {
        diagnostic = "SBM begin-step requested for an undefined level";
        return false;
    }
    if (!std::isfinite(step_old_time)) {
        diagnostic = "SBM begin-step time is nonfinite";
        return false;
    }
    if (state.step_active) {
        diagnostic = "SBM begin-step requested while the prior physical step is active";
        return false;
    }
    if (!state.new_valid) {
        diagnostic = "SBM begin-step has no accepted current/new spectrum";
        return false;
    }
    const bool same_time = dt_min.has_value()
        ? semantic_time_matches(state.new_time, step_old_time, *dt_min)
        : semantic_time_matches_exact(state.new_time, step_old_time);
    if (!same_time) {
        std::ostringstream message;
        message.precision(17);
        message << "SBM begin-step time mismatch at level " << lev
                << ": expected host step-old time " << step_old_time
                << ", accepted spectral time " << state.new_time
                << ", lifecycle=accepted new state, step inactive";
        diagnostic = message.str();
        return false;
    }

    // Match ERF::Advance: rotate storage by swapping ownership, never by
    // copying the full spectral state.  When the host scheduler's current
    // step-start expression was admitted as the same endpoint within the
    // bounded fine-step roundoff test above, the rotated old FAB now owns
    // that step-start semantic label.  The spectrum values are unchanged;
    // this only keeps the owner timestamp in lockstep with S_old/t_old at
    // stage 0.  With no dt_min, the exact-time path above still applies.
    std::swap(state.old_state, state.new_state);
    state.old_time = step_old_time;
    state.old_valid = true;
    state.new_valid = false;
    state.step_active = true;
    return true;
}

bool SBMStateManager::accept_stage_target (const int lev, const double target_time,
                                          const bool physical_step_complete,
                                          std::string& diagnostic)
{
    diagnostic.clear();
    auto& state = level(lev);
    if (!state.old_state || !state.new_state || !state.step_active || !state.old_valid) {
        diagnostic = "SBM stage target cannot be accepted without an active old/new lifecycle";
        return false;
    }
    if (!std::isfinite(target_time)) {
        diagnostic = "SBM stage target semantic time is nonfinite";
        return false;
    }
    state.new_time = target_time;
    state.new_valid = true;
    if (physical_step_complete) state.step_active = false;
    return true;
}

bool SBMStateManager::accept_synchronized_state (
    const int lev, const amrex::MultiFab& candidate,
    const double semantic_time, std::string& diagnostic)
{
    if (!can_accept_synchronized_state(lev, candidate, semantic_time,
                                       diagnostic)) {
        return false;
    }
    publish_synchronized_state(lev, candidate);
    return true;
}

bool SBMStateManager::can_accept_synchronized_state (
    const int lev, const amrex::MultiFab& candidate,
    const double semantic_time, std::string& diagnostic) const
{
    diagnostic.clear();
    const auto& state = level(lev);
    const double scale = std::max({1.0, std::abs(semantic_time),
                                   std::abs(state.new_time)});
    const bool time_matches = std::isfinite(semantic_time) &&
        std::abs(semantic_time - state.new_time) <=
            64.0 * std::numeric_limits<double>::epsilon() * scale;
    const bool local_ok = state.new_state != nullptr && state.new_valid &&
        !state.step_active && time_matches &&
        candidate.nComp() == m_layout.ncomp() &&
        erf_auxiliary::SameCellLayout(candidate, *state.new_state);
    int vote = local_ok ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntMin(vote);
    if (vote == 0) {
        diagnostic = local_ok
            ? "SBM synchronized candidate is invalid on another MPI rank"
            : "SBM synchronized candidate must match the accepted level layout and time";
        return false;
    }
    return true;
}

void SBMStateManager::publish_synchronized_state (
    const int lev, const amrex::MultiFab& candidate)
{
    auto& state = level(lev);
    AMREX_ALWAYS_ASSERT(state.new_state != nullptr && state.new_valid &&
                        !state.step_active &&
                        candidate.nComp() == m_layout.ncomp() &&
                        erf_auxiliary::SameCellLayout(candidate,
                                                      *state.new_state));
    amrex::MultiFab::Copy(*state.new_state, candidate, 0, 0,
                          m_layout.ncomp(), 0);
}

bool SBMStateManager::step_active (const int lev) const
{
    return level(lev).step_active;
}

void SBMStateManager::project_to_core (const int lev, amrex::MultiFab& core,
                                       const int qc_component,
                                       const int qr_component) const
{
    const auto& spectrum = new_state(lev);
    if (qc_component < 0 || qr_component < 0 ||
        qc_component >= core.nComp() || qr_component >= core.nComp()) {
        throw std::invalid_argument("SBM compact projection target is outside core state");
    }
    if (!(spectrum.boxArray() == core.boxArray()) ||
        !(spectrum.DistributionMap() == core.DistributionMap())) {
        throw std::invalid_argument("SBM spectrum and core state must share a level layout");
    }
    for (amrex::MFIter mfi(core); mfi.isValid(); ++mfi) {
        m_projection.apply_to_core(mfi.validbox(), spectrum.const_array(mfi),
                                   core.array(mfi), qc_component, qr_component);
    }
}

} // namespace erf_sbm
