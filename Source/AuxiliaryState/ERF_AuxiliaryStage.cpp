#include "ERF_AuxiliaryStage.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace erf_auxiliary {
namespace {

int stage_count (const HostIntegrator method) noexcept
{
    switch (method) {
    case HostIntegrator::CompressibleRK3: return 3;
    case HostIntegrator::AnelasticHeun: return 2;
    case HostIntegrator::AnelasticMidPoint: return 2;
    }
    return 0;
}

bool finite_nonnegative (const double value) noexcept
{
    return std::isfinite(value) && value >= 0.0;
}

bool valid_const_view (const ConstTimedFieldView& view,
                       const double expected_time) noexcept
{
    return view.field != nullptr && view.component >= 0 &&
           view.component < view.field->nComp() &&
           std::isfinite(view.time) && view.time == expected_time;
}

bool storage_overlaps (const amrex::MultiFab& lhs, const amrex::MultiFab& rhs)
{
    if (&lhs == &rhs) { return true; }
    if (!SameCellLayout(lhs, rhs)) { return false; }
    for (amrex::MFIter mfi(lhs); mfi.isValid(); ++mfi) {
        const auto& lhs_fab = lhs[mfi];
        const auto& rhs_fab = rhs[mfi];
        const auto lhs_begin = reinterpret_cast<std::uintptr_t>(lhs_fab.dataPtr());
        const auto rhs_begin = reinterpret_cast<std::uintptr_t>(rhs_fab.dataPtr());
        const auto lhs_size = static_cast<std::uintptr_t>(lhs_fab.box().numPts()) *
                              static_cast<std::uintptr_t>(lhs_fab.nComp()) * sizeof(amrex::Real);
        const auto rhs_size = static_cast<std::uintptr_t>(rhs_fab.box().numPts()) *
                              static_cast<std::uintptr_t>(rhs_fab.nComp()) * sizeof(amrex::Real);
        if (lhs_begin < rhs_begin + rhs_size && rhs_begin < lhs_begin + lhs_size) {
            return true;
        }
    }
    return false;
}

} // namespace

const char* HostIntegratorName (const HostIntegrator method) noexcept
{
    switch (method) {
    case HostIntegrator::CompressibleRK3: return "CompressibleRK3";
    case HostIntegrator::AnelasticHeun: return "AnelasticHeun";
    case HostIntegrator::AnelasticMidPoint: return "AnelasticMidPoint";
    }
    return "Unknown";
}

bool MakeAuxiliaryStageRecipe (const HostIntegrator method,
                               const int stage,
                               const double host_stage_interval,
                               AuxiliaryStageRecipe& recipe,
                               std::string& diagnostic)
{
    diagnostic.clear();
    recipe = {};
    if (!std::isfinite(host_stage_interval) || !(host_stage_interval > 0.0)) {
        diagnostic = "host stage interval must be finite and positive";
        return false;
    }
    if (method == HostIntegrator::AnelasticMidPoint) {
        diagnostic = "AnelasticMidPoint is represented but is not qualified by the M2 fixture";
        return false;
    }
    if (stage < 0 || stage >= stage_count(method)) {
        diagnostic = "stage index is outside the selected host integrator sequence";
        return false;
    }

    if (method == HostIntegrator::CompressibleRK3) {
        recipe.limiter_trial_base = LimiterTrialBase::Anchor;
        recipe.anchor_weight = 1.0;
        recipe.limiter_trial_interval = host_stage_interval;
        recipe.face_rate_time_coefficient = host_stage_interval;
        if (stage == 2) { recipe.completed_ledger_time = host_stage_interval; }
    } else if (method == HostIntegrator::AnelasticHeun) {
        if (stage == 0) {
            recipe.limiter_trial_base = LimiterTrialBase::Anchor;
            recipe.anchor_weight = 1.0;
            recipe.limiter_trial_interval = host_stage_interval;
            recipe.face_rate_time_coefficient = host_stage_interval;
            recipe.completed_ledger_time = 0.5 * host_stage_interval;
        } else {
            recipe.limiter_trial_base = LimiterTrialBase::Input;
            recipe.anchor_weight = 0.5;
            recipe.input_weight = 0.5;
            recipe.limiter_trial_interval = host_stage_interval;
            recipe.face_rate_time_coefficient = 0.5 * host_stage_interval;
            recipe.completed_ledger_time = 0.5 * host_stage_interval;
        }
    }
    return true;
}

bool BuildAuxiliaryIntensiveState (const ConstTimedFieldView& state_input,
                                   const ConstTimedFieldView& rho_input,
                                   const double input_time,
                                   amrex::MultiFab& intensive,
                                   const AuxiliaryFieldValidationPolicy validation,
                                   std::string& diagnostic)
{
    diagnostic.clear();
    if (!std::isfinite(input_time) ||
        !valid_const_view(state_input, input_time) ||
        !valid_const_view(rho_input, input_time)) {
        diagnostic = "auxiliary intensive state requires state and density views at input_time";
        return false;
    }
    if (!SameCellLayout(intensive, *state_input.field) ||
        !SameCellLayout(intensive, *rho_input.field) ||
        intensive.nComp() < 1) {
        diagnostic = "auxiliary intensive state fields must share a cell layout";
        return false;
    }
    if (validation == AuxiliaryFieldValidationPolicy::Global) {
        if (!ValidatePositiveFiniteComponent(*rho_input.field, rho_input.component, diagnostic)) {
            return false;
        }
    } else if (validation != AuxiliaryFieldValidationPolicy::AssumeValid) {
        diagnostic = "auxiliary intensive state received an unknown validation policy";
        return false;
    }
    for (amrex::MFIter mfi(intensive, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto state = state_input.field->const_array(mfi);
        const auto rho = rho_input.field->const_array(mfi);
        const auto out = intensive.array(mfi);
        const int state_comp = state_input.component;
        const int rho_comp = rho_input.component;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = state(i, j, k, state_comp) / rho(i, j, k, rho_comp);
        });
    }
    if (validation == AuxiliaryFieldValidationPolicy::Global) {
        if (!ValidateFiniteComponent(intensive, 0, diagnostic)) {
            return false;
        }
    }
    return true;
}

bool AuxiliaryStageTargetIsDisjoint (const AuxiliaryStageContext& context)
{
    return context.state_target.field != nullptr &&
           (context.state_anchor.field == nullptr ||
            !storage_overlaps(*context.state_target.field, *context.state_anchor.field)) &&
           (context.state_input.field == nullptr ||
            !storage_overlaps(*context.state_target.field, *context.state_input.field));
}

void ApplyAuxiliaryMappedStage (const AuxiliaryStageContext& context,
                                const MappedFaceFluxRate& rate,
                                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
                                const int rate_component)
{
    AMREX_ALWAYS_ASSERT(context.level >= 0 && context.stage >= 0);
    AMREX_ALWAYS_ASSERT(std::isfinite(context.step_old_time));
    AMREX_ALWAYS_ASSERT(std::isfinite(context.input_time));
    AMREX_ALWAYS_ASSERT(std::isfinite(context.target_time));
    AMREX_ALWAYS_ASSERT(valid_const_view(context.state_anchor, context.step_old_time));
    AMREX_ALWAYS_ASSERT(valid_const_view(context.state_input, context.input_time));
    AMREX_ALWAYS_ASSERT(context.state_target.field != nullptr);
    AMREX_ALWAYS_ASSERT(context.state_target.component >= 0 &&
                        context.state_target.component < context.state_target.field->nComp());
    AMREX_ALWAYS_ASSERT(std::isfinite(context.state_target.time) &&
                        context.state_target.time == context.target_time);
    AMREX_ALWAYS_ASSERT(valid_const_view(context.rho_anchor, context.step_old_time));
    AMREX_ALWAYS_ASSERT(valid_const_view(context.rho_input, context.input_time));
    AMREX_ALWAYS_ASSERT(valid_const_view(context.rho_target, context.target_time));
    AMREX_ALWAYS_ASSERT(valid_const_view(context.measure_anchor, context.step_old_time));
    AMREX_ALWAYS_ASSERT(valid_const_view(context.measure_input, context.input_time));
    AMREX_ALWAYS_ASSERT(valid_const_view(context.measure_target, context.target_time));
    AMREX_ALWAYS_ASSERT(context.carrier.x != nullptr && context.carrier.y != nullptr &&
                        context.carrier.z != nullptr);
    AMREX_ALWAYS_ASSERT(rate.is_defined() && rate_component >= 0 &&
                        rate_component < rate.nComp());
    AMREX_ALWAYS_ASSERT(AuxiliaryStageTargetIsDisjoint(context));

    const auto& target = *context.state_target.field;
    const auto& anchor = *context.state_anchor.field;
    const auto& input = *context.state_input.field;
    const auto& omega_anchor = *context.measure_anchor.field;
    const auto& omega_input = *context.measure_input.field;
    const auto& omega_target = *context.measure_target.field;
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, anchor));
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, input));
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, *context.rho_anchor.field));
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, *context.rho_input.field));
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, *context.rho_target.field));
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, omega_anchor));
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, omega_input));
    AMREX_ALWAYS_ASSERT(SameCellLayout(target, omega_target));
    AMREX_ALWAYS_ASSERT(MappedFaceLayoutMatchesCellLayout(rate, target));

    const amrex::Real anchor_weight = static_cast<amrex::Real>(context.recurrence.anchor_weight);
    const amrex::Real input_weight = static_cast<amrex::Real>(context.recurrence.input_weight);
    const amrex::Real face_weight =
        static_cast<amrex::Real>(context.recurrence.face_rate_time_coefficient);
    const amrex::Real dx = dx_inv[0];
    const amrex::Real dy = dx_inv[1];
    const amrex::Real dz = dx_inv[2];

    for (amrex::MFIter mfi(target, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto ua = anchor.const_array(mfi);
        const auto ui = input.const_array(mfi);
        const auto ut = context.state_target.field->array(mfi);
        const auto oa = omega_anchor.const_array(mfi);
        const auto oi = omega_input.const_array(mfi);
        const auto omega_target_array = omega_target.const_array(mfi);
        const auto fx = rate.dir(0).const_array(mfi);
        const auto fy = rate.dir(1).const_array(mfi);
        const auto fz = rate.dir(2).const_array(mfi);
        const int ua_comp = context.state_anchor.component;
        const int ui_comp = context.state_input.component;
        const int ut_comp = context.state_target.component;
        const int oa_comp = context.measure_anchor.component;
        const int oi_comp = context.measure_input.component;
        const int omega_target_component = context.measure_target.component;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real div = ComputationalMappedDivergence(
                fx(i + 1, j, k, rate_component), fx(i, j, k, rate_component),
                fy(i, j + 1, k, rate_component), fy(i, j, k, rate_component),
                fz(i, j, k + 1, rate_component), fz(i, j, k, rate_component),
                dx, dy, dz);
            const amrex::Real h_target =
                anchor_weight * oa(i, j, k, oa_comp) * ua(i, j, k, ua_comp) +
                input_weight * oi(i, j, k, oi_comp) * ui(i, j, k, ui_comp) -
                face_weight * div;
            ut(i, j, k, ut_comp) = h_target /
                omega_target_array(i, j, k, omega_target_component);
        });
    }
}

void CompletedStepFluxLedger::define (const amrex::BoxArray& cell_ba,
                                    const amrex::DistributionMapping& dm,
                                    const int ncomp)
{
    m_integral.define(cell_ba, dm, ncomp, 0);
    m_integral.setVal(amrex::Real(0.0));
    m_step_active = false;
    m_step_complete = false;
    m_stage_open = false;
    m_next_stage = 0;
    m_open_stage = -1;
    m_open_stage_weight = amrex::Real(0.0);
    m_step_old_time = 0.0;
    m_stage_components_seen.assign(static_cast<std::size_t>(ncomp), 0);
    m_stage_failure_diagnostic.clear();
}

void CompletedStepFluxLedger::latch_stage_failure (
    const std::string& diagnostic)
{
    if (m_stage_failure_diagnostic.empty()) {
        m_stage_failure_diagnostic = diagnostic.empty()
                                         ? "completed-step ledger stage operation failed"
                                         : diagnostic;
    }
}

void CompletedStepFluxLedger::discard_step ()
{
    if (m_integral.is_defined()) {
        m_integral.setVal(amrex::Real(0.0));
    }
    m_step_active = false;
    m_step_complete = false;
    m_stage_open = false;
    m_next_stage = 0;
    m_open_stage = -1;
    m_open_stage_weight = amrex::Real(0.0);
    m_step_old_time = 0.0;
    std::fill(m_stage_components_seen.begin(), m_stage_components_seen.end(), 0);
    m_stage_failure_diagnostic.clear();
}

bool CompletedStepFluxLedger::accept_stage (const HostIntegrator method,
                                            const int stage,
                                            const double step_old_time,
                                            const AuxiliaryStageRecipe& recipe,
                                            const MappedFaceFluxRate& rate,
                                            std::string& diagnostic)
{
    diagnostic.clear();
    bool local_stage_ok = true;
    if (!is_defined() || !rate.is_defined() || rate.nComp() != m_integral.nComp()) {
        diagnostic = "completed-step ledger and face-rate layouts are not defined compatibly";
        latch_stage_failure(diagnostic);
        local_stage_ok = false;
    } else if (!SameMappedFaceLayout(m_integral, rate)) {
        diagnostic = "completed-step ledger and face-rate BoxArray/DistributionMapping do not match";
        latch_stage_failure(diagnostic);
        local_stage_ok = false;
    }
    if (local_stage_ok) {
        local_stage_ok = begin_stage(method, stage, step_old_time, recipe, diagnostic);
    }
    if (local_stage_ok) {
        for (int component = 0; component < rate.nComp(); ++component) {
            if (!accumulate_stage_component(rate, component, component, diagnostic)) {
                local_stage_ok = false;
                break;
            }
        }
    }
    const bool accepted = finish_stage(diagnostic);
    return local_stage_ok && accepted;
}

bool CompletedStepFluxLedger::begin_stage (const HostIntegrator method,
                                           const int stage,
                                           const double step_old_time,
                                           const AuxiliaryStageRecipe& recipe,
                                           std::string& diagnostic)
{
    diagnostic.clear();
    if (!m_stage_failure_diagnostic.empty()) {
        diagnostic = m_stage_failure_diagnostic;
        return false;
    }
    if (!is_defined()) {
        diagnostic = "completed-step ledger storage is not defined";
        latch_stage_failure(diagnostic);
        return false;
    }
    if (method == HostIntegrator::AnelasticMidPoint) {
        diagnostic = "AnelasticMidPoint completed-step ledger is not qualified by M2";
        latch_stage_failure(diagnostic);
        return false;
    }
    if (stage < 0 || stage >= stage_count(method)) {
        diagnostic = "stage index is outside the completed-step ledger sequence";
        latch_stage_failure(diagnostic);
        return false;
    }
    if (!std::isfinite(step_old_time) ||
        !finite_nonnegative(recipe.anchor_weight) ||
        !finite_nonnegative(recipe.input_weight) ||
        !std::isfinite(recipe.limiter_trial_interval) || recipe.limiter_trial_interval <= 0.0 ||
        !std::isfinite(recipe.face_rate_time_coefficient) || recipe.face_rate_time_coefficient <= 0.0 ||
        !finite_nonnegative(recipe.completed_ledger_time)) {
        diagnostic = "completed-step ledger received an invalid temporal coefficient";
        latch_stage_failure(diagnostic);
        return false;
    }
    const amrex::Real stage_weight =
        static_cast<amrex::Real>(recipe.completed_ledger_time);
    if (!amrex::Math::isfinite(stage_weight) || stage_weight < amrex::Real(0.0)) {
        diagnostic = "completed-step ledger weight is not representable in amrex::Real";
        latch_stage_failure(diagnostic);
        return false;
    }
    if (m_stage_open) {
        diagnostic = "a completed-step ledger stage transaction is already open";
        latch_stage_failure(diagnostic);
        return false;
    }

    if (stage == 0) {
        if (m_step_active) {
            diagnostic = "stage 0 arrived while the previous auxiliary step was unfinished";
            latch_stage_failure(diagnostic);
            return false;
        }
        m_integral.setVal(amrex::Real(0.0));
        m_step_active = true;
        m_step_complete = false;
        m_method = method;
        m_step_old_time = step_old_time;
        m_next_stage = 0;
    } else {
        if (!m_step_active) {
            diagnostic = "auxiliary stage arrived before stage 0";
            latch_stage_failure(diagnostic);
            return false;
        }
        if (method != m_method) {
            diagnostic = "host integrator changed during an auxiliary timestep";
            latch_stage_failure(diagnostic);
            return false;
        }
        if (step_old_time != m_step_old_time) {
            diagnostic = "step-old time changed during an auxiliary timestep";
            latch_stage_failure(diagnostic);
            return false;
        }
    }
    if (stage != m_next_stage) {
        diagnostic = "duplicate, skipped, or out-of-order auxiliary stage";
        latch_stage_failure(diagnostic);
        return false;
    }

    std::fill(m_stage_components_seen.begin(), m_stage_components_seen.end(), 0);
    m_stage_open = true;
    m_open_stage = stage;
    m_open_stage_weight = stage_weight;
    return true;
}

bool CompletedStepFluxLedger::accumulate_stage_component (
    const MappedFaceFluxRate& rate, const int source_component,
    const int ledger_component, std::string& diagnostic)
{
    diagnostic.clear();
    const auto fail = [this, &diagnostic] () {
        latch_stage_failure(diagnostic);
        return false;
    };
    if (!m_stage_open || !m_step_active || m_open_stage != m_next_stage) {
        diagnostic = "spectral component accumulation requires an open host stage";
        return fail();
    }
    if (!rate.is_defined() || !SameMappedFaceLayout(m_integral, rate)) {
        diagnostic = "completed-step ledger and chunk face-rate layouts do not match";
        return fail();
    }
    if (source_component < 0 || source_component >= rate.nComp() ||
        ledger_component < 0 || ledger_component >= m_integral.nComp()) {
        diagnostic = "completed-step chunk component is outside its source or destination";
        return fail();
    }
    auto& seen = m_stage_components_seen[static_cast<std::size_t>(ledger_component)];
    if (seen != 0) {
        diagnostic = "completed-step stage contains a duplicate destination component";
        return fail();
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Saxpy(m_integral.dir(dir), m_open_stage_weight,
                               rate.dir(dir), source_component,
                               ledger_component, 1, 0);
    }
    seen = 1;
    return true;
}

bool CompletedStepFluxLedger::finish_stage (std::string& diagnostic)
{
    diagnostic.clear();
    if (!m_stage_failure_diagnostic.empty()) {
        diagnostic = m_stage_failure_diagnostic;
    }
    const bool local_stage_valid =
        m_stage_failure_diagnostic.empty() && m_stage_open && m_step_active &&
        m_open_stage == m_next_stage &&
        m_open_stage >= 0 && m_open_stage < stage_count(m_method);
    if (!local_stage_valid && diagnostic.empty()) {
        diagnostic = "completed-step stage finish requires an open host stage";
    }

    bool all_components_seen = true;
    for (std::size_t component = 0; component < m_stage_components_seen.size(); ++component) {
        if (m_stage_components_seen[component] == 0) {
            all_components_seen = false;
            if (diagnostic.empty()) {
                diagnostic =
                    "completed-step stage omitted destination component " +
                    std::to_string(component);
            }
        }
    }

    // Reserve disjoint ranges for each logical host integrator.  A zero code
    // means that this rank has no locally valid, complete stage transaction.
    int local_stage_code = 0;
    if (local_stage_valid && all_components_seen) {
        switch (m_method) {
        case HostIntegrator::CompressibleRK3:
            local_stage_code = 1 + m_open_stage;
            break;
        case HostIntegrator::AnelasticHeun:
            local_stage_code = 4 + m_open_stage;
            break;
        case HostIntegrator::AnelasticMidPoint:
            local_stage_code = 7 + m_open_stage;
            break;
        }
    }
    const bool completes_step = local_stage_code != 0 &&
                                m_next_stage + 1 == stage_count(m_method);
    bool locally_finite = true;
    if (completes_step) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            locally_finite =
                m_integral.dir(dir).is_finite(0, m_integral.nComp(), 0, true) &&
                locally_finite;
        }
    }
    const int prior_complete_and_idle =
        m_step_complete && !m_step_active && !m_stage_open ? 1 : 0;
    int vote[4] = {local_stage_code, -local_stage_code, locally_finite ? 0 : -1,
                   prior_complete_and_idle};
    amrex::ParallelDescriptor::ReduceIntMin(vote, 4);
    const bool stage_agreed = vote[0] != 0 && vote[0] == -vote[1];
    const bool globally_finite = vote[2] == 0;
    const bool all_prior_complete_and_idle = vote[3] == 1;
    if (!stage_agreed || !globally_finite) {
        if (!stage_agreed && diagnostic.empty()) {
            diagnostic = vote[0] == 0
                             ? "completed-step stage is invalid or incomplete on another MPI rank"
                             : "completed-step stage identity differs across MPI ranks";
        } else if (stage_agreed && !globally_finite && diagnostic.empty()) {
            diagnostic = locally_finite
                             ? "completed-step integrated face transfer is nonfinite on another MPI rank"
                             : "completed-step integrated face transfer is nonfinite";
        }
        if (all_prior_complete_and_idle) {
            // Keep the diagnostic for this rejected call, but allow a later
            // valid stage zero to start a fresh transaction.
            m_stage_failure_diagnostic.clear();
        } else {
            discard_step();
        }
        return false;
    }

    m_stage_open = false;
    m_open_stage = -1;
    m_open_stage_weight = amrex::Real(0.0);
    ++m_next_stage;
    if (completes_step) {
        m_step_active = false;
        m_step_complete = true;
    }
    m_stage_failure_diagnostic.clear();
    return true;
}

} // namespace erf_auxiliary
