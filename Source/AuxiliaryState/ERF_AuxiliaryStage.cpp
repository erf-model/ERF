#include "ERF_AuxiliaryStage.H"

#include <cmath>

namespace erf_auxiliary {
namespace {

int stage_count(const HostIntegrator method) noexcept
{
    switch (method) {
    case HostIntegrator::CompressibleRK3: return 3;
    case HostIntegrator::AnelasticHeun: return 2;
    case HostIntegrator::AnelasticMidPoint: return 2;
    }
    return 0;
}

bool finite_nonnegative(const double value) noexcept
{
    return std::isfinite(value) && value >= 0.0;
}

bool valid_const_view(const ConstTimedFieldView& view,
                      const double expected_time) noexcept
{
    return view.field != nullptr && view.component >= 0 &&
           view.component < view.field->nComp() &&
           std::isfinite(view.time) && view.time == expected_time;
}

bool same_layout(const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    return a.boxArray() == b.boxArray() && a.DistributionMap() == b.DistributionMap();
}

} // namespace

const char* HostIntegratorName(const HostIntegrator method) noexcept
{
    switch (method) {
    case HostIntegrator::CompressibleRK3: return "CompressibleRK3";
    case HostIntegrator::AnelasticHeun: return "AnelasticHeun";
    case HostIntegrator::AnelasticMidPoint: return "AnelasticMidPoint";
    }
    return "Unknown";
}

bool MakeAuxiliaryStageRecipe(const HostIntegrator method,
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
        recipe.anchor_weight = 1.0;
        if (stage == 0) {
            recipe.limiter_trial_interval = host_stage_interval;
            recipe.face_rate_time_coefficient = host_stage_interval;
        } else if (stage == 1) {
            recipe.limiter_trial_interval = host_stage_interval;
            recipe.face_rate_time_coefficient = host_stage_interval;
        } else {
            recipe.limiter_trial_interval = host_stage_interval;
            recipe.face_rate_time_coefficient = host_stage_interval;
            recipe.completed_ledger_time = host_stage_interval;
        }
    } else if (method == HostIntegrator::AnelasticHeun) {
        if (stage == 0) {
            recipe.anchor_weight = 1.0;
            recipe.limiter_trial_interval = host_stage_interval;
            recipe.face_rate_time_coefficient = host_stage_interval;
            recipe.completed_ledger_time = 0.5 * host_stage_interval;
        } else {
            recipe.anchor_weight = 0.5;
            recipe.input_weight = 0.5;
            recipe.limiter_trial_interval = host_stage_interval;
            recipe.face_rate_time_coefficient = 0.5 * host_stage_interval;
            recipe.completed_ledger_time = 0.5 * host_stage_interval;
        }
    }
    return true;
}

bool BuildAuxiliaryIntensiveState(const ConstTimedFieldView& state_input,
                                  const ConstTimedFieldView& rho_input,
                                  const double input_time,
                                  amrex::MultiFab& intensive,
                                  std::string& diagnostic)
{
    diagnostic.clear();
    if (!std::isfinite(input_time) ||
        !valid_const_view(state_input, input_time) ||
        !valid_const_view(rho_input, input_time)) {
        diagnostic = "auxiliary intensive state requires state and density views at input_time";
        return false;
    }
    if (!same_layout(intensive, *state_input.field) ||
        !same_layout(intensive, *rho_input.field) ||
        intensive.nComp() < 1) {
        diagnostic = "auxiliary intensive state fields must share a cell layout";
        return false;
    }
    if (!ValidatePositiveFiniteComponent(*rho_input.field, rho_input.component, diagnostic)) {
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
    if (!ValidateFiniteComponent(intensive, 0, diagnostic)) {
        return false;
    }
    return true;
}

void ApplyAuxiliaryMappedStage(const AuxiliaryStageContext& context,
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

    const auto& target = *context.state_target.field;
    const auto& anchor = *context.state_anchor.field;
    const auto& input = *context.state_input.field;
    const auto& omega_anchor = *context.measure_anchor.field;
    const auto& omega_input = *context.measure_input.field;
    const auto& omega_target = *context.measure_target.field;
    AMREX_ALWAYS_ASSERT(same_layout(target, anchor));
    AMREX_ALWAYS_ASSERT(same_layout(target, input));
    AMREX_ALWAYS_ASSERT(same_layout(target, omega_anchor));
    AMREX_ALWAYS_ASSERT(same_layout(target, omega_input));
    AMREX_ALWAYS_ASSERT(same_layout(target, omega_target));

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

void CompletedStepFluxLedger::define(const amrex::BoxArray& cell_ba,
                                    const amrex::DistributionMapping& dm,
                                    const int ncomp)
{
    m_integral.define(cell_ba, dm, ncomp, 0);
    m_integral.setVal(amrex::Real(0.0));
    m_step_active = false;
    m_step_complete = false;
    m_next_stage = 0;
    m_step_old_time = 0.0;
}

bool CompletedStepFluxLedger::accept_stage(const HostIntegrator method,
                                           const int stage,
                                           const double step_old_time,
                                           const AuxiliaryStageRecipe& recipe,
                                           const MappedFaceFluxRate& rate,
                                           std::string& diagnostic)
{
    diagnostic.clear();
    if (!is_defined() || !rate.is_defined() || rate.nComp() != m_integral.nComp()) {
        diagnostic = "completed-step ledger and face-rate layouts are not defined compatibly";
        return false;
    }
    if (method == HostIntegrator::AnelasticMidPoint) {
        diagnostic = "AnelasticMidPoint completed-step ledger is not qualified by M2";
        return false;
    }
    if (stage < 0 || stage >= stage_count(method)) {
        diagnostic = "stage index is outside the completed-step ledger sequence";
        return false;
    }
    if (!std::isfinite(step_old_time) ||
        !finite_nonnegative(recipe.anchor_weight) ||
        !finite_nonnegative(recipe.input_weight) ||
        !std::isfinite(recipe.limiter_trial_interval) || recipe.limiter_trial_interval <= 0.0 ||
        !std::isfinite(recipe.face_rate_time_coefficient) || recipe.face_rate_time_coefficient <= 0.0 ||
        !finite_nonnegative(recipe.completed_ledger_time)) {
        diagnostic = "completed-step ledger received an invalid temporal coefficient";
        return false;
    }

    if (stage == 0) {
        if (m_step_active) {
            diagnostic = "stage 0 arrived while the previous auxiliary step was unfinished";
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
            return false;
        }
        if (method != m_method) {
            diagnostic = "host integrator changed during an auxiliary timestep";
            return false;
        }
        if (step_old_time != m_step_old_time) {
            diagnostic = "step-old time changed during an auxiliary timestep";
            return false;
        }
    }

    if (stage != m_next_stage) {
        diagnostic = "duplicate, skipped, or out-of-order auxiliary stage";
        return false;
    }

    AccumulateIntegratedFaceFlux(m_integral, rate,
        static_cast<amrex::Real>(recipe.completed_ledger_time));
    ++m_next_stage;
    if (m_next_stage == stage_count(method)) {
        m_step_active = false;
        m_step_complete = true;
    }
    return true;
}

} // namespace erf_auxiliary
