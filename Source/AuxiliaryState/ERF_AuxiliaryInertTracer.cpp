#include "ERF_AuxiliaryInertTracer.H"

#include "ERF_AuxiliaryMappedTransport.H"
#include "ERF_AuxiliaryStage.H"
#include "ERF_IndexDefines.H"
#include "Advection/ERF_AdvectionSrcForScalars.H"

#include <AMReX_Math.H>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

namespace erf_auxiliary {

struct AuxiliaryInertTracer::LevelStorage {
    amrex::MultiFab state;
    amrex::MultiFab target;
    amrex::MultiFab anchor;
    amrex::MultiFab intensive;
    amrex::MultiFab measure;
    MappedFaceFluxRate rate;
    MappedFaceFluxRate previous_rate;
    CompletedStepFluxLedger ledger;
    bool have_previous_rate{false};
    bool measure_ready{false};
    amrex::Real max_stage_rate_delta{0.0};
};

AuxiliaryInertTracer::~AuxiliaryInertTracer () = default;

AuxiliaryInertTracer::AuxiliaryInertTracer (const int number_of_levels)
    : m_levels(static_cast<std::size_t>(number_of_levels))
{
    AMREX_ALWAYS_ASSERT(number_of_levels > 0);
}

const AuxiliaryStateLayout& AuxiliaryInertTracer::layout ()
{
    static const AuxiliaryStateLayout value(
        "erf-auxiliary-inert-tracer-m2-v1",
        {{"inert_tracer", "generic.inert.tracer", "kg m^-3"}});
    return value;
}

void AuxiliaryInertTracer::define (const int level,
                                 const amrex::BoxArray& cell_ba,
                                 const amrex::DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(level >= 0 && level < static_cast<int>(m_levels.size()));
    AMREX_ALWAYS_ASSERT(level == 0);
    AMREX_ALWAYS_ASSERT(m_levels[static_cast<std::size_t>(level)] == nullptr);
    auto data = std::make_unique<LevelStorage>();
    data->state.define(cell_ba, dm, 1, 0);
    data->target.define(cell_ba, dm, 1, 0);
    data->anchor.define(cell_ba, dm, 1, 0);
    data->intensive.define(cell_ba, dm, 1, 1);
    data->measure.define(cell_ba, dm, 1, 0);
    data->rate.define(cell_ba, dm, 1, 0);
    data->previous_rate.define(cell_ba, dm, 1, 0);
    data->ledger.define(cell_ba, dm, 1);
    data->state.setVal(amrex::Real(0.0));
    data->target.setVal(amrex::Real(0.0));
    data->anchor.setVal(amrex::Real(0.0));
    data->intensive.setVal(amrex::Real(0.0));
    data->rate.setVal(amrex::Real(0.0));
    data->previous_rate.setVal(amrex::Real(0.0));

    m_levels[static_cast<std::size_t>(level)] = std::move(data);
}

bool AuxiliaryInertTracer::rebuild_static_measure (const int level,
                                                   const amrex::MultiFab& detJ,
                                                   const amrex::MultiFab& mx,
                                                   const amrex::MultiFab& my,
                                                   std::string& diagnostic)
{
    diagnostic.clear();
    if (level < 0 || level >= static_cast<int>(m_levels.size()) ||
        m_levels[static_cast<std::size_t>(level)] == nullptr) {
        diagnostic = "static mapped measure requires a defined auxiliary level";
        return false;
    }
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    data.measure_ready = false;
    if (!BuildMappedCellMeasure(data.measure, detJ, mx, my, diagnostic)) {
        return false;
    }
    data.measure_ready = true;
    return true;
}

void AuxiliaryInertTracer::initialize (const int level,
                                     const amrex::MultiFab& conserved,
                                     const amrex::Geometry& geometry)
{
    AMREX_ALWAYS_ASSERT(is_defined(level));
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(data.measure_ready,
        "M2 auxiliary inert tracer cannot initialize before its static measure is built");
    std::string diagnostic;
    if (!ValidatePositiveFiniteComponent(conserved, Rho_comp, diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer initialization: " + diagnostic);
    }

    const auto prob_lo = geometry.ProbLoArray();
    const auto prob_hi = geometry.ProbHiArray();
    const auto cell_size = geometry.CellSizeArray();
    for (amrex::MFIter mfi(data.state, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto rho = conserved.const_array(mfi);
        const auto out = data.state.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real x = prob_lo[0] + (i + amrex::Real(0.5)) * cell_size[0];
            const amrex::Real y = prob_lo[1] + (j + amrex::Real(0.5)) * cell_size[1];
            const amrex::Real z = prob_lo[2] + (k + amrex::Real(0.5)) * cell_size[2];
            const amrex::Real lx = prob_hi[0] - prob_lo[0];
            const amrex::Real ly = prob_hi[1] - prob_lo[1];
            const amrex::Real lz = prob_hi[2] - prob_lo[2];
            const amrex::Real ratio = amrex::Real(1.0) + amrex::Real(0.10) *
                std::sin(amrex::Real(2.0) * amrex::Real(PI) * (x - prob_lo[0]) / lx) +
                amrex::Real(0.07) * std::sin(amrex::Real(2.0) * amrex::Real(PI) * (y - prob_lo[1]) / ly) +
                amrex::Real(0.05) * std::sin(amrex::Real(2.0) * amrex::Real(PI) * (z - prob_lo[2]) / lz);
            out(i, j, k, 0) = rho(i, j, k, Rho_comp) * ratio;
        });
    }
    amrex::MultiFab::Copy(data.anchor, data.state, 0, 0, 1, 0);
    if (!ValidateFiniteComponent(data.state, 0, diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer initialization: " + diagnostic);
    }
    amrex::Print() << "AUX_M2_INITIALIZED schema=" << layout().schema_id()
                   << " components=" << layout().ncomp() << std::endl;
}

void AuxiliaryInertTracer::advance_stage (
    const int level,
    const HostIntegrator method,
    const int stage,
    const double step_old_time,
    const double input_time,
    const double target_time,
    const double host_stage_interval,
    const amrex::MultiFab& conserved_anchor,
    const amrex::MultiFab& conserved_input,
    const amrex::MultiFab& conserved_target,
    const amrex::MultiFab& avg_xmom,
    const amrex::MultiFab& avg_ymom,
    const amrex::MultiFab& avg_zmom,
    const AdvChoice& advection,
    const amrex::Geometry& geometry)
{
    AMREX_ALWAYS_ASSERT(is_defined(level));
    auto& data = *m_levels[static_cast<std::size_t>(level)];
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(data.measure_ready,
        "M2 auxiliary inert tracer cannot advance before its static measure is built");
    AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    if (!MakeAuxiliaryStageRecipe(method, stage, host_stage_interval, recipe, diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer stage recipe: " + diagnostic);
    }
    if (!std::isfinite(step_old_time) || !std::isfinite(input_time) ||
        !std::isfinite(target_time)) {
        amrex::Abort("M2 auxiliary inert tracer received a nonfinite semantic time");
    }
    if (!data.ledger.is_defined()) {
        amrex::Abort("M2 auxiliary inert tracer ledger is not defined");
    }
    if (stage == 0 && data.ledger.step_active()) {
        amrex::Abort("M2 auxiliary inert tracer received stage 0 before the prior step completed");
    }

    if (!ValidatePositiveFiniteComponent(conserved_anchor, Rho_comp, diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer anchor density: " + diagnostic);
    }
    if (!ValidatePositiveFiniteComponent(conserved_target, Rho_comp, diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer target density: " + diagnostic);
    }

    if (stage == 0) {
        amrex::MultiFab::Copy(data.anchor, data.state, 0, 0, 1, 0);
        data.have_previous_rate = false;
        data.max_stage_rate_delta = amrex::Real(0.0);
    }

    const ConstTimedFieldView state_input{&data.state, 0, input_time};
    const ConstTimedFieldView rho_input{&conserved_input, Rho_comp, input_time};
    if (!BuildAuxiliaryIntensiveState(state_input, rho_input, input_time,
                                      data.intensive,
                                      AuxiliaryFieldValidationPolicy::Global,
                                      diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer intensive state: " + diagnostic);
    }
    data.intensive.FillBoundary(geometry.periodicity());

    AMREX_ALWAYS_ASSERT(advection.dryscal_horiz_adv_type == AdvType::Centered_2nd);
    AMREX_ALWAYS_ASSERT(advection.dryscal_vert_adv_type == AdvType::Centered_2nd);
    data.rate.setVal(amrex::Real(0.0));
    for (amrex::MFIter mfi(data.intensive, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flux_views{
            data.rate.dir(0).array(mfi), data.rate.dir(1).array(mfi),
            data.rate.dir(2).array(mfi)};
        BuildScalarAdvectionFluxes(
            bx, data.intensive.const_array(mfi), 0, flux_views, 0,
            avg_xmom.const_array(mfi), avg_ymom.const_array(mfi), avg_zmom.const_array(mfi),
            advection.dryscal_horiz_adv_type, advection.dryscal_vert_adv_type,
            advection.dryscal_horiz_upw_frac, advection.dryscal_vert_upw_frac);
    }

    // These full-domain reductions are qualification-only diagnostics for this
    // test fixture, not part of the reusable transport operator path.
    const amrex::Real carrier_max = amrex::max(
        amrex::max(avg_xmom.norm0(0), avg_ymom.norm0(0)), avg_zmom.norm0(0));
    const amrex::Real rate_max = amrex::max(
        amrex::max(data.rate.dir(0).norm0(0), data.rate.dir(1).norm0(0)),
        data.rate.dir(2).norm0(0));
    if (!(carrier_max > amrex::Real(0.0))) {
        amrex::Abort("M2 auxiliary inert tracer fixture requires a nonzero native host carrier");
    }

    const amrex::Real stage_rate_delta = data.have_previous_rate ?
        MaxFaceFieldDifference(data.rate, 0, data.previous_rate, 0) : amrex::Real(0.0);
    data.max_stage_rate_delta = amrex::max(data.max_stage_rate_delta, stage_rate_delta);

    // accept_stage is MPI-collective: all ranks on this level must enter even
    // after a local recoverable error; do not add rank-local early returns.
    if (!data.ledger.accept_stage(method, stage, step_old_time, recipe, data.rate, diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer stage sequence: " + diagnostic);
    }

    AuxiliaryStageContext context;
    context.method = method;
    context.level = level;
    context.stage = stage;
    context.step_old_time = step_old_time;
    context.input_time = input_time;
    context.target_time = target_time;
    context.recurrence = recipe;
    context.state_anchor = {&data.anchor, 0, step_old_time};
    context.state_input = {&data.state, 0, input_time};
    context.state_target = {&data.target, 0, target_time};
    context.rho_anchor = {&conserved_anchor, Rho_comp, step_old_time};
    context.rho_input = {&conserved_input, Rho_comp, input_time};
    context.rho_target = {&conserved_target, Rho_comp, target_time};
    context.measure_anchor = {&data.measure, 0, step_old_time};
    context.measure_input = {&data.measure, 0, input_time};
    context.measure_target = {&data.measure, 0, target_time};
    // This proof consumer is qualified only for a time-invariant mapped
    // measure. Its three semantic measure views intentionally reference the
    // same field. A moving-mesh consumer must supply distinct anchor/input/
    // target measures and is outside this fixture's supported envelope.
    context.carrier = {&avg_xmom, &avg_ymom, &avg_zmom};
    ApplyAuxiliaryMappedStage(context, data.rate, geometry.InvCellSizeArray(), 0);

    if (!ValidateFiniteComponent(data.target, 0, diagnostic)) {
        amrex::Abort("M2 auxiliary inert tracer produced a nonfinite state: " + diagnostic);
    }
    // All reads from the old active state have completed. Promote the result
    // only after the stage kernel and qualification check have succeeded.
    std::swap(data.state, data.target);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Copy(data.previous_rate.dir(dir), data.rate.dir(dir), 0, 0, 1, 0);
    }
    data.have_previous_rate = true;

    amrex::Print() << "AUX_M2_STAGE method=" << HostIntegratorName(method)
                   << " stage=" << stage
                   << " step_time=" << step_old_time
                   << " input_time=" << input_time
                   << " target_time=" << target_time
                   << " carrier_max=" << carrier_max
                   << " rate_max=" << rate_max
                   << " rate_stage_delta=" << stage_rate_delta
                   << " finite_state=true"
                   << " ledger_weight=" << recipe.completed_ledger_time << std::endl;

    if (data.ledger.step_complete()) {
        if (!(data.max_stage_rate_delta > std::numeric_limits<amrex::Real>::epsilon() *
              amrex::max(amrex::Real(1.0), rate_max))) {
            amrex::Abort("M2 auxiliary inert tracer did not produce stage-varying mapped face rates");
        }

        amrex::MultiFab residual(data.state.boxArray(), data.state.DistributionMap(), 1, 0);
        amrex::MultiFab scaled_residual(data.state.boxArray(), data.state.DistributionMap(), 1, 0);
        const auto dx_inv = geometry.InvCellSizeArray();
        const amrex::Real min_scale = std::numeric_limits<amrex::Real>::min();
        for (amrex::MFIter mfi(data.state, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto final_u = data.state.const_array(mfi);
            const auto old_u = data.anchor.const_array(mfi);
            const auto omega = data.measure.const_array(mfi);
            const auto fx = data.ledger.integrated_flux().dir(0).const_array(mfi);
            const auto fy = data.ledger.integrated_flux().dir(1).const_array(mfi);
            const auto fz = data.ledger.integrated_flux().dir(2).const_array(mfi);
            const auto out = residual.array(mfi);
            const auto scaled = scaled_residual.array(mfi);
            const auto dx = dx_inv[0];
            const auto dy = dx_inv[1];
            const auto dz = dx_inv[2];
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real h_final = omega(i, j, k, 0) * final_u(i, j, k, 0);
                const amrex::Real h_old = omega(i, j, k, 0) * old_u(i, j, k, 0);
                const amrex::Real div = ComputationalMappedDivergence(
                    fx(i + 1, j, k, 0), fx(i, j, k, 0),
                    fy(i, j + 1, k, 0), fy(i, j, k, 0),
                    fz(i, j, k + 1, 0), fz(i, j, k, 0), dx, dy, dz);
                const amrex::Real r = h_final - h_old + div;
                const amrex::Real scale = amrex::max(
                    amrex::max(amrex::Math::abs(h_final), amrex::Math::abs(h_old)),
                    amrex::max(amrex::Math::abs(div), min_scale));
                out(i, j, k, 0) = r;
                scaled(i, j, k, 0) = amrex::Math::abs(r) / scale;
            });
        }
        const amrex::Real max_absolute = residual.norm0(0);
        const amrex::Real max_scaled = scaled_residual.norm0(0);
        const int nstages = method == HostIntegrator::CompressibleRK3 ? 3 : 2;
        const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon() *
            amrex::Real(64.0 * (nstages + 1));
        amrex::Print() << "AUX_M2_LEDGER max_abs_residual=" << max_absolute
                       << " max_scaled_residual=" << max_scaled
                       << " roundoff_tolerance=" << tolerance
                       << " max_stage_rate_delta=" << data.max_stage_rate_delta << std::endl;
        if (!(max_scaled <= tolerance)) {
            std::ostringstream message;
            message << "completed-step mapped flux ledger failed local closure: abs="
                    << max_absolute << " scaled=" << max_scaled
                    << " tolerance=" << tolerance;
            amrex::Abort(message.str());
        }
    }
}

void AuxiliaryInertTracer::destroy (const int level)
{
    AMREX_ALWAYS_ASSERT(level >= 0 && level < static_cast<int>(m_levels.size()));
    m_levels[static_cast<std::size_t>(level)].reset();
}

bool AuxiliaryInertTracer::is_defined (const int level) const
{
    return level >= 0 && level < static_cast<int>(m_levels.size()) &&
           m_levels[static_cast<std::size_t>(level)] != nullptr;
}

bool AuxiliaryInertTracer::measure_is_ready (const int level) const
{
    return is_defined(level) && m_levels[static_cast<std::size_t>(level)]->measure_ready;
}

const amrex::MultiFab& AuxiliaryInertTracer::static_measure (const int level) const
{
    AMREX_ALWAYS_ASSERT(measure_is_ready(level));
    return m_levels[static_cast<std::size_t>(level)]->measure;
}

const amrex::MultiFab& AuxiliaryInertTracer::state (const int level) const
{
    AMREX_ALWAYS_ASSERT(is_defined(level));
    return m_levels[static_cast<std::size_t>(level)]->state;
}

} // namespace erf_auxiliary
