#ifndef ERF_AUXILIARY_STAGE_CONTEXT_H_
#define ERF_AUXILIARY_STAGE_CONTEXT_H_

#include <AMReX_MultiFab.H>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace erf_auxiliary {

enum class IntegrationMethod {
    CompressibleRK3,
    AnelasticHeun
};

/** Explicit temporal contract passed to an auxiliary transport kernel. */
struct StageContext {
    IntegrationMethod method{IntegrationMethod::CompressibleRK3};
    int stage_index{0};
    double old_time{0.0};
    double evaluation_time{0.0};
    double output_time{0.0};
    double stage_interval{0.0};
    bool completes_level_step{false};

    // These are intentionally explicit rather than inferred from stage_index.
    const amrex::MultiFab* baseline_state{nullptr};
    const amrex::MultiFab* evaluation_state{nullptr};

    // Stage update coefficient relative to the explicit RHS.  For the ERF
    // callback intervals this is h/3, h/2, h for RK3 and h for both Heun calls.
    [[nodiscard]] double rhs_interval() const noexcept { return stage_interval; }

    // Weight in the accepted full-step integrated face-transfer ledger.
    [[nodiscard]] double accepted_ledger_weight() const noexcept {
        if (method == IntegrationMethod::AnelasticHeun) return 0.5;
        return completes_level_step ? 1.0 : 0.0;
    }
};

inline StageContext make_compressible_stage(const int stage,
                                            const double old_time,
                                            const double evaluation_time,
                                            const double output_time,
                                            const double full_step,
                                            const amrex::MultiFab* baseline,
                                            const amrex::MultiFab* evaluation)
{
    if (stage < 0 || stage > 2 || full_step < 0.0) {
        throw std::invalid_argument("invalid compressible auxiliary stage");
    }
    const double expected[] = {full_step / 3.0, full_step / 2.0, full_step};
    const double interval = output_time - old_time;
    const double scale = 1.0 + std::abs(expected[stage]);
    if (std::abs(interval - expected[stage]) > 64.0 * std::numeric_limits<double>::epsilon() * scale) {
        throw std::invalid_argument("compressible auxiliary stage interval does not match ERF RK3");
    }
    return {IntegrationMethod::CompressibleRK3, stage, old_time, evaluation_time,
            output_time, interval, stage == 2, baseline, evaluation};
}

inline StageContext make_anelastic_stage(const int stage,
                                         const double old_time,
                                         const double evaluation_time,
                                         const double output_time,
                                         const double full_step,
                                         const amrex::MultiFab* baseline,
                                         const amrex::MultiFab* evaluation)
{
    if (stage < 0 || stage > 1 || full_step < 0.0) {
        throw std::invalid_argument("invalid anelastic auxiliary stage");
    }
    const double interval = output_time - old_time;
    const double scale = 1.0 + std::abs(full_step);
    if (std::abs(interval - full_step) > 64.0 * std::numeric_limits<double>::epsilon() * scale) {
        throw std::invalid_argument("anelastic auxiliary stage interval does not match ERF Heun");
    }
    return {IntegrationMethod::AnelasticHeun, stage, old_time, evaluation_time,
            output_time, interval, stage == 1, baseline, evaluation};
}

} // namespace erf_auxiliary

#endif
