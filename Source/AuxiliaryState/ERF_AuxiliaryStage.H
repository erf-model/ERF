#ifndef ERF_AUXILIARY_STAGE_H_
#define ERF_AUXILIARY_STAGE_H_

#include "ERF_AuxiliaryMappedTransport.H"

#include <string>

namespace erf_auxiliary {

enum class HostIntegrator {
    CompressibleRK3,
    AnelasticHeun,
    AnelasticMidPoint
};

/** State view that bounds the low-order limiter trial for this host stage. */
enum class LimiterTrialBase {
    Anchor,
    Input
};

struct ConstTimedFieldView {
    const amrex::MultiFab* field{nullptr};
    int component{-1};
    double time{0.0};
};

struct MutableTimedFieldView {
    amrex::MultiFab* field{nullptr};
    int component{-1};
    double time{0.0};
};

struct HostCarrierFluxView {
    const amrex::MultiFab* x{nullptr};
    const amrex::MultiFab* y{nullptr};
    const amrex::MultiFab* z{nullptr};
};

struct AuxiliaryStageRecipe {
    LimiterTrialBase limiter_trial_base{LimiterTrialBase::Anchor};
    double anchor_weight{0.0};
    double input_weight{0.0};
    double limiter_trial_interval{0.0};
    double face_rate_time_coefficient{0.0};
    double completed_ledger_time{0.0};
};

enum class AuxiliaryFieldValidationPolicy {
    /** Perform global positive-density and finite-result validation. */
    Global,
    /**
     * Skip global checks. The caller guarantees every input state value is
     * finite, density is finite and strictly positive, and each state/density
     * quotient is finite and representable in amrex::Real. No global
     * reduction is performed by this policy.
     */
    AssumeValid
};

/** Semantic host and state views at the pre-overwrite slow-RHS post seam. */
struct AuxiliaryStageContext {
    HostIntegrator method{HostIntegrator::CompressibleRK3};
    int level{-1};
    int stage{-1};

    double step_old_time{0.0};
    double input_time{0.0};
    double target_time{0.0};

    AuxiliaryStageRecipe recurrence;

    ConstTimedFieldView state_anchor;
    ConstTimedFieldView state_input;
    MutableTimedFieldView state_target;

    ConstTimedFieldView rho_anchor;
    ConstTimedFieldView rho_input;
    ConstTimedFieldView rho_target;

    ConstTimedFieldView measure_anchor;
    ConstTimedFieldView measure_input;
    ConstTimedFieldView measure_target;

    HostCarrierFluxView carrier;
};

/** Return the explicit ERF recurrence recipe for one supported host stage. */
bool MakeAuxiliaryStageRecipe (HostIntegrator method,
                               int stage,
                               double host_stage_interval,
                               AuxiliaryStageRecipe& recipe,
                               std::string& diagnostic);

/** Form z=U_input/rho_input from the explicit input-time views. */
bool BuildAuxiliaryIntensiveState (const ConstTimedFieldView& state_input,
                                   const ConstTimedFieldView& rho_input,
                                   double input_time,
                                   amrex::MultiFab& intensive,
                                   AuxiliaryFieldValidationPolicy validation,
                                   std::string& diagnostic);

const char* HostIntegratorName (HostIntegrator method) noexcept;

/** Apply H_target = a H_anchor + b H_input - c D_xi(F_rate). */
void ApplyAuxiliaryMappedStage (const AuxiliaryStageContext& context,
                                const MappedFaceFluxRate& rate,
                                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
                                int rate_component = 0);

/** Whether the mutable stage target is a distinct MultiFab from both inputs. */
bool AuxiliaryStageTargetIsDisjoint (const AuxiliaryStageContext& context);

/** Completed-step time integral of accepted mapped face rates and stage order. */
class CompletedStepFluxLedger {
public:
    CompletedStepFluxLedger () = default;
    CompletedStepFluxLedger (const CompletedStepFluxLedger&) = delete;
    CompletedStepFluxLedger& operator= (const CompletedStepFluxLedger&) = delete;

    void define (const amrex::BoxArray& cell_ba,
                 const amrex::DistributionMapping& dm,
                 int ncomp = 1);

    bool accept_stage (HostIntegrator method,
                       int stage,
                       double step_old_time,
                       const AuxiliaryStageRecipe& recipe,
                       const MappedFaceFluxRate& rate,
                       std::string& diagnostic);

    [[nodiscard]] bool is_defined () const noexcept { return m_integral.is_defined(); }
    [[nodiscard]] bool step_active () const noexcept { return m_step_active; }
    [[nodiscard]] bool step_complete () const noexcept { return m_step_complete; }
    [[nodiscard]] int next_stage () const noexcept { return m_next_stage; }
    [[nodiscard]] const IntegratedMappedFaceFlux& integrated_flux () const noexcept {
        return m_integral;
    }

private:
    IntegratedMappedFaceFlux m_integral;
    bool m_step_active{false};
    bool m_step_complete{false};
    HostIntegrator m_method{HostIntegrator::CompressibleRK3};
    int m_next_stage{0};
    double m_step_old_time{0.0};
};

} // namespace erf_auxiliary

#endif
