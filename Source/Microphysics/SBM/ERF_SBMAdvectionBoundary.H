#ifndef ERF_SBM_ADVECTION_BOUNDARY_H_
#define ERF_SBM_ADVECTION_BOUNDARY_H_

#include <AMReX_REAL.H>

namespace erf_sbm {

/** Physical face categories understood by the bounded M3 advection contract. */
enum class AdvectionFaceKind : int {
    Periodic,
    ImpermeableWall,
    AdvectiveOutflow,
    PrescribedInflow
};

/** Result of applying the M3 spectral state policy at one oriented face. */
enum class AdvectionFaceAction : int {
    PeriodicSharedFace,
    ZeroNormalFlux,
    InteriorDonor,
    ExplicitSpectralInflow,
    RejectInwardOutflow,
    RejectMissingInflowState,
    RejectIncompleteInflowState
};

struct AdvectionFacePolicyInput {
    AdvectionFaceKind kind{AdvectionFaceKind::Periodic};
    // Positive is outward from the physical domain at this particular face.
    amrex::Real outward_carrier{amrex::Real(0.0)};
    bool explicit_spectral_state_present{false};
    // True only when every member of the atomic group, including all attached
    // properties, is present. Compact qc/qr fields are not spectral inflow.
    bool explicit_spectral_group_complete{false};
    bool high_order_stencil_available{false};
};

struct AdvectionFacePolicy {
    AdvectionFaceAction action{AdvectionFaceAction::ZeroNormalFlux};
    bool use_interior_donor{false};
    bool use_explicit_spectral_state{false};
    bool use_high_order_candidate{false};

    [[nodiscard]] constexpr bool accepted () const noexcept
    {
        return action != AdvectionFaceAction::RejectInwardOutflow &&
               action != AdvectionFaceAction::RejectMissingInflowState &&
               action != AdvectionFaceAction::RejectIncompleteInflowState;
    }
};

/**
 * Decide the face-level spectral advection semantics without performing ERF's
 * deferred physical-boundary ghost-fill lifecycle. At a physical face whose
 * high-order stencil is not qualified, the admissible donor proposal remains
 * available and high-order transport falls back to it.
 */
[[nodiscard]] constexpr AdvectionFacePolicy
make_advection_face_policy (const AdvectionFacePolicyInput& input) noexcept
{
    using Action = AdvectionFaceAction;
    if (input.kind == AdvectionFaceKind::Periodic) {
        return {Action::PeriodicSharedFace, false, false, true};
    }
    if (input.kind == AdvectionFaceKind::ImpermeableWall ||
        input.outward_carrier == amrex::Real(0.0)) {
        return {Action::ZeroNormalFlux, false, false, false};
    }

    if (input.outward_carrier > amrex::Real(0.0)) {
        // Outward flow uses the interior spectrum. The donor path is always
        // defined; high order is used only with a qualified physical stencil.
        return {Action::InteriorDonor, true, false,
                input.high_order_stencil_available};
    }

    // Reverse flow contradicts a face declared as advective outflow. External
    // spectral data must not silently change that declared boundary policy.
    if (input.kind == AdvectionFaceKind::AdvectiveOutflow) {
        return {Action::RejectInwardOutflow, false, false, false};
    }

    // The remaining negative-carrier boundary category is prescribed inflow.
    // A complete explicit atomic-group spectrum is required; qc/qr cannot
    // supply missing bin or attached-property composition.
    if (!input.explicit_spectral_state_present) {
        return {Action::RejectMissingInflowState, false, false, false};
    }
    if (!input.explicit_spectral_group_complete) {
        return {Action::RejectIncompleteInflowState, false, false, false};
    }
    return {Action::ExplicitSpectralInflow, false, true,
            input.high_order_stencil_available};
}

} // namespace erf_sbm

#endif // ERF_SBM_ADVECTION_BOUNDARY_H_
