#ifndef ERF_SBM_CONSTRAINT_GROUPS_H_
#define ERF_SBM_CONSTRAINT_GROUPS_H_

#include "ERF_SBMLayout.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace erf_sbm {

/** A fixed coefficient in a device-ready linear admissibility form. */
struct ConstraintTerm {
    int component{0};
    amrex::Real coefficient{0.0};
};

/** Host-side semantic name plus a linear form evaluated as form >= 0. */
struct LinearConstraint {
    std::string semantic_id;
    std::vector<ConstraintTerm> terms;
};

/**
 * Complete atomic constraint group over one population bin.
 */
struct ConstraintGroup {
    int population_id{0};
    int bin{0};
    std::string semantic_id;
    MomentMode moment_mode{MomentMode::OneMoment};
    std::vector<int> members;
    std::vector<int> transport_members;
    std::vector<int> attached_property_indices;
    std::vector<LinearConstraint> constraints;

    [[nodiscard]] bool contains (const int component) const noexcept;
    [[nodiscard]] amrex::Real evaluate (const LinearConstraint& constraint,
                                        const std::vector<amrex::Real>& state) const;
    [[nodiscard]] bool admissible (const std::vector<amrex::Real>& state,
                                   amrex::Real* minimum_margin = nullptr,
                                   std::string* failed_constraint = nullptr) const;
};

/**
 * Device-ready flattened form of one semantic constraint. Terms are stored
 * in one sparse array and addressed by offset/count, with no fixed term limit.
 */
struct ConstraintDescriptor {
    int group_index{0};
    int population_id{0};
    int bin{0};
    int constraint_index{0};
    int term_offset{0};
    int term_count{0};
};

/** Host-owned compact storage copied as contiguous views to production kernels. */
struct FlattenedConstraintSet {
    std::vector<ConstraintDescriptor> constraints;
    std::vector<ConstraintTerm> terms;
};

/** Evaluate a validated flattened constraint in host or device code. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_constraint_descriptor (const ConstraintDescriptor& constraint,
                                            const ConstraintTerm* terms,
                                            const amrex::Real* state) noexcept
{
    amrex::Real result = amrex::Real(0.0);
    for (int i = 0; i < constraint.term_count; ++i) {
        const auto& term = terms[constraint.term_offset + i];
        result = std::fma(term.coefficient, state[term.component], result);
    }
    return result;
}

/** Device-ready metadata for a donor-derived attached-property envelope. */
struct AttachedPropertySupportDescriptor {
    int group_index{0};
    int property_index{0};
    int property_component{-1};
    int mass_component{-1};
    int number_component{-1};
    int carrier_component{-1};
    int two_moment{0};
    amrex::Real pivot{amrex::Real(1.0)};
    amrex::Real hard_min{amrex::Real(0.0)};
    amrex::Real hard_max{amrex::Real(0.0)};
    int has_hard_min{0};
    int has_hard_max{0};
};

/** Build deterministic groups from the immutable runtime layout. */
[[nodiscard]] std::vector<ConstraintGroup>
make_constraint_groups (const SBMLayout& layout);

/** Flatten complete semantic groups into bounded device descriptors. */
[[nodiscard]] FlattenedConstraintSet
make_constraint_descriptors (const SBMLayout& layout);

/** Flatten arbitrary semantic groups through the production sparse format. */
[[nodiscard]] FlattenedConstraintSet
flatten_constraint_groups (const std::vector<ConstraintGroup>& groups,
                           int state_components);

[[nodiscard]] std::vector<AttachedPropertySupportDescriptor>
make_attached_property_support_descriptors (const SBMLayout& layout);

/**
 * A production scratch chunk is a sequence of complete constraint groups.
 * `max_groups` is a policy limit on atomic groups, never a component limit;
 * all mass/number/property members of each group stay together.
 */
struct ConstraintClosureChunk {
    std::vector<int> group_indices;
    std::vector<int> components;
};

[[nodiscard]] std::vector<ConstraintClosureChunk>
make_constraint_closure_chunks (const SBMLayout& layout, int max_groups);

/**
 * Compatibility entry point delegating to the canonical grid-layer transform.
 */
[[nodiscard]] EndpointTransform
transform_two_moment (amrex::Real C, amrex::Real M,
                      amrex::Real lower, amrex::Real upper);

/** Inverse transform with finite-state and nonnegative-endpoint checks. */
[[nodiscard]] std::pair<amrex::Real, amrex::Real>
inverse_two_moment (amrex::Real L, amrex::Real H,
                    amrex::Real lower, amrex::Real upper);

/** Host reference for the finite property-support constraints. */
[[nodiscard]] bool
property_support_is_admissible (amrex::Real property, amrex::Real carrier_number,
                                amrex::Real lower, amrex::Real upper,
                                amrex::Real tolerance = 0.0) noexcept;

} // namespace erf_sbm

#endif
