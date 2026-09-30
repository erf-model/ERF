#ifndef ERF_SBM_CONSTRAINT_GROUPS_H_
#define ERF_SBM_CONSTRAINT_GROUPS_H_

#include "ERF_SBMLayout.H"

#include <AMReX_REAL.H>

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
 * Complete atomic transport group.  The component offsets are physical layout
 * offsets reused as bounded scratch offsets: for a two-moment population the
 * mass slot contains L and the number slot contains H only while transport
 * scratch is active.  Checkpoint/storage remains (M,C).
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

    [[nodiscard]] bool contains(const int component) const noexcept;
    [[nodiscard]] amrex::Real evaluate(const LinearConstraint& constraint,
                                        const std::vector<amrex::Real>& state) const;
    [[nodiscard]] bool admissible(const std::vector<amrex::Real>& state,
                                  amrex::Real* minimum_margin = nullptr,
                                  std::string* failed_constraint = nullptr) const;
};

/**
 * Device-ready flattened form of one semantic constraint.  P2 constraint
 * groups are deliberately linear and currently have at most two terms; this
 * representation lets the production limiter and validation path consume the
 * exact same group definition without capturing host strings or vectors in a
 * GPU kernel.
 */
struct ConstraintDescriptor {
    int group_index{0};
    int population_id{0};
    int bin{0};
    int constraint_index{0};
    int term_count{0};
    int component0{-1};
    int component1{-1};
    amrex::Real coefficient0{amrex::Real(0.0)};
    amrex::Real coefficient1{amrex::Real(0.0)};
};

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
make_constraint_groups(const SBMLayout& layout);

/** Flatten complete semantic groups into bounded device descriptors. */
[[nodiscard]] std::vector<ConstraintDescriptor>
make_constraint_descriptors(const SBMLayout& layout);

[[nodiscard]] std::vector<AttachedPropertySupportDescriptor>
make_attached_property_support_descriptors(const SBMLayout& layout);

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
make_constraint_closure_chunks(const SBMLayout& layout, int max_groups);

struct EndpointTransform {
    amrex::Real L{0.0};
    amrex::Real H{0.0};
    amrex::Real tolerance{0.0};
    bool normalized_roundoff{false};
};

/**
 * Stable two-moment transform.  Only cancellation-sized negative endpoints
 * are normalized to zero; materially negative endpoints throw.
 */
[[nodiscard]] EndpointTransform
transform_two_moment(amrex::Real C, amrex::Real M,
                     amrex::Real lower, amrex::Real upper);

/** Inverse transform with finite-state and nonnegative-endpoint checks. */
[[nodiscard]] std::pair<amrex::Real, amrex::Real>
inverse_two_moment(amrex::Real L, amrex::Real H,
                   amrex::Real lower, amrex::Real upper);

/** Host reference for the finite property-support constraints. */
[[nodiscard]] bool
property_support_is_admissible(amrex::Real property, amrex::Real carrier_number,
                               amrex::Real lower, amrex::Real upper,
                               amrex::Real tolerance = 0.0) noexcept;

} // namespace erf_sbm

#endif
