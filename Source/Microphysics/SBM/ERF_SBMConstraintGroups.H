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

/** Build deterministic groups from the immutable runtime layout. */
[[nodiscard]] std::vector<ConstraintGroup>
make_constraint_groups(const SBMLayout& layout);

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
