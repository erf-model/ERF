#include "ERF_SBMFCT.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace erf_sbm {

namespace {

amrex::Real form_value(const LinearConstraint& constraint,
                       const std::vector<amrex::Real>& values)
{
    amrex::Real result = amrex::Real(0.0);
    for (const auto& term : constraint.terms) {
        if (term.component < 0 || term.component >= static_cast<int>(values.size())) {
            throw std::invalid_argument("FCT constraint component is outside the state");
        }
        result = std::fma(term.coefficient, values[static_cast<std::size_t>(term.component)], result);
    }
    return result;
}

amrex::Real form_delta(const LinearConstraint& constraint,
                       const std::vector<amrex::Real>& delta,
                       const amrex::Real sign,
                       const amrex::Real volume)
{
    amrex::Real result = amrex::Real(0.0);
    for (const auto& term : constraint.terms) {
        result = std::fma(term.coefficient,
                          sign * delta[static_cast<std::size_t>(term.component)] / volume,
                          result);
    }
    return result;
}

} // namespace

FCTResult limit_grouped(const std::vector<amrex::Real>& low_state,
                        const int ncell, const int ncomp,
                        const std::vector<FCTFaceTransfer>& faces,
                        const std::vector<ConstraintGroup>& groups,
                        const int chunk_size)
{
    if (ncell <= 0 || ncomp <= 0 || low_state.size() != static_cast<std::size_t>(ncell * ncomp)) {
        throw std::invalid_argument("invalid flattened FCT state");
    }
    if (chunk_size < 0) throw std::invalid_argument("FCT chunk size must be nonnegative");
    std::set<std::pair<int,int>> owned_faces;
    for (const auto& face : faces) {
        if (face.left_cell < 0 || face.left_cell >= ncell || face.right_cell < 0 ||
            face.right_cell >= ncell || !(face.left_volume > amrex::Real(0.0)) ||
            !(face.right_volume > amrex::Real(0.0)) || face.low.size() != static_cast<std::size_t>(ncomp) ||
            face.high.size() != static_cast<std::size_t>(ncomp)) {
            throw std::invalid_argument("invalid FCT face transfer");
        }
        const auto ownership_key = std::minmax(face.left_cell, face.right_cell);
        if (!owned_faces.emplace(ownership_key.first, ownership_key.second).second) {
            throw std::invalid_argument("duplicate SBM face ownership in grouped FCT input");
        }
    }
    for (const auto& group : groups) {
        for (const auto& constraint : group.constraints) {
            for (const auto& term : constraint.terms) {
                if (term.component < 0 || term.component >= ncomp) {
                    throw std::invalid_argument("FCT group references an unavailable component");
                }
            }
        }
    }

    // The budgets are per cell, per complete group, per atomic constraint.
    // A single budget per group is insufficient: different constraints can
    // have different available margins and adverse face demands.
    std::vector<std::size_t> constraint_offsets(groups.size() + 1, 0);
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
        constraint_offsets[gi + 1] = constraint_offsets[gi] + groups[gi].constraints.size();
    }
    const std::size_t nconstraints = constraint_offsets.back();
    std::vector<amrex::Real> budgets(static_cast<std::size_t>(ncell) * nconstraints, 0.0);
    std::vector<amrex::Real> margins(static_cast<std::size_t>(ncell) * nconstraints,
                                     std::numeric_limits<amrex::Real>::infinity());
    for (int cell = 0; cell < ncell; ++cell) {
        std::vector<amrex::Real> state(low_state.begin() + cell * ncomp,
                                       low_state.begin() + (cell + 1) * ncomp);
        for (std::size_t gi = 0; gi < groups.size(); ++gi) {
            std::string failed;
            if (!groups[gi].admissible(state, nullptr, &failed)) {
                throw std::domain_error("FCT low-order state is inadmissible in " + groups[gi].semantic_id +
                                        " constraint " + failed);
            }
            for (std::size_t ci = 0; ci < groups[gi].constraints.size(); ++ci) {
                margins[(static_cast<std::size_t>(cell) * nconstraints) +
                        constraint_offsets[gi] + ci] =
                    form_value(groups[gi].constraints[ci], state);
            }
        }
    }
    for (const auto& face : faces) {
        std::vector<amrex::Real> delta(static_cast<std::size_t>(ncomp));
        for (int c = 0; c < ncomp; ++c) {
            const amrex::Real difference = face.high[static_cast<std::size_t>(c)] -
                                           face.low[static_cast<std::size_t>(c)];
            if (!std::isfinite(difference)) throw std::domain_error("nonfinite FCT candidate transfer");
            delta[static_cast<std::size_t>(c)] = difference;
        }
        for (std::size_t gi = 0; gi < groups.size(); ++gi) {
            const auto& group = groups[gi];
            const std::size_t first_constraint = constraint_offsets[gi];
            const std::size_t group_chunk = chunk_size > 0 ? static_cast<std::size_t>(chunk_size) : group.constraints.size();
            for (std::size_t chunk_begin = 0; chunk_begin < group.constraints.size(); chunk_begin += group_chunk) {
                const std::size_t chunk_end = std::min(group.constraints.size(), chunk_begin + group_chunk);
                for (std::size_t ci = chunk_begin; ci < chunk_end; ++ci) {
                const auto& constraint = group.constraints[ci];
                const amrex::Real left_change = form_delta(constraint, delta, -1.0, face.left_volume);
                const amrex::Real right_change = form_delta(constraint, delta, 1.0, face.right_volume);
                budgets[static_cast<std::size_t>(face.left_cell) * nconstraints + first_constraint + ci] +=
                    std::max(amrex::Real(0.0), -left_change);
                budgets[static_cast<std::size_t>(face.right_cell) * nconstraints + first_constraint + ci] +=
                    std::max(amrex::Real(0.0), -right_change);
                }
            }
        }
    }

    FCTResult result;
    result.accepted_faces = faces;
    result.limiter.resize(faces.size(), amrex::Real(1.0));
    result.updated_state = low_state;
    result.minimum_margin = std::numeric_limits<amrex::Real>::infinity();

    for (std::size_t fi = 0; fi < faces.size(); ++fi) {
        const auto& face = faces[fi];
        amrex::Real lambda = amrex::Real(1.0);
        for (std::size_t gi = 0; gi < groups.size(); ++gi) {
            const auto& group = groups[gi];
            const std::size_t first_constraint = constraint_offsets[gi];
            const std::size_t group_chunk = chunk_size > 0 ? static_cast<std::size_t>(chunk_size) : group.constraints.size();
            for (std::size_t chunk_begin = 0; chunk_begin < group.constraints.size(); chunk_begin += group_chunk) {
                const std::size_t chunk_end = std::min(group.constraints.size(), chunk_begin + group_chunk);
                for (std::size_t ci = chunk_begin; ci < chunk_end; ++ci) {
                const auto& constraint = group.constraints[ci];
                // The limiter acts on the antidiffusive correction only.
                // `high - low` is the correction because low_state already
                // contains the complete low-order update.
                const amrex::Real left_change = form_delta(constraint, face.high, -1.0, face.left_volume) -
                    form_delta(constraint, face.low, -1.0, face.left_volume);
                const amrex::Real right_change = form_delta(constraint, face.high, 1.0, face.right_volume) -
                    form_delta(constraint, face.low, 1.0, face.right_volume);
                const auto budget_left = budgets[static_cast<std::size_t>(face.left_cell) * nconstraints + first_constraint + ci];
                const auto budget_right = budgets[static_cast<std::size_t>(face.right_cell) * nconstraints + first_constraint + ci];
                const amrex::Real margin_left = margins[static_cast<std::size_t>(face.left_cell) * nconstraints + first_constraint + ci];
                const amrex::Real margin_right = margins[static_cast<std::size_t>(face.right_cell) * nconstraints + first_constraint + ci];
                if (left_change < amrex::Real(0.0) && budget_left > amrex::Real(0.0)) {
                    lambda = std::min(lambda, std::min(amrex::Real(1.0), margin_left / budget_left));
                }
                if (right_change < amrex::Real(0.0) && budget_right > amrex::Real(0.0)) {
                    lambda = std::min(lambda, std::min(amrex::Real(1.0), margin_right / budget_right));
                }
                }
            }
        }
        if (!std::isfinite(lambda) || lambda < amrex::Real(0.0)) {
            throw std::domain_error("FCT produced an invalid face limiter");
        }
        lambda = std::min(lambda, amrex::Real(1.0));
        result.limiter[fi] = lambda;
        for (int c = 0; c < ncomp; ++c) {
            result.accepted_faces[fi].low[static_cast<std::size_t>(c)] =
                face.low[static_cast<std::size_t>(c)] + lambda *
                (face.high[static_cast<std::size_t>(c)] - face.low[static_cast<std::size_t>(c)]);
            // low_state is the already-updated donor/diffusive candidate, so
            // FCT applies only the bounded antidiffusive correction.  The
            // accepted face vector above is still retained for ledger and
            // conservation diagnostics.
            const amrex::Real correction = lambda *
                (face.high[static_cast<std::size_t>(c)] - face.low[static_cast<std::size_t>(c)]);
            result.updated_state[static_cast<std::size_t>(face.left_cell * ncomp + c)] -= correction / face.left_volume;
            result.updated_state[static_cast<std::size_t>(face.right_cell * ncomp + c)] += correction / face.right_volume;
        }
    }

    for (int cell = 0; cell < ncell; ++cell) {
        std::vector<amrex::Real> state(result.updated_state.begin() + cell * ncomp,
                                       result.updated_state.begin() + (cell + 1) * ncomp);
        for (const auto& group : groups) {
            amrex::Real margin = 0.0;
            std::string failed;
            if (!group.admissible(state, &margin, &failed)) {
                throw std::domain_error("grouped FCT accepted state is inadmissible in " + group.semantic_id +
                                        " constraint " + failed);
            }
            result.minimum_margin = std::min(result.minimum_margin, margin);
        }
    }
    if (result.minimum_margin == std::numeric_limits<amrex::Real>::infinity()) {
        result.minimum_margin = amrex::Real(0.0);
    }
    return result;
}

std::vector<amrex::Real> accepted_stage_correction(const std::vector<amrex::Real>& low_transfer,
                                                   const std::vector<amrex::Real>& high_transfer,
                                                   const amrex::Real stage_interval)
{
    if (low_transfer.size() != high_transfer.size() || !std::isfinite(stage_interval) || stage_interval < 0.0) {
        throw std::invalid_argument("invalid FCT stage correction");
    }
    std::vector<amrex::Real> result(low_transfer.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = stage_interval * (high_transfer[i] - low_transfer[i]);
    }
    return result;
}

StageWeightContract stage_weight_contract(const ::erf_auxiliary::IntegrationMethod method, const int stage,
                                          const bool completes_step) noexcept
{
    if (method == ::erf_auxiliary::IntegrationMethod::CompressibleRK3) {
        return {1.0, 1.0, 1.0, completes_step ? 1.0 : 0.0};
    }
    if (stage == 0) return {1.0, 1.0, 1.0, 0.5};
    return {0.5, 0.5, 0.5, 0.5};
}

} // namespace erf_sbm
