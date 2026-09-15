#include "ERF_SBMConstraintGroups.H"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

namespace {

void add_constraint(ConstraintGroup& group, std::string id,
                    std::initializer_list<ConstraintTerm> terms)
{
    group.constraints.push_back({std::move(id), std::vector<ConstraintTerm>(terms)});
}

bool finite(const amrex::Real x) noexcept { return std::isfinite(x); }

} // namespace

bool ConstraintGroup::contains(const int component) const noexcept
{
    for (const int member : members) if (member == component) return true;
    return false;
}

amrex::Real ConstraintGroup::evaluate(const LinearConstraint& constraint,
                                      const std::vector<amrex::Real>& state) const
{
    amrex::Real result = amrex::Real(0.0);
    for (const auto& term : constraint.terms) {
        if (term.component < 0 || term.component >= static_cast<int>(state.size())) {
            throw std::out_of_range("constraint references an unavailable component");
        }
        result = std::fma(term.coefficient, state[static_cast<std::size_t>(term.component)], result);
    }
    return result;
}

bool ConstraintGroup::admissible(const std::vector<amrex::Real>& state,
                                 amrex::Real* minimum_margin,
                                 std::string* failed_constraint) const
{
    amrex::Real minimum = std::numeric_limits<amrex::Real>::infinity();
    for (const auto& constraint : constraints) {
        const amrex::Real value = evaluate(constraint, state);
        if (!finite(value)) {
            if (failed_constraint) *failed_constraint = constraint.semantic_id;
            if (minimum_margin) *minimum_margin = value;
            return false;
        }
        minimum = std::min(minimum, value);
        const amrex::Real scale = [&]() {
            amrex::Real s = amrex::Real(0.0);
            for (const auto& term : constraint.terms) {
                if (term.component >= 0 && term.component < static_cast<int>(state.size())) {
                    s += std::abs(term.coefficient * state[static_cast<std::size_t>(term.component)]);
                }
            }
            return s;
        }();
        const amrex::Real tolerance = amrex::Real(128.0) *
            std::numeric_limits<amrex::Real>::epsilon() * scale;
        if (value < -tolerance) {
            if (failed_constraint) *failed_constraint = constraint.semantic_id;
            if (minimum_margin) *minimum_margin = value;
            return false;
        }
    }
    if (minimum_margin) *minimum_margin = minimum;
    return true;
}

std::vector<ConstraintGroup> make_constraint_groups(const SBMLayout& layout)
{
    std::vector<ConstraintGroup> groups;
    for (const auto& population : layout.populations()) {
        const int nbins = population.grid.nbins();
        for (int bin = 0; bin < nbins; ++bin) {
            ConstraintGroup group;
            group.population_id = population.population_id;
            group.bin = bin;
            group.semantic_id = population.semantic_id + ".bin." + std::to_string(bin);
            group.moment_mode = population.moment_mode;
            const int mass = population.mass_offset + bin;
            group.members.push_back(mass);
            group.transport_members.push_back(mass);
            const int number = population.number_offset >= 0 ? population.number_offset + bin : -1;
            if (number >= 0) {
                group.members.push_back(number);
                group.transport_members.push_back(number);
                const amrex::Real lower = population.grid.edges()[static_cast<std::size_t>(bin)];
                const amrex::Real upper = population.grid.edges()[static_cast<std::size_t>(bin + 1)];
                const amrex::Real denominator = upper - lower;
                add_constraint(group, "number_nonnegative", {{number, 1.0}});
                add_constraint(group, "endpoint_low",
                               {{number, upper / denominator},
                                {mass, -1.0 / denominator}});
                add_constraint(group, "endpoint_high",
                               {{mass, 1.0 / denominator},
                                {number, -lower / denominator}});
            } else {
                add_constraint(group, "mass_nonnegative", {{mass, 1.0}});
            }

            for (std::size_t property_index = 0;
                 property_index < layout.attached_properties().size(); ++property_index) {
                const auto& property = layout.attached_properties()[property_index];
                if (property.carrier_population != population.population_id) continue;
                const int property_component = layout.property_offset(static_cast<int>(property_index)) + bin;
                group.members.push_back(property_component);
                group.transport_members.push_back(property_component);
                group.attached_property_indices.push_back(static_cast<int>(property_index));
                add_constraint(group, property.semantic_id + ".nonnegative", {{property_component, 1.0}});

                const bool has_upper = finite(property.support_max) && property.support_max >= property.support_min;
                if (has_upper) {
                    if (number >= 0) {
                        // Two-moment storage is physical (M,C), so support
                        // bounds are bounds on the carrier number C.
                        add_constraint(group, property.semantic_id + ".support_upper",
                                       {{number, property.support_max},
                                        {property_component, -1.0}});
                        if (property.support_min > amrex::Real(0.0)) {
                            add_constraint(group, property.semantic_id + ".support_lower",
                                           {{number, -property.support_min},
                                            {property_component, 1.0}});
                        }
                    } else {
                        const amrex::Real pivot = population.grid.pivot(bin);
                        add_constraint(group, property.semantic_id + ".support_upper",
                                       {{mass, property.support_max / pivot},
                                        {property_component, -1.0}});
                        if (property.support_min > amrex::Real(0.0)) {
                            add_constraint(group, property.semantic_id + ".support_lower",
                                           {{mass, -property.support_min / pivot},
                                            {property_component, 1.0}});
                        }
                    }
                }

                if (property.kind == PropertyKind::MassBoundedSubset) {
                    if (number >= 0) {
                        // The carrier mass is already the authoritative M
                        // component in two-moment storage.
                        add_constraint(group, property.semantic_id + ".carrier_mass_minus_subset",
                                       {{mass, 1.0}, {property_component, -1.0}});
                    } else {
                        add_constraint(group, property.semantic_id + ".carrier_mass_minus_subset",
                                       {{mass, 1.0}, {property_component, -1.0}});
                    }
                }
            }
            groups.push_back(std::move(group));
        }
    }
    return groups;
}

EndpointTransform transform_two_moment(const amrex::Real C, const amrex::Real M,
                                       const amrex::Real lower, const amrex::Real upper)
{
    if (!finite(C) || !finite(M) || !finite(lower) || !finite(upper) || !(lower < upper)) {
        throw std::invalid_argument("two-moment transform requires finite edges and strict lower < upper");
    }
    if (C < amrex::Real(0.0)) throw std::invalid_argument("two-moment count is negative");
    const amrex::Real scale = std::abs(M) + std::abs(lower * C) + std::abs(upper * C);
    const amrex::Real tolerance = amrex::Real(128.0) *
        std::numeric_limits<amrex::Real>::epsilon() * scale;
    const amrex::Real low_numerator = std::fma(upper, C, -M);
    const amrex::Real high_numerator = std::fma(-lower, C, M);
    if (low_numerator < -tolerance || high_numerator < -tolerance) {
        throw std::invalid_argument("materially non-realizable two-moment state");
    }
    if (C == amrex::Real(0.0)) {
        if (std::abs(M) > tolerance) throw std::invalid_argument("zero-number state carries mass");
        return {amrex::Real(0.0), amrex::Real(0.0), tolerance, std::abs(M) > amrex::Real(0.0)};
    }
    const amrex::Real denominator = upper - lower;
    amrex::Real L = low_numerator / denominator;
    amrex::Real H = high_numerator / denominator;
    bool normalized = false;
    if (L < amrex::Real(0.0)) { L = amrex::Real(0.0); normalized = true; }
    if (H < amrex::Real(0.0)) { H = amrex::Real(0.0); normalized = true; }
    return {L, H, tolerance / denominator, normalized};
}

std::pair<amrex::Real, amrex::Real> inverse_two_moment(const amrex::Real L, const amrex::Real H,
                                                       const amrex::Real lower, const amrex::Real upper)
{
    if (!finite(L) || !finite(H) || !finite(lower) || !finite(upper) || !(lower < upper) ||
        L < amrex::Real(0.0) || H < amrex::Real(0.0)) {
        throw std::invalid_argument("invalid endpoint transport state");
    }
    return {L + H, std::fma(lower, L, upper * H)};
}

bool property_support_is_admissible(const amrex::Real property,
                                    const amrex::Real carrier_number,
                                    const amrex::Real lower,
                                    const amrex::Real upper,
                                    const amrex::Real tolerance) noexcept
{
    if (!finite(property) || !finite(carrier_number) || !finite(lower) || !finite(upper) ||
        carrier_number < amrex::Real(0.0) || lower < amrex::Real(0.0) || upper < lower) return false;
    return property >= lower * carrier_number - tolerance &&
           property <= upper * carrier_number + tolerance;
}

} // namespace erf_sbm
