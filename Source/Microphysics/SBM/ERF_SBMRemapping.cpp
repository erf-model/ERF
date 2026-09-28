#include "ERF_SBMRemapping.H"

#include "ERF_SpectralGrid.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace erf_sbm {

namespace {

const PopulationLayout* find_population(const SBMLayout& layout, const int population_id)
{
    const auto& populations = layout.populations();
    const auto found = std::find_if(populations.begin(), populations.end(),
        [population_id](const PopulationLayout& population) {
            return population.population_id == population_id;
        });
    return found == populations.end() ? nullptr : &*found;
}

std::vector<int> property_indices(const SBMLayout& layout, const int population_id)
{
    std::vector<int> indices;
    for (std::size_t index = 0; index < layout.attached_properties().size(); ++index) {
        if (layout.attached_properties()[index].carrier_population == population_id) {
            indices.push_back(static_cast<int>(index));
        }
    }
    return indices;
}

bool finite_nonnegative(const amrex::Real value) noexcept
{
    return std::isfinite(value) && value >= amrex::Real(0.0);
}

bool finite_product(const amrex::Real left, const amrex::Real right,
                    amrex::Real& result) noexcept
{
    result = left * right;
    return std::isfinite(result) && result >= amrex::Real(0.0) &&
           !(left > amrex::Real(0.0) && right > amrex::Real(0.0) &&
             result == amrex::Real(0.0));
}

} // namespace

PacketApplicationResult
apply_packet_routing(const SBMLayout& layout, const int population_id,
                     const PacketRoutingPlan& plan,
                     const std::vector<amrex::Real>& property_packet_amounts,
                     std::vector<amrex::Real>& candidate_state)
{
    PacketApplicationResult result;
    result.status = plan.status;
    if (plan.status != RemapStatus::Ok && plan.status != RemapStatus::ZeroWaterResidual) {
        return result;
    }
    const auto* population = find_population(layout, population_id);
    if (population == nullptr || candidate_state.size() != static_cast<std::size_t>(layout.ncomp()) ||
        plan.destination_count < 0 || plan.destination_count > 2 ||
        !finite_nonnegative(plan.packet_number) || !finite_nonnegative(plan.packet_particle_mass) ||
        !std::isfinite(plan.residual_number_weight) || plan.residual_number_weight < amrex::Real(0.0) ||
        plan.residual_number_weight > amrex::Real(1.0)) {
        result.status = RemapStatus::Invalid;
        return result;
    }
    const auto properties = property_indices(layout, population_id);
    if (property_packet_amounts.size() != properties.size()) {
        result.status = RemapStatus::Invalid;
        return result;
    }
    for (const auto amount : property_packet_amounts) {
        if (!finite_nonnegative(amount)) {
            result.status = RemapStatus::Invalid;
            return result;
        }
    }
    if (plan.packet_number == amrex::Real(0.0) &&
        std::any_of(property_packet_amounts.begin(), property_packet_amounts.end(),
                    [](const amrex::Real amount) { return amount != amrex::Real(0.0); })) {
        result.status = RemapStatus::Invalid;
        return result;
    }

    amrex::Real weight_sum = plan.residual_number_weight;
    for (int destination = 0; destination < plan.destination_count; ++destination) {
        const auto& route = plan.destinations[destination];
        if (route.bin < 0 || route.bin >= population->grid.nbins() ||
            !std::isfinite(route.number_weight) || route.number_weight < amrex::Real(0.0) ||
            route.number_weight > amrex::Real(1.0) || !finite_nonnegative(route.particle_mass)) {
            result.status = RemapStatus::Invalid;
            return result;
        }
        for (int previous = 0; previous < destination; ++previous) {
            if (plan.destinations[previous].bin == route.bin) {
                result.status = RemapStatus::Invalid;
                return result;
            }
        }
        weight_sum += route.number_weight;
    }
    if (plan.packet_number > amrex::Real(0.0)) {
        const amrex::Real weight_tolerance = amrex::Real(16.0) *
            std::numeric_limits<amrex::Real>::epsilon();
        if (!std::isfinite(weight_sum) || std::abs(weight_sum - amrex::Real(1.0)) > weight_tolerance ||
            ((plan.status == RemapStatus::ZeroWaterResidual) !=
             (plan.residual_number_weight > amrex::Real(0.0)))) {
            result.status = RemapStatus::Invalid;
            return result;
        }
    } else if (plan.destination_count != 0 || plan.residual_number_weight != amrex::Real(0.0)) {
        result.status = RemapStatus::Invalid;
        return result;
    }

    std::vector<amrex::Real> next = candidate_state;
    if (!finite_product(plan.packet_number, plan.residual_number_weight, result.residual_number)) {
        result.status = RemapStatus::Invalid;
        return result;
    }
    result.residual_water_mass = amrex::Real(0.0);
    result.residual_properties.resize(properties.size(), amrex::Real(0.0));
    for (std::size_t local_property = 0; local_property < properties.size(); ++local_property) {
        amrex::Real residual = amrex::Real(0.0);
        if (!finite_product(property_packet_amounts[local_property],
                            plan.residual_number_weight, residual)) {
            result.status = RemapStatus::Invalid;
            return result;
        }
        result.residual_properties[local_property] = residual;
    }

    for (int destination = 0; destination < plan.destination_count; ++destination) {
        const auto& route = plan.destinations[destination];
        amrex::Real packet_number_here = amrex::Real(0.0);
        if (!finite_product(plan.packet_number, route.number_weight, packet_number_here)) {
            result.status = RemapStatus::Invalid;
            return result;
        }
        amrex::Real mass_increment = amrex::Real(0.0);
        const amrex::Real deposited_mass = population->moment_mode == MomentMode::OneMoment ?
            route.particle_mass : plan.packet_particle_mass;
        if (!finite_product(packet_number_here, deposited_mass, mass_increment)) {
            result.status = RemapStatus::Invalid;
            return result;
        }
        const int mass_component = population->mass_offset + route.bin;
        if (!finite_nonnegative(next[static_cast<std::size_t>(mass_component)])) {
            result.status = RemapStatus::Invalid;
            return result;
        }
        const amrex::Real new_mass = next[static_cast<std::size_t>(mass_component)] + mass_increment;
        if (!finite_nonnegative(new_mass)) {
            result.status = RemapStatus::Invalid;
            return result;
        }
        next[static_cast<std::size_t>(mass_component)] = new_mass;

        if (population->number_offset >= 0) {
            const int number_component = population->number_offset + route.bin;
            if (!finite_nonnegative(next[static_cast<std::size_t>(number_component)])) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            const amrex::Real new_number =
                next[static_cast<std::size_t>(number_component)] + packet_number_here;
            if (!finite_nonnegative(new_number)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            next[static_cast<std::size_t>(number_component)] = new_number;
        }

        for (std::size_t local_property = 0; local_property < properties.size(); ++local_property) {
            const int component = layout.property_offset(properties[local_property]) + route.bin;
            amrex::Real increment = amrex::Real(0.0);
            if (!finite_product(property_packet_amounts[local_property], route.number_weight, increment) ||
                !finite_nonnegative(next[static_cast<std::size_t>(component)])) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            const amrex::Real value = next[static_cast<std::size_t>(component)] + increment;
            if (!finite_nonnegative(value)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            next[static_cast<std::size_t>(component)] = value;
        }
    }

    candidate_state.swap(next);
    return result;
}

ReconstructionStatus
reconstruct_bin(const SBMLayout& layout, const int population_id, const int bin,
                const std::vector<amrex::Real>& persisted_state,
                ReconstructionDelta& reconstruction)
{
    const auto* population = find_population(layout, population_id);
    if (population == nullptr || bin < 0 || bin >= population->grid.nbins() ||
        persisted_state.size() != static_cast<std::size_t>(layout.ncomp())) {
        return ReconstructionStatus::Invalid;
    }
    ReconstructionDelta next;
    next.population_id = population_id;
    next.bin = bin;
    const auto properties = property_indices(layout, population_id);
    next.property_per_particle.assign(properties.size(), amrex::Real(0.0));
    const amrex::Real mass = persisted_state[static_cast<std::size_t>(population->mass_offset + bin)];
    if (!finite_nonnegative(mass)) return ReconstructionStatus::Invalid;

    if (population->moment_mode == MomentMode::OneMoment) {
        const amrex::Real pivot = population->grid.pivot(bin);
        next.number = mass / pivot;
        next.particle_mass = pivot;
        if (mass > amrex::Real(0.0) && next.number == amrex::Real(0.0)) {
            return ReconstructionStatus::Invalid;
        }
    } else {
        const amrex::Real number =
            persisted_state[static_cast<std::size_t>(population->number_offset + bin)];
        if (!finite_nonnegative(number)) return ReconstructionStatus::Invalid;
        EndpointTransform endpoints;
        try {
            endpoints = SpectralGrid::two_moment_to_endpoints(
                number, mass, population->grid.edges()[static_cast<std::size_t>(bin)],
                population->grid.edges()[static_cast<std::size_t>(bin + 1)]);
        } catch (const std::invalid_argument&) {
            return ReconstructionStatus::Invalid;
        }
        next.normalized_roundoff = endpoints.normalized_roundoff;
        next.number = number;
        if (number > amrex::Real(0.0)) {
            next.particle_mass = mass / number;
            const amrex::Real lower = population->grid.edges()[static_cast<std::size_t>(bin)];
            const amrex::Real upper = population->grid.edges()[static_cast<std::size_t>(bin + 1)];
            if (next.particle_mass < lower && endpoints.H == amrex::Real(0.0)) {
                next.particle_mass = lower;
                next.normalized_roundoff = true;
            } else if (next.particle_mass > upper && endpoints.L == amrex::Real(0.0)) {
                next.particle_mass = upper;
                next.normalized_roundoff = true;
            }
            if (!std::isfinite(next.particle_mass) ||
                (mass > amrex::Real(0.0) && next.particle_mass == amrex::Real(0.0)) ||
                next.particle_mass < lower || next.particle_mass > upper) {
                return ReconstructionStatus::Invalid;
            }
        }
    }

    if (!std::isfinite(next.number) || next.number < amrex::Real(0.0)) {
        return ReconstructionStatus::Invalid;
    }
    next.empty = next.number == amrex::Real(0.0);
    if (next.empty) next.particle_mass = amrex::Real(0.0);
    for (std::size_t local_property = 0; local_property < properties.size(); ++local_property) {
        const amrex::Real amount = persisted_state[static_cast<std::size_t>(
            layout.property_offset(properties[local_property]) + bin)];
        if (!finite_nonnegative(amount)) return ReconstructionStatus::Invalid;
        if (next.empty) {
            if (amount != amrex::Real(0.0)) return ReconstructionStatus::Invalid;
        } else {
            next.property_per_particle[local_property] = amount / next.number;
            if (!finite_nonnegative(next.property_per_particle[local_property]) ||
                (amount > amrex::Real(0.0) &&
                 next.property_per_particle[local_property] == amrex::Real(0.0))) {
                return ReconstructionStatus::Invalid;
            }
        }
    }
    reconstruction = std::move(next);
    return reconstruction.empty ? ReconstructionStatus::Empty : ReconstructionStatus::Populated;
}

bool integrate_interval(const ReconstructionDelta& reconstruction,
                        const amrex::Real lower, const amrex::Real upper,
                        const bool include_upper, IntegratedMoments& integral)
{
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(lower < upper) ||
        !finite_nonnegative(reconstruction.number) ||
        !finite_nonnegative(reconstruction.particle_mass)) {
        return false;
    }
    IntegratedMoments next;
    next.attached_properties.assign(reconstruction.property_per_particle.size(), amrex::Real(0.0));
    for (const auto amount_per_particle : reconstruction.property_per_particle) {
        if (!finite_nonnegative(amount_per_particle)) return false;
    }
    const bool contains = !reconstruction.empty && reconstruction.particle_mass >= lower &&
        (reconstruction.particle_mass < upper ||
         (include_upper && reconstruction.particle_mass == upper));
    if (contains) {
        next.number = reconstruction.number;
        if (!finite_product(reconstruction.number, reconstruction.particle_mass, next.water_mass)) {
            return false;
        }
        for (std::size_t property = 0; property < reconstruction.property_per_particle.size(); ++property) {
            if (!finite_product(reconstruction.number, reconstruction.property_per_particle[property],
                                next.attached_properties[property])) {
                return false;
            }
        }
    }
    integral = std::move(next);
    return true;
}

bool project_reconstruction_bin(const SBMLayout& layout,
                                const ReconstructionDelta& reconstruction,
                                std::vector<amrex::Real>& candidate_state)
{
    const auto* population = find_population(layout, reconstruction.population_id);
    if (population == nullptr || reconstruction.bin < 0 || reconstruction.bin >= population->grid.nbins() ||
        candidate_state.size() != static_cast<std::size_t>(layout.ncomp()) ||
        !finite_nonnegative(reconstruction.number) || !finite_nonnegative(reconstruction.particle_mass)) {
        return false;
    }
    const auto properties = property_indices(layout, reconstruction.population_id);
    if (properties.size() != reconstruction.property_per_particle.size()) return false;
    if ((!reconstruction.empty && reconstruction.number <= amrex::Real(0.0)) ||
        (reconstruction.empty && (reconstruction.number != amrex::Real(0.0) ||
                                  reconstruction.particle_mass != amrex::Real(0.0)))) {
        return false;
    }
    if (!reconstruction.empty) {
        if (population->moment_mode == MomentMode::OneMoment) {
            if (reconstruction.particle_mass != population->grid.pivot(reconstruction.bin)) return false;
        } else {
            const amrex::Real lower = population->grid.edges()[static_cast<std::size_t>(reconstruction.bin)];
            const amrex::Real upper = population->grid.edges()[static_cast<std::size_t>(reconstruction.bin + 1)];
            if (reconstruction.particle_mass < lower || reconstruction.particle_mass > upper) return false;
        }
    }

    std::vector<amrex::Real> next = candidate_state;
    amrex::Real mass = amrex::Real(0.0);
    if (!finite_product(reconstruction.number, reconstruction.particle_mass, mass)) return false;
    next[static_cast<std::size_t>(population->mass_offset + reconstruction.bin)] = mass;
    if (population->number_offset >= 0) {
        next[static_cast<std::size_t>(population->number_offset + reconstruction.bin)] =
            reconstruction.number;
    }
    for (std::size_t local_property = 0; local_property < properties.size(); ++local_property) {
        amrex::Real amount = amrex::Real(0.0);
        if (!finite_product(reconstruction.number, reconstruction.property_per_particle[local_property], amount)) {
            return false;
        }
        next[static_cast<std::size_t>(layout.property_offset(properties[local_property]) + reconstruction.bin)] =
            amount;
    }
    candidate_state.swap(next);
    return true;
}

} // namespace erf_sbm
