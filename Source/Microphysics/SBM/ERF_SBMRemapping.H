#ifndef ERF_SBM_REMAPPING_H_
#define ERF_SBM_REMAPPING_H_

#include "ERF_SBMLayout.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>
#include <vector>

namespace erf_sbm {

enum class RemapStatus {
    Ok, ZeroWaterResidual, BelowSupportedGrid, Overflow, NumericalUnderflow, Invalid
};
enum class ReconstructionStatus { Populated, Empty, Invalid };

/**
 * Compact, non-owning scientific context for one population. Its referenced
 * arrays must be immutable and accessible to the execution backend for the
 * full lifetime of any plan or kernel work using the view. Device callers
 * construct the same POD with device-accessible metadata arrays.
 */
struct PopulationRemapView {
    int population_id{-1};
    MomentMode moment_mode{MomentMode::OneMoment};
    CoordinateKind coordinate_kind{CoordinateKind::Mass};
    const amrex::Real* edges{nullptr};
    const amrex::Real* pivots{nullptr};
    int nbins{0};
    int mass_offset{-1};
    int number_offset{-1};
    int state_components{0};
    const int* property_component_offsets{nullptr};
    int property_count{0};
};

/** Host adapter building a stable view into the immutable SBMLayout. */
[[nodiscard]] PopulationRemapView
population_remap_view (const SBMLayout& layout, int population_id) noexcept;

/**
 * Constant-time structural/context checks for already qualified immutable
 * metadata. Callers must derive views from a validated SBMLayout or call the
 * full validator below once before using a raw view in the remapping cores.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_population_remap_context (const PopulationRemapView& population) noexcept
{
    if (population.population_id < 0 ||
        population.coordinate_kind != CoordinateKind::Mass ||
        (population.moment_mode != MomentMode::OneMoment &&
         population.moment_mode != MomentMode::TwoMoment) ||
        population.nbins <= 0 || population.edges == nullptr ||
        population.state_components < population.nbins || population.mass_offset < 0 ||
        population.property_count < 0 ||
        population.property_count > population.state_components / population.nbins ||
        (population.property_count > 0 && population.property_component_offsets == nullptr)) {
        return false;
    }
    if (population.mass_offset > population.state_components - population.nbins) return false;
    if (population.moment_mode == MomentMode::OneMoment) {
        if (population.number_offset != -1 || population.pivots == nullptr) return false;
    } else if (population.number_offset < 0 ||
               population.number_offset > population.state_components - population.nbins) {
        return false;
    }
    if (population.moment_mode == MomentMode::TwoMoment &&
        population.number_offset < population.mass_offset + population.nbins &&
        population.mass_offset < population.number_offset + population.nbins) {
        return false;
    }
    return true;
}

/** Full O(nbins + properties^2) validator for untrusted/arbitrary views. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_population_remap_view (const PopulationRemapView& population) noexcept
{
    if (!valid_population_remap_context(population)) return false;
    for (int edge = 0; edge <= population.nbins; ++edge) {
        const amrex::Real value = population.edges[edge];
        if (!amrex::Math::isfinite(value) || value < amrex::Real(0.0) ||
            (edge > 0 && !(population.edges[edge - 1] < value))) {
            return false;
        }
    }
    if (population.moment_mode == MomentMode::OneMoment) {
        for (int bin = 0; bin < population.nbins; ++bin) {
            const amrex::Real pivot = population.pivots[bin];
            if (!amrex::Math::isfinite(pivot) || pivot <= amrex::Real(0.0) ||
                pivot < population.edges[bin] || pivot > population.edges[bin + 1] ||
                (bin > 0 && !(population.pivots[bin - 1] < pivot))) {
                return false;
            }
        }
    }
    for (int property = 0; property < population.property_count; ++property) {
        const int offset = population.property_component_offsets[property];
        if (offset < 0 || offset > population.state_components - population.nbins ||
            (offset < population.mass_offset + population.nbins &&
             population.mass_offset < offset + population.nbins)) {
            return false;
        }
        if (population.moment_mode == MomentMode::TwoMoment &&
            offset < population.number_offset + population.nbins &&
            population.number_offset < offset + population.nbins) {
            return false;
        }
        for (int previous = 0; previous < property; ++previous) {
            const int previous_offset = population.property_component_offsets[previous];
            if (offset < previous_offset + population.nbins &&
                previous_offset < offset + population.nbins) {
                return false;
            }
        }
    }
    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool same_population_context (const PopulationRemapView& left,
                              const PopulationRemapView& right) noexcept
{
    return left.population_id == right.population_id &&
           left.moment_mode == right.moment_mode &&
           left.coordinate_kind == right.coordinate_kind &&
           left.edges == right.edges && left.pivots == right.pivots &&
           left.nbins == right.nbins && left.mass_offset == right.mass_offset &&
           left.number_offset == right.number_offset &&
           left.state_components == right.state_components &&
           left.property_component_offsets == right.property_component_offsets &&
           left.property_count == right.property_count;
}

struct RouteDestination {
    int bin{-1};
    amrex::Real number_weight{0.0};
    amrex::Real particle_mass{0.0};
};

/** Device-safe, fixed-capacity route bound to its originating population view. */
struct PacketRoutingPlan {
    PopulationRemapView population{};
    RemapStatus status{RemapStatus::Invalid};
    amrex::Real packet_number{0.0};
    amrex::Real packet_particle_mass{0.0};
    int destination_count{0};
    RouteDestination destinations[2]{};
    amrex::Real residual_number_weight{0.0};
    bool normalized_roundoff{false};
};

namespace remap_detail {

enum class ProductStatus { Ok, Underflow, Invalid };
enum class QuotientStatus { Ok, Underflow, Invalid };

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool finite_nonnegative (const amrex::Real value) noexcept
{
    return amrex::Math::isfinite(value) && value >= amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ProductStatus checked_product (const amrex::Real left, const amrex::Real right,
                               amrex::Real& result) noexcept
{
    if (!finite_nonnegative(left) || !finite_nonnegative(right)) return ProductStatus::Invalid;
    result = left * right;
    if (!amrex::Math::isfinite(result)) return ProductStatus::Invalid;
    if (left > amrex::Real(0.0) && right > amrex::Real(0.0) &&
        result == amrex::Real(0.0)) return ProductStatus::Underflow;
    return ProductStatus::Ok;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
QuotientStatus checked_quotient (const amrex::Real numerator,
                                 const amrex::Real denominator,
                                 amrex::Real& result) noexcept
{
    if (!finite_nonnegative(numerator) || !amrex::Math::isfinite(denominator) ||
        !(denominator > amrex::Real(0.0))) return QuotientStatus::Invalid;
    result = numerator / denominator;
    if (!amrex::Math::isfinite(result)) return QuotientStatus::Invalid;
    if (numerator > amrex::Real(0.0) && result == amrex::Real(0.0)) {
        return QuotientStatus::Underflow;
    }
    return QuotientStatus::Ok;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool finite_product (const amrex::Real left, const amrex::Real right,
                     amrex::Real& result) noexcept
{
    return checked_product(left, right, result) == ProductStatus::Ok;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RemapStatus remap_product_status (const ProductStatus status) noexcept
{
    return status == ProductStatus::Underflow ? RemapStatus::NumericalUnderflow :
           status == ProductStatus::Ok ? RemapStatus::Ok : RemapStatus::Invalid;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RemapStatus remap_quotient_status (const QuotientStatus status) noexcept
{
    return status == QuotientStatus::Underflow ? RemapStatus::NumericalUnderflow :
           status == QuotientStatus::Ok ? RemapStatus::Ok : RemapStatus::Invalid;
}

/** Accept at most four representable steps outward from a global endpoint. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool within_boundary_roundoff (const amrex::Real value, const amrex::Real endpoint,
                               const bool lower_boundary) noexcept
{
    const amrex::Real direction = lower_boundary ?
        -std::numeric_limits<amrex::Real>::infinity() :
         std::numeric_limits<amrex::Real>::infinity();
    amrex::Real allowance = endpoint;
    for (int step = 0; step < 4; ++step) {
        allowance = std::nextafter(allowance, direction);
    }
    return lower_boundary ? value >= allowance : value <= allowance;
}

/** Validate a route against its local claimed bin(s), without searching. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool route_plan_is_consistent (const PopulationRemapView& population,
                               const PacketRoutingPlan& plan) noexcept
{
    if (!valid_population_remap_context(population) ||
        !same_population_context(population, plan.population) ||
        !amrex::Math::isfinite(plan.packet_number) ||
        !amrex::Math::isfinite(plan.packet_particle_mass) ||
        plan.packet_number < amrex::Real(0.0) || plan.packet_particle_mass < amrex::Real(0.0) ||
        plan.destination_count < 0 || plan.destination_count > 2 ||
        !finite_nonnegative(plan.residual_number_weight) ||
        plan.residual_number_weight > amrex::Real(1.0)) return false;

    if (plan.packet_number == amrex::Real(0.0)) {
        return plan.status == RemapStatus::Ok && plan.destination_count == 0 &&
               plan.residual_number_weight == amrex::Real(0.0) &&
               !plan.normalized_roundoff;
    }

    if (population.moment_mode == MomentMode::TwoMoment) {
        if (plan.packet_particle_mass == amrex::Real(0.0)) {
            return plan.status == RemapStatus::ZeroWaterResidual &&
                   plan.destination_count == 0 &&
                   plan.residual_number_weight == amrex::Real(1.0) &&
                   !plan.normalized_roundoff;
        }
        const amrex::Real lower = population.edges[0];
        const amrex::Real top = population.edges[population.nbins];
        int required_bin = -1;
        amrex::Real expected_mass = plan.packet_particle_mass;
        bool normalized = false;
        if (plan.packet_particle_mass < lower) {
            if (!within_boundary_roundoff(plan.packet_particle_mass, lower, true)) return false;
            required_bin = 0;
            expected_mass = lower;
            normalized = true;
        } else if (plan.packet_particle_mass > top) {
            if (!within_boundary_roundoff(plan.packet_particle_mass, top, false)) return false;
            required_bin = population.nbins - 1;
            expected_mass = top;
            normalized = true;
        }
        if (plan.status != RemapStatus::Ok || plan.destination_count != 1 ||
            plan.residual_number_weight != amrex::Real(0.0) ||
            plan.normalized_roundoff != normalized) return false;
        const auto& route = plan.destinations[0];
        return (required_bin < 0 || route.bin == required_bin) &&
               route.number_weight == amrex::Real(1.0) &&
               route.particle_mass == expected_mass && route.bin >= 0 &&
               route.bin < population.nbins && expected_mass >= population.edges[route.bin] &&
               (expected_mass < population.edges[route.bin + 1] ||
                (route.bin == population.nbins - 1 &&
                 expected_mass <= population.edges[route.bin + 1]));
    }

    const amrex::Real first = population.pivots[0];
    const amrex::Real last = population.pivots[population.nbins - 1];
    if (plan.packet_particle_mass == amrex::Real(0.0)) {
        return plan.status == RemapStatus::ZeroWaterResidual &&
               plan.destination_count == 0 &&
               plan.residual_number_weight == amrex::Real(1.0) &&
               !plan.normalized_roundoff;
    }
    if (plan.packet_particle_mass < first) {
        amrex::Real wet_weight = amrex::Real(0.0);
        if (checked_quotient(plan.packet_particle_mass, first, wet_weight) !=
            QuotientStatus::Ok) return false;
        return plan.status == RemapStatus::ZeroWaterResidual && !plan.normalized_roundoff &&
               plan.destination_count == 1 && plan.residual_number_weight ==
               amrex::Real(1.0) - wet_weight && plan.destinations[0].bin == 0 &&
               plan.destinations[0].number_weight == wet_weight &&
               plan.destinations[0].particle_mass == first;
    }
    if (plan.packet_particle_mass > last || plan.status == RemapStatus::ZeroWaterResidual) {
        return false;
    }
    if (plan.destination_count == 1) {
        const int pivot_bin = plan.destinations[0].bin;
        if (pivot_bin < 0 || pivot_bin >= population.nbins) return false;
        return plan.status == RemapStatus::Ok && !plan.normalized_roundoff &&
               plan.destination_count == 1 && plan.residual_number_weight == amrex::Real(0.0) &&
               plan.packet_particle_mass == population.pivots[pivot_bin] &&
               plan.destinations[0].bin == pivot_bin &&
               plan.destinations[0].number_weight == amrex::Real(1.0) &&
               plan.destinations[0].particle_mass == population.pivots[pivot_bin];
    }
    if (plan.destination_count != 2) return false;
    const int lower_bin = plan.destinations[0].bin;
    if (lower_bin < 0 || lower_bin + 1 >= population.nbins ||
        plan.destinations[1].bin != lower_bin + 1) return false;
    const amrex::Real lower = population.pivots[lower_bin];
    const amrex::Real upper = population.pivots[lower_bin + 1];
    if (!(plan.packet_particle_mass > lower && plan.packet_particle_mass < upper)) return false;
    const amrex::Real width = upper - lower;
    amrex::Real lower_weight = amrex::Real(0.0);
    amrex::Real upper_weight = amrex::Real(0.0);
    if (checked_quotient(upper - plan.packet_particle_mass, width, lower_weight) !=
            QuotientStatus::Ok ||
        checked_quotient(plan.packet_particle_mass - lower, width, upper_weight) !=
            QuotientStatus::Ok) return false;
    return plan.status == RemapStatus::Ok && !plan.normalized_roundoff &&
           plan.destination_count == 2 && plan.residual_number_weight == amrex::Real(0.0) &&
           plan.destinations[0].bin == lower_bin &&
           plan.destinations[0].number_weight == lower_weight &&
           plan.destinations[0].particle_mass == lower &&
           plan.destinations[1].bin == lower_bin + 1 &&
           plan.destinations[1].number_weight == upper_weight &&
           plan.destinations[1].particle_mass == upper;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool canonical_two_moment_moments (const PopulationRemapView& population,
                                   const int bin, const amrex::Real number,
                                   const amrex::Real mass) noexcept
{
    if (!valid_population_remap_context(population) ||
        population.moment_mode != MomentMode::TwoMoment || bin < 0 || bin >= population.nbins ||
        !finite_nonnegative(number) || !finite_nonnegative(mass)) return false;
    if (number == amrex::Real(0.0)) return mass == amrex::Real(0.0);
    if (mass <= amrex::Real(0.0)) return false;
    amrex::Real mean = amrex::Real(0.0);
    if (checked_quotient(mass, number, mean) != QuotientStatus::Ok) return false;
    const amrex::Real lower = population.edges[bin];
    const amrex::Real upper = population.edges[bin + 1];
    return amrex::Math::isfinite(mean) && mean >= lower &&
           (mean < upper || (bin == population.nbins - 1 && mean <= upper));
}

/** Strict persisted-state test shared by reconstruction, restart, and fixtures. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool canonical_two_moment_bin_state (const PopulationRemapView& population,
                                     const int bin, const amrex::Real* state,
                                     const int state_count) noexcept
{
    if (!valid_population_remap_context(population) || state == nullptr ||
        state_count != population.state_components || bin < 0 || bin >= population.nbins ||
        population.moment_mode != MomentMode::TwoMoment) return false;
    const amrex::Real mass = state[population.mass_offset + bin];
    const amrex::Real number = state[population.number_offset + bin];
    if (!canonical_two_moment_moments(population, bin, number, mass)) return false;
    if (number == amrex::Real(0.0)) {
        for (int property = 0; property < population.property_count; ++property) {
            if (state[population.property_component_offsets[property] + bin] != amrex::Real(0.0)) {
                return false;
            }
        }
    }
    for (int property = 0; property < population.property_count; ++property) {
        const amrex::Real amount = state[population.property_component_offsets[property] + bin];
        if (!finite_nonnegative(amount)) return false;
        if (number == amrex::Real(0.0)) continue;
        amrex::Real per_particle = amrex::Real(0.0);
        if (checked_quotient(amount, number, per_particle) != QuotientStatus::Ok) return false;
    }
    return true;
}

/** Strict persisted-state test for a one-moment bin and its attached inventories. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool canonical_one_moment_bin_state (const PopulationRemapView& population,
                                     const int bin, const amrex::Real* state,
                                     const int state_count) noexcept
{
    if (!valid_population_remap_context(population) || state == nullptr ||
        state_count != population.state_components || bin < 0 || bin >= population.nbins ||
        population.moment_mode != MomentMode::OneMoment) return false;
    const amrex::Real mass = state[population.mass_offset + bin];
    if (!finite_nonnegative(mass)) return false;
    amrex::Real number = amrex::Real(0.0);
    if (checked_quotient(mass, population.pivots[bin], number) != QuotientStatus::Ok) {
        return false;
    }
    for (int property = 0; property < population.property_count; ++property) {
        const amrex::Real amount = state[population.property_component_offsets[property] + bin];
        if (!finite_nonnegative(amount)) return false;
        if (number == amrex::Real(0.0)) {
            if (amount != amrex::Real(0.0)) return false;
        } else {
            amrex::Real per_particle = amrex::Real(0.0);
            if (checked_quotient(amount, number, per_particle) != QuotientStatus::Ok) {
                return false;
            }
        }
    }
    return true;
}

/** Shared admission rule for every persisted bin accepted by reconstruction. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool canonical_persisted_bin_state (const PopulationRemapView& population,
                                    const int bin, const amrex::Real* state,
                                    const int state_count) noexcept
{
    return population.moment_mode == MomentMode::OneMoment ?
        canonical_one_moment_bin_state(population, bin, state, state_count) :
        canonical_two_moment_bin_state(population, bin, state, state_count);
}

} // namespace remap_detail

/**
 * Classify one packet without allocation or state mutation. Two-moment routing
 * uses half-open intervals, includes the global top edge, and admits only four
 * representable steps of outer-boundary roundoff. The plan keeps the original
 * packet mass while its destination records the exact boundary mass. The
 * immutable metadata must have been fully qualified once before this call;
 * the per-packet path performs only constant-time context checks and a
 * logarithmic interval search.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PacketRoutingPlan plan_packet_routing (const PopulationRemapView population,
                                       const amrex::Real packet_number,
                                       const amrex::Real particle_mass) noexcept
{
    PacketRoutingPlan plan;
    plan.population = population;
    plan.packet_number = packet_number;
    plan.packet_particle_mass = particle_mass;
    if (!valid_population_remap_context(population) ||
        !amrex::Math::isfinite(packet_number) || !amrex::Math::isfinite(particle_mass) ||
        packet_number < amrex::Real(0.0) || particle_mass < amrex::Real(0.0)) {
        return plan;
    }
    if (packet_number == amrex::Real(0.0)) {
        plan.status = RemapStatus::Ok;
        return plan;
    }

    if (population.moment_mode == MomentMode::TwoMoment) {
        // A dry positive-number packet is residual material, not a populated
        // zero-water liquid interval, even when the first edge is exactly zero.
        if (particle_mass == amrex::Real(0.0)) {
            plan.status = RemapStatus::ZeroWaterResidual;
            plan.residual_number_weight = amrex::Real(1.0);
            return plan;
        }
        if (particle_mass < population.edges[0]) {
            if (!remap_detail::within_boundary_roundoff(
                    particle_mass, population.edges[0], true)) {
                plan.status = RemapStatus::BelowSupportedGrid;
                return plan;
            }
            plan.status = RemapStatus::Ok;
            plan.destination_count = 1;
            plan.destinations[0] = {0, amrex::Real(1.0), population.edges[0]};
            plan.normalized_roundoff = true;
            return plan;
        }
        if (particle_mass > population.edges[population.nbins]) {
            if (!remap_detail::within_boundary_roundoff(
                    particle_mass, population.edges[population.nbins], false)) {
                plan.status = RemapStatus::Overflow;
                return plan;
            }
            plan.status = RemapStatus::Ok;
            plan.destination_count = 1;
            plan.destinations[0] = {population.nbins - 1, amrex::Real(1.0),
                                    population.edges[population.nbins]};
            plan.normalized_roundoff = true;
            return plan;
        }
        const int bin = particle_mass == population.edges[population.nbins] ?
            population.nbins - 1 :
            amrex::bisect(population.edges, 0, population.nbins, particle_mass);
        plan.status = RemapStatus::Ok;
        plan.destination_count = 1;
        plan.destinations[0] = {bin, amrex::Real(1.0), particle_mass};
        return plan;
    }

    if (particle_mass == amrex::Real(0.0)) {
        plan.status = RemapStatus::ZeroWaterResidual;
        plan.residual_number_weight = amrex::Real(1.0);
        return plan;
    }
    if (particle_mass < population.pivots[0]) {
        amrex::Real wet_weight = amrex::Real(0.0);
        const auto quotient_status = remap_detail::checked_quotient(
            particle_mass, population.pivots[0], wet_weight);
        if (quotient_status != remap_detail::QuotientStatus::Ok) {
            plan.status = remap_detail::remap_quotient_status(quotient_status);
            return plan;
        }
        plan.status = RemapStatus::ZeroWaterResidual;
        plan.destination_count = 1;
        plan.destinations[0] = {0, wet_weight, population.pivots[0]};
        plan.residual_number_weight = amrex::Real(1.0) - wet_weight;
        return plan;
    }
    if (particle_mass > population.pivots[population.nbins - 1]) {
        plan.status = RemapStatus::Overflow;
        return plan;
    }
    const int lower_bin = particle_mass == population.pivots[population.nbins - 1] ?
        population.nbins - 1 :
        amrex::bisect(population.pivots, 0, population.nbins - 1, particle_mass);
    if (particle_mass == population.pivots[lower_bin]) {
        plan.status = RemapStatus::Ok;
        plan.destination_count = 1;
        plan.destinations[0] = {lower_bin, amrex::Real(1.0), population.pivots[lower_bin]};
        return plan;
    }
    if (lower_bin + 1 < population.nbins) {
        const amrex::Real lower = population.pivots[lower_bin];
        const amrex::Real upper = population.pivots[lower_bin + 1];
        const amrex::Real width = upper - lower;
        amrex::Real lower_weight = amrex::Real(0.0);
        amrex::Real upper_weight = amrex::Real(0.0);
        const auto lower_status = remap_detail::checked_quotient(
            upper - particle_mass, width, lower_weight);
        const auto upper_status = remap_detail::checked_quotient(
            particle_mass - lower, width, upper_weight);
        if (lower_status != remap_detail::QuotientStatus::Ok ||
            upper_status != remap_detail::QuotientStatus::Ok) {
            const auto quotient_status =
                lower_status != remap_detail::QuotientStatus::Ok ? lower_status : upper_status;
            plan.status = remap_detail::remap_quotient_status(quotient_status);
            return plan;
        }
        plan.status = RemapStatus::Ok;
        plan.destination_count = 2;
        plan.destinations[0] = {lower_bin, lower_weight, lower};
        plan.destinations[1] = {lower_bin + 1, upper_weight, upper};
        return plan;
    }
    return plan;
}

struct PacketApplicationCoreResult {
    RemapStatus status{RemapStatus::Invalid};
    amrex::Real residual_number{0.0};
    amrex::Real residual_water_mass{0.0};
    amrex::Real roundoff_water_mass_correction{0.0};
    bool normalized_roundoff{false};
};

/**
 * Allocation-free packet application. The first pass validates the complete
 * plan, residual outputs, and every candidate increment. The second pass
 * commits only after all failure paths have been ruled out. Input/output and
 * property/state storage supplied by the caller must not overlap.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PacketApplicationCoreResult apply_packet_routing_core (
    const PopulationRemapView population, const PacketRoutingPlan plan,
    const amrex::Real* property_packet_amounts, const int property_count,
    amrex::Real* residual_properties, const int residual_capacity,
    amrex::Real* candidate_state, const int state_count) noexcept
{
    PacketApplicationCoreResult result;
    result.status = plan.status;
    if (plan.status != RemapStatus::Ok && plan.status != RemapStatus::ZeroWaterResidual) {
        return result;
    }
    if (!valid_population_remap_context(population) ||
        !same_population_context(population, plan.population) ||
        !amrex::Math::isfinite(plan.packet_number) ||
        !amrex::Math::isfinite(plan.packet_particle_mass) ||
        plan.packet_number < amrex::Real(0.0) || plan.packet_particle_mass < amrex::Real(0.0) ||
        property_count != population.property_count ||
        (property_count > 0 && (property_packet_amounts == nullptr || residual_properties == nullptr)) ||
        residual_capacity < property_count || candidate_state == nullptr ||
        state_count != population.state_components || plan.destination_count < 0 ||
        plan.destination_count > 2 || !amrex::Math::isfinite(plan.residual_number_weight) ||
        plan.residual_number_weight < amrex::Real(0.0) ||
        plan.residual_number_weight > amrex::Real(1.0)) {
        result.status = RemapStatus::Invalid;
        return result;
    }
    if (!remap_detail::route_plan_is_consistent(population, plan)) {
        result.status = RemapStatus::Invalid;
        return result;
    }
    for (int property = 0; property < property_count; ++property) {
        if (!remap_detail::finite_nonnegative(property_packet_amounts[property])) {
            result.status = RemapStatus::Invalid;
            return result;
        }
    }
    if (plan.packet_number == amrex::Real(0.0)) {
        for (int property = 0; property < property_count; ++property) {
            if (property_packet_amounts[property] != amrex::Real(0.0)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
        }
    }
    amrex::Real weight_sum = plan.residual_number_weight;
    for (int destination = 0; destination < plan.destination_count; ++destination) {
        const auto& route = plan.destinations[destination];
        if (route.bin < 0 || route.bin >= population.nbins ||
            !amrex::Math::isfinite(route.number_weight) || route.number_weight < amrex::Real(0.0) ||
            route.number_weight > amrex::Real(1.0) ||
            !remap_detail::finite_nonnegative(route.particle_mass)) {
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
        if (!amrex::Math::isfinite(weight_sum) ||
            std::abs(weight_sum - amrex::Real(1.0)) > weight_tolerance ||
            ((plan.status == RemapStatus::ZeroWaterResidual) !=
             (plan.residual_number_weight > amrex::Real(0.0)))) {
            result.status = RemapStatus::Invalid;
            return result;
        }
    } else if (plan.destination_count != 0 ||
               plan.residual_number_weight != amrex::Real(0.0)) {
        result.status = RemapStatus::Invalid;
        return result;
    }

    amrex::Real residual_number = amrex::Real(0.0);
    auto product_status = remap_detail::checked_product(
        plan.packet_number, plan.residual_number_weight, residual_number);
    if (product_status != remap_detail::ProductStatus::Ok) {
        result.status = remap_detail::remap_product_status(product_status);
        return result;
    }
    amrex::Real roundoff_water_mass_correction = amrex::Real(0.0);
    // Validate all outputs and candidate increments before changing either.
    for (int property = 0; property < property_count; ++property) {
        amrex::Real residual = amrex::Real(0.0);
        product_status = remap_detail::checked_product(
            property_packet_amounts[property], plan.residual_number_weight, residual);
        if (product_status != remap_detail::ProductStatus::Ok) {
            result.status = remap_detail::remap_product_status(product_status);
            return result;
        }
    }
    for (int destination = 0; destination < plan.destination_count; ++destination) {
        const auto& route = plan.destinations[destination];
        amrex::Real number_increment = amrex::Real(0.0);
        amrex::Real mass_increment = amrex::Real(0.0);
        const amrex::Real deposited_mass = route.particle_mass;
        product_status = remap_detail::checked_product(
            plan.packet_number, route.number_weight, number_increment);
        if (product_status == remap_detail::ProductStatus::Ok) {
            product_status = remap_detail::checked_product(
                number_increment, deposited_mass, mass_increment);
        }
        if (product_status != remap_detail::ProductStatus::Ok) {
            result.status = remap_detail::remap_product_status(product_status);
            return result;
        }
        if (population.moment_mode == MomentMode::TwoMoment && plan.normalized_roundoff) {
            const amrex::Real mass_delta =
                route.particle_mass - plan.packet_particle_mass;
            if (!amrex::Math::isfinite(mass_delta) ||
                mass_delta == amrex::Real(0.0)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            const amrex::Real mass_delta_magnitude = mass_delta < amrex::Real(0.0) ?
                -mass_delta : mass_delta;
            amrex::Real correction_magnitude = amrex::Real(0.0);
            product_status = remap_detail::checked_product(
                plan.packet_number, mass_delta_magnitude, correction_magnitude);
            if (product_status != remap_detail::ProductStatus::Ok) {
                result.status = remap_detail::remap_product_status(product_status);
                return result;
            }
            const amrex::Real correction = mass_delta < amrex::Real(0.0) ?
                -correction_magnitude : correction_magnitude;
            if (!amrex::Math::isfinite(correction) ||
                !amrex::Math::isfinite(roundoff_water_mass_correction + correction)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            roundoff_water_mass_correction += correction;
        }
        const int mass_component = population.mass_offset + route.bin;
        const amrex::Real old_mass = candidate_state[mass_component];
        if (!remap_detail::finite_nonnegative(old_mass) ||
            !remap_detail::finite_nonnegative(old_mass + mass_increment)) {
            result.status = RemapStatus::Invalid;
            return result;
        }
        const amrex::Real next_mass = old_mass + mass_increment;
        amrex::Real next_number = amrex::Real(0.0);
        if (population.number_offset >= 0) {
            const int number_component = population.number_offset + route.bin;
            const amrex::Real old_number = candidate_state[number_component];
            if (!remap_detail::finite_nonnegative(old_number)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            next_number = old_number + number_increment;
            if (!remap_detail::finite_nonnegative(next_number)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            if (next_number == amrex::Real(0.0)) {
                if (next_mass != amrex::Real(0.0)) {
                    result.status = RemapStatus::Invalid;
                    return result;
                }
            } else {
                amrex::Real mean = amrex::Real(0.0);
                const auto quotient_status = remap_detail::checked_quotient(
                    next_mass, next_number, mean);
                if (quotient_status != remap_detail::QuotientStatus::Ok) {
                    result.status = remap_detail::remap_quotient_status(quotient_status);
                    return result;
                }
                if (!remap_detail::canonical_two_moment_moments(
                        population, route.bin, next_number, next_mass)) {
                    result.status = RemapStatus::Invalid;
                    return result;
                }
            }
        } else {
            const auto quotient_status = remap_detail::checked_quotient(
                next_mass, population.pivots[route.bin], next_number);
            if (quotient_status != remap_detail::QuotientStatus::Ok) {
                result.status = remap_detail::remap_quotient_status(quotient_status);
                return result;
            }
        }
        for (int property = 0; property < property_count; ++property) {
            amrex::Real increment = amrex::Real(0.0);
            const int component = population.property_component_offsets[property] + route.bin;
            const amrex::Real old_value = candidate_state[component];
            product_status = remap_detail::checked_product(
                property_packet_amounts[property], route.number_weight, increment);
            if (product_status != remap_detail::ProductStatus::Ok) {
                result.status = remap_detail::remap_product_status(product_status);
                return result;
            }
            if (!remap_detail::finite_nonnegative(old_value)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            const amrex::Real next_value = old_value + increment;
            if (!remap_detail::finite_nonnegative(next_value)) {
                result.status = RemapStatus::Invalid;
                return result;
            }
            if (next_number == amrex::Real(0.0)) {
                if (next_value != amrex::Real(0.0)) {
                    result.status = RemapStatus::Invalid;
                    return result;
                }
            } else {
                amrex::Real per_particle = amrex::Real(0.0);
                const auto quotient_status = remap_detail::checked_quotient(
                    next_value, next_number, per_particle);
                if (quotient_status != remap_detail::QuotientStatus::Ok) {
                    result.status = remap_detail::remap_quotient_status(quotient_status);
                    return result;
                }
            }
        }
    }

    for (int destination = 0; destination < plan.destination_count; ++destination) {
        const auto& route = plan.destinations[destination];
        const amrex::Real number_increment = plan.packet_number * route.number_weight;
        const amrex::Real deposited_mass = route.particle_mass;
        candidate_state[population.mass_offset + route.bin] += number_increment * deposited_mass;
        if (population.number_offset >= 0) {
            candidate_state[population.number_offset + route.bin] += number_increment;
        }
        for (int property = 0; property < property_count; ++property) {
            candidate_state[population.property_component_offsets[property] + route.bin] +=
                property_packet_amounts[property] * route.number_weight;
        }
    }
    for (int property = 0; property < property_count; ++property) {
        residual_properties[property] =
            property_packet_amounts[property] * plan.residual_number_weight;
    }
    result.residual_number = residual_number;
    result.residual_water_mass = amrex::Real(0.0);
    result.roundoff_water_mass_correction = roundoff_water_mass_correction;
    result.normalized_roundoff = plan.normalized_roundoff;
    return result;
}

/** Non-owning result from allocation-free one-bin reconstruction. */
struct ReconstructionDeltaView {
    PopulationRemapView population{};
    int bin{-1};
    bool empty{true};
    bool normalized_roundoff{false};
    amrex::Real number{0.0};
    amrex::Real particle_mass{0.0};
    const amrex::Real* property_per_particle{nullptr};
    int property_count{0};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ReconstructionStatus reconstruct_bin_core (
    const PopulationRemapView population, const int bin,
    const amrex::Real* persisted_state, const int state_count,
    amrex::Real* property_per_particle, const int property_capacity,
    ReconstructionDeltaView& reconstruction) noexcept
{
    if (!valid_population_remap_context(population) || bin < 0 || bin >= population.nbins ||
        persisted_state == nullptr || state_count != population.state_components ||
        property_capacity < population.property_count ||
        (population.property_count > 0 && property_per_particle == nullptr)) {
        return ReconstructionStatus::Invalid;
    }
    if (!remap_detail::canonical_persisted_bin_state(
            population, bin, persisted_state, state_count)) {
        return ReconstructionStatus::Invalid;
    }
    const amrex::Real mass = persisted_state[population.mass_offset + bin];
    if (!remap_detail::finite_nonnegative(mass)) return ReconstructionStatus::Invalid;
    ReconstructionDeltaView next;
    next.population = population;
    next.bin = bin;
    next.property_per_particle = property_per_particle;
    next.property_count = population.property_count;
    if (population.moment_mode == MomentMode::OneMoment) {
        next.particle_mass = population.pivots[bin];
        if (remap_detail::checked_quotient(
                mass, next.particle_mass, next.number) != remap_detail::QuotientStatus::Ok) {
            return ReconstructionStatus::Invalid;
        }
    } else {
        const amrex::Real number = persisted_state[population.number_offset + bin];
        if (!remap_detail::finite_nonnegative(number)) return ReconstructionStatus::Invalid;
        if (!remap_detail::canonical_two_moment_bin_state(
                population, bin, persisted_state, state_count)) {
            return ReconstructionStatus::Invalid;
        }
        const amrex::Real lower = population.edges[bin];
        const amrex::Real upper = population.edges[bin + 1];
        EndpointTransform endpoints;
        if (!try_two_moment_to_endpoints(number, mass, lower, upper, endpoints)) {
            return ReconstructionStatus::Invalid;
        }
        next.number = number;
        // Endpoint tolerance is only a realizability aid for the transform.
        // Persisted state has already passed strict canonical ownership and is
        // never normalized by no-process reconstruction.
        next.normalized_roundoff = false;
        if (number > amrex::Real(0.0)) {
            if (remap_detail::checked_quotient(
                    mass, number, next.particle_mass) != remap_detail::QuotientStatus::Ok) {
                return ReconstructionStatus::Invalid;
            }
        }
    }
    next.empty = next.number == amrex::Real(0.0);
    if (next.empty) next.particle_mass = amrex::Real(0.0);
    for (int property = 0; property < population.property_count; ++property) {
        const amrex::Real amount = persisted_state[
            population.property_component_offsets[property] + bin];
        if (!remap_detail::finite_nonnegative(amount)) return ReconstructionStatus::Invalid;
        if (next.empty) {
            if (amount != amrex::Real(0.0)) return ReconstructionStatus::Invalid;
        } else {
            amrex::Real per_particle = amrex::Real(0.0);
            if (remap_detail::checked_quotient(
                    amount, next.number, per_particle) != remap_detail::QuotientStatus::Ok) {
                return ReconstructionStatus::Invalid;
            }
        }
    }
    for (int property = 0; property < population.property_count; ++property) {
        const amrex::Real amount = persisted_state[
            population.property_component_offsets[property] + bin];
        if (next.empty) {
            property_per_particle[property] = amrex::Real(0.0);
        } else if (remap_detail::checked_quotient(
                       amount, next.number, property_per_particle[property]) !=
                   remap_detail::QuotientStatus::Ok) {
            return ReconstructionStatus::Invalid;
        }
    }
    reconstruction = next;
    return next.empty ? ReconstructionStatus::Empty : ReconstructionStatus::Populated;
}

struct IntegratedMomentsCoreResult {
    amrex::Real number{0.0};
    amrex::Real water_mass{0.0};
};

/** Integrate one reconstructed delta over [lower, upper), optionally including upper. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool integrate_interval_core (const ReconstructionDeltaView reconstruction,
                              const amrex::Real lower, const amrex::Real upper,
                              const bool include_upper, amrex::Real* attached_properties,
                              const int property_capacity,
                              IntegratedMomentsCoreResult& integral) noexcept
{
    if (!amrex::Math::isfinite(lower) || !amrex::Math::isfinite(upper) || !(lower < upper) ||
        !remap_detail::finite_nonnegative(reconstruction.number) ||
        !remap_detail::finite_nonnegative(reconstruction.particle_mass) ||
        reconstruction.property_count < 0 || property_capacity < reconstruction.property_count ||
        (reconstruction.property_count > 0 &&
         (reconstruction.property_per_particle == nullptr || attached_properties == nullptr)) ||
        (!reconstruction.empty && reconstruction.number <= amrex::Real(0.0)) ||
        (reconstruction.empty && (reconstruction.number != amrex::Real(0.0) ||
                                  reconstruction.particle_mass != amrex::Real(0.0)))) {
        return false;
    }
    for (int property = 0; property < reconstruction.property_count; ++property) {
        const amrex::Real value = reconstruction.property_per_particle[property];
        if (!remap_detail::finite_nonnegative(value) ||
            (reconstruction.empty && value != amrex::Real(0.0))) {
            return false;
        }
    }
    const bool contains = !reconstruction.empty && reconstruction.particle_mass >= lower &&
        (reconstruction.particle_mass < upper ||
         (include_upper && reconstruction.particle_mass == upper));
    IntegratedMomentsCoreResult next;
    if (contains) {
        next.number = reconstruction.number;
        if (!remap_detail::finite_product(reconstruction.number,
                                           reconstruction.particle_mass, next.water_mass)) {
            return false;
        }
        for (int property = 0; property < reconstruction.property_count; ++property) {
            amrex::Real amount = amrex::Real(0.0);
            if (!remap_detail::finite_product(reconstruction.number,
                    reconstruction.property_per_particle[property], amount)) {
                return false;
            }
        }
    }
    for (int property = 0; property < reconstruction.property_count; ++property) {
        amrex::Real amount = amrex::Real(0.0);
        if (contains && !remap_detail::finite_product(reconstruction.number,
                reconstruction.property_per_particle[property], amount)) {
            return false;
        }
        attached_properties[property] = amount;
    }
    integral = next;
    return true;
}

/** Project one reconstructed bin into caller-owned state after full validation. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool project_reconstruction_bin_core (const PopulationRemapView population,
                                      const ReconstructionDeltaView reconstruction,
                                      amrex::Real* candidate_state,
                                      const int state_count) noexcept
{
    if (!valid_population_remap_context(population) ||
        !same_population_context(population, reconstruction.population) ||
        reconstruction.bin < 0 || reconstruction.bin >= population.nbins ||
        candidate_state == nullptr || state_count != population.state_components ||
        reconstruction.property_count != population.property_count ||
        (reconstruction.property_count > 0 && reconstruction.property_per_particle == nullptr) ||
        !remap_detail::finite_nonnegative(reconstruction.number) ||
        !remap_detail::finite_nonnegative(reconstruction.particle_mass) ||
        (!reconstruction.empty && reconstruction.number <= amrex::Real(0.0)) ||
        (reconstruction.empty && (reconstruction.number != amrex::Real(0.0) ||
                                  reconstruction.particle_mass != amrex::Real(0.0)))) {
        return false;
    }
    if (!reconstruction.empty) {
        if (population.moment_mode == MomentMode::OneMoment) {
            if (reconstruction.particle_mass != population.pivots[reconstruction.bin]) return false;
        }
    }
    amrex::Real mass = amrex::Real(0.0);
    if (!remap_detail::finite_product(reconstruction.number,
                                      reconstruction.particle_mass, mass)) {
        return false;
    }
    amrex::Real projected_number = reconstruction.number;
    if (population.moment_mode == MomentMode::OneMoment) {
        if (remap_detail::checked_quotient(
                mass, population.pivots[reconstruction.bin], projected_number) !=
            remap_detail::QuotientStatus::Ok) return false;
    } else if (!remap_detail::canonical_two_moment_moments(
                   population, reconstruction.bin, projected_number, mass)) {
        return false;
    }
    for (int property = 0; property < population.property_count; ++property) {
        const amrex::Real value = reconstruction.property_per_particle[property];
        amrex::Real amount = amrex::Real(0.0);
        if (!remap_detail::finite_nonnegative(value) ||
            (reconstruction.empty && value != amrex::Real(0.0)) ||
            !remap_detail::finite_product(reconstruction.number, value, amount)) {
            return false;
        }
        if (projected_number == amrex::Real(0.0)) {
            if (amount != amrex::Real(0.0)) return false;
        } else {
            amrex::Real per_particle = amrex::Real(0.0);
            if (remap_detail::checked_quotient(
                    amount, projected_number, per_particle) !=
                remap_detail::QuotientStatus::Ok) return false;
        }
    }
    const int bin = reconstruction.bin;
    candidate_state[population.mass_offset + bin] = mass;
    if (population.number_offset >= 0) {
        candidate_state[population.number_offset + bin] = reconstruction.number;
    }
    for (int property = 0; property < population.property_count; ++property) {
        candidate_state[population.property_component_offsets[property] + bin] =
            reconstruction.number * reconstruction.property_per_particle[property];
    }
    return true;
}

// Host-owning reference adapters. These allocate diagnostic vectors as
// convenient for tests, then delegate every numerical operation to the POD
// routines above.
struct PacketApplicationResult {
    RemapStatus status{RemapStatus::Invalid};
    amrex::Real residual_number{0.0};
    amrex::Real residual_water_mass{0.0};
    amrex::Real roundoff_water_mass_correction{0.0};
    bool normalized_roundoff{false};
    std::vector<amrex::Real> residual_properties;
};

[[nodiscard]] PacketApplicationResult
apply_packet_routing (const SBMLayout& layout, const PacketRoutingPlan& plan,
                      const std::vector<amrex::Real>& property_packet_amounts,
                      std::vector<amrex::Real>& candidate_state);

struct ReconstructionDelta {
    int population_id{-1};
    int bin{-1};
    bool empty{true};
    bool normalized_roundoff{false};
    amrex::Real number{0.0};
    amrex::Real particle_mass{0.0};
    std::vector<amrex::Real> property_per_particle;
};

[[nodiscard]] ReconstructionStatus
reconstruct_bin (const SBMLayout& layout, int population_id, int bin,
                 const std::vector<amrex::Real>& persisted_state,
                 ReconstructionDelta& reconstruction);

struct IntegratedMoments {
    amrex::Real number{0.0};
    amrex::Real water_mass{0.0};
    std::vector<amrex::Real> attached_properties;
};

[[nodiscard]] bool
integrate_interval (const ReconstructionDelta& reconstruction,
                    amrex::Real lower, amrex::Real upper, bool include_upper,
                    IntegratedMoments& integral);

[[nodiscard]] bool
project_reconstruction_bin (const SBMLayout& layout,
                            const ReconstructionDelta& reconstruction,
                            std::vector<amrex::Real>& candidate_state);

} // namespace erf_sbm

#endif
