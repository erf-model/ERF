#ifndef ERF_SBM_REMAPPING_H_
#define ERF_SBM_REMAPPING_H_

#include "ERF_SBMLayout.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <vector>

namespace erf_sbm {

enum class RemapStatus { Ok, ZeroWaterResidual, BelowSupportedGrid, Overflow, Invalid };
enum class ReconstructionStatus { Populated, Empty, Invalid };

/** Non-owning grid view suitable for passing to host/device routing code. */
struct SpectralGridView {
    const amrex::Real* edges{nullptr};
    const amrex::Real* pivots{nullptr};
    int nbins{0};
};

struct RouteDestination {
    int bin{-1};
    amrex::Real number_weight{0.0};
    amrex::Real particle_mass{0.0};
};

/** Pure fixed-capacity routing plan; a 1M packet has at most two destinations. */
struct PacketRoutingPlan {
    RemapStatus status{RemapStatus::Invalid};
    amrex::Real packet_number{0.0};
    amrex::Real packet_particle_mass{0.0};
    int destination_count{0};
    RouteDestination destinations[2]{};
    amrex::Real residual_number_weight{0.0};
};

/**
 * Classify a packet without changing scientific state or allocating memory.
 * Interval2M uses [edge_j, edge_(j+1)); the global upper edge is included.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PacketRoutingPlan plan_packet_routing(const SpectralGridView grid,
                                      const MomentMode mode,
                                      const amrex::Real packet_number,
                                      const amrex::Real particle_mass) noexcept
{
    PacketRoutingPlan plan;
    plan.packet_number = packet_number;
    plan.packet_particle_mass = particle_mass;
    if (grid.nbins <= 0 || grid.edges == nullptr ||
        !std::isfinite(packet_number) || !std::isfinite(particle_mass) ||
        packet_number < amrex::Real(0.0) || particle_mass < amrex::Real(0.0) ||
        (mode != MomentMode::OneMoment && mode != MomentMode::TwoMoment)) {
        return plan;
    }
    for (int edge = 0; edge <= grid.nbins; ++edge) {
        if (!std::isfinite(grid.edges[edge]) || grid.edges[edge] < amrex::Real(0.0) ||
            (edge > 0 && !(grid.edges[edge - 1] < grid.edges[edge]))) {
            return plan;
        }
    }
    if (mode == MomentMode::OneMoment) {
        if (grid.pivots == nullptr) return plan;
        for (int bin = 0; bin < grid.nbins; ++bin) {
            if (!std::isfinite(grid.pivots[bin]) || !(grid.edges[bin] <= grid.pivots[bin]) ||
                !(grid.pivots[bin] <= grid.edges[bin + 1]) ||
                (bin > 0 && !(grid.pivots[bin - 1] < grid.pivots[bin]))) {
                return plan;
            }
        }
    }

    if (packet_number == amrex::Real(0.0)) {
        plan.status = RemapStatus::Ok;
        return plan;
    }

    if (mode == MomentMode::TwoMoment) {
        if (particle_mass < grid.edges[0]) {
            plan.status = RemapStatus::BelowSupportedGrid;
            return plan;
        }
        for (int bin = 0; bin < grid.nbins; ++bin) {
            const bool in_interval = particle_mass >= grid.edges[bin] &&
                (particle_mass < grid.edges[bin + 1] ||
                 (bin == grid.nbins - 1 && particle_mass <= grid.edges[bin + 1]));
            if (in_interval) {
                plan.status = RemapStatus::Ok;
                plan.destination_count = 1;
                plan.destinations[0] = {bin, amrex::Real(1.0), particle_mass};
                return plan;
            }
        }
        plan.status = RemapStatus::Overflow;
        return plan;
    }

    if (particle_mass == amrex::Real(0.0)) {
        plan.status = RemapStatus::ZeroWaterResidual;
        plan.residual_number_weight = amrex::Real(1.0);
        return plan;
    }
    if (particle_mass < grid.pivots[0]) {
        const amrex::Real wet_weight = particle_mass / grid.pivots[0];
        plan.status = RemapStatus::ZeroWaterResidual;
        plan.destination_count = 1;
        plan.destinations[0] = {0, wet_weight, grid.pivots[0]};
        plan.residual_number_weight = amrex::Real(1.0) - wet_weight;
        return plan;
    }
    if (particle_mass > grid.pivots[grid.nbins - 1]) {
        plan.status = RemapStatus::Overflow;
        return plan;
    }
    for (int bin = 0; bin < grid.nbins; ++bin) {
        const amrex::Real lower = grid.pivots[bin];
        if (particle_mass == lower) {
            plan.status = RemapStatus::Ok;
            plan.destination_count = 1;
            plan.destinations[0] = {bin, amrex::Real(1.0), lower};
            return plan;
        }
        if (bin + 1 < grid.nbins) {
            const amrex::Real upper = grid.pivots[bin + 1];
            if (particle_mass == upper) {
                plan.status = RemapStatus::Ok;
                plan.destination_count = 1;
                plan.destinations[0] = {bin + 1, amrex::Real(1.0), upper};
                return plan;
            }
            if (particle_mass > lower && particle_mass < upper) {
                const amrex::Real width = upper - lower;
                plan.status = RemapStatus::Ok;
                plan.destination_count = 2;
                plan.destinations[0] = {bin, (upper - particle_mass) / width, lower};
                plan.destinations[1] = {bin + 1, (particle_mass - lower) / width, upper};
                return plan;
            }
        }
    }
    return plan;
}

struct PacketApplicationResult {
    RemapStatus status{RemapStatus::Invalid};
    amrex::Real residual_number{0.0};
    amrex::Real residual_water_mass{0.0};
    // In the order of attached properties carried by the selected population.
    std::vector<amrex::Real> residual_properties;
};

/**
 * Apply a valid plan transactionally to caller-owned candidate state. Property
 * packet amounts must follow layout declaration order for this population;
 * every attached property receives the same particle-number weights.
 */
[[nodiscard]] PacketApplicationResult
apply_packet_routing(const SBMLayout& layout, int population_id,
                     const PacketRoutingPlan& plan,
                     const std::vector<amrex::Real>& property_packet_amounts,
                     std::vector<amrex::Real>& candidate_state);

/** A temporary one-node delta view of one persisted spectral bin. */
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
reconstruct_bin(const SBMLayout& layout, int population_id, int bin,
                const std::vector<amrex::Real>& persisted_state,
                ReconstructionDelta& reconstruction);

struct IntegratedMoments {
    amrex::Real number{0.0};
    amrex::Real water_mass{0.0};
    std::vector<amrex::Real> attached_properties;
};

/** Integrate a delta over [lower, upper), optionally including upper. */
[[nodiscard]] bool
integrate_interval(const ReconstructionDelta& reconstruction,
                   amrex::Real lower, amrex::Real upper, bool include_upper,
                   IntegratedMoments& integral);

/** Project one reconstructed bin into a candidate persisted state atomically. */
[[nodiscard]] bool
project_reconstruction_bin(const SBMLayout& layout,
                           const ReconstructionDelta& reconstruction,
                           std::vector<amrex::Real>& candidate_state);

} // namespace erf_sbm

#endif
