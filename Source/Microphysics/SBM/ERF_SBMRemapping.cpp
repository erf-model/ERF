#include "ERF_SBMRemapping.H"

#include <algorithm>
#include <utility>

namespace erf_sbm {

PopulationRemapView population_remap_view (const SBMLayout& layout,
                                         const int population_id) noexcept
{
    const auto& populations = layout.populations();
    const auto found = std::find_if(populations.begin(), populations.end(),
        [population_id](const PopulationLayout& population) {
            return population.population_id == population_id;
        });
    if (found == populations.end()) return {};

    PopulationRemapView view;
    view.population_id = found->population_id;
    view.moment_mode = found->moment_mode;
    view.coordinate_kind = found->grid.coordinate_kind();
    view.edges = found->grid.edges().data();
    view.pivots = found->grid.pivots().data();
    view.nbins = found->grid.nbins();
    view.mass_offset = found->mass_offset;
    view.number_offset = found->number_offset;
    view.state_components = layout.ncomp();
    view.property_component_offsets = found->property_component_offsets.data();
    view.property_count = static_cast<int>(found->property_component_offsets.size());
    // SBMLayout validates these immutable arrays and constructs the component
    // offsets. This adapter is used on per-packet paths, so do not repeat the
    // O(nbins + properties^2) full qualification here. Keep the full validator
    // for arbitrary raw views; this layout-derived view needs only a constant-
    // time context check before entering the device/core routines.
    if (!valid_population_remap_context(view)) return {};
    return view;
}

PacketApplicationResult
apply_packet_routing (const SBMLayout& layout, const PacketRoutingPlan& plan,
                      const std::vector<amrex::Real>& property_packet_amounts,
                      std::vector<amrex::Real>& candidate_state)
{
    PacketApplicationResult result;
    result.status = plan.status;
    if (plan.status != RemapStatus::Ok && plan.status != RemapStatus::ZeroWaterResidual) {
        return result;
    }
    if (&property_packet_amounts == &candidate_state) {
        result.status = RemapStatus::Invalid;
        return result;
    }
    const auto population = population_remap_view(layout, plan.population.population_id);
    if (population.property_count >= 0) {
        result.residual_properties.resize(static_cast<std::size_t>(population.property_count),
                                          amrex::Real(0.0));
    }
    const auto core = apply_packet_routing_core(
        population, plan,
        property_packet_amounts.data(), static_cast<int>(property_packet_amounts.size()),
        result.residual_properties.data(), static_cast<int>(result.residual_properties.size()),
        candidate_state.data(), static_cast<int>(candidate_state.size()));
    result.status = core.status;
    result.residual_number = core.residual_number;
    result.residual_water_mass = core.residual_water_mass;
    result.roundoff_water_mass_correction = core.roundoff_water_mass_correction;
    result.normalized_roundoff = core.normalized_roundoff;
    return result;
}

ReconstructionStatus
reconstruct_bin (const SBMLayout& layout, const int population_id, const int bin,
                 const std::vector<amrex::Real>& persisted_state,
                 ReconstructionDelta& reconstruction)
{
    const auto population = population_remap_view(layout, population_id);
    std::vector<amrex::Real> properties(
        population.property_count >= 0 ? static_cast<std::size_t>(population.property_count) : 0u,
        amrex::Real(0.0));
    ReconstructionDeltaView view;
    const auto status = reconstruct_bin_core(
        population, bin, persisted_state.data(), static_cast<int>(persisted_state.size()),
        properties.data(), static_cast<int>(properties.size()), view);
    if (status == ReconstructionStatus::Invalid) return status;

    ReconstructionDelta next;
    next.population_id = population_id;
    next.bin = bin;
    next.empty = view.empty;
    next.normalized_roundoff = view.normalized_roundoff;
    next.number = view.number;
    next.particle_mass = view.particle_mass;
    next.property_per_particle = std::move(properties);
    reconstruction = std::move(next);
    return status;
}

bool integrate_interval (const ReconstructionDelta& reconstruction,
                         const amrex::Real lower, const amrex::Real upper,
                         const bool include_upper, IntegratedMoments& integral)
{
    std::vector<amrex::Real> properties(reconstruction.property_per_particle.size(),
                                        amrex::Real(0.0));
    ReconstructionDeltaView view;
    view.bin = reconstruction.bin;
    view.empty = reconstruction.empty;
    view.normalized_roundoff = reconstruction.normalized_roundoff;
    view.number = reconstruction.number;
    view.particle_mass = reconstruction.particle_mass;
    view.property_per_particle = reconstruction.property_per_particle.data();
    view.property_count = static_cast<int>(reconstruction.property_per_particle.size());
    IntegratedMomentsCoreResult core;
    if (!integrate_interval_core(view, lower, upper, include_upper,
                                 properties.data(), static_cast<int>(properties.size()), core)) {
        return false;
    }
    IntegratedMoments next;
    next.number = core.number;
    next.water_mass = core.water_mass;
    next.attached_properties = std::move(properties);
    integral = std::move(next);
    return true;
}

bool project_reconstruction_bin (const SBMLayout& layout,
                                 const ReconstructionDelta& reconstruction,
                                 std::vector<amrex::Real>& candidate_state)
{
    const auto population = population_remap_view(layout, reconstruction.population_id);
    ReconstructionDeltaView view;
    view.population = population;
    view.bin = reconstruction.bin;
    view.empty = reconstruction.empty;
    view.normalized_roundoff = reconstruction.normalized_roundoff;
    view.number = reconstruction.number;
    view.particle_mass = reconstruction.particle_mass;
    view.property_per_particle = reconstruction.property_per_particle.data();
    view.property_count = static_cast<int>(reconstruction.property_per_particle.size());
    return project_reconstruction_bin_core(population, view, candidate_state.data(),
                                            static_cast<int>(candidate_state.size()));
}

} // namespace erf_sbm
