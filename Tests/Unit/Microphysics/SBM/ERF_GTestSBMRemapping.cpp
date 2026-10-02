#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuLaunch.H>

#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMRemapping.H"
#include "ERF_SBMRepresentation.H"
#include "ERF_SBMRestart.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using amrex::Real;

erf_sbm::SpectralGridSpec make_grid(std::vector<Real> edges,
                                   std::vector<Real> pivots)
{
    erf_sbm::SpectralGridSpec grid;
    grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    grid.coordinate_units = "kg particle^-1";
    grid.edges = std::move(edges);
    grid.pivots = std::move(pivots);
    return grid;
}

erf_sbm::SpectralPopulationSpec make_population (
    const int id, const erf_sbm::MomentMode mode,
    std::vector<Real> edges = {Real(0.5), Real(1.5), Real(2.5)},
    std::vector<Real> pivots = {Real(1.0), Real(2.0)},
    const erf_sbm::PopulationPhase phase = erf_sbm::PopulationPhase::Liquid)
{
    erf_sbm::SpectralPopulationSpec population;
    population.population_id = id;
    population.semantic_id = "population_" + std::to_string(id);
    population.phase = phase;
    population.grid = make_grid(std::move(edges), std::move(pivots));
    population.moment_mode = mode;
    return population;
}

erf_sbm::AttachedPropertyDescriptor make_property (const std::string& name,
                                                   const int population_id)
{
    erf_sbm::AttachedPropertyDescriptor property;
    property.name = name;
    property.semantic_id = name + ".semantic";
    property.units = "kg m^-3";
    property.carrier_population = population_id;
    property.kind = erf_sbm::PropertyKind::ExtensiveMass;
    return property;
}

erf_sbm::SBMLayout make_layout (
    const erf_sbm::MomentMode mode,
    std::vector<erf_sbm::AttachedPropertyDescriptor> properties = {})
{
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(make_population(0, mode));
    spec.liquid_projection = {0, 1};
    spec.attached_properties = std::move(properties);
    return erf_sbm::SBMLayout(std::move(spec));
}

const erf_sbm::PopulationLayout& population (const erf_sbm::SBMLayout& layout,
                                           const int id)
{
    const auto found = std::find_if(layout.populations().begin(), layout.populations().end(),
        [id](const erf_sbm::PopulationLayout& value) { return value.population_id == id; });
    if (found == layout.populations().end()) throw std::runtime_error("missing test population");
    return *found;
}

void expect_close (const Real actual, const Real expected, const Real scale = Real(1.0))
{
    EXPECT_NEAR(actual, expected,
                Real(32.0) * std::numeric_limits<Real>::epsilon() *
                std::max({std::abs(actual), std::abs(expected), scale}));
}

TEST(SBMRemapping, FixedPivotOneMomentConservesNumberWaterPropertiesAndVariance)
{
    const auto layout = make_layout(erf_sbm::MomentMode::OneMoment,
        {make_property("solute_a", 0), make_property("solute_b", 0)});
    const auto& liquid = population(layout, 0);
    constexpr Real packet_number = Real(3.0);
    constexpr Real packet_mass = Real(1.4);
    const auto plan = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), packet_number, packet_mass);
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(plan.destination_count, 2);
    EXPECT_EQ(plan.destinations[0].bin, 0);
    EXPECT_EQ(plan.destinations[1].bin, 1);
    expect_close(plan.destinations[0].number_weight, Real(0.6));
    expect_close(plan.destinations[1].number_weight, Real(0.4));

    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto result = erf_sbm::apply_packet_routing(
        layout, plan, {Real(6.0), Real(15.0)}, state);
    ASSERT_EQ(result.status, erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(result.residual_number, Real(0.0));
    expect_close(state[0], Real(1.8));
    expect_close(state[1], Real(2.4));
    expect_close(state[0] / liquid.grid.pivot(0) + state[1] / liquid.grid.pivot(1),
                 packet_number);
    expect_close(state[0] + state[1], packet_number * packet_mass);

    const int solute_a = layout.property_offset(0);
    const int solute_b = layout.property_offset(1);
    expect_close(state[static_cast<std::size_t>(solute_a)], Real(3.6));
    expect_close(state[static_cast<std::size_t>(solute_a + 1)], Real(2.4));
    expect_close(state[static_cast<std::size_t>(solute_b)], Real(9.0));
    expect_close(state[static_cast<std::size_t>(solute_b + 1)], Real(6.0));
    expect_close(state[static_cast<std::size_t>(solute_a)] + state[static_cast<std::size_t>(solute_a + 1)],
                 Real(6.0));
    expect_close(state[static_cast<std::size_t>(solute_b)] + state[static_cast<std::size_t>(solute_b + 1)],
                 Real(15.0));

    const Real remapped_second_moment = state[0] * liquid.grid.pivot(0) +
                                         state[1] * liquid.grid.pivot(1);
    const Real added_variance = remapped_second_moment - packet_number * packet_mass * packet_mass;
    expect_close(added_variance,
                 packet_number * (packet_mass - liquid.grid.pivot(0)) *
                 (liquid.grid.pivot(1) - packet_mass));

    std::vector<Real> projected(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    for (int bin = 0; bin < liquid.grid.nbins(); ++bin) {
        erf_sbm::ReconstructionDelta reconstruction;
        ASSERT_EQ(erf_sbm::reconstruct_bin(layout, 0, bin, state, reconstruction),
                  erf_sbm::ReconstructionStatus::Populated);
        expect_close(reconstruction.number, state[static_cast<std::size_t>(bin)] /
                     liquid.grid.pivot(bin));
        EXPECT_DOUBLE_EQ(reconstruction.particle_mass, liquid.grid.pivot(bin));
        ASSERT_EQ(reconstruction.property_per_particle.size(), 2u);
        ASSERT_TRUE(erf_sbm::project_reconstruction_bin(layout, reconstruction, projected));
    }
    ASSERT_EQ(projected.size(), state.size());
    for (std::size_t component = 0; component < state.size(); ++component) {
        expect_close(projected[component], state[component]);
    }
}

TEST(SBMRemapping, FixedPivotExactPivotResidualOverflowAndAtomicFailure)
{
    const auto layout = make_layout(erf_sbm::MomentMode::OneMoment,
        {make_property("material_a", 0), make_property("material_b", 0)});
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));

    const auto exact = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(4.0), Real(2.0));
    ASSERT_EQ(exact.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(exact.destination_count, 1);
    EXPECT_EQ(exact.destinations[0].bin, 1);
    EXPECT_DOUBLE_EQ(exact.destinations[0].number_weight, Real(1.0));
    ASSERT_EQ(erf_sbm::apply_packet_routing(layout, exact,
        {Real(8.0), Real(12.0)}, state).status, erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(state[0], Real(0.0));
    EXPECT_DOUBLE_EQ(state[1], Real(8.0));

    const auto subpivot = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(5.0), Real(0.4));
    ASSERT_EQ(subpivot.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    ASSERT_EQ(subpivot.destination_count, 1);
    expect_close(subpivot.destinations[0].number_weight, Real(0.4));
    expect_close(subpivot.residual_number_weight, Real(0.6));
    const auto before_subpivot = state;
    const auto residual = erf_sbm::apply_packet_routing(
        layout, subpivot, {Real(10.0), Real(20.0)}, state);
    ASSERT_EQ(residual.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    expect_close(residual.residual_number, Real(3.0));
    EXPECT_DOUBLE_EQ(residual.residual_water_mass, Real(0.0));
    expect_close(residual.residual_properties[0], Real(6.0));
    expect_close(residual.residual_properties[1], Real(12.0));
    expect_close(state[0] - before_subpivot[0], Real(2.0));
    expect_close(state[1] - before_subpivot[1], Real(0.0));
    expect_close(state[static_cast<std::size_t>(layout.property_offset(0))], Real(4.0));
    expect_close(state[static_cast<std::size_t>(layout.property_offset(1))], Real(8.0));

    const auto zero_water = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(2.0), Real(0.0));
    ASSERT_EQ(zero_water.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    EXPECT_EQ(zero_water.destination_count, 0);
    const auto all_residual = erf_sbm::apply_packet_routing(
        layout, zero_water, {Real(3.0), Real(7.0)}, state);
    expect_close(all_residual.residual_number, Real(2.0));
    EXPECT_EQ(all_residual.residual_properties, (std::vector<Real>{Real(3.0), Real(7.0)}));
    EXPECT_EQ(all_residual.residual_water_mass, Real(0.0));

    const auto overflow = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), Real(2.1));
    ASSERT_EQ(overflow.status, erf_sbm::RemapStatus::Overflow);
    const auto before_overflow = state;
    EXPECT_EQ(erf_sbm::apply_packet_routing(
        layout, overflow, {Real(1.0), Real(1.0)}, state).status,
        erf_sbm::RemapStatus::Overflow);
    EXPECT_EQ(state, before_overflow);

    const auto invalid = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0),
        std::numeric_limits<Real>::quiet_NaN(), Real(1.0));
    EXPECT_EQ(invalid.status, erf_sbm::RemapStatus::Invalid);
    const auto atomic_plan = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), Real(1.4));
    auto atomic_state = std::vector<Real>(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    atomic_state[static_cast<std::size_t>(layout.property_offset(0) + 1)] =
        std::numeric_limits<Real>::max();
    const auto unchanged = atomic_state;
    const auto failed_apply = erf_sbm::apply_packet_routing(
        layout, atomic_plan, {std::numeric_limits<Real>::max(), Real(0.0)}, atomic_state);
    EXPECT_EQ(failed_apply.status, erf_sbm::RemapStatus::Invalid);
    EXPECT_EQ(atomic_state, unchanged);
}

TEST(SBMRemapping, IntervalTwoMomentPreservesActualMassAndUsesHalfOpenEdges)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment,
        {make_property("coating", 0), make_property("rime", 0)});
    const auto& liquid = population(layout, 0);

    const auto packet = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), Real(1.4));
    ASSERT_EQ(packet.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(packet.destination_count, 1);
    EXPECT_EQ(packet.destinations[0].bin, 0);
    EXPECT_DOUBLE_EQ(packet.destinations[0].particle_mass, Real(1.4));
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto packet_application = erf_sbm::apply_packet_routing(
        layout, packet, {Real(2.0), Real(5.0)}, state);
    ASSERT_EQ(packet_application.status, erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(packet_application.roundoff_water_mass_correction, Real(0.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.mass_offset)], Real(1.4));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.number_offset)], Real(1.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.mass_offset + 1)], Real(0.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.number_offset + 1)], Real(0.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(layout.property_offset(0))], Real(2.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(layout.property_offset(1))], Real(5.0));

    const Real shared = Real(1.5);
    const auto at_shared_edge = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), shared);
    ASSERT_EQ(at_shared_edge.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(at_shared_edge.destinations[0].bin, 1);
    const auto below_shared = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), std::nextafter(shared, Real(0.0)));
    const auto above_shared = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0),
        std::nextafter(shared, std::numeric_limits<Real>::infinity()));
    EXPECT_EQ(below_shared.destinations[0].bin, 0);
    EXPECT_EQ(above_shared.destinations[0].bin, 1);

    const auto at_lower = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), Real(0.5));
    const auto at_upper = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), Real(2.5));
    const auto just_inside_lower = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0),
        std::nextafter(Real(0.5), std::numeric_limits<Real>::infinity()));
    const auto just_inside_upper = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0),
        std::nextafter(Real(2.5), Real(0.0)));
    EXPECT_EQ(at_lower.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(at_lower.destinations[0].bin, 0);
    EXPECT_EQ(at_upper.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(at_upper.destinations[0].bin, 1);
    EXPECT_EQ(just_inside_lower.status, erf_sbm::RemapStatus::Ok);
    EXPECT_FALSE(just_inside_lower.normalized_roundoff);
    EXPECT_EQ(just_inside_upper.status, erf_sbm::RemapStatus::Ok);
    EXPECT_FALSE(just_inside_upper.normalized_roundoff);
    const auto one_step_below = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0),
        std::nextafter(Real(0.5), Real(0.0)));
    EXPECT_EQ(one_step_below.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(one_step_below.destinations[0].bin, 0);
    EXPECT_DOUBLE_EQ(one_step_below.destinations[0].particle_mass, Real(0.5));
    EXPECT_DOUBLE_EQ(one_step_below.packet_particle_mass,
                     std::nextafter(Real(0.5), Real(0.0)));
    EXPECT_TRUE(one_step_below.normalized_roundoff);
    const auto one_step_above = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0),
        std::nextafter(Real(2.5), std::numeric_limits<Real>::infinity()));
    EXPECT_EQ(one_step_above.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(one_step_above.destinations[0].bin, 1);
    EXPECT_DOUBLE_EQ(one_step_above.destinations[0].particle_mass, Real(2.5));
    EXPECT_DOUBLE_EQ(one_step_above.packet_particle_mass,
                     std::nextafter(Real(2.5), std::numeric_limits<Real>::infinity()));
    EXPECT_TRUE(one_step_above.normalized_roundoff);
    auto four_steps_below_mass = Real(0.5);
    for (int step = 0; step < 4; ++step) {
        four_steps_below_mass = std::nextafter(four_steps_below_mass, Real(0.0));
    }
    const auto four_steps_below = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), four_steps_below_mass);
    EXPECT_EQ(four_steps_below.status, erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(four_steps_below.destinations[0].particle_mass, Real(0.5));
    EXPECT_TRUE(four_steps_below.normalized_roundoff);
    const auto five_steps_below_mass = std::nextafter(four_steps_below_mass, Real(0.0));
    EXPECT_EQ(erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0),
        five_steps_below_mass).status, erf_sbm::RemapStatus::BelowSupportedGrid);
    auto four_steps_above_mass = Real(2.5);
    for (int step = 0; step < 4; ++step) {
        four_steps_above_mass = std::nextafter(
            four_steps_above_mass, std::numeric_limits<Real>::infinity());
    }
    const auto four_steps_above = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), four_steps_above_mass);
    EXPECT_EQ(four_steps_above.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(four_steps_above.destinations[0].bin, 1);
    EXPECT_DOUBLE_EQ(four_steps_above.destinations[0].particle_mass, Real(2.5));
    EXPECT_TRUE(four_steps_above.normalized_roundoff);
    const auto five_steps_above_mass = std::nextafter(
        four_steps_above_mass, std::numeric_limits<Real>::infinity());
    EXPECT_EQ(erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0),
        five_steps_above_mass).status, erf_sbm::RemapStatus::Overflow);
    EXPECT_EQ(erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), Real(3.0)).status,
        erf_sbm::RemapStatus::Overflow);
    EXPECT_EQ(erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), Real(0.1)).status,
        erf_sbm::RemapStatus::BelowSupportedGrid);
}

TEST(SBMRemapping, TwoMomentZeroWaterPacketReturnsCompleteResidual)
{
    auto wet_population = make_population(0, erf_sbm::MomentMode::TwoMoment,
        {Real(0.0), Real(1.0), Real(2.0)}, {Real(0.5), Real(1.5)});
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(wet_population);
    spec.liquid_projection = {0, 1};
    spec.attached_properties = {make_property("solute_a", 0), make_property("solute_b", 0)};
    const erf_sbm::SBMLayout layout(std::move(spec));
    const auto view = erf_sbm::population_remap_view(layout, 0);
    const auto plan = erf_sbm::plan_packet_routing(view, Real(2.0), Real(0.0));
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    EXPECT_EQ(plan.destination_count, 0);
    EXPECT_DOUBLE_EQ(plan.residual_number_weight, Real(1.0));

    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto before = state;
    const auto applied = erf_sbm::apply_packet_routing(
        layout, plan, {Real(3.0), Real(8.0)}, state);
    ASSERT_EQ(applied.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    EXPECT_EQ(state, before);
    EXPECT_DOUBLE_EQ(applied.residual_number, Real(2.0));
    EXPECT_DOUBLE_EQ(applied.residual_water_mass, Real(0.0));
    ASSERT_EQ(applied.residual_properties.size(), 2u);
    EXPECT_DOUBLE_EQ(applied.residual_properties[0], Real(3.0));
    EXPECT_DOUBLE_EQ(applied.residual_properties[1], Real(8.0));

    erf_sbm::SBMLayoutSpec no_property_spec;
    no_property_spec.populations.push_back(wet_population);
    no_property_spec.liquid_projection = {0, 1};
    const erf_sbm::SBMLayout no_property_layout(std::move(no_property_spec));
    const auto no_property_plan = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(no_property_layout, 0), Real(1.0), Real(0.0));
    std::vector<Real> no_property_state(static_cast<std::size_t>(no_property_layout.ncomp()), Real(0.0));
    const auto no_property_result = erf_sbm::apply_packet_routing(
        no_property_layout, no_property_plan, {}, no_property_state);
    EXPECT_EQ(no_property_result.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    EXPECT_DOUBLE_EQ(no_property_result.residual_number, Real(1.0));
    EXPECT_TRUE(no_property_result.residual_properties.empty());
    EXPECT_EQ(no_property_state,
              (std::vector<Real>(static_cast<std::size_t>(no_property_layout.ncomp()), Real(0.0))));
}

TEST(SBMRemapping, StrictTwoMomentOwnershipAndNoProcessFixedPoint)
{
    auto layout_spec = erf_sbm::SBMLayoutSpec{};
    layout_spec.populations.push_back(make_population(
        0, erf_sbm::MomentMode::TwoMoment, {Real(0.0), Real(1.0), Real(2.0)},
        {Real(0.5), Real(1.5)}));
    layout_spec.liquid_projection = {0, 1};
    layout_spec.attached_properties = {make_property("coating", 0)};
    const erf_sbm::SBMLayout layout(std::move(layout_spec));
    const auto& pop = population(layout, 0);
    const auto view = erf_sbm::population_remap_view(layout, 0);
    ASSERT_TRUE(erf_sbm::valid_population_remap_view(view));

    // A dry packet still uses the residual path, but the corresponding liquid
    // persisted state is noncanonical even though this grid starts at zero.
    const auto dry_packet = erf_sbm::plan_packet_routing(view, Real(1.0), Real(0.0));
    ASSERT_EQ(dry_packet.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    EXPECT_EQ(dry_packet.destination_count, 0);
    std::vector<Real> dry_state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    dry_state[static_cast<std::size_t>(pop.number_offset)] = Real(1.0);
    EXPECT_FALSE(erf_sbm::remap_detail::canonical_two_moment_bin_state(
        view, 0, dry_state.data(), static_cast<int>(dry_state.size())));
    erf_sbm::ReconstructionDelta rejected;
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, dry_state, rejected),
              erf_sbm::ReconstructionStatus::Invalid);
    auto bad_projection = erf_sbm::ReconstructionDeltaView{};
    bad_projection.population = view;
    bad_projection.bin = 0;
    bad_projection.empty = false;
    bad_projection.number = Real(1.0);
    bad_projection.particle_mass = Real(0.0);
    const Real zero_property = Real(0.0);
    bad_projection.property_per_particle = &zero_property;
    bad_projection.property_count = 1;
    const auto unchanged = dry_state;
    EXPECT_FALSE(erf_sbm::project_reconstruction_bin_core(
        view, bad_projection, dry_state.data(), static_cast<int>(dry_state.size())));
    EXPECT_EQ(dry_state, unchanged);

    // Shared edges are upper-owned. Persisting one in the lower bin is
    // rejected, and the same moment pair in the upper bin is canonical.
    auto shared_plan = erf_sbm::plan_packet_routing(view, Real(1.0), Real(1.0));
    ASSERT_EQ(shared_plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(shared_plan.destination_count, 1);
    EXPECT_EQ(shared_plan.destinations[0].bin, 1);
    std::vector<Real> edge_state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    edge_state[static_cast<std::size_t>(pop.mass_offset)] = Real(1.0);
    edge_state[static_cast<std::size_t>(pop.number_offset)] = Real(1.0);
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, edge_state, rejected),
              erf_sbm::ReconstructionStatus::Invalid);
    edge_state[static_cast<std::size_t>(pop.mass_offset)] = Real(0.0);
    edge_state[static_cast<std::size_t>(pop.number_offset)] = Real(0.0);
    edge_state[static_cast<std::size_t>(pop.mass_offset + 1)] = Real(1.0);
    edge_state[static_cast<std::size_t>(pop.number_offset + 1)] = Real(1.0);
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 1, edge_state, rejected),
              erf_sbm::ReconstructionStatus::Populated);

    // The tolerant endpoint transform remains tolerant, while persisted-state
    // admission rejects the out-of-owned-interval mean without normalizing it.
    const Real upper = Real(1.0);
    const Real tolerance_probe = upper + Real(64.0) * std::numeric_limits<Real>::epsilon();
    erf_sbm::EndpointTransform endpoint_probe;
    ASSERT_TRUE(erf_sbm::try_two_moment_to_endpoints(
        Real(1.0), tolerance_probe, Real(0.0), upper, endpoint_probe));
    EXPECT_TRUE(endpoint_probe.normalized_roundoff);
    std::vector<Real> noncanonical(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    noncanonical[static_cast<std::size_t>(pop.mass_offset)] = tolerance_probe;
    noncanonical[static_cast<std::size_t>(pop.number_offset)] = Real(1.0);
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, noncanonical, rejected),
              erf_sbm::ReconstructionStatus::Invalid);
    EXPECT_DOUBLE_EQ(noncanonical[static_cast<std::size_t>(pop.mass_offset)], tolerance_probe);
    bad_projection.number = Real(1.0);
    bad_projection.particle_mass = tolerance_probe;
    const auto noncanonical_before = noncanonical;
    EXPECT_FALSE(erf_sbm::project_reconstruction_bin_core(
        view, bad_projection, noncanonical.data(), static_cast<int>(noncanonical.size())));
    EXPECT_EQ(noncanonical, noncanonical_before);

    // A positive tiny mean in the zero-based first interval remains valid.
    const Real tiny_positive_mass = std::numeric_limits<Real>::epsilon();
    std::vector<Real> tiny_state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    tiny_state[static_cast<std::size_t>(pop.mass_offset)] = tiny_positive_mass;
    tiny_state[static_cast<std::size_t>(pop.number_offset)] = Real(1.0);
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, tiny_state, rejected),
              erf_sbm::ReconstructionStatus::Populated);

    // Every bin reconstructs and redeposits with its original moments and
    // attached extensive inventory, including the global top endpoint.
    std::vector<Real> original(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    original[static_cast<std::size_t>(pop.mass_offset)] = Real(1.0);
    original[static_cast<std::size_t>(pop.number_offset)] = Real(2.0);
    original[static_cast<std::size_t>(layout.property_offset(0))] = Real(6.0);
    original[static_cast<std::size_t>(pop.mass_offset + 1)] = Real(2.0);
    original[static_cast<std::size_t>(pop.number_offset + 1)] = Real(1.0);
    original[static_cast<std::size_t>(layout.property_offset(0) + 1)] = Real(4.0);
    std::vector<Real> redeposited(original.size(), Real(0.0));
    for (int bin = 0; bin < pop.grid.nbins(); ++bin) {
        erf_sbm::ReconstructionDelta reconstructed;
        ASSERT_EQ(erf_sbm::reconstruct_bin(layout, 0, bin, original, reconstructed),
                  erf_sbm::ReconstructionStatus::Populated);
        const Real property_amount = reconstructed.number *
            reconstructed.property_per_particle[0];
        const auto route = erf_sbm::plan_packet_routing(
            view, reconstructed.number, reconstructed.particle_mass);
        ASSERT_EQ(route.status, erf_sbm::RemapStatus::Ok);
        ASSERT_EQ(route.destination_count, 1);
        EXPECT_EQ(route.destinations[0].bin, bin);
        const auto applied = erf_sbm::apply_packet_routing(
            layout, route, {property_amount}, redeposited);
        ASSERT_EQ(applied.status, erf_sbm::RemapStatus::Ok);
    }
    for (std::size_t component = 0; component < original.size(); ++component) {
        expect_close(redeposited[component], original[component]);
    }
}

TEST(SBMRemapping, PositiveProductUnderflowHasDistinctAtomicStatus)
{
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    ASSERT_GT(tiny, Real(0.0));
    Real product = Real(0.0);
    EXPECT_EQ(erf_sbm::remap_detail::checked_product(tiny, Real(0.5), product),
              erf_sbm::remap_detail::ProductStatus::Underflow);
    EXPECT_DOUBLE_EQ(product, Real(0.0));

    const auto layout = make_layout(erf_sbm::MomentMode::OneMoment,
        {make_property("coating", 0)});
    const auto view = erf_sbm::population_remap_view(layout, 0);
    const auto plan = erf_sbm::plan_packet_routing(view, Real(1.0), Real(0.5));
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    std::vector<Real> residual_properties{Real(77.0)};
    const auto unchanged = state;
    const auto applied = erf_sbm::apply_packet_routing_core(
        view, plan, &tiny, 1, residual_properties.data(), 1,
        state.data(), static_cast<int>(state.size()));
    EXPECT_EQ(applied.status, erf_sbm::RemapStatus::NumericalUnderflow);
    EXPECT_DOUBLE_EQ(applied.residual_number, Real(0.0));
    EXPECT_DOUBLE_EQ(applied.residual_water_mass, Real(0.0));
    EXPECT_DOUBLE_EQ(applied.roundoff_water_mass_correction, Real(0.0));
    EXPECT_FALSE(applied.normalized_roundoff);
    EXPECT_EQ(state, unchanged);
    ASSERT_EQ(residual_properties.size(), 1u);
    EXPECT_DOUBLE_EQ(residual_properties[0], Real(77.0));

    const auto two_moment = make_layout(erf_sbm::MomentMode::TwoMoment);
    const auto two_view = erf_sbm::population_remap_view(two_moment, 0);
    const auto two_plan = erf_sbm::plan_packet_routing(two_view, tiny, Real(0.5));
    ASSERT_EQ(two_plan.status, erf_sbm::RemapStatus::Ok);
    std::vector<Real> two_state(static_cast<std::size_t>(two_moment.ncomp()), Real(0.0));
    const auto two_unchanged = two_state;
    const auto two_applied = erf_sbm::apply_packet_routing_core(
        two_view, two_plan, nullptr, 0, nullptr, 0,
        two_state.data(), static_cast<int>(two_state.size()));
    EXPECT_EQ(two_applied.status, erf_sbm::RemapStatus::NumericalUnderflow);
    EXPECT_EQ(two_state, two_unchanged);
}

TEST(SBMRemapping, PositiveQuotientUnderflowRejectsPersistedAndDepositedState)
{
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    ASSERT_GT(tiny, Real(0.0));

    erf_sbm::SBMLayoutSpec tiny_route_spec;
    tiny_route_spec.populations.push_back(make_population(
        0, erf_sbm::MomentMode::OneMoment,
        {Real(1.0), Real(3.0), Real(5.0)}, {Real(2.0), Real(4.0)}));
    tiny_route_spec.liquid_projection = {0, 1};
    const erf_sbm::SBMLayout tiny_route_layout(std::move(tiny_route_spec));
    const Real tiny_pivot = tiny_route_layout.populations()[0].grid.pivot(0);
    ASSERT_GT(tiny_pivot, Real(0.0));
    EXPECT_EQ(tiny / tiny_pivot, Real(0.0));
    const auto tiny_route = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(tiny_route_layout, 0), Real(1.0), tiny);
    EXPECT_EQ(tiny_route.status, erf_sbm::RemapStatus::NumericalUnderflow);

   const auto two_moment = make_layout(erf_sbm::MomentMode::TwoMoment,
        {make_property("coating", 0)});
    const auto& two_pop = population(two_moment, 0);
    const auto two_view = erf_sbm::population_remap_view(two_moment, 0);
    std::vector<Real> persisted(static_cast<std::size_t>(two_moment.ncomp()), Real(0.0));
    persisted[static_cast<std::size_t>(two_pop.mass_offset)] = Real(2.0);
    persisted[static_cast<std::size_t>(two_pop.number_offset)] = Real(2.0);
    persisted[static_cast<std::size_t>(two_moment.property_offset(0))] = tiny;
    EXPECT_FALSE(erf_sbm::remap_detail::canonical_two_moment_bin_state(
        two_view, 0, persisted.data(), static_cast<int>(persisted.size())));
    EXPECT_FALSE(erf_sbm::remap_detail::canonical_persisted_bin_state(
        two_view, 0, persisted.data(), static_cast<int>(persisted.size())));
    erf_sbm::ReconstructionDelta reconstruction;
    EXPECT_EQ(erf_sbm::reconstruct_bin(two_moment, 0, 0, persisted, reconstruction),
              erf_sbm::ReconstructionStatus::Invalid);

    const auto two_plan = erf_sbm::plan_packet_routing(two_view, Real(2.0), Real(1.0));
    ASSERT_EQ(two_plan.status, erf_sbm::RemapStatus::Ok);
    std::vector<Real> two_candidate(static_cast<std::size_t>(two_moment.ncomp()), Real(0.0));
    const auto two_before = two_candidate;
    Real residual_property = Real(77.0);
    const auto two_applied = erf_sbm::apply_packet_routing_core(
        two_view, two_plan, &tiny, 1, &residual_property, 1, two_candidate.data(),
        static_cast<int>(two_candidate.size()));
    EXPECT_EQ(two_applied.status, erf_sbm::RemapStatus::NumericalUnderflow);
    EXPECT_DOUBLE_EQ(two_applied.residual_number, Real(0.0));
    EXPECT_DOUBLE_EQ(two_applied.residual_water_mass, Real(0.0));
    EXPECT_DOUBLE_EQ(two_applied.roundoff_water_mass_correction, Real(0.0));
    EXPECT_FALSE(two_applied.normalized_roundoff);
    EXPECT_EQ(two_candidate, two_before);
    EXPECT_DOUBLE_EQ(residual_property, Real(77.0));

    const auto one_moment = make_layout(erf_sbm::MomentMode::OneMoment,
        {make_property("coating", 0)});
    const auto& one_pop = population(one_moment, 0);
    const auto one_view = erf_sbm::population_remap_view(one_moment, 0);
    std::vector<Real> one_persisted(static_cast<std::size_t>(one_moment.ncomp()), Real(0.0));
    one_persisted[static_cast<std::size_t>(one_pop.mass_offset)] = Real(4.0);
    one_persisted[static_cast<std::size_t>(one_moment.property_offset(0))] = tiny;
    EXPECT_FALSE(erf_sbm::remap_detail::canonical_persisted_bin_state(
        one_view, 0, one_persisted.data(), static_cast<int>(one_persisted.size())));
    EXPECT_EQ(erf_sbm::reconstruct_bin(one_moment, 0, 0, one_persisted, reconstruction),
              erf_sbm::ReconstructionStatus::Invalid);

    const auto one_plan = erf_sbm::plan_packet_routing(one_view, Real(4.0), Real(1.0));
    ASSERT_EQ(one_plan.status, erf_sbm::RemapStatus::Ok);
    std::vector<Real> one_candidate(static_cast<std::size_t>(one_moment.ncomp()), Real(0.0));
    const auto one_before = one_candidate;
    residual_property = Real(77.0);
    const auto one_applied = erf_sbm::apply_packet_routing_core(
        one_view, one_plan, &tiny, 1, &residual_property, 1, one_candidate.data(),
        static_cast<int>(one_candidate.size()));
    EXPECT_EQ(one_applied.status, erf_sbm::RemapStatus::NumericalUnderflow);
    EXPECT_DOUBLE_EQ(one_applied.residual_number, Real(0.0));
    EXPECT_DOUBLE_EQ(one_applied.residual_water_mass, Real(0.0));
    EXPECT_DOUBLE_EQ(one_applied.roundoff_water_mass_correction, Real(0.0));
    EXPECT_FALSE(one_applied.normalized_roundoff);
    EXPECT_EQ(one_candidate, one_before);
    EXPECT_DOUBLE_EQ(residual_property, Real(77.0));
}

TEST(SBMRemapping, NormalizedTwoMomentBoundariesPreserveZeroProcessIdentity)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    const auto& pop = population(layout, 0);
    const auto groups = erf_sbm::make_constraint_groups(layout);

    for (const bool upper_boundary : {false, true}) {
        SCOPED_TRACE(upper_boundary ? "upper global edge" : "lower global edge");
        const int bin = upper_boundary ? pop.grid.nbins() - 1 : 0;
        const Real edge = upper_boundary ? pop.grid.edges().back() : pop.grid.edges().front();
        const Real input_mass = std::nextafter(edge, upper_boundary ?
            std::numeric_limits<Real>::infinity() : Real(0.0));
        const auto plan = erf_sbm::plan_packet_routing(
            erf_sbm::population_remap_view(layout, 0), Real(1.0), input_mass);
        ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
        ASSERT_TRUE(plan.normalized_roundoff);
        ASSERT_EQ(plan.destination_count, 1);
        EXPECT_EQ(plan.destinations[0].bin, bin);
        EXPECT_DOUBLE_EQ(plan.packet_particle_mass, input_mass);
        EXPECT_DOUBLE_EQ(plan.destinations[0].particle_mass, edge);

        std::vector<Real> persisted_state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
        const auto applied = erf_sbm::apply_packet_routing(
            layout, plan, {}, persisted_state);
        ASSERT_EQ(applied.status, erf_sbm::RemapStatus::Ok);
        EXPECT_TRUE(applied.normalized_roundoff);
        const Real expected_correction = edge - input_mass;
        EXPECT_DOUBLE_EQ(applied.roundoff_water_mass_correction, expected_correction);
        if (upper_boundary) {
            EXPECT_LT(applied.roundoff_water_mass_correction, Real(0.0));
        } else {
            EXPECT_GT(applied.roundoff_water_mass_correction, Real(0.0));
        }
        const Real stored_mass = persisted_state[static_cast<std::size_t>(pop.mass_offset + bin)];
        const Real stored_number = persisted_state[static_cast<std::size_t>(pop.number_offset + bin)];
        EXPECT_DOUBLE_EQ(stored_number, Real(1.0));
        EXPECT_DOUBLE_EQ(stored_mass, edge);
        EXPECT_DOUBLE_EQ(input_mass + applied.roundoff_water_mass_correction, stored_mass);

        erf_sbm::EndpointTransform endpoints;
        ASSERT_TRUE(erf_sbm::try_two_moment_to_endpoints(
            stored_number, stored_mass, pop.grid.edges()[static_cast<std::size_t>(bin)],
            pop.grid.edges()[static_cast<std::size_t>(bin + 1)], endpoints));
        EXPECT_FALSE(endpoints.normalized_roundoff);
        const auto group = std::find_if(groups.begin(), groups.end(),
            [bin](const erf_sbm::ConstraintGroup& candidate) {
                return candidate.population_id == 0 && candidate.bin == bin;
            });
        ASSERT_NE(group, groups.end());
        EXPECT_TRUE(group->admissible(persisted_state));

        erf_sbm::ReconstructionDelta reconstruction;
        ASSERT_EQ(erf_sbm::reconstruct_bin(layout, 0, bin, persisted_state, reconstruction),
                  erf_sbm::ReconstructionStatus::Populated);
        EXPECT_FALSE(reconstruction.normalized_roundoff);
        std::vector<Real> projected_state(persisted_state.size(), Real(0.0));
        ASSERT_TRUE(erf_sbm::project_reconstruction_bin(layout, reconstruction, projected_state));
        EXPECT_EQ(projected_state, persisted_state);
    }
}

TEST(SBMRemapping, BoundaryCorrectionUnderflowFailsAtomically)
{
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(make_population(0, erf_sbm::MomentMode::TwoMoment,
        {Real(1.0), Real(2.0), Real(3.0)}, {Real(1.5), Real(2.5)}));
    spec.attached_properties.push_back(make_property("coating", 0));
    spec.liquid_projection = {0, 1};
    const erf_sbm::SBMLayout layout(std::move(spec));
    const auto& pop = population(layout, 0);
    const Real boundary_mass = pop.grid.edges().front();
    const Real input_mass = std::nextafter(boundary_mass, Real(0.0));
    const Real packet_number = std::numeric_limits<Real>::denorm_min();

    ASSERT_GT(packet_number, Real(0.0));
    ASSERT_NE(input_mass, boundary_mass);
    const Real mass_delta = boundary_mass - input_mass;
    ASSERT_NE(mass_delta, Real(0.0));
    const Real deposited_water_increment = packet_number * boundary_mass;
    ASSERT_GT(deposited_water_increment, Real(0.0));
    const Real naive_correction = packet_number * std::abs(mass_delta);
    ASSERT_EQ(naive_correction, Real(0.0));

    const auto view = erf_sbm::population_remap_view(layout, 0);
    const auto plan = erf_sbm::plan_packet_routing(view, packet_number, input_mass);
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_TRUE(plan.normalized_roundoff);
    ASSERT_EQ(plan.destination_count, 1);
    EXPECT_DOUBLE_EQ(plan.destinations[0].particle_mass, boundary_mass);

    std::vector<Real> candidate(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto before = candidate;
    Real residual_property = Real(77.0);
    const auto applied = erf_sbm::apply_packet_routing_core(
        view, plan, &packet_number, 1, &residual_property, 1,
        candidate.data(), static_cast<int>(candidate.size()));
    EXPECT_EQ(applied.status, erf_sbm::RemapStatus::NumericalUnderflow);
    EXPECT_DOUBLE_EQ(applied.residual_number, Real(0.0));
    EXPECT_DOUBLE_EQ(applied.residual_water_mass, Real(0.0));
    EXPECT_DOUBLE_EQ(applied.roundoff_water_mass_correction, Real(0.0));
    EXPECT_FALSE(applied.normalized_roundoff);
    EXPECT_EQ(candidate, before);
    EXPECT_DOUBLE_EQ(residual_property, Real(77.0));
}

TEST(SBMRemapping, ZeroNumberPacketIsNoOpAndRejectsUnattachedMaterial)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment,
        {make_property("solute_a", 0), make_property("solute_b", 0)});
    const auto plan = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(0.0), Real(0.0));
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(plan.destination_count, 0);
    EXPECT_DOUBLE_EQ(plan.residual_number_weight, Real(0.0));
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    state[0] = Real(1.25);
    const auto before = state;
    const auto no_op = erf_sbm::apply_packet_routing(
        layout, plan, {Real(0.0), Real(0.0)}, state);
    EXPECT_EQ(no_op.status, erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(no_op.residual_number, Real(0.0));
    EXPECT_EQ(no_op.residual_properties, (std::vector<Real>{Real(0.0), Real(0.0)}));
    EXPECT_EQ(state, before);
    EXPECT_EQ(erf_sbm::apply_packet_routing(layout, plan, {Real(1.0), Real(0.0)}, state).status,
              erf_sbm::RemapStatus::Invalid);
    EXPECT_EQ(state, before);
}

TEST(SBMRemapping, ReferenceLayoutRejectsNonMassCoordinates)
{
    auto radius_population = make_population(0, erf_sbm::MomentMode::TwoMoment);
    radius_population.grid.coordinate_kind = erf_sbm::CoordinateKind::Radius;
    radius_population.grid.coordinate_units = "m";
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(std::move(radius_population));
    spec.liquid_projection = {0, 1};
    EXPECT_TRUE(erf_sbm::SpectralGrid::validate(spec.populations[0].grid).valid);
    const auto validation = erf_sbm::SBMLayout::validate(spec);
    EXPECT_FALSE(validation.valid);
    EXPECT_NE(validation.message.find("mass"), std::string::npos);

    const auto valid_layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    auto non_mass_view = erf_sbm::population_remap_view(valid_layout, 0);
    non_mass_view.coordinate_kind = erf_sbm::CoordinateKind::Radius;
    EXPECT_EQ(erf_sbm::plan_packet_routing(non_mass_view, Real(1.0), Real(1.4)).status,
              erf_sbm::RemapStatus::Invalid);
}

TEST(SBMRemapping, RoutingPlanCannotCrossPopulationOrGridContexts)
{
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(make_population(0, erf_sbm::MomentMode::OneMoment));
    spec.populations.push_back(make_population(1, erf_sbm::MomentMode::TwoMoment,
        {Real(0.5), Real(1.5), Real(2.5)}, {Real(1.0), Real(2.0)},
        erf_sbm::PopulationPhase::Aerosol));
    spec.liquid_projection = {0, 1};
    const erf_sbm::SBMLayout layout(std::move(spec));
    const auto first_view = erf_sbm::population_remap_view(layout, 0);
    const auto second_view = erf_sbm::population_remap_view(layout, 1);
    const auto first_plan = erf_sbm::plan_packet_routing(first_view, Real(1.0), Real(1.4));
    std::vector<Real> candidate(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto unchanged = candidate;
    const auto mismatch = erf_sbm::apply_packet_routing_core(
        second_view, first_plan, nullptr, 0, nullptr, 0,
        candidate.data(), static_cast<int>(candidate.size()));
    EXPECT_EQ(mismatch.status, erf_sbm::RemapStatus::Invalid);
    EXPECT_EQ(candidate, unchanged);

    auto forged_plan = first_plan;
    forged_plan.destinations[0].bin = 1;
    const auto forged = erf_sbm::apply_packet_routing_core(
        first_view, forged_plan, nullptr, 0, nullptr, 0,
        candidate.data(), static_cast<int>(candidate.size()));
    EXPECT_EQ(forged.status, erf_sbm::RemapStatus::Invalid);
    EXPECT_EQ(candidate, unchanged);

    const auto second_plan = erf_sbm::plan_packet_routing(second_view, Real(1.0), Real(1.4));
    ASSERT_EQ(second_plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(second_plan.destination_count, 1);
    EXPECT_EQ(second_plan.destinations[0].bin, 0);
    EXPECT_DOUBLE_EQ(second_plan.destinations[0].particle_mass, Real(1.4));
    const auto accepted = erf_sbm::apply_packet_routing_core(
        second_view, second_plan, nullptr, 0, nullptr, 0,
        candidate.data(), static_cast<int>(candidate.size()));
    ASSERT_EQ(accepted.status, erf_sbm::RemapStatus::Ok);
    const auto& second = population(layout, 1);
    EXPECT_DOUBLE_EQ(candidate[static_cast<std::size_t>(second.mass_offset)], Real(1.4));
    EXPECT_DOUBLE_EQ(candidate[static_cast<std::size_t>(second.number_offset)], Real(1.0));

    erf_sbm::SBMLayoutSpec other_spec;
    other_spec.populations.push_back(make_population(0, erf_sbm::MomentMode::OneMoment,
        {Real(0.5), Real(1.25), Real(2.5)}, {Real(0.9), Real(1.75)}));
    other_spec.liquid_projection = {0, 1};
    const erf_sbm::SBMLayout other_layout(std::move(other_spec));
    const auto other_view = erf_sbm::population_remap_view(other_layout, 0);
    const auto other_plan = erf_sbm::plan_packet_routing(first_view, Real(1.0), Real(1.4));
    std::vector<Real> other_state(static_cast<std::size_t>(other_layout.ncomp()), Real(0.0));
    const auto other_before = other_state;
    const auto stale_grid = erf_sbm::apply_packet_routing_core(
        other_view, other_plan, nullptr, 0, nullptr, 0,
        other_state.data(), static_cast<int>(other_state.size()));
    EXPECT_EQ(stale_grid.status, erf_sbm::RemapStatus::Invalid);
    EXPECT_EQ(other_state, other_before);
}

TEST(SBMRemapping, MeanDeltaReconstructionProjectsExactlyAndIntegratesIntervals)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment,
        {make_property("coating", 0), make_property("rime", 0)});
    const auto& liquid = population(layout, 0);
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    state[0] = Real(1.4);
    state[1] = Real(4.4);
    state[static_cast<std::size_t>(liquid.number_offset)] = Real(1.0);
    state[static_cast<std::size_t>(liquid.number_offset + 1)] = Real(2.0);
    state[static_cast<std::size_t>(layout.property_offset(0))] = Real(2.0);
    state[static_cast<std::size_t>(layout.property_offset(0) + 1)] = Real(6.0);
    state[static_cast<std::size_t>(layout.property_offset(1))] = Real(5.0);
    state[static_cast<std::size_t>(layout.property_offset(1) + 1)] = Real(8.0);

    std::vector<Real> projected(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    erf_sbm::ReconstructionDelta first;
    ASSERT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, state, first),
              erf_sbm::ReconstructionStatus::Populated);
    EXPECT_DOUBLE_EQ(first.number, Real(1.0));
    EXPECT_DOUBLE_EQ(first.particle_mass, Real(1.4));
    EXPECT_DOUBLE_EQ(first.property_per_particle[0], Real(2.0));
    EXPECT_DOUBLE_EQ(first.property_per_particle[1], Real(5.0));
    ASSERT_TRUE(erf_sbm::project_reconstruction_bin(layout, first, projected));

    erf_sbm::ReconstructionDelta second;
    ASSERT_EQ(erf_sbm::reconstruct_bin(layout, 0, 1, state, second),
              erf_sbm::ReconstructionStatus::Populated);
    EXPECT_DOUBLE_EQ(second.number, Real(2.0));
    EXPECT_DOUBLE_EQ(second.particle_mass, Real(2.2));
    EXPECT_DOUBLE_EQ(second.property_per_particle[0], Real(3.0));
    EXPECT_DOUBLE_EQ(second.property_per_particle[1], Real(4.0));
    ASSERT_TRUE(erf_sbm::project_reconstruction_bin(layout, second, projected));
    EXPECT_EQ(projected, state);

    erf_sbm::IntegratedMoments integral;
    ASSERT_TRUE(erf_sbm::integrate_interval(first, Real(0.5), Real(1.0), false, integral));
    EXPECT_DOUBLE_EQ(integral.number, Real(0.0));
    EXPECT_EQ(integral.attached_properties, (std::vector<Real>{Real(0.0), Real(0.0)}));
    ASSERT_TRUE(erf_sbm::integrate_interval(first, Real(1.3), Real(1.4), false, integral));
    EXPECT_DOUBLE_EQ(integral.number, Real(0.0));
    ASSERT_TRUE(erf_sbm::integrate_interval(first, Real(1.4), Real(1.5), false, integral));
    EXPECT_DOUBLE_EQ(integral.number, Real(1.0));
    EXPECT_DOUBLE_EQ(integral.water_mass, Real(1.4));
    EXPECT_EQ(integral.attached_properties, (std::vector<Real>{Real(2.0), Real(5.0)}));
    ASSERT_TRUE(erf_sbm::integrate_interval(first, Real(1.3), Real(1.4), true, integral));
    EXPECT_DOUBLE_EQ(integral.number, Real(1.0));

    std::vector<Real> empty_state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    erf_sbm::ReconstructionDelta empty;
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, empty_state, empty),
              erf_sbm::ReconstructionStatus::Empty);
    EXPECT_TRUE(empty.empty);
    EXPECT_DOUBLE_EQ(empty.particle_mass, Real(0.0));
    ASSERT_TRUE(erf_sbm::integrate_interval(empty, Real(0.5), Real(1.5), false, integral));
    EXPECT_DOUBLE_EQ(integral.number, Real(0.0));
    empty_state[0] = Real(1.0);
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, empty_state, empty),
              erf_sbm::ReconstructionStatus::Invalid);
    empty_state[0] = Real(0.0);
    empty_state[static_cast<std::size_t>(layout.property_offset(0))] = Real(1.0);
    EXPECT_EQ(erf_sbm::reconstruct_bin(layout, 0, 0, empty_state, empty),
              erf_sbm::ReconstructionStatus::Invalid);
}

TEST(SBMRemapping, CanonicalEndpointTransformIsUnitInvariantAndDimensioned)
{
    Real reference_low = Real(0.0);
    Real reference_high = Real(0.0);
    for (const Real scale : {Real(1.0e-12), Real(1.0), Real(1.0e12)}) {
        const Real lower = Real(0.5) * scale;
        const Real upper = Real(1.5) * scale;
        const Real count = Real(2.0);
        const Real mass = Real(2.8) * scale;
        const auto endpoints = erf_sbm::SpectralGrid::two_moment_to_endpoints(
            count, mass, lower, upper);
        expect_close(endpoints.L, Real(0.2));
        expect_close(endpoints.H, Real(1.8));
        if (scale == Real(1.0e-12)) {
            reference_low = endpoints.L;
            reference_high = endpoints.H;
        } else {
            expect_close(endpoints.L, reference_low);
            expect_close(endpoints.H, reference_high);
        }

        const Real upper_mass = upper * count;
        const Real just_above = std::nextafter(upper_mass,
            std::numeric_limits<Real>::infinity());
        const auto roundoff = erf_sbm::transform_two_moment(
            count, just_above, lower, upper);
        EXPECT_TRUE(roundoff.normalized_roundoff);
        EXPECT_DOUBLE_EQ(roundoff.L, Real(0.0));
        expect_close(roundoff.H, count);
        EXPECT_NEAR(roundoff.endpoint_tolerance,
                    roundoff.moment_tolerance / (upper - lower),
                    Real(8.0) * std::numeric_limits<Real>::epsilon() *
                    std::abs(roundoff.endpoint_tolerance));
    }

    // A fixed negative endpoint-number violation is rejected at either mass
    // scale. Comparing it directly with tau_M would incorrectly accept only
    // the large-unit case because tau_M carries mass units.
    for (const Real scale : {Real(1.0e-6), Real(1.0e6)}) {
        const Real lower = Real(0.5) * scale;
        const Real upper = Real(1.5) * scale;
        const Real count = Real(1.0);
        const Real negative_endpoint_number = Real(1024.0) *
            std::numeric_limits<Real>::epsilon();
        const Real mass = upper * count + negative_endpoint_number * (upper - lower);
        EXPECT_FALSE(erf_sbm::SpectralGrid::two_moment_realizable(count, mass, lower, upper));
        EXPECT_THROW(static_cast<void>(erf_sbm::SpectralGrid::two_moment_to_endpoints(
                         count, mass, lower, upper)), std::invalid_argument);
    }
    const auto exact = erf_sbm::SpectralGrid::two_moment_to_endpoints(
        Real(1.0), Real(1.5), Real(0.5), Real(1.5));
    EXPECT_DOUBLE_EQ(exact.L, Real(0.0));
    EXPECT_DOUBLE_EQ(exact.H, Real(1.0));
}

TEST(SBMRemapping, CanonicalInverseAndCompatibilityInverseAreIdentical)
{
    for (const auto& endpoints : {
             std::pair<Real, Real>{Real(0.0), Real(3.0)},
             std::pair<Real, Real>{Real(2.0), Real(0.0)},
             std::pair<Real, Real>{Real(0.25), Real(0.75)},
             std::pair<Real, Real>{Real(0.0), Real(0.0)}}) {
        const auto canonical = erf_sbm::SpectralGrid::endpoints_to_two_moment(
            endpoints.first, endpoints.second, Real(0.5), Real(1.5));
        const auto compatibility = erf_sbm::inverse_two_moment(
            endpoints.first, endpoints.second, Real(0.5), Real(1.5));
        EXPECT_DOUBLE_EQ(canonical.first, compatibility.first);
        EXPECT_DOUBLE_EQ(canonical.second, compatibility.second);
        EXPECT_DOUBLE_EQ(canonical.first, endpoints.first + endpoints.second);
        EXPECT_DOUBLE_EQ(canonical.second,
            std::fma(Real(0.5), endpoints.first, Real(1.5) * endpoints.second));
    }
    Real count = Real(0.0);
    Real mass = Real(0.0);
    EXPECT_FALSE(erf_sbm::try_endpoints_to_two_moment(
        std::numeric_limits<Real>::infinity(), Real(0.0), Real(0.5), Real(1.5), count, mass));
    EXPECT_FALSE(erf_sbm::try_endpoints_to_two_moment(
        Real(-1.0), Real(0.0), Real(0.5), Real(1.5), count, mass));
    EXPECT_THROW(static_cast<void>(erf_sbm::inverse_two_moment(
        std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(),
        Real(0.5), Real(1.5))), std::invalid_argument);
}

TEST(SBMRemapping, AllocationFreeCoreMatchesHostReferenceAdapters)
{
    static_assert(std::is_trivially_copyable_v<erf_sbm::PopulationRemapView>);
    static_assert(std::is_trivially_copyable_v<erf_sbm::RouteDestination>);
    static_assert(std::is_trivially_copyable_v<erf_sbm::PacketRoutingPlan>);
    static_assert(std::is_trivially_copyable_v<erf_sbm::PacketApplicationCoreResult>);
    static_assert(std::is_trivially_copyable_v<erf_sbm::ReconstructionDeltaView>);
    static_assert(std::is_trivially_copyable_v<erf_sbm::IntegratedMomentsCoreResult>);

    for (const auto mode : {erf_sbm::MomentMode::OneMoment,
                            erf_sbm::MomentMode::TwoMoment}) {
        const auto layout = make_layout(mode,
            {make_property("solute_a", 0), make_property("solute_b", 0)});
        const auto population_view = erf_sbm::population_remap_view(layout, 0);
        const auto plan = erf_sbm::plan_packet_routing(population_view, Real(3.0), Real(1.4));
        ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
        if (mode == erf_sbm::MomentMode::OneMoment) {
            ASSERT_EQ(plan.destination_count, 2);
            EXPECT_EQ(plan.destinations[0].bin, 0);
            EXPECT_EQ(plan.destinations[1].bin, 1);
            expect_close(plan.destinations[0].number_weight, Real(0.6));
            expect_close(plan.destinations[1].number_weight, Real(0.4));
        } else {
            ASSERT_EQ(plan.destination_count, 1);
            EXPECT_EQ(plan.destinations[0].bin, 0);
            EXPECT_DOUBLE_EQ(plan.destinations[0].particle_mass, Real(1.4));
        }
        const std::vector<Real> property_amounts{Real(6.0), Real(15.0)};
        std::vector<Real> core_state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
        std::vector<Real> host_state(core_state.size(), Real(0.0));
        Real residual_properties[2] = {Real(-1.0), Real(-1.0)};
        const auto core_application = erf_sbm::apply_packet_routing_core(
            population_view, plan, property_amounts.data(), 2,
            residual_properties, 2, core_state.data(), static_cast<int>(core_state.size()));
        const auto host_application = erf_sbm::apply_packet_routing(
            layout, plan, property_amounts, host_state);
        EXPECT_EQ(core_application.status, host_application.status);
        EXPECT_DOUBLE_EQ(core_application.residual_number, host_application.residual_number);
        EXPECT_DOUBLE_EQ(core_application.residual_water_mass, host_application.residual_water_mass);
        EXPECT_DOUBLE_EQ(core_application.roundoff_water_mass_correction,
                         host_application.roundoff_water_mass_correction);
        EXPECT_DOUBLE_EQ(core_application.roundoff_water_mass_correction, Real(0.0));
        EXPECT_EQ(core_application.normalized_roundoff, host_application.normalized_roundoff);
        ASSERT_EQ(host_application.residual_properties.size(), 2u);
        for (int property = 0; property < 2; ++property) {
            EXPECT_DOUBLE_EQ(residual_properties[property],
                             host_application.residual_properties[static_cast<std::size_t>(property)]);
        }
        ASSERT_EQ(core_state.size(), host_state.size());
        for (std::size_t component = 0; component < core_state.size(); ++component) {
            EXPECT_DOUBLE_EQ(core_state[component], host_state[component]);
        }

        std::vector<Real> core_projected(core_state.size(), Real(0.0));
        std::vector<Real> host_projected(host_state.size(), Real(0.0));
        const auto& pop = population(layout, 0);
        for (int bin = 0; bin < pop.grid.nbins(); ++bin) {
            Real core_properties[2] = {Real(-1.0), Real(-1.0)};
            erf_sbm::ReconstructionDeltaView core_reconstruction;
            const auto core_status = erf_sbm::reconstruct_bin_core(
                population_view, bin, core_state.data(), static_cast<int>(core_state.size()),
                core_properties, 2, core_reconstruction);
            erf_sbm::ReconstructionDelta host_reconstruction;
            const auto host_status = erf_sbm::reconstruct_bin(
                layout, 0, bin, host_state, host_reconstruction);
            ASSERT_EQ(core_status, host_status);
            EXPECT_EQ(core_reconstruction.empty, host_reconstruction.empty);
            EXPECT_EQ(core_reconstruction.normalized_roundoff, host_reconstruction.normalized_roundoff);
            EXPECT_DOUBLE_EQ(core_reconstruction.number, host_reconstruction.number);
            EXPECT_DOUBLE_EQ(core_reconstruction.particle_mass, host_reconstruction.particle_mass);
            ASSERT_EQ(host_reconstruction.property_per_particle.size(), 2u);
            for (int property = 0; property < 2; ++property) {
                EXPECT_DOUBLE_EQ(core_properties[property],
                    host_reconstruction.property_per_particle[static_cast<std::size_t>(property)]);
            }

            Real core_integrated_properties[2] = {Real(-1.0), Real(-1.0)};
            erf_sbm::IntegratedMomentsCoreResult core_integral;
            ASSERT_TRUE(erf_sbm::integrate_interval_core(
                core_reconstruction, Real(0.0), Real(3.0), true,
                core_integrated_properties, 2, core_integral));
            erf_sbm::IntegratedMoments host_integral;
            ASSERT_TRUE(erf_sbm::integrate_interval(
                host_reconstruction, Real(0.0), Real(3.0), true, host_integral));
            EXPECT_DOUBLE_EQ(core_integral.number, host_integral.number);
            EXPECT_DOUBLE_EQ(core_integral.water_mass, host_integral.water_mass);
            for (int property = 0; property < 2; ++property) {
                EXPECT_DOUBLE_EQ(core_integrated_properties[property],
                    host_integral.attached_properties[static_cast<std::size_t>(property)]);
            }

            ASSERT_TRUE(erf_sbm::project_reconstruction_bin_core(
                population_view, core_reconstruction, core_projected.data(),
                static_cast<int>(core_projected.size())));
            ASSERT_TRUE(erf_sbm::project_reconstruction_bin(
                layout, host_reconstruction, host_projected));
        }
        EXPECT_EQ(core_projected, core_state);
        EXPECT_EQ(host_projected, host_state);
    }
}

#if defined(AMREX_USE_GPU)
void launch_device_remap_kernel(const amrex::Box& box,
                               const erf_sbm::PopulationRemapView population_view,
                               const Real packet_mass,
                               const Real* edges,
                               Real* state,
                               Real* output)
{
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int, int, int) noexcept {
        const auto plan = erf_sbm::plan_packet_routing(
            population_view, Real(1.0), packet_mass);
        const auto applied = erf_sbm::apply_packet_routing_core(
            population_view, plan, nullptr, 0, nullptr, 0, state, 4);
        erf_sbm::EndpointTransform endpoints;
        const bool endpoint_valid = erf_sbm::try_two_moment_to_endpoints(
            state[3], state[1], edges[1], edges[2], endpoints);
        output[0] = static_cast<Real>(static_cast<int>(plan.status));
        output[1] = static_cast<Real>(plan.destination_count);
        output[2] = static_cast<Real>(plan.destinations[0].bin);
        output[3] = plan.normalized_roundoff ? Real(1.0) : Real(0.0);
        output[4] = static_cast<Real>(static_cast<int>(applied.status));
        output[5] = state[1];
        output[6] = state[3];
        output[7] = endpoint_valid ? Real(1.0) : Real(0.0);
        output[8] = endpoints.normalized_roundoff ? Real(1.0) : Real(0.0);
        output[9] = applied.roundoff_water_mass_correction;
    });
}

TEST(SBMRemapping, AllocationFreeTwoMomentCoreRunsInDeviceKernel)
{
    const std::vector<Real> host_edges{Real(0.5), Real(1.5), Real(2.5)};
    const Real packet_mass = std::nextafter(
        Real(2.5), std::numeric_limits<Real>::infinity());
    amrex::Gpu::DeviceVector<Real> device_edges(host_edges.size());
    amrex::Gpu::DeviceVector<Real> device_state(4, Real(0.0));
    amrex::Gpu::DeviceVector<Real> device_output(10, Real(-1.0));
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_edges.begin(), host_edges.end(),
                     device_edges.begin());
    erf_sbm::PopulationRemapView population_view;
    population_view.population_id = 7;
    population_view.moment_mode = erf_sbm::MomentMode::TwoMoment;
    population_view.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    population_view.edges = device_edges.data();
    population_view.nbins = 2;
    population_view.mass_offset = 0;
    population_view.number_offset = 2;
    population_view.state_components = 4;
    const auto* edges = device_edges.data();
    auto* state = device_state.data();
    auto* output = device_output.data();
    const amrex::Box box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    launch_device_remap_kernel(box, population_view, packet_mass, edges, state, output);
    amrex::Gpu::streamSynchronize();
    std::vector<Real> result(10, Real(0.0));
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_output.begin(), device_output.end(),
                     result.begin());
    EXPECT_EQ(static_cast<int>(result[0]), static_cast<int>(erf_sbm::RemapStatus::Ok));
    EXPECT_DOUBLE_EQ(result[1], Real(1.0));
    EXPECT_DOUBLE_EQ(result[2], Real(1.0));
    EXPECT_DOUBLE_EQ(result[3], Real(1.0));
    EXPECT_EQ(static_cast<int>(result[4]), static_cast<int>(erf_sbm::RemapStatus::Ok));
    EXPECT_DOUBLE_EQ(result[5], Real(2.5));
    EXPECT_DOUBLE_EQ(result[6], Real(1.0));
    EXPECT_DOUBLE_EQ(result[7], Real(1.0));
    EXPECT_DOUBLE_EQ(result[8], Real(0.0));
    EXPECT_DOUBLE_EQ(result[9], Real(2.5) - packet_mass);
    EXPECT_LT(result[9], Real(0.0));
}
#endif

TEST(SBMRemapping, EndpointAndConstraintAdmissibilityAgreeForRoundoffStates)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    const auto group = std::find_if(groups.begin(), groups.end(),
        [](const erf_sbm::ConstraintGroup& candidate) {
            return candidate.population_id == 0 && candidate.bin == 1;
        });
    ASSERT_NE(group, groups.end());
    const auto& pop = population(layout, 0);
    const Real lower = pop.grid.edges()[1];
    const Real upper = pop.grid.edges()[2];
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    state[static_cast<std::size_t>(pop.number_offset + 1)] = Real(1.0);

    for (const Real edge_mass : {lower, upper}) {
        state[static_cast<std::size_t>(pop.mass_offset + 1)] = edge_mass;
        erf_sbm::EndpointTransform endpoints;
        EXPECT_TRUE(erf_sbm::try_two_moment_to_endpoints(
            Real(1.0), edge_mass, lower, upper, endpoints));
        EXPECT_TRUE(group->admissible(state));
    }

    const Real just_above_upper = std::nextafter(
        upper, std::numeric_limits<Real>::infinity());
    state[static_cast<std::size_t>(pop.mass_offset + 1)] = Real(0.0);
    state[static_cast<std::size_t>(pop.number_offset + 1)] = Real(0.0);
    const auto plan = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 0), Real(1.0), just_above_upper);
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_TRUE(plan.normalized_roundoff);
    ASSERT_EQ(plan.destinations[0].bin, 1);
    EXPECT_DOUBLE_EQ(plan.packet_particle_mass, just_above_upper);
    EXPECT_DOUBLE_EQ(plan.destinations[0].particle_mass, upper);
    const auto applied = erf_sbm::apply_packet_routing_core(
        erf_sbm::population_remap_view(layout, 0), plan, nullptr, 0, nullptr, 0,
        state.data(), static_cast<int>(state.size()));
    ASSERT_EQ(applied.status, erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(applied.roundoff_water_mass_correction, upper - just_above_upper);
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(pop.mass_offset + 1)], upper);
    erf_sbm::EndpointTransform normalized;
    EXPECT_TRUE(erf_sbm::try_two_moment_to_endpoints(
        state[static_cast<std::size_t>(pop.number_offset + 1)],
        state[static_cast<std::size_t>(pop.mass_offset + 1)], lower, upper, normalized));
    EXPECT_FALSE(normalized.normalized_roundoff);
    EXPECT_TRUE(group->admissible(state));

    state[static_cast<std::size_t>(pop.mass_offset + 1)] = Real(3.0);
    EXPECT_FALSE(erf_sbm::try_two_moment_to_endpoints(
        Real(1.0), Real(3.0), lower, upper, normalized));
    EXPECT_FALSE(group->admissible(state));
}

TEST(SBMRemapping, ProductionSparseConstraintFlattenerSupportsThreeTerms)
{
    erf_sbm::ConstraintGroup synthetic;
    synthetic.population_id = 12;
    synthetic.bin = 3;
    synthetic.members = {0, 1, 2};
    synthetic.constraints.push_back({"total_minus_coating_minus_rime",
        {{0, Real(1.0)}, {1, Real(-1.0)}, {2, Real(-1.0)}}});
    const auto flattened = erf_sbm::flatten_constraint_groups({synthetic}, 3);
    ASSERT_EQ(flattened.constraints.size(), 1u);
    ASSERT_EQ(flattened.terms.size(), 3u);
    const auto& descriptor = flattened.constraints[0];
    EXPECT_EQ(descriptor.term_offset, 0);
    EXPECT_EQ(descriptor.term_count, 3);
    EXPECT_EQ(descriptor.population_id, 12);
    EXPECT_EQ(descriptor.bin, 3);
    const Real state[] = {Real(10.0), Real(3.0), Real(2.0)};
    const Real independent = state[0] - state[1] - state[2];
    EXPECT_DOUBLE_EQ(erf_sbm::evaluate_constraint_descriptor(
        descriptor, flattened.terms.data(), state), independent);

    auto truncated = descriptor;
    truncated.term_count = 2;
    EXPECT_NE(erf_sbm::evaluate_constraint_descriptor(
                  truncated, flattened.terms.data(), state), independent);

    auto out_of_bounds = synthetic;
    out_of_bounds.constraints[0].terms[2].component = 3;
    EXPECT_THROW(static_cast<void>(erf_sbm::flatten_constraint_groups(
                     {out_of_bounds}, 3)), std::invalid_argument);

    const auto real_layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    const auto production = erf_sbm::make_constraint_descriptors(real_layout);
    ASSERT_FALSE(production.constraints.empty());
    ASSERT_FALSE(production.terms.empty());
    for (const auto& constraint : production.constraints) {
        EXPECT_GT(constraint.term_count, 0);
        EXPECT_GE(constraint.term_offset, 0);
        EXPECT_LE(constraint.term_offset + constraint.term_count,
                  static_cast<int>(production.terms.size()));
    }
}

TEST(SBMRemapping, UnequalSecondPopulationUsesItsGridOffsetsAndProperties)
{
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(make_population(0, erf_sbm::MomentMode::OneMoment));
    spec.populations.push_back(make_population(1, erf_sbm::MomentMode::TwoMoment,
        {Real(0.1), Real(0.4), Real(1.0), Real(2.0)},
        {Real(0.2), Real(0.7), Real(1.5)}, erf_sbm::PopulationPhase::Aerosol));
    spec.liquid_projection = {0, 1};
    spec.attached_properties.push_back(make_property("second_population_material", 1));
    const erf_sbm::SBMLayout layout(std::move(spec));
    const auto& first = population(layout, 0);
    const auto& second = population(layout, 1);
    ASSERT_EQ(first.grid.nbins(), 2);
    ASSERT_EQ(second.grid.nbins(), 3);
    ASSERT_EQ(second.mass_offset, 2);
    ASSERT_EQ(second.number_offset, 5);
    ASSERT_EQ(layout.property_offset(0), 8);

    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto plan = erf_sbm::plan_packet_routing(
        erf_sbm::population_remap_view(layout, 1), Real(4.0), Real(0.6));
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(plan.destination_count, 1);
    EXPECT_EQ(plan.destinations[0].bin, 1);
    ASSERT_EQ(erf_sbm::apply_packet_routing(layout, plan, {Real(7.0)}, state).status,
              erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(state[0], Real(0.0));
    EXPECT_DOUBLE_EQ(state[1], Real(0.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(second.mass_offset + 1)], Real(2.4));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(second.number_offset + 1)], Real(4.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(layout.property_offset(0) + 1)], Real(7.0));

    erf_sbm::ReconstructionDelta reconstructed;
    ASSERT_EQ(erf_sbm::reconstruct_bin(layout, 1, 1, state, reconstructed),
              erf_sbm::ReconstructionStatus::Populated);
    EXPECT_DOUBLE_EQ(reconstructed.particle_mass, Real(0.6));
    EXPECT_DOUBLE_EQ(reconstructed.number, Real(4.0));
    EXPECT_DOUBLE_EQ(reconstructed.property_per_particle[0], Real(1.75));
}

TEST(SBMRemapping, ScientificIdentityAndRestartAreExactAndVersioned)
{
    const auto one_moment = make_layout(erf_sbm::MomentMode::OneMoment);
    const auto two_moment = make_layout(erf_sbm::MomentMode::TwoMoment);
    EXPECT_EQ(one_moment.schema_identity(), one_moment.schema_identity());
    EXPECT_NE(one_moment.schema_identity(), two_moment.schema_identity());
    EXPECT_NE(one_moment.schema_identity().find("fixed-pivot-1m-v1"), std::string::npos);
    EXPECT_NE(one_moment.schema_identity().find("fixed-pivot-delta-v1"), std::string::npos);
    EXPECT_NE(one_moment.schema_identity().find("fixed-pivot-two-center-v1"), std::string::npos);
    EXPECT_NE(one_moment.schema_identity().find("sbm-layout-v1"), std::string::npos);
    EXPECT_EQ(one_moment.schema_identity().find("M2R"), std::string::npos);
    EXPECT_EQ(one_moment.schema_identity().find("m2r"), std::string::npos);
    EXPECT_EQ(one_moment.schema_identity().find("M2-R"), std::string::npos);
    EXPECT_NE(two_moment.schema_identity().find("interval-2m-v1"), std::string::npos);
    EXPECT_NE(two_moment.schema_identity().find("mean-delta-2m-v1"), std::string::npos);
    EXPECT_NE(two_moment.schema_identity().find("actual-mass-interval-v1"), std::string::npos);
    EXPECT_NE(one_moment.inspection().find("reconstruction=fixed-pivot-delta-v1"), std::string::npos);

    const std::string persisted = erf_sbm::restart_schema(one_moment);
    EXPECT_NE(persisted.find("ERF-SBM-RESTART-v1"), std::string::npos);
    EXPECT_EQ(persisted.find("M2R"), std::string::npos);
    EXPECT_EQ(persisted.find("m2r"), std::string::npos);
    EXPECT_EQ(persisted.find("M2-R"), std::string::npos);
    EXPECT_EQ(persisted.find("population=0:representation="), std::string::npos);
    const auto layout_representation = persisted.find(":representation=");
    ASSERT_NE(layout_representation, std::string::npos);
    EXPECT_EQ(persisted.find(":representation=", layout_representation + 1), std::string::npos);
    EXPECT_TRUE(erf_sbm::restart_schema_matches(one_moment, persisted));
    std::string changed = persisted;
    const auto identity = changed.find("fixed-pivot-delta-v1");
    ASSERT_NE(identity, std::string::npos);
    changed.replace(identity, std::string("fixed-pivot-delta-v1").size(), "different-delta-v9");
    EXPECT_FALSE(erf_sbm::restart_schema_matches(one_moment, changed));
    EXPECT_EQ(one_moment.ncomp(), 2);
    EXPECT_EQ(two_moment.ncomp(), 4);
}

TEST(SBMRemapping, OneMomentRejectsCoincidentPivotsWithoutRestrictingTwoMoment)
{
    const std::vector<Real> edges{Real(0.5), Real(1.5), Real(2.5)};
    const std::vector<Real> coincident{Real(1.5), Real(1.5)};
    auto one_moment = make_population(0, erf_sbm::MomentMode::OneMoment, edges, coincident);
    erf_sbm::SBMLayoutSpec one_spec;
    one_spec.populations.push_back(one_moment);
    one_spec.liquid_projection = {0, 1};
    EXPECT_FALSE(erf_sbm::SBMLayout::validate(one_spec).valid);

    auto two_moment = make_population(0, erf_sbm::MomentMode::TwoMoment, edges, coincident);
    erf_sbm::SBMLayoutSpec two_spec;
    two_spec.populations.push_back(two_moment);
    two_spec.liquid_projection = {0, 1};
    EXPECT_TRUE(erf_sbm::SBMLayout::validate(two_spec).valid);
}

} // namespace
