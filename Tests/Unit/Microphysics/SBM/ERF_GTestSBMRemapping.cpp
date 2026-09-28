#include <gtest/gtest.h>

#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMRemapping.H"
#include "ERF_SBMRepresentation.H"
#include "ERF_SBMRestart.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
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

erf_sbm::SpectralPopulationSpec make_population(
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

erf_sbm::AttachedPropertyDescriptor make_property(const std::string& name,
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

erf_sbm::SBMLayout make_layout(
    const erf_sbm::MomentMode mode,
    std::vector<erf_sbm::AttachedPropertyDescriptor> properties = {})
{
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(make_population(0, mode));
    spec.liquid_projection = {0, 1};
    spec.attached_properties = std::move(properties);
    return erf_sbm::SBMLayout(std::move(spec));
}

erf_sbm::SpectralGridView grid_view(const erf_sbm::PopulationLayout& population)
{
    return {population.grid.edges().data(), population.grid.pivots().data(),
            population.grid.nbins()};
}

const erf_sbm::PopulationLayout& population(const erf_sbm::SBMLayout& layout,
                                           const int id)
{
    const auto found = std::find_if(layout.populations().begin(), layout.populations().end(),
        [id](const erf_sbm::PopulationLayout& value) { return value.population_id == id; });
    if (found == layout.populations().end()) throw std::runtime_error("missing test population");
    return *found;
}

void expect_close(const Real actual, const Real expected, const Real scale = Real(1.0))
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
        grid_view(liquid), liquid.moment_mode, packet_number, packet_mass);
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(plan.destination_count, 2);
    EXPECT_EQ(plan.destinations[0].bin, 0);
    EXPECT_EQ(plan.destinations[1].bin, 1);
    expect_close(plan.destinations[0].number_weight, Real(0.6));
    expect_close(plan.destinations[1].number_weight, Real(0.4));

    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto result = erf_sbm::apply_packet_routing(
        layout, 0, plan, {Real(6.0), Real(15.0)}, state);
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
    const auto& liquid = population(layout, 0);
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));

    const auto exact = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(4.0), Real(2.0));
    ASSERT_EQ(exact.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(exact.destination_count, 1);
    EXPECT_EQ(exact.destinations[0].bin, 1);
    EXPECT_DOUBLE_EQ(exact.destinations[0].number_weight, Real(1.0));
    ASSERT_EQ(erf_sbm::apply_packet_routing(layout, 0, exact,
        {Real(8.0), Real(12.0)}, state).status, erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(state[0], Real(0.0));
    EXPECT_DOUBLE_EQ(state[1], Real(8.0));

    const auto subpivot = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(5.0), Real(0.4));
    ASSERT_EQ(subpivot.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    ASSERT_EQ(subpivot.destination_count, 1);
    expect_close(subpivot.destinations[0].number_weight, Real(0.4));
    expect_close(subpivot.residual_number_weight, Real(0.6));
    const auto before_subpivot = state;
    const auto residual = erf_sbm::apply_packet_routing(
        layout, 0, subpivot, {Real(10.0), Real(20.0)}, state);
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
        grid_view(liquid), liquid.moment_mode, Real(2.0), Real(0.0));
    ASSERT_EQ(zero_water.status, erf_sbm::RemapStatus::ZeroWaterResidual);
    EXPECT_EQ(zero_water.destination_count, 0);
    const auto all_residual = erf_sbm::apply_packet_routing(
        layout, 0, zero_water, {Real(3.0), Real(7.0)}, state);
    expect_close(all_residual.residual_number, Real(2.0));
    EXPECT_EQ(all_residual.residual_properties, (std::vector<Real>{Real(3.0), Real(7.0)}));
    EXPECT_EQ(all_residual.residual_water_mass, Real(0.0));

    const auto overflow = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0), Real(2.1));
    ASSERT_EQ(overflow.status, erf_sbm::RemapStatus::Overflow);
    const auto before_overflow = state;
    EXPECT_EQ(erf_sbm::apply_packet_routing(
        layout, 0, overflow, {Real(1.0), Real(1.0)}, state).status,
        erf_sbm::RemapStatus::Overflow);
    EXPECT_EQ(state, before_overflow);

    const auto invalid = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode,
        std::numeric_limits<Real>::quiet_NaN(), Real(1.0));
    EXPECT_EQ(invalid.status, erf_sbm::RemapStatus::Invalid);
    const auto atomic_plan = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0), Real(1.4));
    auto atomic_state = std::vector<Real>(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    atomic_state[static_cast<std::size_t>(layout.property_offset(0) + 1)] =
        std::numeric_limits<Real>::max();
    const auto unchanged = atomic_state;
    const auto failed_apply = erf_sbm::apply_packet_routing(
        layout, 0, atomic_plan, {std::numeric_limits<Real>::max(), Real(0.0)}, atomic_state);
    EXPECT_EQ(failed_apply.status, erf_sbm::RemapStatus::Invalid);
    EXPECT_EQ(atomic_state, unchanged);
}

TEST(SBMRemapping, IntervalTwoMomentPreservesActualMassAndUsesHalfOpenEdges)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment,
        {make_property("coating", 0), make_property("rime", 0)});
    const auto& liquid = population(layout, 0);

    const auto packet = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0), Real(1.4));
    ASSERT_EQ(packet.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(packet.destination_count, 1);
    EXPECT_EQ(packet.destinations[0].bin, 0);
    EXPECT_DOUBLE_EQ(packet.destinations[0].particle_mass, Real(1.4));
    std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    ASSERT_EQ(erf_sbm::apply_packet_routing(
        layout, 0, packet, {Real(2.0), Real(5.0)}, state).status,
        erf_sbm::RemapStatus::Ok);
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.mass_offset)], Real(1.4));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.number_offset)], Real(1.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.mass_offset + 1)], Real(0.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(liquid.number_offset + 1)], Real(0.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(layout.property_offset(0))], Real(2.0));
    EXPECT_DOUBLE_EQ(state[static_cast<std::size_t>(layout.property_offset(1))], Real(5.0));

    const Real shared = Real(1.5);
    const auto at_shared_edge = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0), shared);
    ASSERT_EQ(at_shared_edge.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(at_shared_edge.destinations[0].bin, 1);
    const auto below_shared = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0), std::nextafter(shared, Real(0.0)));
    const auto above_shared = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0),
        std::nextafter(shared, std::numeric_limits<Real>::infinity()));
    EXPECT_EQ(below_shared.destinations[0].bin, 0);
    EXPECT_EQ(above_shared.destinations[0].bin, 1);

    const auto at_lower = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0), Real(0.5));
    const auto at_upper = erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0), Real(2.5));
    EXPECT_EQ(at_lower.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(at_lower.destinations[0].bin, 0);
    EXPECT_EQ(at_upper.status, erf_sbm::RemapStatus::Ok);
    EXPECT_EQ(at_upper.destinations[0].bin, 1);
    EXPECT_EQ(erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0),
        std::nextafter(Real(0.5), Real(0.0))).status,
        erf_sbm::RemapStatus::BelowSupportedGrid);
    EXPECT_EQ(erf_sbm::plan_packet_routing(
        grid_view(liquid), liquid.moment_mode, Real(1.0),
        std::nextafter(Real(2.5), std::numeric_limits<Real>::infinity())).status,
        erf_sbm::RemapStatus::Overflow);
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
        ASSERT_NEAR(endpoints.L, Real(0.2), Real(2.0e-14));
        ASSERT_NEAR(endpoints.H, Real(1.8), Real(2.0e-14));
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
        EXPECT_NEAR(roundoff.H, count, Real(2.0e-14));
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
        const Real negative_endpoint_number = Real(1.0e-12);
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

TEST(SBMRemapping, ProductionSparseConstraintFlattenerSupportsThreeTerms)
{
    erf_sbm::ConstraintGroup synthetic;
    synthetic.population_id = 12;
    synthetic.bin = 3;
    synthetic.members = {0, 1, 2};
    synthetic.constraints.push_back({"total_minus_coating_minus_rime",
        {{0, Real(1.0)}, {1, Real(-1.0)}, {2, Real(-1.0)}}});
    const auto flattened = erf_sbm::flatten_constraint_groups({synthetic});
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
        grid_view(second), second.moment_mode, Real(4.0), Real(0.6));
    ASSERT_EQ(plan.status, erf_sbm::RemapStatus::Ok);
    ASSERT_EQ(plan.destination_count, 1);
    EXPECT_EQ(plan.destinations[0].bin, 1);
    ASSERT_EQ(erf_sbm::apply_packet_routing(layout, 1, plan, {Real(7.0)}, state).status,
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
    EXPECT_NE(two_moment.schema_identity().find("interval-2m-v1"), std::string::npos);
    EXPECT_NE(two_moment.schema_identity().find("mean-delta-2m-v1"), std::string::npos);
    EXPECT_NE(two_moment.schema_identity().find("actual-mass-interval-v1"), std::string::npos);
    EXPECT_NE(one_moment.inspection().find("reconstruction=fixed-pivot-delta-v1"), std::string::npos);

    const std::string persisted = erf_sbm::restart_schema(one_moment);
    EXPECT_NE(persisted.find("ERF-SBM-RESTART-M2R-v1"), std::string::npos);
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
