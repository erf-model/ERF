#include <gtest/gtest.h>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_MultiFab.H>
#include <AMReX_RealBox.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "ERF_AuxiliaryProjection.H"
#include "ERF_AuxiliaryStateLayout.H"
#include "ERF_AuxiliaryStateManager.H"
#include "ERF_AuxiliaryStageContext.H"
#include "ERF_IndexDefines.H"
#include "ERF_SBMContracts.H"
#include "ERF_SBMBulkProjection.H"
#include "ERF_SBMLayout.H"
#include "ERF_SBMTransferClosure.H"
#include "ERF_SBMTransportPrototype.H"
#include "ERF_SpectralGrid.H"

namespace {

using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::Real;

erf_sbm::SpectralGridSpec make_grid (const int nbins)
{
    erf_sbm::SpectralGridSpec spec;
    spec.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    spec.coordinate_units = "kg";
    spec.edges.resize(static_cast<std::size_t>(nbins + 1));
    spec.pivots.resize(static_cast<std::size_t>(nbins));
    for (int n = 0; n <= nbins; ++n) spec.edges[static_cast<std::size_t>(n)] = Real(n);
    for (int n = 0; n < nbins; ++n) {
        spec.pivots[static_cast<std::size_t>(n)] = Real(n) + Real(0.5);
    }
    return spec;
}

erf_sbm::SpectralPopulationSpec make_population (
    const int nbins, const int population = 0,
    const erf_sbm::PopulationPhase phase = erf_sbm::PopulationPhase::Liquid)
{
    erf_sbm::SpectralPopulationSpec spec;
    spec.population_id = population;
    spec.semantic_id = population == 0 ? "liquid_mass" : "aerosol_mass";
    spec.phase = phase;
    spec.grid = make_grid(nbins);
    spec.moment_mode = erf_sbm::MomentMode::OneMoment;
    spec.mass_state_units = "kg m^-3";
    spec.number_state_units = "m^-3";
    return spec;
}

erf_sbm::SBMLayout make_layout (const int nbins,
                                const erf_sbm::MomentMode mode = erf_sbm::MomentMode::OneMoment)
{
    erf_sbm::SBMLayoutSpec spec;
    auto population = make_population(nbins);
    population.moment_mode = mode;
    spec.populations.push_back(std::move(population));
    spec.liquid_projection = {0, nbins / 2};
    return erf_sbm::SBMLayout(std::move(spec));
}

Geometry make_geometry (const Box& domain)
{
    const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                   {AMREX_D_DECL(4.0, 2.0, 2.0)});
    amrex::Array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    return Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodicity.data());
}

TEST (SBMP0, SpectralGridValidationAndRuntimeSizes)
{
    const auto valid = erf_sbm::SpectralGrid::validate(make_grid(4));
    EXPECT_TRUE(valid.valid) << valid.message;

    auto negative = make_grid(4);
    negative.edges[1] = Real(-1.0);
    EXPECT_FALSE(erf_sbm::SpectralGrid::validate(negative).valid);

    auto repeated = make_grid(4);
    repeated.edges[2] = repeated.edges[1];
    EXPECT_FALSE(erf_sbm::SpectralGrid::validate(repeated).valid);

    auto nonfinite = make_grid(4);
    nonfinite.edges[2] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_FALSE(erf_sbm::SpectralGrid::validate(nonfinite).valid);

    auto radius = make_grid(4);
    radius.coordinate_kind = erf_sbm::CoordinateKind::Radius;
    radius.coordinate_units = "m";
    EXPECT_TRUE(erf_sbm::SpectralGrid::validate(radius).valid);

    auto bad_split = erf_sbm::SBMLayoutSpec{};
    bad_split.populations.push_back(make_population(4));
    bad_split.liquid_projection = {0, 0};
    EXPECT_FALSE(erf_sbm::SBMLayout::validate(bad_split).valid);
    bad_split.liquid_projection.cloud_rain_split = 4;
    EXPECT_FALSE(erf_sbm::SBMLayout::validate(bad_split).valid);
    EXPECT_FALSE(erf_sbm::validate_runtime_bin_count(-1).empty());
    EXPECT_FALSE(erf_sbm::validate_runtime_bin_count(1000001).empty());
    EXPECT_TRUE(erf_sbm::validate_runtime_bin_count(4).empty());

    for (const int nbins : {4, 16, 64}) {
        const auto one_moment = make_layout(nbins);
        const auto two_moment = make_layout(nbins, erf_sbm::MomentMode::TwoMoment);
        EXPECT_EQ(one_moment.ncomp(), nbins);
        EXPECT_EQ(two_moment.ncomp(), 2 * nbins);
        EXPECT_EQ(one_moment.auxiliary_layout().ncomp(), one_moment.ncomp());
        EXPECT_NE(one_moment.schema_identity(), two_moment.schema_identity());
        EXPECT_EQ(make_layout(nbins).schema_identity(), one_moment.schema_identity());
    }
}

TEST (SBMP0, TwoMomentEndpointAlgebraAndDistinctSemantics)
{
    constexpr Real lower = Real(2.0);
    constexpr Real upper = Real(6.0);
    for (const auto pair : {std::pair<Real, Real>{Real(3.0), Real(12.0)},
                            std::pair<Real, Real>{Real(2.0), Real(4.0)},
                            std::pair<Real, Real>{Real(2.0), Real(12.0)}}) {
        const auto endpoints = erf_sbm::SpectralGrid::two_moment_to_endpoints(
            pair.first, pair.second, lower, upper);
        const auto inverse = erf_sbm::SpectralGrid::endpoints_to_two_moment(
            endpoints.first, endpoints.second, lower, upper);
        EXPECT_NEAR(inverse.first, pair.first, 1.e-14);
        EXPECT_NEAR(inverse.second, pair.second, 1.e-14);
    }

    const auto empty = erf_sbm::SpectralGrid::two_moment_to_endpoints(Real(0.0), Real(0.0), lower, upper);
    EXPECT_EQ(empty.first, Real(0.0));
    EXPECT_EQ(empty.second, Real(0.0));
    EXPECT_TRUE(erf_sbm::SpectralGrid::two_moment_realizable(Real(1.0), lower, lower, upper));
    EXPECT_TRUE(erf_sbm::SpectralGrid::two_moment_realizable(Real(1.0), upper, lower, upper));
    EXPECT_FALSE(erf_sbm::SpectralGrid::two_moment_realizable(Real(1.0), Real(7.0), lower, upper));
    EXPECT_THROW((void)erf_sbm::SpectralGrid::two_moment_to_endpoints(Real(1.0), Real(7.0), lower, upper),
                 std::invalid_argument);

    const auto one_moment = make_layout(4, erf_sbm::MomentMode::OneMoment);
    const auto two_moment = make_layout(4, erf_sbm::MomentMode::TwoMoment);
    EXPECT_EQ(one_moment.populations().front().number_offset, -1);
    EXPECT_GE(two_moment.populations().front().number_offset, 0);
}

TEST (SBMP0, ProjectionPartitionAndInvalidComponentClasses)
{
    const auto layout = make_layout(4);
    std::vector<Real> bins{Real(1.0), Real(2.0), Real(4.0), Real(8.0)};
    const auto bulk = erf_sbm::SBMBulkProjection(layout).apply(bins);
    EXPECT_EQ(bulk.qc, Real(3.0));
    EXPECT_EQ(bulk.qr, Real(12.0));

    const auto& projection = layout.bulk_projection();
    ASSERT_EQ(projection.rules().size(), 2U);
    EXPECT_EQ(projection.rules()[0].source_begin, 0);
    EXPECT_EQ(projection.rules()[0].source_count, 2);
    EXPECT_EQ(projection.rules()[1].source_begin, 2);
    EXPECT_EQ(projection.rules()[1].source_count, 2);
    EXPECT_TRUE(projection.validate(layout.ncomp()).valid);

    erf_auxiliary::AuxiliaryFaceTransfer spectral_faces;
    erf_auxiliary::AuxiliaryFaceTransfer bulk_faces;
    const Box face_domain(IntVect(0, 0, 0), IntVect(1, 1, 1));
    const BoxArray face_boxes(face_domain);
    const DistributionMapping face_dm(face_boxes);
    spectral_faces.define(face_boxes, face_dm, 4, 0);
    bulk_faces.define(face_boxes, face_dm, 2, 0);
    spectral_faces.setVal(Real(1.0));
    erf_sbm::SBMBulkProjection(layout).apply_to_face_transfer(spectral_faces, bulk_faces);
    EXPECT_DOUBLE_EQ(erf_sbm::SBMBulkProjection(layout).max_face_projection_error(
                         spectral_faces, bulk_faces), Real(0.0));
    // Negative control: an independently altered compact face flux no longer
    // satisfies the accepted projection invariant.
    bulk_faces.x().setVal(Real(0.0), 0, 1, 0);
    EXPECT_GT(erf_sbm::SBMBulkProjection(layout).max_face_projection_error(
                  spectral_faces, bulk_faces), Real(0.0));

    erf_auxiliary::AuxiliaryProjection duplicate({{"qc", 0, 1}, {"qc", 1, 1}});
    EXPECT_FALSE(duplicate.validate(4).valid);
    erf_auxiliary::AuxiliaryProjection overlap({{"qc", 0, 2}, {"qr", 1, 2}});
    EXPECT_FALSE(overlap.validate(4).valid);
    erf_auxiliary::AuxiliaryProjection vapor({{"qv", 0, 1,
        erf_auxiliary::ProjectionTargetKind::Vapor}});
    EXPECT_FALSE(vapor.validate(4).valid);
    erf_auxiliary::AuxiliaryProjection number({{"number", 0, 1,
        erf_auxiliary::ProjectionTargetKind::BulkCoupling,
        erf_auxiliary::ProjectionSourceKind::AuxiliaryNumber}});
    EXPECT_FALSE(number.validate(4).valid);
}

TEST (SBMP0, AttachedPropertiesAndSecondPopulationAreExtensible)
{
    erf_sbm::SBMLayoutSpec spec;
    spec.populations = {make_population(4, 0, erf_sbm::PopulationPhase::Liquid),
                        make_population(3, 1, erf_sbm::PopulationPhase::Aerosol)};
    spec.liquid_projection = {0, 2};
    spec.attached_properties.push_back({"dry_solute", "attached_dry_solute", "kg m^-3", 0,
                                        erf_sbm::PropertyKind::ExtensiveMass,
                                        erf_sbm::SupportRequirement::PositiveMass,
                                        erf_sbm::PropertyRemapPolicy::CarrierBinConservative,
                                        true, true});
    spec.attached_properties.push_back({"core_number", "aerosol_core_number", "m^-3", 1,
                                        erf_sbm::PropertyKind::NumberCarried,
                                        erf_sbm::SupportRequirement::PositiveNumber,
                                        erf_sbm::PropertyRemapPolicy::CarrierBinConservative});
    const erf_sbm::SBMLayout layout(std::move(spec));
    EXPECT_EQ(layout.populations()[0].mass_offset, 0);
    EXPECT_EQ(layout.populations()[1].mass_offset, 4);
    EXPECT_EQ(layout.property_offset(0), 7);
    EXPECT_EQ(layout.property_offset(1), 11);
    EXPECT_EQ(layout.ncomp(), 14);
    ASSERT_EQ(layout.auxiliary_layout().components().size(), 14U);
    EXPECT_EQ(layout.attached_properties()[0].kind, erf_sbm::PropertyKind::ExtensiveMass);
    EXPECT_EQ(layout.attached_properties()[0].support, erf_sbm::SupportRequirement::PositiveMass);
    EXPECT_TRUE(layout.attached_properties()[0].may_overlap_mass);
    EXPECT_EQ(layout.attached_properties()[0].remap_policy,
              erf_sbm::PropertyRemapPolicy::CarrierBinConservative);
    EXPECT_EQ(layout.attached_properties()[1].kind, erf_sbm::PropertyKind::NumberCarried);
    EXPECT_EQ(layout.attached_properties()[1].units, "m^-3");
    EXPECT_EQ(layout.populations()[1].phase, erf_sbm::PopulationPhase::Aerosol);
    EXPECT_EQ(layout.auxiliary_layout().components()[4].units, "kg m^-3");
    EXPECT_EQ(layout.auxiliary_layout().components()[7].units, "kg m^-3");
    EXPECT_EQ(layout.auxiliary_layout().components()[11].units, "m^-3");
    EXPECT_NE(layout.inspection().find("property dry_solute"), std::string::npos);
    EXPECT_NE(layout.inspection().find("remap="), std::string::npos);

    std::vector<Real> with_properties(14, Real(0.0));
    with_properties[0] = Real(1.0);
    with_properties[1] = Real(2.0);
    with_properties[2] = Real(4.0);
    with_properties[3] = Real(8.0);
    with_properties[7] = Real(100.0); // extensive dry mass is not water
    with_properties[11] = Real(200.0); // carried number is not water
    const auto projected = erf_sbm::SBMBulkProjection(layout).apply(with_properties);
    EXPECT_EQ(projected.qc, Real(3.0));
    EXPECT_EQ(projected.qr, Real(12.0));
}

TEST (SBMP0, GenericManagerPreservesBaselineAndPublishesAcceptedState)
{
    const erf_auxiliary::AuxiliaryStateLayout state_layout(
        "manager-test", {{"x", "test", "1"}, {"y", "test", "1"}});
    erf_auxiliary::AuxiliaryStateManager manager(state_layout);
    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    manager.define_level(0, boxes, dm, 2);
    manager.output(0).setVal(Real(2.0));
    manager.begin_step(0);
    manager.output(0).setVal(Real(5.0));
    manager.accept_stage(0);
    EXPECT_EQ(manager.old(0).min(0), Real(2.0));
    EXPECT_EQ(manager.evaluation(0).min(0), Real(5.0));
    EXPECT_EQ(manager.output(0).min(0), Real(5.0));
    const auto one_state_bytes = erf_auxiliary::allocated_payload_bytes(manager.old(0));
    EXPECT_EQ(manager.state_resident_bytes(), 4U * one_state_bytes);
    EXPECT_GT(one_state_bytes,
              static_cast<std::size_t>(boxes.numPts()) * 2U * sizeof(Real));
    erf_auxiliary::AuxiliaryFaceTransfer face_transfer;
    face_transfer.define(boxes, dm, 2, 0);
    EXPECT_EQ(manager.face_transfer_resident_bytes(), 2U * face_transfer.resident_bytes());
    EXPECT_GT(manager.resident_bytes(), 0U);
}

TEST (SBMP0, CapabilityReportFailsClosedAndIsStable)
{
    const erf_sbm::CapabilityInput supported_input;
    const auto supported = erf_sbm::evaluate_p1_capabilities(supported_input);
    EXPECT_TRUE(supported.supported);
    EXPECT_EQ(supported.stable_description(),
              erf_sbm::evaluate_p1_capabilities(supported_input).stable_description());

    auto rejected_input = supported_input;
    rejected_input.max_level = 1;
    rejected_input.diffusion = true;
    rejected_input.shoc_or_macrophysics = true;
    rejected_input.moving_terrain = true;
    rejected_input.embedded_boundary = true;
    rejected_input.sedimentation = true;
    rejected_input.two_moment_transport = true;
    rejected_input.condensation = true;
    rejected_input.activation = true;
    rejected_input.collision = true;
    rejected_input.custom_moisture_forcing = true;
    rejected_input.large_scale_forcing = true;
    rejected_input.sounding_nudging = true;
    rejected_input.sponge_or_wall_modification = true;
    rejected_input.periodic_cartesian = false;
    const auto rejected = erf_sbm::evaluate_p1_capabilities(rejected_input);
    EXPECT_FALSE(rejected.supported);
    EXPECT_GE(rejected.rejected_reasons.size(), 10U);

    const auto inspection = erf_sbm::stable_inspection(
        make_layout(4), supported_input, "erf-base", "amrex-base", "design-sha");
    EXPECT_EQ(inspection,
              erf_sbm::stable_inspection(
                  make_layout(4), supported_input, "erf-base", "amrex-base", "design-sha"));
    EXPECT_NE(inspection.find("format=erf-sbm-inspection-v1"), std::string::npos);
    EXPECT_NE(inspection.find("SBM-BULK-PROJECTION"), std::string::npos);
    EXPECT_NE(inspection.find("identity.design=design-sha"), std::string::npos);
}

TEST (SBMP1, ExactReducedRecurrencesAndNegativeControls)
{
    const std::vector<Real> old{Real(10.0), Real(20.0)};
    const std::vector<Real> predictor{Real(12.0), Real(24.0)};
    const std::vector<Real> rhs{Real(3.0), Real(-4.0)};
    const auto comp0 = erf_auxiliary::make_compressible_stage(0, 0.0, 0.0, 1.0/3.0, 1.0, nullptr, nullptr);
    const auto comp1 = erf_auxiliary::make_compressible_stage(1, 0.0, 0.0, 0.5, 1.0, nullptr, nullptr);
    const auto comp2 = erf_auxiliary::make_compressible_stage(2, 0.0, 0.0, 1.0, 1.0, nullptr, nullptr);
    std::vector<Real> output;
    erf_sbm::update_stage(comp0, old, old, rhs, output);
    EXPECT_DOUBLE_EQ(output[0], 11.0);
    erf_sbm::update_stage(comp1, old, predictor, rhs, output);
    EXPECT_DOUBLE_EQ(output[0], 11.5);
    erf_sbm::update_stage(comp2, old, predictor, rhs, output);
    EXPECT_DOUBLE_EQ(output[0], 13.0);
    // Negative control: treating stage 1 as a successive update from the
    // predictor would produce 13.5, which this contract must detect.
    EXPECT_NE(output[0], predictor[0] + Real(0.5) * rhs[0]);

    const auto anelastic0 = erf_auxiliary::make_anelastic_stage(0, 0.0, 0.0, 1.0, 1.0, nullptr, nullptr);
    const auto anelastic1 = erf_auxiliary::make_anelastic_stage(1, 0.0, 1.0, 1.0, 1.0, nullptr, nullptr);
    erf_sbm::update_stage(anelastic0, old, old, rhs, output);
    EXPECT_DOUBLE_EQ(output[0], 13.0);
    erf_sbm::update_stage(anelastic1, old, output, std::vector<Real>{Real(5.0), Real(7.0)}, output);
    EXPECT_DOUBLE_EQ(output[0], 14.0);

    const std::vector<erf_sbm::HostState> stage_fluxes{
        {Real(1.0), Real(2.0)}, {Real(3.0), Real(5.0)}};
    const auto accepted = erf_sbm::accepted_ledger(anelastic1, stage_fluxes, 2.0);
    EXPECT_DOUBLE_EQ(accepted[0], 4.0);
    EXPECT_DOUBLE_EQ(accepted[1], 7.0);
    // Negative control: final-stage-only would be (6,10), not the Heun ledger.
    EXPECT_NE(accepted[0], 2.0 * stage_fluxes[1][0]);

    const auto comp_ledger = erf_sbm::accepted_ledger(comp2,
        {{Real(1.0)}, {Real(3.0)}, {Real(7.0)}}, 2.0);
    EXPECT_DOUBLE_EQ(comp_ledger[0], 14.0);
}

TEST (SBMP1, OwnershipGuardDetectsIndependentProjectedWrite)
{
    const erf_sbm::OwnershipRegistry inactive(false);
    EXPECT_FALSE(inactive.owns_cloud_or_rain(RhoQ2_comp));
    const erf_sbm::OwnershipRegistry active(true);
    for (const auto path : {erf_sbm::NativeWritePath::Advection,
                            erf_sbm::NativeWritePath::Diffusion,
                            erf_sbm::NativeWritePath::Source,
                            erf_sbm::NativeWritePath::PositivityClip,
                            erf_sbm::NativeWritePath::Microphysics,
                            erf_sbm::NativeWritePath::Wall}) {
        EXPECT_TRUE(active.owns(RhoQ2_comp, path));
        EXPECT_TRUE(active.owns(RhoQ3_comp, path));
    }
    EXPECT_FALSE(active.owns(RhoQ1_comp, erf_sbm::NativeWritePath::Advection));

    const Real spectral_projection = Real(3.0);
    const Real independently_written_bulk = spectral_projection + Real(1.0);
    // Negative control: any independent qc/qr write leaves the projection
    // invariant false and is therefore observable rather than clipped away.
    EXPECT_NE(independently_written_bulk, spectral_projection);
}

TEST (SBMP1, AcceptedFaceTransferLedgerUsesTheActualStageFlux)
{
    const Box domain(IntVect(0, 0, 0), IntVect(1, 1, 1));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    erf_auxiliary::AuxiliaryFaceTransferLedger ledger;
    ledger.define(boxes, dm, 1, 3, 0);
    erf_auxiliary::AuxiliaryFaceTransfer stage0;
    erf_auxiliary::AuxiliaryFaceTransfer stage1;
    erf_auxiliary::AuxiliaryFaceTransfer stage2;
    stage0.define(boxes, dm, 1, 0);
    stage1.define(boxes, dm, 1, 0);
    stage2.define(boxes, dm, 1, 0);

    const auto c0 = erf_auxiliary::make_compressible_stage(
        0, 0.0, 0.0, 2.0/3.0, 2.0, nullptr, nullptr);
    const auto c1 = erf_auxiliary::make_compressible_stage(
        1, 0.0, 2.0/3.0, 1.0, 2.0, nullptr, nullptr);
    const auto c2 = erf_auxiliary::make_compressible_stage(
        2, 0.0, 1.0, 2.0, 2.0, nullptr, nullptr);
    stage0.setVal(1.0);
    stage1.setVal(2.0);
    stage2.setVal(3.0);
    ledger.record_stage(c0, stage0);
    ledger.record_stage(c1, stage1);
    ledger.record_stage(c2, stage2);
    EXPECT_DOUBLE_EQ(ledger.accepted().x().min(0), 6.0);
    EXPECT_DOUBLE_EQ(ledger.stage().x().min(0), 3.0);
    EXPECT_EQ(ledger.resident_bytes(), 2U * stage0.resident_bytes());
    EXPECT_EQ(ledger.stage_resident_bytes(), stage0.resident_bytes());
    EXPECT_EQ(ledger.accepted_resident_bytes(), stage0.resident_bytes());

    ledger.begin_step();
    const auto a0 = erf_auxiliary::make_anelastic_stage(
        0, 0.0, 0.0, 2.0, 2.0, nullptr, nullptr);
    const auto a1 = erf_auxiliary::make_anelastic_stage(
        1, 0.0, 2.0, 2.0, 2.0, nullptr, nullptr);
    stage0.setVal(1.0);
    stage1.setVal(3.0);
    ledger.record_stage(a0, stage0);
    ledger.record_stage(a1, stage1);
    EXPECT_DOUBLE_EQ(ledger.accepted().y().min(0), 4.0);
}

TEST (SBMP1, AcceptedTransferClosesAgainstAcceptedState)
{
    const auto layout = make_layout(4);
    const Box domain(IntVect(0, 0, 0), IntVect(3, 0, 0));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    const Geometry geometry = make_geometry(domain);
    MultiFab old_spectral(boxes, dm, layout.ncomp(), 0);
    MultiFab new_spectral(boxes, dm, layout.ncomp(), 0);
    MultiFab old_compact(boxes, dm, 2, 0);
    MultiFab new_compact(boxes, dm, 2, 0);
    erf_auxiliary::AuxiliaryFaceTransfer accepted_spectral;
    erf_auxiliary::AuxiliaryFaceTransfer accepted_compact;
    accepted_spectral.define(boxes, dm, layout.ncomp(), 0);
    accepted_compact.define(boxes, dm, 2, 0);
    accepted_spectral.setVal(Real(0.0));
    for (MFIter mfi(accepted_spectral.x()); mfi.isValid(); ++mfi) {
        const auto flux = accepted_spectral.x().array(mfi);
        const Real known_flux[] = {Real(0.0), Real(0.25), Real(0.5), Real(0.25), Real(0.0)};
        for (int i = 0; i <= 4; ++i) flux(i, 0, 0, 0) = known_flux[i];
    }

    for (MFIter mfi(old_spectral); mfi.isValid(); ++mfi) {
        const auto old_arr = old_spectral.array(mfi);
        const auto new_arr = new_spectral.array(mfi);
        const auto old_bulk = old_compact.array(mfi);
        const auto new_bulk = new_compact.array(mfi);
        for (int i = 0; i <= 3; ++i) {
            const Real divergence = (i == 0 || i == 1) ? Real(0.25) : Real(-0.25);
            for (int b = 0; b < layout.ncomp(); ++b) {
                old_arr(i, 0, 0, b) = Real(10.0);
                new_arr(i, 0, 0, b) = Real(10.0) - (b == 0 ? divergence : Real(0.0));
            }
            old_bulk(i, 0, 0, 0) = Real(20.0);
            old_bulk(i, 0, 0, 1) = Real(20.0);
            new_bulk(i, 0, 0, 0) = Real(20.0) - divergence;
            new_bulk(i, 0, 0, 1) = Real(20.0);
        }
    }
    const erf_sbm::SBMBulkProjection projection(layout);
    projection.apply_to_face_transfer(accepted_spectral, accepted_compact);
    const auto closure = erf_sbm::evaluate_accepted_transfer_closure(
        layout, old_spectral, new_spectral, accepted_spectral,
        old_compact, 0, 1, new_compact, 0, 1, accepted_compact, geometry);
    EXPECT_LE(closure.spectral_max, closure.spectral_tolerance);
    EXPECT_LE(closure.qc_max, closure.qc_tolerance);
    EXPECT_LE(closure.qr_max, closure.qr_tolerance);
    EXPECT_TRUE(closure.passes());

    // Exact zero states use an exactly zero scale and therefore an exactly
    // zero closure tolerance.  Small physical magnitudes still get a
    // magnitude-scaled tolerance without an order-one floor.
    MultiFab zero_spectral(boxes, dm, layout.ncomp(), 0);
    MultiFab zero_compact(boxes, dm, 2, 0);
    MultiFab small_spectral(boxes, dm, layout.ncomp(), 0);
    MultiFab small_compact(boxes, dm, 2, 0);
    zero_spectral.setVal(Real(0.0));
    zero_compact.setVal(Real(0.0));
    small_spectral.setVal(Real(1.0e-14));
    small_compact.setVal(Real(1.0e-14));
    erf_auxiliary::AuxiliaryFaceTransfer zero_transfer_spectral;
    erf_auxiliary::AuxiliaryFaceTransfer zero_transfer_compact;
    zero_transfer_spectral.define(boxes, dm, layout.ncomp(), 0);
    zero_transfer_compact.define(boxes, dm, 2, 0);
    zero_transfer_spectral.setVal(Real(0.0));
    zero_transfer_compact.setVal(Real(0.0));
    const auto zero_closure = erf_sbm::evaluate_accepted_transfer_closure(
        layout, zero_spectral, zero_spectral, zero_transfer_spectral,
        zero_compact, 0, 1, zero_compact, 0, 1, zero_transfer_compact, geometry);
    EXPECT_DOUBLE_EQ(zero_closure.spectral_max, Real(0.0));
    EXPECT_DOUBLE_EQ(zero_closure.qc_max, Real(0.0));
    EXPECT_DOUBLE_EQ(zero_closure.qr_max, Real(0.0));
    EXPECT_DOUBLE_EQ(zero_closure.spectral_tolerance, Real(0.0));
    EXPECT_DOUBLE_EQ(zero_closure.qc_tolerance, Real(0.0));
    EXPECT_DOUBLE_EQ(zero_closure.qr_tolerance, Real(0.0));
    EXPECT_TRUE(zero_closure.passes());
    const auto small_closure = erf_sbm::evaluate_accepted_transfer_closure(
        layout, small_spectral, small_spectral, zero_transfer_spectral,
        small_compact, 0, 1, small_compact, 0, 1, zero_transfer_compact, geometry);
    EXPECT_DOUBLE_EQ(small_closure.spectral_max, Real(0.0));
    EXPECT_GT(small_closure.spectral_tolerance, Real(0.0));
    EXPECT_TRUE(small_closure.passes());

    // Negative control for the production second-step bug: using the
    // destination compact state as the old baseline leaves the material
    // accepted-transfer divergence in the residual.
    MultiFab destination_compact(boxes, dm, 2, 0);
    amrex::MultiFab::Copy(destination_compact, new_compact, 0, 0, 2, 0);
    const auto stale_baseline = erf_sbm::evaluate_accepted_transfer_closure(
        layout, old_spectral, new_spectral, accepted_spectral,
        destination_compact, 0, 1, new_compact, 0, 1, accepted_compact, geometry);
    EXPECT_GT(stale_baseline.qc_max, stale_baseline.qc_tolerance);

    // Negative control: using final-stage-only anelastic flux when stage 0 is
    // zero and stage 1 is twice the accepted average leaves a material local
    // closure residual.
    erf_auxiliary::AuxiliaryFaceTransfer final_stage_only;
    erf_auxiliary::AuxiliaryFaceTransfer wrong_compact;
    final_stage_only.define(boxes, dm, layout.ncomp(), 0);
    wrong_compact.define(boxes, dm, 2, 0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Copy(final_stage_only.direction(dir), accepted_spectral.direction(dir),
                              0, 0, layout.ncomp(), 0);
        amrex::MultiFab::Saxpy(final_stage_only.direction(dir), Real(1.0),
                               accepted_spectral.direction(dir), 0, 0, layout.ncomp(), 0);
    }
    projection.apply_to_face_transfer(final_stage_only, wrong_compact);
    const auto wrong = erf_sbm::evaluate_accepted_transfer_closure(
        layout, old_spectral, new_spectral, final_stage_only,
        old_compact, 0, 1, new_compact, 0, 1, wrong_compact, geometry);
    EXPECT_GT(wrong.spectral_max, wrong.spectral_tolerance);

    // Negative control: perturbing one accepted spectral face is also visible
    // directly in the local residual and cannot be hidden by a global sum.
    final_stage_only.x().setVal(Real(0.0));
    for (MFIter mfi(final_stage_only.x()); mfi.isValid(); ++mfi) {
        final_stage_only.x().array(mfi)(2, 0, 0, 0) = Real(0.75);
    }
    const auto perturbed = erf_sbm::evaluate_accepted_transfer_closure(
        layout, old_spectral, new_spectral, final_stage_only,
        old_compact, 0, 1, new_compact, 0, 1, wrong_compact, geometry);
    EXPECT_GT(perturbed.spectral_max, perturbed.spectral_tolerance);
}

TEST (SBMP1, FiniteAndMaterialNegativityChecksAreFailClosed)
{
    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    MultiFab state(boxes, dm, 2, 0);
    const auto set_state = [&](const Real background, const Real component_one) {
        state.setVal(background);
        for (MFIter mfi(state); mfi.isValid(); ++mfi) {
            const auto arr = state.array(mfi);
            arr(domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2), 1) = component_one;
        }
    };

    set_state(Real(0.0), Real(0.0));
    EXPECT_NO_THROW(erf_sbm::validate_nonnegative_state(state, 2, "exact zero"));
    set_state(Real(1.0e-14), Real(5.0e-15));
    EXPECT_NO_THROW(erf_sbm::validate_nonnegative_state(state, 2, "small positive state"));
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real small_scale = Real(1.0e-14);
    set_state(small_scale, -Real(64.0) * eps * small_scale);
    EXPECT_NO_THROW(erf_sbm::validate_nonnegative_state(state, 2, "roundoff negative"));
    set_state(small_scale, -Real(256.0) * eps * small_scale);
    EXPECT_THROW(erf_sbm::validate_nonnegative_state(state, 2, "scaled negative"), std::runtime_error);
    set_state(small_scale, -Real(1.0e-2) * small_scale);
    EXPECT_THROW(erf_sbm::validate_nonnegative_state(state, 2, "material negative"), std::runtime_error);
    set_state(Real(1.0), std::numeric_limits<Real>::quiet_NaN());
    EXPECT_THROW(erf_sbm::validate_nonnegative_state(state, 2, "non-min NaN"), std::runtime_error);
    set_state(Real(1.0), std::numeric_limits<Real>::infinity());
    EXPECT_THROW(erf_sbm::validate_nonnegative_state(state, 2, "+Inf"), std::runtime_error);
    set_state(Real(1.0), -std::numeric_limits<Real>::infinity());
    EXPECT_THROW(erf_sbm::validate_nonnegative_state(state, 2, "-Inf"), std::runtime_error);
}

void run_manufactured_transport (const int nbins, const bool anelastic,
                                 const erf_sbm::TransportMethod method = erf_sbm::TransportMethod::DonorCell,
                                 const erf_sbm::MomentMode mode = erf_sbm::MomentMode::OneMoment)
{
    const auto layout = make_layout(nbins, mode);
    const Box domain(IntVect(0, 0, 0), IntVect(3, 1, 1));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    const Geometry geometry = make_geometry(domain);
    MultiFab rho(boxes, dm, 1, 1);
    MultiFab core(boxes, dm, RhoQ3_comp + 1, 1);
    MultiFab carrier_x(amrex::convert(boxes, IntVect(1, 0, 0)), dm, 1, 1);
    MultiFab carrier_y(amrex::convert(boxes, IntVect(0, 1, 0)), dm, 1, 1);
    MultiFab carrier_z(amrex::convert(boxes, IntVect(0, 0, 1)), dm, 1, 1);
    rho.setVal(Real(0.0));
    core.setVal(Real(0.0));
    carrier_x.setVal(Real(0.125));
    carrier_y.setVal(Real(0.0));
    carrier_z.setVal(Real(0.0));
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
        const auto rho_arr = rho.array(mfi);
        const auto core_arr = core.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real density = Real(1.0) + Real(0.05) * Real(i + 1);
            rho_arr(i,j,k,0) = density;
            core_arr(i,j,k,Rho_comp) = density;
        });
    }
    rho.FillBoundary(geometry.periodicity());
    carrier_x.FillBoundary(geometry.periodicity());
    carrier_y.FillBoundary(geometry.periodicity());
    carrier_z.FillBoundary(geometry.periodicity());

    erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
    manager.define_level(0, boxes, dm, 2);
    auto& initial = manager.output(0);
    const auto& population = layout.populations().front();
    const int mass_offset = population.mass_offset;
    const int number_offset = population.number_offset;
    for (MFIter mfi(initial); mfi.isValid(); ++mfi) {
        const auto aux = initial.array(mfi);
        const auto density = rho.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real variation = Real(1.0) + Real(0.25) *
                std::sin(Real(6.2831853071795864769) * (Real(i) + Real(0.5)) / Real(4.0));
            for (int b = 0; b < nbins; ++b) {
                const Real mass = density(i,j,k,0) * Real(1.e-3) * Real(b + 1) * variation;
                aux(i,j,k,mass_offset+b) = mass;
                if (mode == erf_sbm::MomentMode::TwoMoment) {
                    const Real pivot = Real(b) + Real(0.5);
                    aux(i,j,k,number_offset+b) = mass / pivot;
                }
            }
        });
    }
    initial.FillBoundary(geometry.periodicity());
    manager.begin_step(0);

    const auto check_stage_invariants = [&]() {
        for (MFIter mfi(manager.output(0)); mfi.isValid(); ++mfi) {
            const auto aux = manager.output(0).const_array(mfi);
            const auto density = rho.const_array(mfi);
            const auto compact = core.const_array(mfi);
            for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
                for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
                    for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
                        const Real ratio = aux(i,j,k,0) / density(i,j,k,0);
                        if (mode == erf_sbm::MomentMode::OneMoment) {
                            for (int b = 1; b < nbins; ++b) {
                                EXPECT_NEAR(aux(i,j,k,b) / density(i,j,k,0),
                                            ratio * Real(b + 1), 2.e-14);
                            }
                        } else {
                            for (int b = 0; b < nbins; ++b) {
                                const Real mass = aux(i,j,k,population.mass_offset+b);
                                const Real number = aux(i,j,k,population.number_offset+b);
                                const Real lower = population.grid.edges()[static_cast<std::size_t>(b)];
                                const Real upper = population.grid.edges()[static_cast<std::size_t>(b+1)];
                                EXPECT_GE(number, -2.e-14);
                                EXPECT_GE(mass - lower*number, -2.e-14);
                                EXPECT_GE(upper*number - mass, -2.e-14);
                            }
                        }
                        Real qc = 0.0;
                        Real qr = 0.0;
                        for (int b = 0; b < nbins / 2; ++b) {
                            qc += aux(i,j,k,population.mass_offset+b);
                        }
                        for (int b = nbins / 2; b < nbins; ++b) {
                            qr += aux(i,j,k,population.mass_offset+b);
                        }
                        EXPECT_NEAR(compact(i,j,k,RhoQ2_comp), qc, 2.e-14);
                        EXPECT_NEAR(compact(i,j,k,RhoQ3_comp), qr, 2.e-14);
                    }
                }
            }
        }
    };

    const auto c0 = anelastic ?
        erf_auxiliary::make_anelastic_stage(0, 0.0, 0.0, 1.0, 1.0, nullptr, nullptr) :
        erf_auxiliary::make_compressible_stage(0, 0.0, 0.0, 1.0/3.0, 1.0, nullptr, nullptr);
    const Real initial_mass = [&]() {
        Real total = 0.0;
        for (int b = 0; b < nbins; ++b) total += manager.old(0).sum(b);
        return total;
    }();
    erf_sbm::advance_stage(manager, layout, c0, rho, core,
                           carrier_x, carrier_y, carrier_z, geometry,
                           manager.face_transfer_ledger(0).stage(), method);
    check_stage_invariants();
    if (anelastic) {
        const auto c1 = erf_auxiliary::make_anelastic_stage(1, 0.0, 1.0, 1.0, 1.0, nullptr, nullptr);
        erf_sbm::advance_stage(manager, layout, c1, rho, core,
                               carrier_x, carrier_y, carrier_z, geometry,
                               manager.face_transfer_ledger(0).stage(), method);
        check_stage_invariants();
    } else {
        const auto c1 = erf_auxiliary::make_compressible_stage(1, 0.0, 1.0/3.0, 0.5, 1.0, nullptr, nullptr);
        erf_sbm::advance_stage(manager, layout, c1, rho, core,
                               carrier_x, carrier_y, carrier_z, geometry,
                               manager.face_transfer_ledger(0).stage(), method);
        check_stage_invariants();
        const auto c2 = erf_auxiliary::make_compressible_stage(2, 0.0, 0.5, 1.0, 1.0, nullptr, nullptr);
        erf_sbm::advance_stage(manager, layout, c2, rho, core,
                               carrier_x, carrier_y, carrier_z, geometry,
                               manager.face_transfer_ledger(0).stage(), method);
        check_stage_invariants();
    }

    const Real final_mass = [&]() {
        Real total = 0.0;
        for (int b = 0; b < nbins; ++b) total += manager.output(0).sum(b);
        return total;
    }();
    EXPECT_NEAR(final_mass, initial_mass, 2.e-12 * std::max(Real(1.0), std::abs(initial_mass)));
    if (mode == erf_sbm::MomentMode::TwoMoment) {
        const Real initial_number = [&]() {
            Real total = 0.0;
            for (int b = 0; b < nbins; ++b) total += manager.old(0).sum(population.number_offset+b);
            return total;
        }();
        const Real final_number = [&]() {
            Real total = 0.0;
            for (int b = 0; b < nbins; ++b) total += manager.output(0).sum(population.number_offset+b);
            return total;
        }();
        EXPECT_NEAR(final_number, initial_number,
                    2.e-12 * std::max(Real(1.0), std::abs(initial_number)));
    }
    EXPECT_GT(manager.output(0).max(0) - manager.output(0).min(0), Real(0.0));

    erf_auxiliary::AuxiliaryFaceTransfer bulk_transfer;
    // Use the generic container with a temporary one-box decomposition so the
    // projection test exercises all three face directions.
    bulk_transfer.define(boxes, dm, 2, 0);
    const erf_sbm::SBMBulkProjection projection(layout);
    projection.apply_to_face_transfer(manager.face_transfer_ledger(0).accepted(), bulk_transfer);
    EXPECT_NEAR(projection.max_face_projection_error(
                    manager.face_transfer_ledger(0).accepted(), bulk_transfer),
                Real(0.0), 2.e-14);
}

TEST (SBMP1, ManufacturedVariableDensityFreeStreamCompressibleRuntimeBins)
{
    for (const int nbins : {4, 16, 64}) run_manufactured_transport(nbins, false);
}

TEST (SBMP1, ManufacturedVariableDensityFreeStreamAnelasticRuntimeBins)
{
    for (const int nbins : {4, 16, 64}) run_manufactured_transport(nbins, true);
}

TEST (SBMP2, ProductionGroupedFCTWENOZ3RuntimeBins)
{
    for (const int nbins : {4, 16, 64}) {
        run_manufactured_transport(nbins, false, erf_sbm::TransportMethod::GroupedFCT_WENOZ3);
    }
}

TEST (SBMP2, ProductionGroupedFCTWENOZ3TwoMomentRuntimeBins)
{
    for (const int nbins : {4, 16, 64}) {
        run_manufactured_transport(nbins, false, erf_sbm::TransportMethod::GroupedFCT_WENOZ3,
                                   erf_sbm::MomentMode::TwoMoment);
    }
}

} // namespace
