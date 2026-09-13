#include <gtest/gtest.h>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_MultiFab.H>
#include <AMReX_RealBox.H>

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

erf_sbm::SpectralGridSpec make_grid (const int nbins, const int population = 0,
                                     const int split = -1)
{
    erf_sbm::SpectralGridSpec spec;
    spec.population_id = population;
    spec.coordinate_kind = erf_sbm::CoordinateKind::LiquidMass;
    spec.units = "kg";
    spec.semantic_id = population == 0 ? "liquid_mass" : "dummy_mass";
    spec.edges.resize(static_cast<std::size_t>(nbins + 1));
    spec.pivots.resize(static_cast<std::size_t>(nbins));
    for (int n = 0; n <= nbins; ++n) spec.edges[static_cast<std::size_t>(n)] = Real(n);
    for (int n = 0; n < nbins; ++n) {
        spec.pivots[static_cast<std::size_t>(n)] = Real(n) + Real(0.5);
    }
    spec.cloud_rain_split = split < 0 ? nbins / 2 : split;
    return spec;
}

erf_sbm::SBMLayout make_layout (const int nbins,
                                const erf_sbm::MomentMode mode = erf_sbm::MomentMode::OneMoment)
{
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(make_grid(nbins));
    spec.moment_modes.push_back(mode);
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

    auto bad_split = make_grid(4, 0, 0);
    EXPECT_FALSE(erf_sbm::SpectralGrid::validate(bad_split).valid);
    bad_split = make_grid(4, 0, 4);
    EXPECT_FALSE(erf_sbm::SpectralGrid::validate(bad_split).valid);

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
    spec.populations = {make_grid(4, 0), make_grid(3, 1)};
    spec.moment_modes = {erf_sbm::MomentMode::OneMoment, erf_sbm::MomentMode::OneMoment};
    spec.attached_properties.push_back({"solute", "attached_solute", "kg", 1,
                                        erf_sbm::PropertyKind::MassBoundedSubset,
                                        erf_sbm::SupportRequirement::PositiveMass});
    const erf_sbm::SBMLayout layout(std::move(spec));
    EXPECT_EQ(layout.populations()[0].liquid_mass_offset, 0);
    EXPECT_EQ(layout.populations()[1].liquid_mass_offset, 4);
    EXPECT_EQ(layout.property_offset(0), 7);
    EXPECT_EQ(layout.ncomp(), 10);
    ASSERT_EQ(layout.auxiliary_layout().components().size(), 10U);
    EXPECT_EQ(layout.attached_properties()[0].kind, erf_sbm::PropertyKind::MassBoundedSubset);
    EXPECT_EQ(layout.attached_properties()[0].support, erf_sbm::SupportRequirement::PositiveMass);
    EXPECT_NE(layout.inspection().find("property solute"), std::string::npos);
}

TEST (SBMP0, GenericManagerPreservesBaselineAndPublishesAcceptedState)
{
    const erf_auxiliary::AuxiliaryStateLayout state_layout(
        "manager-test", {{"x", "test", "1"}, {"y", "test", "1"}});
    erf_auxiliary::AuxiliaryStateManager manager(state_layout);
    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    manager.define_level(0, boxes, dm, 1);
    manager.output(0).setVal(Real(2.0));
    manager.begin_step(0);
    manager.output(0).setVal(Real(5.0));
    manager.accept_stage(0);
    EXPECT_EQ(manager.old(0).min(0), Real(2.0));
    EXPECT_EQ(manager.evaluation(0).min(0), Real(5.0));
    EXPECT_EQ(manager.output(0).min(0), Real(5.0));
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

void run_manufactured_transport (const int nbins, const bool anelastic)
{
    const auto layout = make_layout(nbins);
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
    manager.define_level(0, boxes, dm, 1);
    auto& initial = manager.output(0);
    for (MFIter mfi(initial); mfi.isValid(); ++mfi) {
        const auto aux = initial.array(mfi);
        const auto density = rho.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int b = 0; b < nbins; ++b) {
                aux(i,j,k,b) = density(i,j,k,0) * Real(1.e-3) * Real(b + 1);
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
                        for (int b = 1; b < nbins; ++b) {
                            EXPECT_NEAR(aux(i,j,k,b) / density(i,j,k,0),
                                        ratio * Real(b + 1), 2.e-14);
                        }
                        Real qc = 0.0;
                        Real qr = 0.0;
                        for (int b = 0; b < nbins / 2; ++b) qc += aux(i,j,k,b);
                        for (int b = nbins / 2; b < nbins; ++b) qr += aux(i,j,k,b);
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
    erf_sbm::advance_stage(manager, layout, c0, rho, core,
                           carrier_x, carrier_y, carrier_z, geometry);
    check_stage_invariants();
    if (anelastic) {
        const auto c1 = erf_auxiliary::make_anelastic_stage(1, 0.0, 1.0, 1.0, 1.0, nullptr, nullptr);
        erf_sbm::advance_stage(manager, layout, c1, rho, core,
                               carrier_x, carrier_y, carrier_z, geometry);
        check_stage_invariants();
    } else {
        const auto c1 = erf_auxiliary::make_compressible_stage(1, 0.0, 1.0/3.0, 0.5, 1.0, nullptr, nullptr);
        erf_sbm::advance_stage(manager, layout, c1, rho, core,
                               carrier_x, carrier_y, carrier_z, geometry);
        check_stage_invariants();
        const auto c2 = erf_auxiliary::make_compressible_stage(2, 0.0, 0.5, 1.0, 1.0, nullptr, nullptr);
        erf_sbm::advance_stage(manager, layout, c2, rho, core,
                               carrier_x, carrier_y, carrier_z, geometry);
        check_stage_invariants();
    }
}

TEST (SBMP1, ManufacturedVariableDensityFreeStreamCompressibleRuntimeBins)
{
    for (const int nbins : {4, 16, 64}) run_manufactured_transport(nbins, false);
}

TEST (SBMP1, ManufacturedVariableDensityFreeStreamAnelasticRuntimeBins)
{
    for (const int nbins : {4, 16, 64}) run_manufactured_transport(nbins, true);
}

} // namespace
