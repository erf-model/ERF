#include <gtest/gtest.h>

#include "ERF_SBMAMR.H"
#include "ERF_AuxiliaryStateManager.H"
#include "ERF_SBMBoundary.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMContracts.H"
#include "ERF_SBMDiffusion.H"
#include "ERF_SBMFCT.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMLayout.H"
#include "ERF_SBMTransportPrototype.H"
#include "ERF_Interpolation_WENO_Z.H"
#include "ERF_IndexDefines.H"

#include <AMReX_MFParallelFor.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>

#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

namespace {

using erf_sbm::MomentMode;
using Real = amrex::Real;
using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;

erf_sbm::SBMLayout make_layout(const int nbins, const MomentMode mode,
                               const bool with_property = false)
{
    erf_sbm::SpectralPopulationSpec population;
    population.population_id = 0;
    population.semantic_id = "liquid";
    population.phase = erf_sbm::PopulationPhase::Liquid;
    population.moment_mode = mode;
    population.grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    population.grid.coordinate_units = "kg";
    for (int b = 0; b <= nbins; ++b) population.grid.edges.push_back(static_cast<Real>(b));
    for (int b = 0; b < nbins; ++b) population.grid.pivots.push_back(static_cast<Real>(b) + 0.5);
    population.mass_state_units = "kg m^-3";
    population.number_state_units = "m^-3";
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(population);
    spec.liquid_projection = {0, nbins/2};
    if (with_property) {
        spec.attached_properties.push_back({"solute", "solute", "kg m^-3", 0,
            erf_sbm::PropertyKind::MassBoundedSubset,
            erf_sbm::SupportRequirement::PositiveMass,
            erf_sbm::PropertyRemapPolicy::CarrierBinConservative,
            true, false, 0.0, 1.0});
    }
    return erf_sbm::SBMLayout(std::move(spec));
}

struct ProductionChunkRun {
    std::unique_ptr<amrex::MultiFab> spectral;
    std::unique_ptr<amrex::MultiFab> compact;
    std::unique_ptr<erf_auxiliary::AuxiliaryFaceTransfer> accepted;
    std::size_t temporary_bytes{0};
};

ProductionChunkRun run_production_chunk_case(const int nbins, const MomentMode mode,
                                             const int chunk_size)
{
    const auto layout = make_layout(nbins, mode);
    const Box domain(IntVect(0, 0, 0), IntVect(15, 1, 1));
    BoxArray boxes(domain);
    boxes.maxSize(IntVect(8, 2, 2));
    const DistributionMapping dm(boxes);
    const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                  {AMREX_D_DECL(1.0, 1.0, 1.0)});
    const std::array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    const Geometry geometry(domain, &real_box, amrex::CoordSys::cartesian,
                            periodicity.data());

    amrex::MultiFab rho(boxes, dm, 1, 2);
    amrex::MultiFab core(boxes, dm, RhoQ3_comp + 1, 0);
    amrex::MultiFab xflux(amrex::convert(boxes, IntVect(1, 0, 0)), dm, 1, 0);
    amrex::MultiFab yflux(amrex::convert(boxes, IntVect(0, 1, 0)), dm, 1, 0);
    amrex::MultiFab zflux(amrex::convert(boxes, IntVect(0, 0, 1)), dm, 1, 0);
    rho.setVal(Real(0.0));
    core.setVal(Real(0.0));
    xflux.setVal(Real(0.125));
    yflux.setVal(Real(0.0));
    zflux.setVal(Real(0.0));
    for (amrex::MFIter mfi(rho); mfi.isValid(); ++mfi) {
        const auto density = rho.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            density(i,j,k,0) = Real(1.0) + Real(0.05) * Real(i + 1);
        });
    }
    rho.FillBoundary(geometry.periodicity());
    xflux.FillBoundary(geometry.periodicity());
    yflux.FillBoundary(geometry.periodicity());
    zflux.FillBoundary(geometry.periodicity());

    erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
    manager.define_level(0, boxes, dm, 2);
    const auto& population = layout.populations().front();
    for (amrex::MFIter mfi(manager.output(0)); mfi.isValid(); ++mfi) {
        const auto state = manager.output(0).array(mfi);
        const auto density = rho.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real variation = Real(1.0) + Real(0.20) *
                std::sin(Real(6.2831853071795864769) *
                         (Real(i) + Real(0.5)) / Real(16.0));
            for (int b = 0; b < nbins; ++b) {
                const Real mass = density(i,j,k,0) * Real(1.e-3) *
                    Real(b + 1) * variation;
                state(i,j,k,population.mass_offset+b) = mass;
                if (mode == MomentMode::TwoMoment) {
                    state(i,j,k,population.number_offset+b) = mass / (Real(b) + Real(0.5));
                }
            }
        });
    }
    manager.output(0).FillBoundary(geometry.periodicity());
    manager.begin_step(0, 0.0);

    erf_auxiliary::AuxiliaryFaceTransfer stage_flux;
    stage_flux.define(boxes, dm, layout.ncomp(), 0);
    const Real dt = Real(0.5) / Real(16.0);
    const auto context = erf_auxiliary::make_compressible_stage(
        2, 0.0, 0.0, dt, dt, nullptr, nullptr);
    erf_sbm::advance_stage(manager, layout, context, rho, core,
                           xflux, yflux, zflux, geometry, stage_flux,
                           erf_sbm::TransportMethod::GroupedFCT_WENOZ3,
                           0, Real(0.0), chunk_size);

    ProductionChunkRun result;
    result.spectral = std::make_unique<amrex::MultiFab>(
        boxes, dm, layout.ncomp(), manager.output(0).nGrowVect());
    amrex::MultiFab::Copy(*result.spectral, manager.output(0), 0, 0,
                          layout.ncomp(), manager.output(0).nGrowVect());
    result.compact = std::make_unique<amrex::MultiFab>(boxes, dm, RhoQ3_comp + 1, 0);
    amrex::MultiFab::Copy(*result.compact, core, 0, 0, RhoQ3_comp + 1, 0);
    result.accepted = std::make_unique<erf_auxiliary::AuxiliaryFaceTransfer>();
    result.accepted->define(boxes, dm, layout.ncomp(), 0);
    const auto& accepted = manager.face_transfer_ledger(0).accepted();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Copy(result.accepted->direction(dir), accepted.direction(dir),
                              0, 0, layout.ncomp(), 0);
    }
    result.temporary_bytes = erf_sbm::grouped_fct_peak_working_bytes(
        layout, boxes, dm, chunk_size);
    return result;
}

Real max_multifab_difference(const amrex::MultiFab& left,
                             const amrex::MultiFab& right, const int ncomp)
{
    amrex::MultiFab difference(left.boxArray(), left.DistributionMap(), ncomp, 0);
    amrex::MultiFab::Copy(difference, left, 0, 0, ncomp, 0);
    amrex::MultiFab::Subtract(difference, right, 0, 0, ncomp, 0);
    // Keep this comparison local.  The helper is called once per rank and a
    // valid one-box decomposition leaves non-owning MPI ranks with no FABs;
    // an empty-rank collective norm is not a useful comparison and is not
    // required because the same DistributionMapping is used for both runs.
    if (difference.local_size() == 0) return Real(0.0);
    return difference.norm0(0, ncomp, amrex::IntVect(0), true);
}

Real max_face_transfer_difference(const erf_auxiliary::AuxiliaryFaceTransfer& left,
                                  const erf_auxiliary::AuxiliaryFaceTransfer& right,
                                  const int ncomp)
{
    Real result = Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        result = amrex::max(result, max_multifab_difference(
            left.direction(dir), right.direction(dir), ncomp));
    }
    return result;
}

TEST(SBMP2, ConstraintGroupsCoverOneAndTwoMomentLayouts)
{
    for (const int nbins : {4, 16, 64}) {
        for (const auto mode : {MomentMode::OneMoment, MomentMode::TwoMoment}) {
            const auto layout = make_layout(nbins, mode, true);
            const auto groups = erf_sbm::make_constraint_groups(layout);
            ASSERT_EQ(groups.size(), static_cast<std::size_t>(nbins));
            for (const auto& group : groups) {
                EXPECT_FALSE(group.constraints.empty());
                EXPECT_TRUE(group.contains(group.members.front()));
                std::vector<Real> state(static_cast<std::size_t>(layout.ncomp()), 0.0);
                for (const int component : group.members) state[static_cast<std::size_t>(component)] =
                    (component >= layout.property_offset(0)) ? 0.001 : 1.0;
                if (mode == MomentMode::TwoMoment) {
                    const auto& population = layout.populations().front();
                    state[static_cast<std::size_t>(population.mass_offset + group.bin)] =
                        0.5 * (population.grid.edges()[static_cast<std::size_t>(group.bin)] +
                               population.grid.edges()[static_cast<std::size_t>(group.bin + 1)]);
                    state[static_cast<std::size_t>(population.number_offset + group.bin)] = 1.0;
                }
                Real margin = 0.0;
                std::string failed;
                EXPECT_TRUE(group.admissible(state, &margin, &failed))
                    << group.semantic_id << " failed " << failed << " margin=" << margin;
            }
        }
    }
}

TEST(SBMP2, ConstraintDescriptorsFlattenCompleteGroups)
{
    const auto layout = make_layout(4, MomentMode::TwoMoment, true);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    const auto descriptors = erf_sbm::make_constraint_descriptors(layout);
    std::size_t expected = 0;
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
        expected += groups[gi].constraints.size();
        for (std::size_t ci = 0; ci < groups[gi].constraints.size(); ++ci) {
            const auto& d = descriptors[expected - groups[gi].constraints.size() + ci];
            EXPECT_EQ(d.group_index, static_cast<int>(gi));
            EXPECT_EQ(d.bin, groups[gi].bin);
            EXPECT_EQ(d.constraint_index, static_cast<int>(ci));
            EXPECT_GE(d.term_count, 1);
            EXPECT_LE(d.term_count, 2);
        }
    }
    EXPECT_EQ(descriptors.size(), expected);
}

TEST(SBMP2, ProductionChunkPlannerKeepsCompleteConstraintGroupsAtomic)
{
    const auto layout = make_layout(4, MomentMode::TwoMoment, true);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    const auto one_group = erf_sbm::make_constraint_closure_chunks(layout, 1);
    ASSERT_EQ(one_group.size(), groups.size());
    for (std::size_t i = 0; i < one_group.size(); ++i) {
        ASSERT_EQ(one_group[i].group_indices.size(), 1U);
        EXPECT_EQ(one_group[i].group_indices.front(), static_cast<int>(i));
        EXPECT_EQ(one_group[i].components.size(), groups[i].members.size());
    }
    const auto two_groups = erf_sbm::make_constraint_closure_chunks(layout, 2);
    ASSERT_EQ(two_groups.size(), 2U);
    EXPECT_EQ(two_groups[0].group_indices.size(), 2U);
    EXPECT_EQ(two_groups[0].components.size(), groups[0].members.size() + groups[1].members.size());
    EXPECT_THROW((void)erf_sbm::make_constraint_closure_chunks(layout, 0), std::invalid_argument);
}

TEST(SBMP2, GroupedFCTUsesCellWideConstraintBudgets)
{
    const auto layout = make_layout(2, MomentMode::OneMoment);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    // Two distinct outgoing faces share cell 0.  Each high-order correction
    // consumes 0.75 of its mass margin; a face-local limiter would accept both
    // at lambda=1 and produce -0.5 in cell 0.  The cell-wide budget requires
    // lambda=2/3 for both and lands exactly on the admissible boundary.
    const std::vector<Real> low_state{1.0, 0.0, 10.0, 0.0, 10.0, 0.0};
    erf_sbm::FCTFaceTransfer face0;
    face0.left_cell = 0; face0.right_cell = 1;
    face0.low = {0.0, 0.0}; face0.high = {0.75, 0.0};
    erf_sbm::FCTFaceTransfer face1;
    face1.left_cell = 0; face1.right_cell = 2;
    face1.low = {0.0, 0.0}; face1.high = {0.75, 0.0};
    const auto result = erf_sbm::limit_grouped(low_state, 3, 2,
                                               {face0, face1}, groups, 1);
    EXPECT_NEAR(result.limiter[0], 2.0/3.0, 1.e-14);
    EXPECT_NEAR(result.limiter[1], 2.0/3.0, 1.e-14);
    EXPECT_NEAR(result.updated_state[0], 0.0, 1.e-14);
    for (int cell = 0; cell < 3; ++cell) {
        const std::vector<Real> state(result.updated_state.begin() + cell*2,
                                      result.updated_state.begin() + (cell+1)*2);
        EXPECT_TRUE(groups.front().admissible(state));
    }
}

TEST(SBMP2, TwoMomentEndpointTransformIsStableAtBounds)
{
    const auto lower = erf_sbm::transform_two_moment(2.0, 2.0, 1.0, 3.0);
    EXPECT_EQ(lower.L, 2.0); EXPECT_EQ(lower.H, 0.0);
    const auto upper = erf_sbm::transform_two_moment(2.0, 6.0, 1.0, 3.0);
    EXPECT_EQ(upper.L, 0.0); EXPECT_EQ(upper.H, 2.0);
    const auto center = erf_sbm::transform_two_moment(2.0, 4.0, 1.0, 3.0);
    EXPECT_NEAR(center.L, 1.0, 32*std::numeric_limits<Real>::epsilon());
    EXPECT_NEAR(center.H, 1.0, 32*std::numeric_limits<Real>::epsilon());
    EXPECT_THROW((void)erf_sbm::transform_two_moment(1.0, 4.1, 1.0, 3.0), std::invalid_argument);
    EXPECT_THROW((void)erf_sbm::transform_two_moment(-1.0, 0.0, 1.0, 3.0), std::invalid_argument);
    const auto recovered = erf_sbm::inverse_two_moment(center.L, center.H, 1.0, 3.0);
    EXPECT_NEAR(recovered.first, 2.0, 1.e-14);
    EXPECT_NEAR(recovered.second, 4.0, 1.e-14);

    // Cancellation-sized endpoint roundoff is normalized only in the
    // derived (L,H) representation; the authoritative (C,M) inputs are not
    // modified.  Violations outside the scale-aware tolerance fail closed.
    const Real eps = std::numeric_limits<Real>::epsilon();
    const auto near_lower = erf_sbm::transform_two_moment(2.0, 2.0 - 32*eps, 1.0, 3.0);
    const auto near_upper = erf_sbm::transform_two_moment(2.0, 6.0 + 32*eps, 1.0, 3.0);
    EXPECT_EQ(near_lower.H, 0.0);
    EXPECT_EQ(near_upper.L, 0.0);
    EXPECT_THROW((void)erf_sbm::transform_two_moment(2.0, 2.0 - 1.e-8, 1.0, 3.0), std::invalid_argument);
    EXPECT_THROW((void)erf_sbm::transform_two_moment(2.0, 6.0 + 1.e-8, 1.0, 3.0), std::invalid_argument);

    const Real tiny_c = Real(1.e-100);
    const Real tiny_m = Real(2.0) * tiny_c;
    const auto tiny = erf_sbm::transform_two_moment(tiny_c, tiny_m, 1.0, 3.0);
    const auto tiny_back = erf_sbm::inverse_two_moment(tiny.L, tiny.H, 1.0, 3.0);
    EXPECT_NEAR(tiny_back.first, tiny_c, 64*eps*tiny_c);
    EXPECT_NEAR(tiny_back.second, tiny_m, 64*eps*tiny_m);

    const Real huge_c = Real(1.e100);
    const Real huge_m = Real(2.0) * huge_c;
    const auto huge = erf_sbm::transform_two_moment(huge_c, huge_m, 1.0, 3.0);
    const auto huge_back = erf_sbm::inverse_two_moment(huge.L, huge.H, 1.0, 3.0);
    EXPECT_NEAR(huge_back.first, huge_c, 64*eps*huge_c);
    EXPECT_NEAR(huge_back.second, huge_m, 64*eps*huge_m);

    const Real lower_edge = Real(2.5), upper_edge = Real(9.25), count = Real(3.75);
    const Real mass = Real(6.125) * count;
    const auto nontrivial = erf_sbm::transform_two_moment(count, mass, lower_edge, upper_edge);
    const auto nontrivial_back = erf_sbm::inverse_two_moment(
        nontrivial.L, nontrivial.H, lower_edge, upper_edge);
    EXPECT_NEAR(nontrivial_back.first, count, 64*eps*count);
    EXPECT_NEAR(nontrivial_back.second, mass, 64*eps*std::abs(mass));
}

TEST(SBMP2, ProductionChunkWorkingMemoryScalesWithAtomicChunkPolicy)
{
    const Box domain(IntVect(0, 0, 0), IntVect(7, 7, 7));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    for (const int nbins : {4, 16, 64}) {
        for (const auto mode : {MomentMode::OneMoment, MomentMode::TwoMoment}) {
            const auto layout = make_layout(nbins, mode);
            const auto one_group = erf_sbm::grouped_fct_peak_working_bytes(layout, boxes, dm, 1);
            const auto all_groups = erf_sbm::grouped_fct_peak_working_bytes(layout, boxes, dm, nbins);
            EXPECT_GT(one_group, std::size_t(0));
            EXPECT_GE(all_groups, one_group);
            EXPECT_EQ(one_group, erf_sbm::grouped_fct_peak_working_bytes(layout, boxes, dm, 1));
            if (nbins > 1) EXPECT_GT(all_groups, one_group);
        }
    }
    EXPECT_THROW((void)erf_sbm::grouped_fct_peak_working_bytes(
        make_layout(4, MomentMode::OneMoment), boxes, dm, 0), std::invalid_argument);
}

TEST(SBMP2, ProductionGroupedTransportIsChunkEquivalentAcrossRuntimeLayouts)
{
    // The all-group run is the production reference.  Every other run still
    // executes the grouped WENO/FCT FAB path; only the complete-group chunk
    // policy changes.  This exercises both physical moments, accepted face
    // transfers, compact projections, and the runtime 4/16/64-bin layouts.
    std::size_t one_group_bytes[2][3]{{0, 0, 0}, {0, 0, 0}};
    std::size_t all_group_bytes[2][3]{{0, 0, 0}, {0, 0, 0}};
    std::ofstream evidence;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        evidence.open("/private/tmp/erf_sbm_p2_chunk_equivalence_memory.csv");
        ASSERT_TRUE(evidence.good());
        evidence << "nbins,moment_mode,chunk_size,temporary_p2_working_bytes\n";
    }
    for (const int nbins : {4, 16, 64}) {
        for (const auto mode : {MomentMode::OneMoment, MomentMode::TwoMoment}) {
            SCOPED_TRACE(std::string("nbins=") + std::to_string(nbins) +
                         " mode=" + std::to_string(static_cast<int>(mode)));
            const auto reference = run_production_chunk_case(nbins, mode, nbins);
            const auto one_group = run_production_chunk_case(nbins, mode, 1);
            const auto two_groups = run_production_chunk_case(nbins, mode, 2);
            const auto four_groups = run_production_chunk_case(nbins, mode, 4);
            const Real tolerance = Real(5.e-13);
            EXPECT_LE(max_multifab_difference(*one_group.spectral,
                                               *reference.spectral, nbins *
                                               (mode == MomentMode::TwoMoment ? 2 : 1)),
                      tolerance);
            EXPECT_LE(max_multifab_difference(*two_groups.spectral,
                                               *reference.spectral, nbins *
                                               (mode == MomentMode::TwoMoment ? 2 : 1)),
                      tolerance);
            EXPECT_LE(max_multifab_difference(*four_groups.spectral,
                                               *reference.spectral, nbins *
                                               (mode == MomentMode::TwoMoment ? 2 : 1)),
                      tolerance);
            EXPECT_LE(max_multifab_difference(*one_group.compact,
                                               *reference.compact, RhoQ3_comp + 1),
                      tolerance);
            EXPECT_LE(max_multifab_difference(*two_groups.compact,
                                               *reference.compact, RhoQ3_comp + 1),
                      tolerance);
            EXPECT_LE(max_multifab_difference(*four_groups.compact,
                                               *reference.compact, RhoQ3_comp + 1),
                      tolerance);
            EXPECT_LE(max_face_transfer_difference(*one_group.accepted,
                                                    *reference.accepted,
                                                    nbins *
                                                    (mode == MomentMode::TwoMoment ? 2 : 1)),
                      tolerance);
            EXPECT_LE(max_face_transfer_difference(*two_groups.accepted,
                                                    *reference.accepted,
                                                    nbins *
                                                    (mode == MomentMode::TwoMoment ? 2 : 1)),
                      tolerance);
            EXPECT_LE(max_face_transfer_difference(*four_groups.accepted,
                                                    *reference.accepted,
                                                    nbins *
                                                    (mode == MomentMode::TwoMoment ? 2 : 1)),
                      tolerance);

            EXPECT_GE(one_group.temporary_bytes, std::size_t(1));
            EXPECT_LE(one_group.temporary_bytes, two_groups.temporary_bytes);
            EXPECT_LE(two_groups.temporary_bytes, four_groups.temporary_bytes);
            EXPECT_LE(four_groups.temporary_bytes, reference.temporary_bytes);

            const int mode_index = mode == MomentMode::TwoMoment ? 1 : 0;
            const int bin_index = nbins == 4 ? 0 : (nbins == 16 ? 1 : 2);
            one_group_bytes[mode_index][bin_index] = one_group.temporary_bytes;
            all_group_bytes[mode_index][bin_index] = reference.temporary_bytes;
            if (evidence) {
                evidence << nbins << ',' << mode_index << ",1," << one_group.temporary_bytes << '\n'
                         << nbins << ',' << mode_index << ",2," << two_groups.temporary_bytes << '\n'
                         << nbins << ',' << mode_index << ",4," << four_groups.temporary_bytes << '\n'
                         << nbins << ',' << mode_index << ',' << nbins << ','
                         << reference.temporary_bytes << '\n';
            }
        }
    }
    for (int mode_index = 0; mode_index < 2; ++mode_index) {
        EXPECT_EQ(one_group_bytes[mode_index][0], one_group_bytes[mode_index][1]);
        EXPECT_EQ(one_group_bytes[mode_index][1], one_group_bytes[mode_index][2]);
        EXPECT_LT(all_group_bytes[mode_index][0], all_group_bytes[mode_index][1]);
        EXPECT_LT(all_group_bytes[mode_index][1], all_group_bytes[mode_index][2]);
    }
}

TEST(SBMP2, GroupedFCTUsesOneFaceLimiterAndConservesEveryComponent)
{
    const auto layout = make_layout(2, MomentMode::OneMoment);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    std::vector<Real> low_state{1.0, 0.0, 2.0, 0.0};
    erf_sbm::FCTFaceTransfer face;
    face.left_cell = 0; face.right_cell = 1;
    face.left_volume = 2.0; face.right_volume = 1.0;
    face.low = {0.0, 0.0}; face.high = {3.0, 0.0};
    auto result = erf_sbm::limit_grouped(low_state, 2, 2, {face}, groups);
    ASSERT_EQ(result.accepted_faces.size(), 1U);
    EXPECT_EQ(result.limiter[0], 2.0/3.0);
    EXPECT_NEAR(result.updated_state[0], 0.0, 1.e-14);
    EXPECT_NEAR(result.updated_state[2], 4.0, 1.e-14);
    EXPECT_NEAR(result.updated_state[0]*2.0 + result.updated_state[2], 4.0, 1.e-14);
    EXPECT_EQ(result.accepted_faces[0].low[0], 2.0);
    EXPECT_THROW((void)erf_sbm::limit_grouped(low_state, 2, 2, {face, face}, groups),
                 std::invalid_argument);
    const auto chunked = erf_sbm::limit_grouped(low_state, 2, 2, {face}, groups, 1);
    EXPECT_NEAR(chunked.updated_state[0], result.updated_state[0], 1.e-14);
    EXPECT_NEAR(chunked.updated_state[2], result.updated_state[2], 1.e-14);
}

TEST(SBMP2, GroupedFCTHonorsCompleteTwoMomentGroupAndSubset)
{
    const auto layout = make_layout(2, MomentMode::TwoMoment, true);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    ASSERT_EQ(groups.size(), 2U);
    // Physical storage is (M,C); the host FCT reference applies the linear
    // cone constraints in that storage.  The candidate would drive both the
    // lower endpoint and the subset negative, so one common lambda must limit
    // every member of the group.
    std::vector<Real> low_state{0.5, 0.0, 0.5, 0.0, 0.25, 0.0,
                                0.5, 0.0, 0.5, 0.0, 0.25, 0.0};
    erf_sbm::FCTFaceTransfer face;
    face.left_cell = 0; face.right_cell = 1;
    face.left_volume = 1.0; face.right_volume = 1.0;
    face.low = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    face.high = {1.0, 0.0, 0.0, 0.0, 0.6, 0.0};
    auto result = erf_sbm::limit_grouped(low_state, 2, 6, {face}, groups);
    EXPECT_GE(result.limiter[0], 0.0);
    EXPECT_LE(result.limiter[0], 0.5 + 1.e-14);
    for (const auto& group : groups) {
        std::vector<Real> state(result.updated_state.begin(), result.updated_state.begin()+6);
        EXPECT_TRUE(group.admissible(state));
    }
}

TEST(SBMP2, CapabilityGateKeepsP2NegativeControlsFailClosed)
{
    erf_sbm::CapabilityInput input;
    input.p2_requested = true;
    const auto supported = erf_sbm::evaluate_p2_capabilities(input);
    EXPECT_TRUE(supported.supported);
    EXPECT_EQ(supported.qualification_status, "qualified");
    EXPECT_EQ(supported.qualified_flags, supported.flags);
    EXPECT_TRUE(supported.qualification_limitations.empty());

    input.moving_terrain = true;
    EXPECT_FALSE(erf_sbm::evaluate_p2_capabilities(input).supported);
    input.moving_terrain = false;
    input.diffusion = true;
    input.explicit_sbm_diffusion = false;
    EXPECT_FALSE(erf_sbm::evaluate_p2_capabilities(input).supported);
    input.explicit_sbm_diffusion = true;
    input.chunk_size = 0;
    EXPECT_FALSE(erf_sbm::evaluate_p2_capabilities(input).supported);
}

TEST(SBMP2, SyntheticSubsetPopulationSurvivesGroupedLimitAndAMRViews)
{
    erf_sbm::SpectralPopulationSpec liquid;
    liquid.population_id = 0;
    liquid.semantic_id = "liquid";
    liquid.phase = erf_sbm::PopulationPhase::Liquid;
    liquid.grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    liquid.grid.coordinate_units = "kg";
    liquid.grid.edges = {0.0, 1.0, 2.0};
    liquid.grid.pivots = {0.5, 1.5};
    liquid.mass_state_units = "kg m^-3";
    liquid.number_state_units = "m^-3";

    erf_sbm::SpectralPopulationSpec synthetic = liquid;
    synthetic.population_id = 1;
    synthetic.semantic_id = "synthetic_ice";
    synthetic.phase = erf_sbm::PopulationPhase::Ice;
    erf_sbm::AttachedPropertyDescriptor rime{
        "rime_like_mass", "rime_like_mass", "kg m^-3", 1,
        erf_sbm::PropertyKind::MassBoundedSubset,
        erf_sbm::SupportRequirement::PositiveMass,
        erf_sbm::PropertyRemapPolicy::CarrierBinConservative,
        true, false, 0.0, std::numeric_limits<Real>::quiet_NaN()};
    erf_sbm::SBMLayoutSpec spec;
    spec.populations = {liquid, synthetic};
    spec.liquid_projection = {0, 1};
    spec.attached_properties = {rime};
    const erf_sbm::SBMLayout layout(std::move(spec));
    const auto groups = erf_sbm::make_constraint_groups(layout);
    ASSERT_EQ(groups.size(), 4U);

    // ncomp = liquid mass(2) + synthetic mass(2) + rime subset(2).
    const std::vector<Real> low_state{
        1.0, 0.0, 2.0, 0.0, 0.5, 0.0,
        1.0, 0.0, 2.0, 0.0, 0.5, 0.0};
    erf_sbm::FCTFaceTransfer face;
    face.left_cell = 0; face.right_cell = 1;
    face.left_volume = 1.0; face.right_volume = 2.0;
    face.low.assign(6, 0.0);
    face.high = {0.0, 0.0, 1.0, 0.0, 0.75, 0.0};
    const auto limited = erf_sbm::limit_grouped(low_state, 2, 6, {face}, groups);
    for (int cell = 0; cell < 2; ++cell) {
        const std::vector<Real> state(limited.updated_state.begin() + cell*6,
                                      limited.updated_state.begin() + (cell+1)*6);
        for (const auto& group : groups) EXPECT_TRUE(group.admissible(state));
    }

    const auto restricted = erf_sbm::volume_weighted_restrict(
        limited.updated_state, 2, 1, 6, {0, 0}, {1.0, 2.0}, {3.0});
    const auto prolonged = erf_sbm::piecewise_constant_prolong(
        restricted, 1, 2, 6, {0, 0});
    for (const auto& state_vector : {restricted, prolonged}) {
        const int count = static_cast<int>(state_vector.size()) / 6;
        for (int cell = 0; cell < count; ++cell) {
            const std::vector<Real> state(state_vector.begin() + cell*6,
                                          state_vector.begin() + (cell+1)*6);
            for (const auto& group : groups) EXPECT_TRUE(group.admissible(state));
        }
    }
    // The synthetic property remains bounded by the synthetic carrier mass;
    // liquid projection is still defined exclusively by population 0.
    EXPECT_GE(restricted[2] - restricted[4], 0.0);
    EXPECT_EQ(layout.bulk_projection().rules().size(), 2U);
}

TEST(SBMP2, FCTStageWeightsAndAcceptedCorrectionAreExact)
{
    EXPECT_EQ(erf_sbm::stage_weight_contract(erf_auxiliary::IntegrationMethod::CompressibleRK3, 1, true).completed_ledger_weight, 1.0);
    EXPECT_EQ(erf_sbm::stage_weight_contract(erf_auxiliary::IntegrationMethod::CompressibleRK3, 1, false).completed_ledger_weight, 0.0);
    EXPECT_EQ(erf_sbm::stage_weight_contract(erf_auxiliary::IntegrationMethod::AnelasticHeun, 0, false).completed_ledger_weight, 0.5);
    EXPECT_EQ(erf_sbm::stage_weight_contract(erf_auxiliary::IntegrationMethod::AnelasticHeun, 1, true).completed_ledger_weight, 0.5);
    const auto correction = erf_sbm::accepted_stage_correction({1.0, -2.0}, {3.0, 2.0}, 0.25);
    EXPECT_DOUBLE_EQ(correction[0], 0.5);
    EXPECT_DOUBLE_EQ(correction[1], 1.0);
}

TEST(SBMP2, DensityWeightedDiffusionUsesIntensiveRatioAndPhysicalGeometry)
{
    erf_sbm::DiffusionFace face;
    face.left_cell = 0; face.right_cell = 1; face.area = 3.0; face.distance = 2.0;
    face.rho_left = 2.0; face.rho_right = 4.0; face.rho_face = 3.0; face.coefficient = 0.5;
    const auto result = erf_sbm::explicit_two_point_diffusion({2.0, 8.0}, 2, 1,
                                                               {4.0, 1.0}, {face}, 0.2);
    // X/rho is 1 and 2; I=-A dt rho_f K grad(X/rho)=-0.45.
    ASSERT_EQ(result.integrated_transfers.size(), 1U);
    EXPECT_NEAR(result.integrated_transfers[0], -0.45, 1.e-14);
    EXPECT_NEAR(result.updated_state[0], 2.1125, 1.e-14);
    EXPECT_NEAR(result.updated_state[1], 7.55, 1.e-14);
    const auto uniform = erf_sbm::explicit_two_point_diffusion({2.0, 4.0}, 2, 1,
                                                                {1.0, 1.0}, {face}, 1.0);
    EXPECT_NEAR(uniform.integrated_transfers[0], 0.0, 1.e-14);
    EXPECT_GT(erf_sbm::admissible_explicit_timestep({2.0, 8.0}, 2, 1,
                                                    {4.0, 1.0}, {face}), 0.0);
}

TEST(SBMP2, ProductionCombinedAdvectionDiffusionFailsClosed)
{
    const auto layout = make_layout(2, MomentMode::OneMoment);
    erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                  {AMREX_D_DECL(2.0, 1.0, 1.0)});
    const std::array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    const Geometry geometry(domain, &real_box, amrex::CoordSys::cartesian,
                            periodicity.data());
    manager.define_level(0, boxes, dm, 2);
    for (amrex::MFIter mfi(manager.output(0)); mfi.isValid(); ++mfi) {
        const auto state = manager.output(0).array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            state(i,j,k,0) = i == 0 ? Real(1.0) : Real(0.0);
            state(i,j,k,1) = Real(0.5);
        });
    }
    manager.begin_step(0, 0.0);

    amrex::MultiFab rho(boxes, dm, 1, 2);
    rho.setVal(1.0);
    rho.FillBoundary(geometry.periodicity());
    amrex::MultiFab xflux(amrex::convert(boxes, IntVect(1, 0, 0)), dm, 1, 0);
    amrex::MultiFab yflux(amrex::convert(boxes, IntVect(0, 1, 0)), dm, 1, 0);
    amrex::MultiFab zflux(amrex::convert(boxes, IntVect(0, 0, 1)), dm, 1, 0);
    xflux.setVal(1.0);
    yflux.setVal(0.0);
    zflux.setVal(0.0);
    amrex::MultiFab core(boxes, dm, RhoQ3_comp + 1, 0);
    erf_auxiliary::AuxiliaryFaceTransfer stage_flux;
    stage_flux.define(boxes, dm, layout.ncomp(), 0);
    const auto context = erf_auxiliary::make_compressible_stage(
        2, 0.0, 0.4, 0.8, 0.8, nullptr, nullptr);

    std::string diagnostic;
    try {
        erf_sbm::advance_stage(manager, layout, context, rho, core,
                               xflux, yflux, zflux, geometry, stage_flux,
                               erf_sbm::TransportMethod::GroupedFCT_WENOZ3,
                               0, Real(0.2), 1);
    } catch (const std::exception& error) {
        diagnostic = error.what();
    }
    EXPECT_NE(diagnostic.find("combined low-order advection+diffusion"), std::string::npos);
}

TEST(SBMP2, ProductionCombinedDemandUsesPreStageBaselineForAllStageContracts)
{
    const auto run_case = [](const Real advection_demand, const Real diffusion,
                             const bool anelastic, const int stage,
                             const MomentMode mode = MomentMode::OneMoment,
                             const bool with_property = false,
                             const int chunk_size = 1) {
        const auto layout = make_layout(2, mode, with_property);
        erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
        const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
        const BoxArray boxes(domain);
        const DistributionMapping dm(boxes);
        const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                      {AMREX_D_DECL(2.0, 1.0, 1.0)});
        const std::array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
        const Geometry geometry(domain, &real_box, amrex::CoordSys::cartesian,
                                periodicity.data());
        manager.define_level(0, boxes, dm, 2);
        manager.output(0).setVal(Real(0.0));
        for (amrex::MFIter mfi(manager.output(0)); mfi.isValid(); ++mfi) {
            const auto state = manager.output(0).array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (mode == MomentMode::TwoMoment) {
                    state(i,j,k,0) = i == 0 ? Real(0.5) : Real(0.0);
                    state(i,j,k,2) = i == 0 ? Real(1.0) : Real(0.0);
                    if (with_property) state(i,j,k,4) = i == 0 ? Real(0.25) : Real(0.0);
                } else {
                    state(i,j,k,0) = i == 0 ? Real(1.0) : Real(0.0);
                }
            });
        }
        manager.begin_step(0, 0.0);
        if (anelastic && stage > 0) manager.evaluation(0).setVal(Real(0.0));
        if (anelastic && stage > 0) {
            for (amrex::MFIter mfi(manager.evaluation(0)); mfi.isValid(); ++mfi) {
                const auto state = manager.evaluation(0).array(mfi);
                amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    if (mode == MomentMode::TwoMoment) {
                        state(i,j,k,0) = i == 0 ? Real(0.5) : Real(0.0);
                        state(i,j,k,2) = i == 0 ? Real(1.0) : Real(0.0);
                        if (with_property) state(i,j,k,4) = i == 0 ? Real(0.25) : Real(0.0);
                    } else {
                        state(i,j,k,0) = i == 0 ? Real(1.0) : Real(0.0);
                    }
                });
            }
        }

        amrex::MultiFab rho(boxes, dm, 1, 2);
        rho.setVal(Real(1.0));
        rho.FillBoundary(geometry.periodicity());
        amrex::MultiFab xflux(amrex::convert(boxes, IntVect(1, 0, 0)), dm, 1, 0);
        amrex::MultiFab yflux(amrex::convert(boxes, IntVect(0, 1, 0)), dm, 1, 0);
        amrex::MultiFab zflux(amrex::convert(boxes, IntVect(0, 0, 1)), dm, 1, 0);
        xflux.setVal(Real(0.0));
        yflux.setVal(Real(0.0));
        zflux.setVal(Real(0.0));
        for (amrex::MFIter mfi(xflux); mfi.isValid(); ++mfi) {
            const auto flux = xflux.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                flux(i,j,k) = i == 1 ? advection_demand : Real(0.0);
            });
        }
        amrex::MultiFab core(boxes, dm, RhoQ3_comp + 1, 0);
        erf_auxiliary::AuxiliaryFaceTransfer stage_flux;
        stage_flux.define(boxes, dm, layout.ncomp(), 0);
        const auto context = anelastic ?
            erf_auxiliary::make_anelastic_stage(stage, 0.0, stage == 0 ? 0.0 : 1.0,
                                                1.0, 1.0, nullptr, nullptr) :
            erf_auxiliary::make_compressible_stage(2, 0.0, 0.0, 1.0, 1.0, nullptr, nullptr);
        std::string diagnostic;
        try {
            erf_sbm::advance_stage(manager, layout, context, rho, core,
                                   xflux, yflux, zflux, geometry, stage_flux,
                                   erf_sbm::TransportMethod::GroupedFCT_WENOZ3,
                                   0, diffusion, chunk_size);
        } catch (const std::exception& error) {
            diagnostic = error.what();
        }
        return diagnostic;
    };

    // B=1, D^out=.6 is valid even though the post-low-order state is .4.
    for (const auto contract : {std::pair<bool,int>{false,0}, {true,0}, {true,1}}) {
        SCOPED_TRACE(std::string(contract.first ? "anelastic" : "compressible") +
                     " stage=" + std::to_string(contract.second));
        const Real stage_scale = contract.first && contract.second > 0 ? Real(2.0) : Real(1.0);
        EXPECT_TRUE(run_case(Real(0.6)*stage_scale, Real(0.0), contract.first, contract.second).empty());
        const auto invalid = run_case(Real(1.2)*stage_scale, Real(0.0), contract.first, contract.second);
        EXPECT_NE(invalid.find("combined low-order advection+diffusion"), std::string::npos) << invalid;
    }

    // Advection and diffusion are each admissible, but their combined demand
    // is not.  The second case deliberately exceeds the final low state while
    // remaining below the pre-stage baseline and must therefore pass.
    for (const auto contract : {std::pair<bool,int>{false,0}, {true,0}, {true,1}}) {
        const Real stage_scale = contract.first && contract.second > 0 ? Real(2.0) : Real(1.0);
        const Real combined_fail_advection = contract.first && contract.second > 0 ? Real(1.9) : Real(0.8);
        const Real combined_pass_advection = Real(0.4) * stage_scale;
        const auto invalid = run_case(combined_fail_advection, Real(0.15), contract.first, contract.second);
        EXPECT_NE(invalid.find("combined low-order advection+diffusion"), std::string::npos) << invalid;
        EXPECT_TRUE(run_case(combined_pass_advection, Real(0.15), contract.first, contract.second).empty());
    }

    for (const auto mode : {MomentMode::TwoMoment}) {
        const auto valid = run_case(Real(0.6), Real(0.0), false, 0, mode, true, 1);
        EXPECT_TRUE(valid.empty()) << valid;
        const auto invalid = run_case(Real(1.2), Real(0.0), false, 0, mode, true, 1);
        EXPECT_NE(invalid.find("combined low-order advection+diffusion"), std::string::npos)
            << invalid;
    }
}

TEST(SBMP2, BoundaryBudgetsHaveNoWallSinkAndValidateInflow)
{
    const auto layout = make_layout(2, MomentMode::OneMoment);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    erf_sbm::BoundaryDescriptor wall;
    wall.kind = erf_sbm::BoundaryKind::ImpermeableWall;
    const auto wall_budget = erf_sbm::make_boundary_budget(wall, {2.0, 3.0}, groups, {0.0, 0.0}, {0.0, 0.0}, 5.0, 0.1);
    EXPECT_EQ(wall_budget.auxiliary[0], 0.0);
    erf_sbm::BoundaryDescriptor inflow;
    inflow.kind = erf_sbm::BoundaryKind::PrescribedSpectralInflow;
    inflow.prescribed_state = {1.0, 2.0};
    EXPECT_TRUE(erf_sbm::validate_prescribed_inflow(inflow, groups));
    inflow.prescribed_state[0] = -1.0;
    EXPECT_FALSE(erf_sbm::validate_prescribed_inflow(inflow, groups));
}

TEST(SBMP2, AMRRestrictionProlongationAndRegisterUnitsAreConservative)
{
    const auto restricted = erf_sbm::volume_weighted_restrict({1.0, 3.0, 5.0, 7.0}, 4, 2, 1,
                                                               {0, 0, 1, 1}, {1.0, 1.0, 2.0, 2.0}, {2.0, 4.0});
    EXPECT_DOUBLE_EQ(restricted[0], 2.0);
    EXPECT_DOUBLE_EQ(restricted[1], 6.0);
    const auto injected = erf_sbm::piecewise_constant_prolong(restricted, 2, 4, 1, {0, 0, 1, 1});
    EXPECT_EQ(injected[0], 2.0); EXPECT_EQ(injected[3], 6.0);
    EXPECT_DOUBLE_EQ(erf_sbm::register_flux_from_integrated_transfer(12.0, 3.0, 2.0), 2.0);
    EXPECT_THROW((void)erf_sbm::register_flux_from_integrated_transfer(1.0, 0.0, 1.0), std::invalid_argument);
}

TEST(SBMP2, IndependentInterfaceOracleCatchesRefluxSignAreaAndTimeErrors)
{
    const std::vector<Real> coarse{2.0};
    const std::vector<Real> fine{1.0, 1.0, 1.0, 1.0};
    const auto good = erf_sbm::check_interface_transfer(
        coarse, fine, 4.0, 1.0, 0.5, 0.5, -2.0,
        64*std::numeric_limits<Real>::epsilon());
    EXPECT_TRUE(good.passes);
    EXPECT_DOUBLE_EQ(good.coarse_transfer, 4.0);
    EXPECT_DOUBLE_EQ(good.fine_transfer, 2.0);
    EXPECT_DOUBLE_EQ(good.mismatch, -2.0);

    const auto missing_reflux = erf_sbm::check_interface_transfer(
        coarse, fine, 4.0, 1.0, 0.5, 0.5, 0.0,
        64*std::numeric_limits<Real>::epsilon());
    EXPECT_FALSE(missing_reflux.passes);
    const auto wrong_sign = erf_sbm::check_interface_transfer(
        coarse, fine, 4.0, 1.0, 0.5, 0.5, 2.0,
        64*std::numeric_limits<Real>::epsilon());
    EXPECT_FALSE(wrong_sign.passes);
    const auto wrong_area = erf_sbm::check_interface_transfer(
        coarse, fine, 4.0, 2.0, 0.5, 0.5, -6.0,
        64*std::numeric_limits<Real>::epsilon());
    EXPECT_FALSE(wrong_area.passes);
    const auto wrong_time = erf_sbm::check_interface_transfer(
        coarse, fine, 4.0, 1.0, 0.5, 1.0, -2.0,
        64*std::numeric_limits<Real>::epsilon());
    EXPECT_FALSE(wrong_time.passes);
    EXPECT_THROW((void)erf_sbm::check_interface_transfer(
        coarse, fine, 4.0, 1.0, 0.5, 0.0, -2.0, 1.e-12), std::invalid_argument);
}

TEST(SBMP2, PostRefluxFailsClosedWithDiagnosticsAndNeverClips)
{
    const auto layout = make_layout(2, MomentMode::OneMoment);
    const auto groups = erf_sbm::make_constraint_groups(layout);
    const auto good = erf_sbm::validate_post_reflux({1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, 2, 1, 2, groups);
    EXPECT_TRUE(good.admissible);
    const auto bad = erf_sbm::validate_post_reflux({1.0, 0.0}, {-2.0, 0.0}, {-1.0, 0.0}, 2, 1, 2, groups);
    EXPECT_FALSE(bad.admissible);
    EXPECT_EQ(bad.level, 2); EXPECT_EQ(bad.cell, 0);
    EXPECT_FALSE(bad.group.empty()); EXPECT_FALSE(bad.constraint.empty());
}

TEST(SBMP2, RestartSchemaAndProjectionComparisonAreStrict)
{
    const auto layout = make_layout(4, MomentMode::TwoMoment);
    const auto schema = erf_sbm::make_checkpoint_schema(layout, "complete-groups-v1", "WENO_Z3+FCT-v1", "gamma-k-v1");
    EXPECT_TRUE(erf_sbm::compare_checkpoint_schema(schema, schema).empty());
    auto altered = schema;
    altered.moment_modes += "changed";
    EXPECT_NE(erf_sbm::compare_checkpoint_schema(schema, altered).find("moment_modes"), std::string::npos);
    EXPECT_TRUE(erf_sbm::compare_projection(1.0, 1.0 + 1.e-14, 1.0, 8));
    EXPECT_FALSE(erf_sbm::compare_projection(1.0, 1.0 + 1.e-4, 1.0, 8));
}

TEST(SBMP2, AuxiliaryManagerOwnsAMRCreateProlongAverageDownRemakeAndDestroy)
{
    const auto layout = make_layout(2, MomentMode::OneMoment);
    erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
    const Box coarse_domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray coarse_boxes(coarse_domain);
    const DistributionMapping coarse_dm(coarse_boxes);
    const IntVect ref_ratio(2, 2, 2);
    const Box fine_domain = amrex::refine(coarse_domain, ref_ratio);
    BoxArray fine_boxes(fine_domain);
    const DistributionMapping fine_dm(fine_boxes);
    const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                  {AMREX_D_DECL(4.0, 2.0, 2.0)});
    const std::array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    const Geometry coarse_geometry(coarse_domain, &real_box, amrex::CoordSys::cartesian,
                                   periodicity.data());
    const Geometry fine_geometry(fine_domain, &real_box, amrex::CoordSys::cartesian,
                                 periodicity.data());
    manager.define_level(0, coarse_boxes, coarse_dm, 2);
    manager.define_level(1, fine_boxes, fine_dm, 2);
    manager.output(0).setVal(3.0);
    manager.old(0).setVal(3.0);
    manager.evaluation(0).setVal(3.0);
    manager.prolong_from_coarse(0, 1, coarse_geometry, fine_geometry, ref_ratio);
    EXPECT_DOUBLE_EQ(manager.output(1).min(0), 3.0);
    manager.average_down_to(0, 1, ref_ratio);
    EXPECT_DOUBLE_EQ(manager.output(0).min(0), 3.0);

    fine_boxes.maxSize(1);
    const DistributionMapping remade_dm(fine_boxes);
    manager.remake_level(1, fine_boxes, remade_dm, 2, fine_geometry.periodicity());
    EXPECT_TRUE(manager.has_level(1));
    EXPECT_GT(manager.resident_bytes(), std::size_t(0));
    manager.destroy_level(1);
    EXPECT_FALSE(manager.has_level(1));
}

TEST(SBMP2, AuxiliaryStageFillPatchAndRemakeUseAuthoritativeCoarseSpectrum)
{
    const auto layout = make_layout(2, MomentMode::OneMoment);
    erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
    const Box coarse_domain(IntVect(0, 0, 0), IntVect(1, 1, 1));
    const BoxArray coarse_boxes(coarse_domain);
    const DistributionMapping coarse_dm(coarse_boxes);
    const IntVect ref_ratio(2, 2, 2);
    const Box fine_domain = amrex::refine(coarse_domain, ref_ratio);
    const BoxArray initial_fine_boxes(Box(IntVect(0, 0, 0), IntVect(1, 3, 3)));
    const BoxArray remade_fine_boxes(fine_domain);
    const DistributionMapping initial_fine_dm(initial_fine_boxes);
    const DistributionMapping remade_fine_dm(remade_fine_boxes);
    const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                  {AMREX_D_DECL(4.0, 4.0, 4.0)});
    const std::array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    const Geometry coarse_geometry(coarse_domain, &real_box, amrex::CoordSys::cartesian,
                                   periodicity.data());
    const Geometry fine_geometry(fine_domain, &real_box, amrex::CoordSys::cartesian,
                                 periodicity.data());
    manager.define_level(0, coarse_boxes, coarse_dm, 2);
    manager.define_level(1, initial_fine_boxes, initial_fine_dm, 2);
    manager.output(0).setVal(2.0);
    manager.begin_step(0, 0.0);
    manager.output(0).setVal(4.0);
    manager.accept_stage(0, 1.0);
    manager.evaluation(1).setVal(0.0);

    manager.fill_stage_from_coarse(0, 1, 0.5, coarse_geometry, fine_geometry, ref_ratio);
    Real stage_ghost_value = 0.0;
    for (amrex::MFIter mfi(manager.evaluation(1)); mfi.isValid(); ++mfi) {
        stage_ghost_value = manager.evaluation(1).const_array(mfi)(2, 1, 1, 0);
    }
    amrex::ParallelDescriptor::ReduceRealMax(stage_ghost_value);
    EXPECT_DOUBLE_EQ(stage_ghost_value, 3.0);

    manager.old(1).setVal(7.0);
    manager.evaluation(1).setVal(7.0);
    manager.output(1).setVal(7.0);
    manager.remake_level_from_coarse(1, remade_fine_boxes, remade_fine_dm, 2,
                                     fine_geometry.periodicity(), 0, coarse_geometry,
                                     fine_geometry, ref_ratio, 0.5);
    EXPECT_DOUBLE_EQ(manager.output(1).min(0), 3.0);
    EXPECT_DOUBLE_EQ(manager.output(1).max(0), 7.0);
}

TEST(SBMP2, CompleteConstraintValidationCoversAttachedProperties)
{
    const auto layout = make_layout(2, MomentMode::TwoMoment, true);
    erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
    const Box domain(IntVect(0, 0, 0), IntVect(0, 0, 0));
    const BoxArray boxes(domain);
    const DistributionMapping dm(boxes);
    manager.define_level(0, boxes, dm, 2);
    auto& state = manager.output(0);
    state.setVal(0.0);
    const auto& population = layout.populations().front();
    for (int b = 0; b < 2; ++b) {
        state.setVal(static_cast<Real>(b) + 0.5, population.mass_offset + b, 1);
        state.setVal(1.0, population.number_offset + b, 1);
        state.setVal(0.25, layout.property_offset(0) + b, 1);
    }
    erf_sbm::validate_admissible_state(manager, layout, 0);
    state.setVal(-1.0, layout.property_offset(0), 1);
    EXPECT_THROW(erf_sbm::validate_admissible_state(manager, layout, 0), std::exception);
}

TEST(SBMP2, WENOZ3ConvergenceBeatsDonorOnPeriodicSmoothOperator)
{
    std::ofstream evidence("/private/tmp/erf_sbm_p2_weno_convergence.csv");
    ASSERT_TRUE(evidence.good());
    evidence << "N,weno_max_error,donor_max_error,weno_order,donor_order\n";
    std::vector<Real> weno_errors;
    std::vector<Real> donor_errors;
    for (const int ncell : {16, 32, 64, 128}) {
        const Box domain(IntVect(0, 0, 0), IntVect(ncell-1, 0, 0));
        const BoxArray boxes(domain);
        const DistributionMapping dm(boxes);
        const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                      {AMREX_D_DECL(1.0, 1.0, 1.0)});
        const std::array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
        const Geometry geometry(domain, &real_box, amrex::CoordSys::cartesian,
                                periodicity.data());
        amrex::MultiFab ratio(boxes, dm, 1, 2);
        for (amrex::MFIter mfi(ratio); mfi.isValid(); ++mfi) {
            const auto values = ratio.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real x = (static_cast<Real>(i) + Real(0.5)) / static_cast<Real>(ncell);
                values(i,j,k) = Real(2.0) + std::sin(Real(6.2831853071795864769) * x);
            });
        }
        ratio.FillBoundary(geometry.periodicity());
        const BoxArray face_boxes = amrex::convert(boxes, IntVect(1, 0, 0));
        amrex::MultiFab errors(face_boxes, dm, 2, 0);
        for (amrex::MFIter mfi(errors); mfi.isValid(); ++mfi) {
            const auto input = ratio.const_array(mfi);
            const auto result = errors.array(mfi);
            const Box face_box = mfi.validbox();
            amrex::ParallelFor(face_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                WENO_Z3 weno(input, Real(0.0));
                Real reconstructed = Real(0.0);
                weno.InterpolateInX(i, j, k, 0, reconstructed, Real(1.0));
                const Real donor = input(i-1, j, k);
                const Real exact = Real(2.0) + std::sin(
                    Real(6.2831853071795864769) * static_cast<Real>(i) / static_cast<Real>(ncell));
                result(i,j,k,0) = std::abs(reconstructed - exact);
                result(i,j,k,1) = std::abs(donor - exact);
            });
        }
        weno_errors.push_back(errors.norm0(0));
        donor_errors.push_back(errors.norm0(1));
    }
    for (std::size_t i = 0; i < weno_errors.size(); ++i) {
        const Real weno_order = i == 0 ? Real(0.0) :
            std::log(weno_errors[i-1] / weno_errors[i]) / std::log(Real(2.0));
        const Real donor_order = i == 0 ? Real(0.0) :
            std::log(donor_errors[i-1] / donor_errors[i]) / std::log(Real(2.0));
        evidence << (16 << i) << ',' << std::setprecision(17) << weno_errors[i] << ','
                 << donor_errors[i] << ',' << weno_order << ',' << donor_order << '\n';
        // Source archaeology shows that this helper consumes pointwise cell
        // center values: its two substencils reconstruct the face between
        // neighboring centers.  This oracle therefore uses point values and
        // compares against the analytic face point value; it does not make a
        // finite-volume cell-average claim.
        if (i > 1) EXPECT_GT(weno_order, Real(1.8));
        if (i > 0) EXPECT_LT(weno_errors[i], donor_errors[i]);
    }
}

TEST(SBMP2, FullGroupedTransportHasSmoothManufacturedConvergence)
{
    std::ofstream evidence("/private/tmp/erf_sbm_p2_full_transport_convergence.csv");
    ASSERT_TRUE(evidence.good());
    evidence << "N,weno_error,weno_order,donor_error,donor_order\n";
    const auto run = [](const int ncell, const erf_sbm::TransportMethod method) {
        const auto layout = make_layout(2, MomentMode::OneMoment);
        erf_auxiliary::AuxiliaryStateManager manager(layout.auxiliary_layout());
        const Box domain(IntVect(0, 0, 0), IntVect(ncell-1, 0, 0));
        const BoxArray boxes(domain);
        const DistributionMapping dm(boxes);
        const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                      {AMREX_D_DECL(1.0, 1.0, 1.0)});
        const std::array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
        const Geometry geometry(domain, &real_box, amrex::CoordSys::cartesian,
                                periodicity.data());
        manager.define_level(0, boxes, dm, 2);
        for (amrex::MFIter mfi(manager.output(0)); mfi.isValid(); ++mfi) {
            const auto state = manager.output(0).array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real x = (static_cast<Real>(i) + Real(0.5)) / static_cast<Real>(ncell);
                state(i,j,k,0) = Real(2.0) + std::sin(Real(6.2831853071795864769) * x);
                state(i,j,k,1) = Real(1.0);
            });
        }
        manager.begin_step(0, 0.0);
        const Real dt = Real(1.0) / static_cast<Real>(ncell);
        amrex::MultiFab rho(boxes, dm, 1, 2);
        rho.setVal(Real(1.0));
        rho.FillBoundary(geometry.periodicity());
        amrex::MultiFab xflux(amrex::convert(boxes, IntVect(1, 0, 0)), dm, 1, 0);
        amrex::MultiFab yflux(amrex::convert(boxes, IntVect(0, 1, 0)), dm, 1, 0);
        amrex::MultiFab zflux(amrex::convert(boxes, IntVect(0, 0, 1)), dm, 1, 0);
        xflux.setVal(Real(0.25));
        yflux.setVal(Real(0.0));
        zflux.setVal(Real(0.0));
        amrex::MultiFab core(boxes, dm, RhoQ3_comp + 1, 0);
        erf_auxiliary::AuxiliaryFaceTransfer stage_flux;
        stage_flux.define(boxes, dm, layout.ncomp(), 0);
        const auto context = erf_auxiliary::make_compressible_stage(
            2, 0.0, 0.0, dt, dt, nullptr, nullptr);
        erf_sbm::advance_stage(manager, layout, context, rho, core,
                               xflux, yflux, zflux, geometry, stage_flux,
                               method, 0, Real(0.0), 1);
        amrex::MultiFab error(boxes, dm, 1, 0);
        for (amrex::MFIter mfi(error); mfi.isValid(); ++mfi) {
            const auto result = error.array(mfi);
            const auto state = manager.output(0).const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real x = (static_cast<Real>(i) + Real(0.5)) / static_cast<Real>(ncell);
                const Real exact = Real(2.0) + std::sin(
                    Real(6.2831853071795864769) * (x - Real(0.25)*dt));
                result(i,j,k) = std::abs(state(i,j,k,0) - exact);
            });
        }
        return error.norm0(0);
    };

    std::vector<Real> weno_errors;
    std::vector<Real> donor_errors;
    for (const int ncell : {16, 32, 64, 128}) {
        weno_errors.push_back(run(ncell, erf_sbm::TransportMethod::GroupedFCT_WENOZ3));
        donor_errors.push_back(run(ncell, erf_sbm::TransportMethod::DonorCell));
    }
    for (std::size_t i = 0; i < weno_errors.size(); ++i) {
        const Real weno_order = i == 0 ? Real(0.0) :
            std::log(weno_errors[i-1] / weno_errors[i]) / std::log(Real(2.0));
        const Real donor_order = i == 0 ? Real(0.0) :
            std::log(donor_errors[i-1] / donor_errors[i]) / std::log(Real(2.0));
        evidence << (16 << i) << ',' << std::setprecision(17)
                 << weno_errors[i] << ',' << weno_order << ','
                 << donor_errors[i] << ',' << donor_order << '\n';
        if (i > 0) {
            EXPECT_GT(weno_order, Real(1.5));
            EXPECT_LT(donor_order, Real(1.5));
            EXPECT_LT(weno_errors[i], donor_errors[i]);
        }
    }
}

} // namespace
