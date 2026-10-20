#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>

#include "AuxiliaryState/ERF_AuxiliaryStage.H"
#include "ERF_IndexDefines.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;
using amrex::MultiFab;
using amrex::Real;

BoxArray
project_to_xy (const BoxArray& cell_ba)
{
    amrex::BoxList boxes = cell_ba.boxList();
    for (auto& box : boxes)
        box.setRange(2, 0);
    return BoxArray(std::move(boxes));
}

erf_sbm::SBMLayout
make_parallel_layout (const erf_sbm::MomentMode mode)
{
    erf_sbm::SpectralGridSpec grid;
    grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    grid.coordinate_units = "kg particle^-1";
    grid.edges = {Real(0.1), Real(0.5), Real(1.0)};
    grid.pivots = {Real(0.3), Real(0.75)};

    erf_sbm::SpectralPopulationSpec population;
    population.population_id = 0;
    population.semantic_id = "liquid";
    population.phase = erf_sbm::PopulationPhase::Liquid;
    population.grid = std::move(grid);
    population.moment_mode = mode;

    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(std::move(population));
    spec.liquid_projection = {0, 1};
    return erf_sbm::SBMLayout(std::move(spec));
}

std::vector<Real>
run_decomposition (const int max_grid_size,
                   const erf_sbm::MomentMode mode)
{
    constexpr int nx = 16;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const Box domain(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    const Geometry geom(domain, &physical, amrex::CoordSys::cartesian,
                        periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    const DistributionMapping dm(ba);
    auto layout = make_parallel_layout(mode);

    erf_sbm::SBMStateManager manager(layout, 1);
    manager.define(0, ba, dm, 0.0);
    erf_sbm::SBMTransport transport(manager.layout(), 1, 1);
    transport.define(0, ba, dm);

    MultiFab detj(ba, dm, 1, 0);
    const BoxArray map_ba = project_to_xy(ba);
    MultiFab mx(map_ba, dm, 1, 0);
    MultiFab my(map_ba, dm, 1, 0);
    detj.setVal(Real(1.0));
    mx.setVal(Real(1.0));
    my.setVal(Real(1.0));
    std::string diagnostic;
    if (!transport.rebuild_static_measure(0, detj, mx, my, diagnostic)) {
        ADD_FAILURE() << diagnostic;
        return {};
    }

    MultiFab conserved_anchor(ba, dm, 3, 0);
    MultiFab conserved_input(ba, dm, 3, 0);
    MultiFab conserved_target(ba, dm, 3, 0);
    conserved_anchor.setVal(Real(0.0));
    conserved_anchor.setVal(Real(1.0), Rho_comp, 1, 0);
    MultiFab::Copy(conserved_input, conserved_anchor, 0, 0, 3, 0);
    MultiFab::Copy(conserved_target, conserved_anchor, 0, 0, 3, 0);

    auto& spectrum = manager.new_state_for_initialization(0);
    spectrum.setVal(Real(0.0));
    const auto& population = manager.layout().populations().front();
    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = spectrum.array(mfi);
        const int mass0 = population.mass_offset;
        const int mass1 = mass0 + 1;
        const int number0 = population.number_offset;
        const int number1 = number0 < 0 ? -1 : number0 + 1;
        const bool two_moment = mode == erf_sbm::MomentMode::TwoMoment;
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (two_moment) {
                    // Deliberately sharp per-bin number and mean-mass changes
                    // exercise endpoint reconstruction and grouped limiting.
                    const int pattern = (17 * i + 31 * j + 13 * k) % 7;
                    const Real n0 = pattern == 0 ? Real(0.01) : Real(0.00001);
                    const Real mean0 = pattern % 2 == 0 ? Real(0.11) : Real(0.45);
                    const Real n1 = pattern == 1 ? Real(0.02) : Real(0.00002);
                    const Real mean1 = pattern % 3 == 0 ? Real(0.51) : Real(0.99);
                    state(i, j, k, mass0) = n0 * mean0;
                    state(i, j, k, mass1) = n1 * mean1;
                    state(i, j, k, number0) = n0;
                    state(i, j, k, number1) = n1;
                } else {
                    const Real phase = static_cast<Real>(i % 8);
                    state(i, j, k, mass0) =
                        Real(0.001) + Real(0.0002) * phase;
                    state(i, j, k, mass1) =
                        Real(0.002) +
                        Real(0.0001) * static_cast<Real>((3 * i + j + k) % 7);
                }
            });
    }

    MultiFab avg_xmom(amrex::convert(ba, IntVect::TheDimensionVector(0)), dm, 1,
                      0);
    MultiFab avg_ymom(amrex::convert(ba, IntVect::TheDimensionVector(1)), dm, 1,
                      0);
    MultiFab avg_zmom(amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1,
                      0);
    avg_xmom.setVal(mode == erf_sbm::MomentMode::TwoMoment ? Real(0.13)
                                                          : Real(0.05));
    avg_ymom.setVal(mode == erf_sbm::MomentMode::TwoMoment ? Real(0.11)
                                                          : Real(0.0));
    avg_zmom.setVal(mode == erf_sbm::MomentMode::TwoMoment ? Real(0.09)
                                                          : Real(0.0));

    const double dt = mode == erf_sbm::MomentMode::TwoMoment ? 0.25 : 0.01;
    if (!manager.begin_step(0, 0.0, diagnostic)) {
        ADD_FAILURE() << diagnostic;
        return {};
    }
    transport.advance_stage_from_host(
        0, erf_auxiliary::HostIntegrator::CompressibleRK3,
                            0, 0.0, 0.0, dt, dt, manager, conserved_anchor,
                            conserved_input, conserved_target, avg_xmom,
                            avg_ymom, avg_zmom, geom, 1, 2);

    const auto& accepted = manager.new_state(0);
    MultiFab signatures(ba, dm, 2 * layout.ncomp(), 0);
    for (amrex::MFIter mfi(accepted); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = accepted.const_array(mfi);
        const auto out = signatures.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                const Real value = state(i, j, k, n);
                const Real weight = Real(1.0) + Real(0.001) * i +
                                    Real(0.0001) * j + Real(0.00001) * k;
                out(i, j, k, n) = weight * value;
                out(i, j, k, ncomp + n) = value * value;
            });
    }

    std::vector<Real> result;
    for (int component = 0; component < layout.ncomp(); ++component) {
        result.push_back(accepted.sum(component));
        result.push_back(signatures.sum(component));
        result.push_back(signatures.sum(layout.ncomp() + component));
    }
    return result;
}

BoxArray
completed_ledger_boxes ()
{
    BoxArray boxes(Box(IntVect(0), IntVect(7)));
    boxes.maxSize(2);
    return boxes;
}

DistributionMapping
completed_ledger_mapping (const BoxArray& boxes)
{
    amrex::Vector<int> owners(static_cast<std::size_t>(boxes.size()));
    for (int box = 0; box < boxes.size(); ++box) {
        owners[box] = box % amrex::ParallelDescriptor::NProcs();
    }
    return DistributionMapping(std::move(owners));
}

struct CompletedLedgerParallelFixture
{
    BoxArray boxes;
    DistributionMapping mapping;
    erf_auxiliary::MappedFaceFluxRate rate;
    erf_auxiliary::CompletedStepFluxLedger ledger;

    CompletedLedgerParallelFixture ()
        : boxes(completed_ledger_boxes()),
          mapping(completed_ledger_mapping(boxes))
    {
        rate.define(boxes, mapping, 2, 0);
        rate.setVal(Real(1.0));
        ledger.define(boxes, mapping, 2);
    }
};

void expect_discarded_ledger (
    const erf_auxiliary::CompletedStepFluxLedger& ledger)
{
    EXPECT_FALSE(ledger.step_complete());
    EXPECT_FALSE(ledger.step_active());
    EXPECT_EQ(ledger.next_stage(), 0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int component = 0; component < ledger.integrated_flux().nComp();
             ++component) {
            EXPECT_EQ(ledger.integrated_flux().dir(dir).norm0(component),
                      Real(0.0))
                << "direction=" << dir << " component=" << component;
        }
    }
}

bool accumulate_ledger_components (CompletedLedgerParallelFixture& fixture,
                                   std::string& diagnostic);

bool restart_heun_stage_zero (CompletedLedgerParallelFixture& fixture,
                              const double dt, std::string& diagnostic)
{
    erf_auxiliary::AuxiliaryStageRecipe recipe;
    const auto method = erf_auxiliary::HostIntegrator::AnelasticHeun;
    if (!erf_auxiliary::MakeAuxiliaryStageRecipe(method, 0, dt, recipe,
                                                diagnostic) ||
        !fixture.ledger.begin_stage(method, 0, 0.0, recipe, diagnostic) ||
        !accumulate_ledger_components(fixture, diagnostic)) {
        return false;
    }
    return fixture.ledger.finish_stage(diagnostic);
}

bool
accumulate_ledger_components (CompletedLedgerParallelFixture& fixture,
                              std::string& diagnostic)
{
    bool complete = true;
    for (int component = 0; component < 2; ++component) {
        const bool accumulated = fixture.ledger.accumulate_stage_component(
            fixture.rate, component, component, diagnostic);
        complete = accumulated && complete;
    }
    return complete;
}

// Keep the device lambda out of the private TestBody generated by TEST for NVCC.
void
set_max_finite_face_value (MultiFab& face_data,
                           const int box_index,
                           const Box& face_cell)
{
    const Real value = std::numeric_limits<Real>::max();
    for (amrex::MFIter mfi(face_data, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        if (mfi.index() != box_index)
            continue;
        const auto face = face_data.array(mfi);
        amrex::ParallelFor(
            face_cell, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                face(i, j, k, 0) = value;
            });
    }
}

TEST(SBMTransportParallel, DecompositionInvariant)
{
    const auto one_box = run_decomposition(16, erf_sbm::MomentMode::OneMoment);
    const auto many_boxes = run_decomposition(4, erf_sbm::MomentMode::OneMoment);
    ASSERT_EQ(one_box.size(), many_boxes.size());
    for (std::size_t i = 0; i < one_box.size(); ++i) {
        EXPECT_NEAR(one_box[i], many_boxes[i],
                    Real(512.0) * std::numeric_limits<Real>::epsilon() *
                        std::max(std::abs(one_box[i]),
                                 std::numeric_limits<Real>::min()));
    }
}

TEST(SBMTransportParallel, TwoMomentEndpointGroupLimiterDecompositionInvariant)
{
    const auto one_box = run_decomposition(16, erf_sbm::MomentMode::TwoMoment);
    const auto many_boxes = run_decomposition(4, erf_sbm::MomentMode::TwoMoment);
    ASSERT_EQ(one_box.size(), many_boxes.size());
    for (std::size_t i = 0; i < one_box.size(); ++i) {
        EXPECT_NEAR(one_box[i], many_boxes[i],
                    Real(512.0) * std::numeric_limits<Real>::epsilon() *
                        std::max(std::abs(one_box[i]),
                                 std::numeric_limits<Real>::min()));
    }
}

TEST(SBMTransportParallel, RankLocalIncompleteLedgerStageFailsCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "collective completed-ledger qualification requires at "
                        "least two MPI ranks";
    }

    CompletedLedgerParallelFixture fixture;
    std::string diagnostic;
    constexpr double dt = 2.0;
    erf_auxiliary::AuxiliaryStageRecipe recipe;
    const auto method = erf_auxiliary::HostIntegrator::AnelasticHeun;

    const bool stage0_recipe_ok = erf_auxiliary::MakeAuxiliaryStageRecipe(
        method, 0, dt, recipe, diagnostic);
    const bool stage0_begun =
        stage0_recipe_ok &&
        fixture.ledger.begin_stage(method, 0, 0.0, recipe, diagnostic);
    const bool stage0_components_ok =
        accumulate_ledger_components(fixture, diagnostic);
    const bool stage0_finished = fixture.ledger.finish_stage(diagnostic);
    EXPECT_TRUE(stage0_recipe_ok) << diagnostic;
    EXPECT_TRUE(stage0_begun) << diagnostic;
    EXPECT_TRUE(stage0_components_ok) << diagnostic;
    EXPECT_TRUE(stage0_finished) << diagnostic;
    EXPECT_GT(fixture.ledger.integrated_flux().dir(0).norm0(0), Real(0.0));

    diagnostic.clear();
    const bool final_recipe_ok = erf_auxiliary::MakeAuxiliaryStageRecipe(
        method, 1, dt, recipe, diagnostic);
    const bool final_stage_begun =
        final_recipe_ok &&
        fixture.ledger.begin_stage(method, 1, 0.0, recipe, diagnostic);
    const bool first_component_ok = fixture.ledger.accumulate_stage_component(
        fixture.rate, 0, 0, diagnostic);
    const bool second_component_ok =
        amrex::ParallelDescriptor::MyProc() == 1
            ? true
            : fixture.ledger.accumulate_stage_component(fixture.rate, 1, 1,
                                                        diagnostic);
    const bool accepted = fixture.ledger.finish_stage(diagnostic);

    EXPECT_TRUE(final_recipe_ok) << diagnostic;
    EXPECT_TRUE(final_stage_begun) << diagnostic;
    EXPECT_TRUE(first_component_ok) << diagnostic;
    EXPECT_TRUE(second_component_ok) << diagnostic;
    EXPECT_FALSE(accepted);
    if (amrex::ParallelDescriptor::MyProc() == 1) {
        EXPECT_NE(diagnostic.find("omitted destination component 1"),
                  std::string::npos)
            << diagnostic;
    } else {
        EXPECT_NE(diagnostic.find("another MPI rank"), std::string::npos)
            << diagnostic;
    }
    expect_discarded_ledger(fixture.ledger);

    int minimum_next_stage = fixture.ledger.next_stage();
    int maximum_next_stage = minimum_next_stage;
    amrex::ParallelDescriptor::ReduceIntMin(minimum_next_stage);
    amrex::ParallelDescriptor::ReduceIntMax(maximum_next_stage);
    EXPECT_EQ(minimum_next_stage, 0);
    EXPECT_EQ(maximum_next_stage, 0);
    EXPECT_TRUE(restart_heun_stage_zero(fixture, dt, diagnostic)) << diagnostic;
    EXPECT_TRUE(fixture.ledger.step_active());
    EXPECT_FALSE(fixture.ledger.step_complete());
    EXPECT_EQ(fixture.ledger.next_stage(), 1);
    EXPECT_GT(fixture.ledger.integrated_flux().dir(0).norm0(0), Real(0.0));
}

TEST(SBMTransportParallel, RankLocalClosedLedgerStageFailsCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "collective completed-ledger qualification requires at "
                        "least two MPI ranks";
    }

    CompletedLedgerParallelFixture fixture;
    std::string diagnostic;
    constexpr double dt = 2.0;
    const auto method = erf_auxiliary::HostIntegrator::AnelasticHeun;
    erf_auxiliary::AuxiliaryStageRecipe recipe;

    bool recipe_ok = erf_auxiliary::MakeAuxiliaryStageRecipe(
        method, 0, dt, recipe, diagnostic);
    bool begun = recipe_ok &&
                 fixture.ledger.begin_stage(method, 0, 0.0, recipe, diagnostic);
    bool components_ok = accumulate_ledger_components(fixture, diagnostic);
    const bool stage0_finished = fixture.ledger.finish_stage(diagnostic);
    EXPECT_TRUE(recipe_ok) << diagnostic;
    EXPECT_TRUE(begun) << diagnostic;
    EXPECT_TRUE(components_ok) << diagnostic;
    EXPECT_TRUE(stage0_finished) << diagnostic;
    EXPECT_GT(fixture.ledger.integrated_flux().dir(0).norm0(0), Real(0.0));

    diagnostic.clear();
    recipe_ok = erf_auxiliary::MakeAuxiliaryStageRecipe(method, 1, dt, recipe,
                                                        diagnostic);
    begun = recipe_ok && amrex::ParallelDescriptor::MyProc() != 1 &&
            fixture.ledger.begin_stage(method, 1, 0.0, recipe, diagnostic);
    if (amrex::ParallelDescriptor::MyProc() != 1) {
        components_ok = accumulate_ledger_components(fixture, diagnostic);
    } else {
        components_ok = false;
    }
    const bool accepted = fixture.ledger.finish_stage(diagnostic);

    EXPECT_TRUE(recipe_ok) << diagnostic;
    EXPECT_EQ(begun, amrex::ParallelDescriptor::MyProc() != 1) << diagnostic;
    EXPECT_EQ(components_ok, amrex::ParallelDescriptor::MyProc() != 1)
        << diagnostic;
    EXPECT_FALSE(accepted);
    if (amrex::ParallelDescriptor::MyProc() == 1) {
        EXPECT_NE(diagnostic.find("open host stage"), std::string::npos)
            << diagnostic;
    } else {
        EXPECT_NE(diagnostic.find("another MPI rank"), std::string::npos)
            << diagnostic;
    }
    expect_discarded_ledger(fixture.ledger);
    int minimum_next_stage = fixture.ledger.next_stage();
    int maximum_next_stage = minimum_next_stage;
    amrex::ParallelDescriptor::ReduceIntMin(minimum_next_stage);
    amrex::ParallelDescriptor::ReduceIntMax(maximum_next_stage);
    EXPECT_EQ(minimum_next_stage, 0);
    EXPECT_EQ(maximum_next_stage, 0);
    EXPECT_TRUE(restart_heun_stage_zero(fixture, dt, diagnostic)) << diagnostic;
    EXPECT_TRUE(fixture.ledger.step_active());
    EXPECT_EQ(fixture.ledger.next_stage(), 1);
}

TEST(SBMTransportParallel, DifferentLedgerStageIdentitiesFailCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "collective completed-ledger qualification requires at "
                        "least two MPI ranks";
    }

    CompletedLedgerParallelFixture fixture;
    std::string diagnostic;
    constexpr double dt = 2.0;
    const auto method = amrex::ParallelDescriptor::MyProc() == 1
                            ? erf_auxiliary::HostIntegrator::AnelasticHeun
                            : erf_auxiliary::HostIntegrator::CompressibleRK3;
    erf_auxiliary::AuxiliaryStageRecipe recipe;
    const bool recipe_ok = erf_auxiliary::MakeAuxiliaryStageRecipe(
        method, 0, dt, recipe, diagnostic);
    const bool begun = recipe_ok && fixture.ledger.begin_stage(
                                        method, 0, 0.0, recipe, diagnostic);
    const bool components_ok =
        accumulate_ledger_components(fixture, diagnostic);
    EXPECT_GT(fixture.ledger.integrated_flux().dir(0).norm0(0), Real(0.0));
    const bool accepted = fixture.ledger.finish_stage(diagnostic);

    EXPECT_TRUE(recipe_ok) << diagnostic;
    EXPECT_TRUE(begun) << diagnostic;
    EXPECT_TRUE(components_ok) << diagnostic;
    EXPECT_FALSE(accepted);
    EXPECT_NE(diagnostic.find("stage identity differs across MPI ranks"),
              std::string::npos)
        << diagnostic;
    expect_discarded_ledger(fixture.ledger);
    int minimum_next_stage = fixture.ledger.next_stage();
    int maximum_next_stage = minimum_next_stage;
    amrex::ParallelDescriptor::ReduceIntMin(minimum_next_stage);
    amrex::ParallelDescriptor::ReduceIntMax(maximum_next_stage);
    EXPECT_EQ(minimum_next_stage, 0);
    EXPECT_EQ(maximum_next_stage, 0);
    EXPECT_TRUE(restart_heun_stage_zero(fixture, dt, diagnostic)) << diagnostic;
    EXPECT_TRUE(fixture.ledger.step_active());
    EXPECT_EQ(fixture.ledger.next_stage(), 1);
}

TEST(SBMTransportParallel,
     AcceptStageRankLocalUndefinedRateRejectsCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "collective accept_stage qualification requires at least "
                        "two MPI ranks";
    }

    CompletedLedgerParallelFixture fixture;
    erf_auxiliary::MappedFaceFluxRate undefined_rate;
    erf_auxiliary::AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    constexpr double dt = 2.0;
    const auto method = erf_auxiliary::HostIntegrator::AnelasticHeun;
    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(
        method, 0, dt, recipe, diagnostic)) << diagnostic;
    const auto& local_rate = amrex::ParallelDescriptor::MyProc() == 1
                                 ? undefined_rate
                                 : fixture.rate;
    const bool accepted = fixture.ledger.accept_stage(
        method, 0, 0.0, recipe, local_rate, diagnostic);
    EXPECT_FALSE(accepted);
    if (amrex::ParallelDescriptor::MyProc() == 1) {
        EXPECT_NE(diagnostic.find("not defined compatibly"), std::string::npos)
            << diagnostic;
    } else {
        EXPECT_NE(diagnostic.find("another MPI rank"), std::string::npos)
            << diagnostic;
    }
    expect_discarded_ledger(fixture.ledger);
    EXPECT_TRUE(restart_heun_stage_zero(fixture, dt, diagnostic)) << diagnostic;
}

TEST(SBMTransportParallel,
     AcceptStageRankLocalBeginErrorDiscardsStepCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "collective accept_stage qualification requires at least "
                        "two MPI ranks";
    }

    CompletedLedgerParallelFixture fixture;
    erf_auxiliary::AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    constexpr double dt = 2.0;
    const auto method = erf_auxiliary::HostIntegrator::AnelasticHeun;
    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(
        method, 0, dt, recipe, diagnostic)) << diagnostic;
    ASSERT_TRUE(fixture.ledger.accept_stage(method, 0, 0.0, recipe,
                                            fixture.rate, diagnostic))
        << diagnostic;

    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(
        method, 1, dt, recipe, diagnostic)) << diagnostic;
    const double step_old_time =
        amrex::ParallelDescriptor::MyProc() == 1 ? 1.0 : 0.0;
    const bool accepted = fixture.ledger.accept_stage(
        method, 1, step_old_time, recipe, fixture.rate, diagnostic);
    EXPECT_FALSE(accepted);
    if (amrex::ParallelDescriptor::MyProc() == 1) {
        EXPECT_NE(diagnostic.find("step-old time changed"), std::string::npos)
            << diagnostic;
    } else {
        EXPECT_NE(diagnostic.find("another MPI rank"), std::string::npos)
            << diagnostic;
    }
    expect_discarded_ledger(fixture.ledger);
    EXPECT_TRUE(restart_heun_stage_zero(fixture, dt, diagnostic)) << diagnostic;
}

TEST(SBMTransportParallel,
     CompletedRecordIsDiscardedWhenOnlySomeRanksBeginNextStep)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "mixed-rank next-step qualification requires at least "
                        "two MPI ranks";
    }

    CompletedLedgerParallelFixture fixture;
    erf_auxiliary::AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    constexpr double dt = 2.0;
    const auto method = erf_auxiliary::HostIntegrator::AnelasticHeun;
    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(method, 0, dt, recipe,
                                                        diagnostic))
        << diagnostic;
    ASSERT_TRUE(fixture.ledger.accept_stage(method, 0, 0.0, recipe,
                                            fixture.rate, diagnostic))
        << diagnostic;
    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(method, 1, dt, recipe,
                                                        diagnostic))
        << diagnostic;
    ASSERT_TRUE(fixture.ledger.accept_stage(method, 1, 0.0, recipe,
                                            fixture.rate, diagnostic))
        << diagnostic;
    ASSERT_TRUE(fixture.ledger.step_complete());
    ASSERT_GT(fixture.ledger.integrated_flux().dir(0).norm0(0), Real(0.0));

    // Rank 0 and any ranks above 1 start the new stage zero. Rank 1 rejects
    // locally before begin_stage; every rank still enters finish_stage.
    erf_auxiliary::MappedFaceFluxRate undefined_rate;
    const auto& local_rate = amrex::ParallelDescriptor::MyProc() == 1
                                 ? undefined_rate
                                 : fixture.rate;
    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(method, 0, dt, recipe,
                                                        diagnostic))
        << diagnostic;
    const bool accepted = fixture.ledger.accept_stage(method, 0, 0.0, recipe,
                                                      local_rate, diagnostic);
    EXPECT_FALSE(accepted);
    if (amrex::ParallelDescriptor::MyProc() == 1) {
        EXPECT_NE(diagnostic.find("not defined compatibly"), std::string::npos)
            << diagnostic;
    } else {
        EXPECT_NE(diagnostic.find("another MPI rank"), std::string::npos)
            << diagnostic;
    }
    expect_discarded_ledger(fixture.ledger);
    EXPECT_TRUE(restart_heun_stage_zero(fixture, dt, diagnostic)) << diagnostic;
    EXPECT_TRUE(fixture.ledger.step_active());
    EXPECT_FALSE(fixture.ledger.step_complete());
    EXPECT_EQ(fixture.ledger.next_stage(), 1);
}

TEST(SBMTransportParallel, NonfiniteCompletedLedgerRejectedCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "collective completed-ledger qualification requires at "
                        "least two MPI ranks";
    }

    Box domain(IntVect(0), IntVect(7));
    BoxArray boxes(domain);
    boxes.maxSize(2);
    amrex::Vector<int> owners(static_cast<std::size_t>(boxes.size()));
    for (int box = 0; box < boxes.size(); ++box) {
        owners[box] = box % amrex::ParallelDescriptor::NProcs();
    }
    DistributionMapping mapping(std::move(owners));

    erf_auxiliary::MappedFaceFluxRate rate;
    rate.define(boxes, mapping, 1, 0);
    erf_auxiliary::CompletedStepFluxLedger ledger;
    ledger.define(boxes, mapping);
    erf_auxiliary::AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    constexpr double dt = 4.0;

    rate.setVal(Real(0.0));
    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(
        erf_auxiliary::HostIntegrator::AnelasticHeun, 0, dt, recipe,
        diagnostic))
        << diagnostic;
    ASSERT_TRUE(
        ledger.accept_stage(erf_auxiliary::HostIntegrator::AnelasticHeun, 0,
                            0.0, recipe, rate, diagnostic))
        << diagnostic;

    int overflow_box = -1;
    for (int box = 0; box < boxes.size(); ++box) {
        if (box % amrex::ParallelDescriptor::NProcs() == 1) {
            overflow_box = box;
            break;
        }
    }
    ASSERT_GE(overflow_box, 0);
    const Box overflow_cells = boxes[overflow_box];
    const IntVect overflow_face = overflow_cells.smallEnd();
    const Box face_cell(overflow_face, overflow_face);
    set_max_finite_face_value(rate.dir(0), overflow_box, face_cell);
    ASSERT_TRUE(rate.dir(0).is_finite(0, 1, 0));
    ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(
        erf_auxiliary::HostIntegrator::AnelasticHeun, 1, dt, recipe,
        diagnostic))
        << diagnostic;
    const bool accepted =
        ledger.accept_stage(erf_auxiliary::HostIntegrator::AnelasticHeun, 1,
                            0.0, recipe, rate, diagnostic);

    int minimum_accepted = accepted ? 1 : 0;
    int maximum_accepted = minimum_accepted;
    amrex::ParallelDescriptor::ReduceIntMin(minimum_accepted);
    amrex::ParallelDescriptor::ReduceIntMax(maximum_accepted);
    EXPECT_EQ(minimum_accepted, 0) << diagnostic;
    EXPECT_EQ(maximum_accepted, 0) << diagnostic;
    EXPECT_NE(diagnostic.find("nonfinite"), std::string::npos) << diagnostic;
    EXPECT_FALSE(ledger.step_complete());
    EXPECT_FALSE(ledger.step_active());
    EXPECT_EQ(ledger.next_stage(), 0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        EXPECT_TRUE(ledger.integrated_flux().dir(dir).is_finite(0, 1, 0, true));
        EXPECT_EQ(ledger.integrated_flux().dir(dir).norm0(0), Real(0.0));
    }
}

} // namespace
