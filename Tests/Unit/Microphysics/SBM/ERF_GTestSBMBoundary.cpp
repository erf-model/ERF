#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include "ERF_IndexDefines.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMPhysicalBoundary.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;
using amrex::MultiFab;
using amrex::Real;
using erf_sbm::AdvectionFaceKind;
using erf_auxiliary::HostIntegrator;

constexpr int nx = 8;
constexpr int ny = 4;
constexpr int nz = 4;
constexpr int xlo = 0;
constexpr int xhi = 3;

erf_sbm::SBMLayout
make_layout (const erf_sbm::MomentMode mode, const bool with_property = false)
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
    if (with_property) {
        erf_sbm::AttachedPropertyDescriptor property;
        property.name = "attached_mass";
        property.semantic_id = "attached_mass";
        property.units = "kg m^-3";
        property.carrier_population = 0;
        property.kind = erf_sbm::PropertyKind::NumberCarried;
        property.support = erf_sbm::SupportRequirement::None;
        property.transported = true;
        property.support_min = Real(0.0);
        property.support_max = Real(0.5);
        spec.attached_properties.push_back(std::move(property));
    }
    return erf_sbm::SBMLayout(std::move(spec));
}

// Dry-air-relative values of one admissible spectrum: 1M (m0, m1) or, for 2M,
// (m0, m1, n0, n1) with bin means 0.3 and 0.75, plus attached material.
std::vector<Real>
spectrum_values (const erf_sbm::SBMLayout& layout, const Real scale)
{
    std::vector<Real> z(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    const auto& population = layout.populations().front();
    const Real n0 = Real(0.01) * scale;
    const Real n1 = Real(0.02) * scale;
    z[static_cast<std::size_t>(population.mass_offset)] = Real(0.3) * n0;
    z[static_cast<std::size_t>(population.mass_offset + 1)] = Real(0.75) * n1;
    if (population.number_offset >= 0) {
        z[static_cast<std::size_t>(population.number_offset)] = n0;
        z[static_cast<std::size_t>(population.number_offset + 1)] = n1;
    }
    if (layout.attached_properties().size() == 1) {
        const int property = layout.property_offset(0);
        z[static_cast<std::size_t>(property)] = Real(0.25) * n0;
        z[static_cast<std::size_t>(property + 1)] = Real(0.25) * n1;
    }
    return z;
}

struct BoundaryCase
{
    erf_sbm::MomentMode mode{erf_sbm::MomentMode::OneMoment};
    bool with_property{false};
    AdvectionFaceKind low{AdvectionFaceKind::PrescribedInflow};
    AdvectionFaceKind high{AdvectionFaceKind::AdvectiveOutflow};
    Real inflow_scale{Real(0.5)};
    Real carrier{Real(0.3)};
    Real low_face_carrier{Real(0.3)};
    Real high_face_carrier{Real(0.3)};
    // Interior spectrum scale varies along x unless uniform.
    bool uniform{false};
    HostIntegrator method{HostIntegrator::CompressibleRK3};
    int stages{-1};
};

struct BoundaryResult
{
    std::vector<Real> initial_inventory;
    std::vector<Real> final_inventory;
    erf_sbm::SBMBoundaryLedger ledger;
    std::unique_ptr<MultiFab> initial;
    std::unique_ptr<MultiFab> final_state;
};

Geometry
x_bounded_geometry ()
{
    const Box domain(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {0, 1, 1};
    return Geometry(domain, &physical, amrex::CoordSys::cartesian, periodic);
}

void
fill_spectrum (MultiFab& spectrum, const std::vector<Real>& base, const bool uniform)
{
    amrex::Gpu::DeviceVector<Real> values(base.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, base.begin(), base.end(), values.begin());
    const Real* z = values.data();
    const int ncomp = spectrum.nComp();
    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const auto state = spectrum.array(mfi);
        amrex::ParallelFor(mfi.validbox(), ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                const Real phase = (static_cast<Real>(i) + Real(0.5)) / Real(nx);
                const Real factor = uniform ? Real(1.0)
                    : Real(1.0) + Real(0.4) * amrex::Math::sinpi(Real(2.0) * phase);
                state(i, j, k, n) = factor * z[n];
            });
    }
    amrex::Gpu::streamSynchronize();
}

void
fill_x_carrier (MultiFab& avg_xmom, const BoundaryCase& c)
{
    const Real interior = c.carrier;
    const Real low = c.low_face_carrier;
    const Real high = c.high_face_carrier;
    for (amrex::MFIter mfi(avg_xmom); mfi.isValid(); ++mfi) {
        const auto flux = avg_xmom.array(mfi);
        amrex::ParallelFor(mfi.validbox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                flux(i, j, k, 0) = i == 0 ? low : (i == nx ? high : interior);
            });
    }
}

std::vector<Real>
inventory (const MultiFab& state, const Geometry& geom)
{
    const Real volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
    std::vector<Real> totals(static_cast<std::size_t>(state.nComp()));
    for (int n = 0; n < state.nComp(); ++n) {
        totals[static_cast<std::size_t>(n)] = volume * state.sum(n);
    }
    return totals;
}

BoundaryResult
run_boundary_step (const BoundaryCase& c, const std::function<void(erf_sbm::SBMTransport&)>& before_step = {})
{
    const Geometry geom = x_bounded_geometry();
    const BoxArray ba(geom.Domain());
    const DistributionMapping dm(ba);
    const auto layout = make_layout(c.mode, c.with_property);

    erf_sbm::SBMStateManager state_manager(layout, 1);
    state_manager.define(0, ba, dm, 0.0);
    erf_sbm::SBMTransport transport(state_manager.layout(), 1);
    transport.define(0, ba, dm);
    erf_sbm::SBMPhysicalBoundaryConfig config;
    config.kind[xlo] = c.low;
    config.kind[xhi] = c.high;
    if (c.low == AdvectionFaceKind::PrescribedInflow) {
        config.inflow_z[xlo] = spectrum_values(layout, c.inflow_scale);
    }
    transport.set_physical_boundaries(config);

    MultiFab detj(ba, dm, 1, 0);
    amrex::BoxList map_boxes = ba.boxList();
    for (auto& box : map_boxes) box.setRange(2, 0);
    const BoxArray map_ba(std::move(map_boxes));
    MultiFab mx(map_ba, dm, 1, 0);
    MultiFab my(map_ba, dm, 1, 0);
    detj.setVal(Real(1.0));
    mx.setVal(Real(1.0));
    my.setVal(Real(1.0));
    std::string diagnostic;
    EXPECT_TRUE(transport.rebuild_static_measure(0, detj, mx, my, diagnostic)) << diagnostic;

    MultiFab conserved_anchor(ba, dm, 3, 0);
    MultiFab conserved_input(ba, dm, 3, 0);
    MultiFab conserved_target(ba, dm, 3, 0);
    for (MultiFab* conserved : {&conserved_anchor, &conserved_input, &conserved_target}) {
        conserved->setVal(Real(0.0));
        conserved->setVal(Real(1.0), Rho_comp, 1);
    }

    auto& spectrum = state_manager.new_state_for_initialization(0);
    fill_spectrum(spectrum, spectrum_values(layout, Real(1.0)), c.uniform);

    BoundaryResult result;
    result.initial = std::make_unique<MultiFab>(ba, dm, layout.ncomp(), 0);
    MultiFab::Copy(*result.initial, spectrum, 0, 0, layout.ncomp(), 0);
    result.initial_inventory = inventory(spectrum, geom);

    MultiFab avg_xmom(amrex::convert(ba, IntVect::TheDimensionVector(0)), dm, 1, 0);
    MultiFab avg_ymom(amrex::convert(ba, IntVect::TheDimensionVector(1)), dm, 1, 0);
    MultiFab avg_zmom(amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1, 0);
    fill_x_carrier(avg_xmom, c);
    avg_ymom.setVal(Real(0.05));
    avg_zmom.setVal(Real(-0.04));

    if (before_step) before_step(transport);
    const double dt = 0.1;
    EXPECT_TRUE(state_manager.begin_step(0, 0.0, diagnostic)) << diagnostic;
    struct Stage { double input; double target; double interval; };
    std::vector<Stage> stages;
    if (c.method == HostIntegrator::CompressibleRK3) {
        stages = {{0.0, dt / 3.0, dt / 3.0}, {dt / 3.0, dt / 2.0, dt / 2.0}, {dt / 2.0, dt, dt}};
    } else if (c.method == HostIntegrator::AnelasticHeun) {
        stages = {{0.0, dt, dt}, {dt, dt, dt}};
    } else {
        stages = {{0.0, dt / 2.0, dt / 2.0}, {dt / 2.0, dt, dt}};
    }
    const int stage_count = c.stages < 0 ? static_cast<int>(stages.size()) : c.stages;
    for (int stage = 0; stage < stage_count; ++stage) {
        const auto& s = stages[static_cast<std::size_t>(stage)];
        transport.advance_stage_from_host(
            0, c.method, stage, 0.0, s.input, s.target, s.interval, state_manager,
            conserved_anchor, conserved_input, conserved_target,
            avg_xmom, avg_ymom, avg_zmom, geom, 1, 2);
    }
    const auto& final_state = state_manager.new_state(0);
    result.final_state = std::make_unique<MultiFab>(ba, dm, layout.ncomp(), 0);
    MultiFab::Copy(*result.final_state, final_state, 0, 0, layout.ncomp(), 0);
    result.final_inventory = inventory(final_state, geom);
    result.ledger = transport.boundary_ledger();
    return result;
}

// Forked death tests hang under MPI on Linux CI, so the expected
// amrex::Abort is raised as an exception for the duration of the check.
template <typename Action>
void expect_amrex_abort (Action&& action, const std::string& expected)
{
    const bool saved = amrex::system::throw_exception;
    amrex::system::throw_exception = true;
    std::string message;
    try {
        action();
    } catch (const std::runtime_error& error) {
        message = error.what();
    }
    amrex::system::throw_exception = saved;
    EXPECT_NE(message.find(expected), std::string::npos)
        << "expected amrex::Abort containing \"" << expected << "\", got \""
        << message << '"';
}

void
expect_budget_closes (const BoundaryResult& result)
{
    for (int n = 0; n < result.ledger.ncomp; ++n) {
        const auto c = static_cast<std::size_t>(n);
        const double change = static_cast<double>(result.final_inventory[c]) -
                              static_cast<double>(result.initial_inventory[c]);
        const double scale = std::abs(result.initial_inventory[c]) +
                             std::abs(result.final_inventory[c]) +
                             result.ledger.total_gross(n);
        EXPECT_LE(std::abs(change + result.ledger.total_outward(n)),
                  1024.0 * std::numeric_limits<Real>::epsilon() * scale)
            << "component " << n;
    }
}

} // namespace

TEST(SBMBoundary, MapsEveryErfBoundaryType)
{
    AdvectionFaceKind kind{};
    std::string reason;
    const auto maps_to = [&] (ERF_BC bc, AdvectionFaceKind expected) {
        EXPECT_TRUE(erf_sbm::map_erf_boundary(bc, kind, reason)) << reason;
        EXPECT_EQ(kind, expected);
    };
    maps_to(ERF_BC::periodic, AdvectionFaceKind::Periodic);
    maps_to(ERF_BC::slip_wall, AdvectionFaceKind::ImpermeableWall);
    maps_to(ERF_BC::no_slip_wall, AdvectionFaceKind::ImpermeableWall);
    maps_to(ERF_BC::outflow, AdvectionFaceKind::AdvectiveOutflow);
    maps_to(ERF_BC::ho_outflow, AdvectionFaceKind::AdvectiveOutflow);
    maps_to(ERF_BC::inflow, AdvectionFaceKind::PrescribedInflow);
    maps_to(ERF_BC::inflow_outflow, AdvectionFaceKind::PrescribedInflow);
    for (const ERF_BC rejected : {ERF_BC::symmetry, ERF_BC::surface_layer,
                                  ERF_BC::open, ERF_BC::undefined}) {
        EXPECT_FALSE(erf_sbm::map_erf_boundary(rejected, kind, reason));
        EXPECT_FALSE(reason.empty());
    }
}

TEST(SBMBoundary, InflowStateMustBeCompleteCanonicalAndAdmissible)
{
    std::string reason;
    const auto one_moment = make_layout(erf_sbm::MomentMode::OneMoment);
    EXPECT_TRUE(erf_sbm::validate_inflow_state(one_moment, {Real(1.0e-3), Real(2.0e-3)}, reason))
        << reason;
    EXPECT_TRUE(erf_sbm::validate_inflow_state(one_moment, {Real(0.0), Real(0.0)}, reason));
    EXPECT_FALSE(erf_sbm::validate_inflow_state(one_moment, {Real(1.0e-3)}, reason));
    EXPECT_NE(reason.find("1 values"), std::string::npos) << reason;
    EXPECT_FALSE(erf_sbm::validate_inflow_state(one_moment, {Real(1.0e-3), Real(-1.0e-9)}, reason));
    EXPECT_FALSE(erf_sbm::validate_inflow_state(
        one_moment, {Real(1.0e-3), std::numeric_limits<Real>::quiet_NaN()}, reason));

    const auto two_moment = make_layout(erf_sbm::MomentMode::TwoMoment);
    const auto complete = spectrum_values(two_moment, Real(1.0));
    EXPECT_TRUE(erf_sbm::validate_inflow_state(two_moment, complete, reason)) << reason;
    // Mass without number, and a bin mean outside its edges, are incomplete groups.
    auto orphan_mass = complete;
    orphan_mass[static_cast<std::size_t>(two_moment.populations().front().number_offset)] = Real(0.0);
    EXPECT_FALSE(erf_sbm::validate_inflow_state(two_moment, orphan_mass, reason));
    auto outside_bin = complete;
    outside_bin[0] = Real(0.9) * outside_bin[static_cast<std::size_t>(
        two_moment.populations().front().number_offset)];
    EXPECT_FALSE(erf_sbm::validate_inflow_state(two_moment, outside_bin, reason));

    const auto attached = make_layout(erf_sbm::MomentMode::TwoMoment, true);
    const auto with_material = spectrum_values(attached, Real(1.0));
    EXPECT_TRUE(erf_sbm::validate_inflow_state(attached, with_material, reason)) << reason;
    // Attached material needs its carrier bin.
    auto empty_carrier = with_material;
    const int number0 = attached.populations().front().number_offset;
    empty_carrier[0] = Real(0.0);
    empty_carrier[static_cast<std::size_t>(number0)] = Real(0.0);
    EXPECT_FALSE(erf_sbm::validate_inflow_state(attached, empty_carrier, reason));
}

TEST(SBMBoundary, RestartIdentityKeepsThePeriodicSchemaAndFingerprintsInflow)
{
    const auto layout = make_layout(erf_sbm::MomentMode::OneMoment);
    erf_sbm::SBMPhysicalBoundaryConfig periodic;
    EXPECT_EQ(periodic.identity(), "boundary=triply-periodic-v1");
    EXPECT_EQ(erf_sbm::restart_schema(layout, periodic.identity()),
              erf_sbm::restart_schema(layout));

    erf_sbm::SBMPhysicalBoundaryConfig bounded;
    bounded.kind[xlo] = AdvectionFaceKind::PrescribedInflow;
    bounded.kind[xhi] = AdvectionFaceKind::AdvectiveOutflow;
    bounded.kind[2] = AdvectionFaceKind::ImpermeableWall;
    bounded.kind[5] = AdvectionFaceKind::ImpermeableWall;
    bounded.inflow_z[xlo] = {Real(1.0e-3), Real(2.0e-3)};
    const std::string identity = bounded.identity();
    EXPECT_EQ(identity.rfind("boundary=physical-faces-v1;xlo=inflow[", 0), 0U) << identity;
    EXPECT_NE(identity.find(";xhi=outflow;ylo=periodic;yhi=periodic;zlo=wall;zhi=wall"),
              std::string::npos) << identity;
    const std::string schema = erf_sbm::restart_schema(layout, identity);
    EXPECT_TRUE(erf_sbm::restart_schema_matches(layout, schema, identity));
    EXPECT_FALSE(erf_sbm::restart_schema_matches(layout, schema));

    auto changed = bounded;
    changed.inflow_z[xlo][1] = std::nextafter(changed.inflow_z[xlo][1], Real(1.0));
    EXPECT_NE(changed.identity(), identity);
    EXPECT_FALSE(erf_sbm::restart_schema_matches(layout, schema, changed.identity()));
}

TEST(SBMBoundary, LedgerSerializationRoundTripsExactly)
{
    erf_sbm::SBMBoundaryLedger ledger;
    ledger.reset(3);
    ledger.face_stage_terms = 123456789012345ULL;
    for (std::size_t slot = 0; slot < ledger.outward.size(); ++slot) {
        ledger.outward[slot] = (slot % 2 == 0 ? -1.0 : 1.0) / (3.0 + static_cast<double>(slot));
        ledger.gross[slot] = 1.0e-300 * static_cast<double>(slot + 1);
    }
    const std::string text = ledger.serialize();
    erf_sbm::SBMBoundaryLedger parsed;
    std::string reason;
    ASSERT_TRUE(erf_sbm::SBMBoundaryLedger::parse(text, 3, parsed, reason)) << reason;
    EXPECT_EQ(parsed.face_stage_terms, ledger.face_stage_terms);
    EXPECT_EQ(parsed.outward, ledger.outward);
    EXPECT_EQ(parsed.gross, ledger.gross);
    EXPECT_EQ(parsed.serialize(), text);

    EXPECT_FALSE(erf_sbm::SBMBoundaryLedger::parse(text, 4, parsed, reason));
    std::string corrupted = text;
    corrupted.replace(corrupted.find("xhi 1 ") + 6, 1, "nan ");
    EXPECT_FALSE(erf_sbm::SBMBoundaryLedger::parse(corrupted, 3, parsed, reason));
    EXPECT_FALSE(erf_sbm::SBMBoundaryLedger::parse(text + "extra\n", 3, parsed, reason));
}

TEST(SBMBoundary, WallsConserveTheSpectrumAndRecordNoTransfer)
{
    for (const auto mode : {erf_sbm::MomentMode::OneMoment, erf_sbm::MomentMode::TwoMoment}) {
        BoundaryCase c;
        c.mode = mode;
        c.low = AdvectionFaceKind::ImpermeableWall;
        c.high = AdvectionFaceKind::ImpermeableWall;
        c.low_face_carrier = Real(0.0);
        c.high_face_carrier = Real(0.0);
        const auto result = run_boundary_step(c);
        expect_budget_closes(result);
        for (std::size_t slot = 0; slot < result.ledger.outward.size(); ++slot) {
            EXPECT_EQ(result.ledger.outward[slot], 0.0);
            EXPECT_EQ(result.ledger.gross[slot], 0.0);
        }
        EXPECT_GT(result.ledger.face_stage_terms, 0U);
    }
}

TEST(SBMBoundary, NonzeroWallCarrierAbortsTheStage)
{
    BoundaryCase c;
    c.low = AdvectionFaceKind::ImpermeableWall;
    c.high = AdvectionFaceKind::ImpermeableWall;
    c.high_face_carrier = Real(0.0);
    expect_amrex_abort([&] { static_cast<void>(run_boundary_step(c)); },
                       "SBM impermeable wall xlo has a nonzero dry-air carrier");
}

TEST(SBMBoundary, ReverseFlowAtOutflowAbortsTheStage)
{
    BoundaryCase c;
    c.high_face_carrier = Real(-0.3);
    expect_amrex_abort([&] { static_cast<void>(run_boundary_step(c)); },
                       "SBM advective outflow xhi received an inward dry-air carrier");
}

TEST(SBMBoundary, OutflowLossEqualsTheLedgerForEveryHostIntegrator)
{
    for (const auto method : {HostIntegrator::CompressibleRK3, HostIntegrator::AnelasticHeun,
                              HostIntegrator::AnelasticMidPoint}) {
        for (const auto mode : {erf_sbm::MomentMode::OneMoment, erf_sbm::MomentMode::TwoMoment}) {
            BoundaryCase c;
            c.mode = mode;
            c.method = method;
            c.inflow_scale = Real(0.0);
            const auto result = run_boundary_step(c);
            expect_budget_closes(result);
            for (int n = 0; n < result.ledger.ncomp; ++n) {
                const auto slot = [&] (int face) {
                    return static_cast<std::size_t>(face * result.ledger.ncomp + n);
                };
                EXPECT_EQ(result.ledger.outward[slot(xlo)], 0.0);
                EXPECT_GT(result.ledger.outward[slot(xhi)], 0.0);
                EXPECT_LT(result.final_inventory[static_cast<std::size_t>(n)],
                          result.initial_inventory[static_cast<std::size_t>(n)]);
            }
        }
    }
}

TEST(SBMBoundary, InflowEqualToTheInteriorKeepsAUniformSpectrum)
{
    for (const auto mode : {erf_sbm::MomentMode::OneMoment, erf_sbm::MomentMode::TwoMoment}) {
        BoundaryCase c;
        c.mode = mode;
        c.uniform = true;
        c.inflow_scale = Real(1.0);
        const auto result = run_boundary_step(c);
        expect_budget_closes(result);
        for (int n = 0; n < result.final_state->nComp(); ++n) {
            const Real expected = result.initial->max(n);
            const Real tolerance = Real(8.0) * std::numeric_limits<Real>::epsilon() * expected;
            EXPECT_NEAR(result.final_state->min(n), expected, tolerance) << "component " << n;
            EXPECT_NEAR(result.final_state->max(n), expected, tolerance) << "component " << n;
        }
    }
}

TEST(SBMBoundary, AttachedMaterialEntersWithItsCarrier)
{
    BoundaryCase c;
    c.mode = erf_sbm::MomentMode::TwoMoment;
    c.with_property = true;
    c.inflow_scale = Real(2.0);
    const auto result = run_boundary_step(c);
    expect_budget_closes(result);
    const auto layout = make_layout(c.mode, true);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::authoritative_state_admissible(*result.final_state, layout, 0,
                                                        &diagnostic)) << diagnostic;
    const int property = layout.property_offset(0);
    const auto slot = static_cast<std::size_t>(xlo * layout.ncomp() + property);
    EXPECT_LT(result.ledger.outward[slot], 0.0);
}

// Faces 0, 1, nx-1 and nx read a ghost cell in the WENO-Z3 stencil, so the
// first stage moves the edge cells by the donor rates alone.
TEST(SBMBoundary, DonorBandSetsTheEdgeCellRates)
{
    BoundaryCase c;
    c.stages = 1;
    const auto result = run_boundary_step(c);
    const auto layout = make_layout(c.mode);
    const auto z_in = spectrum_values(layout, c.inflow_scale);
    const Real tau = Real(0.1) / Real(3.0);
    const Real dx_inv = Real(nx);
    MultiFab initial_host(result.initial->boxArray(), result.initial->DistributionMap(),
                          result.initial->nComp(), 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    MultiFab final_host(result.final_state->boxArray(), result.final_state->DistributionMap(),
                        result.final_state->nComp(), 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    MultiFab::Copy(initial_host, *result.initial, 0, 0, initial_host.nComp(), 0);
    MultiFab::Copy(final_host, *result.final_state, 0, 0, final_host.nComp(), 0);
    amrex::Gpu::streamSynchronize();
    for (amrex::MFIter mfi(final_host); mfi.isValid(); ++mfi) {
        const auto u0 = initial_host.const_array(mfi);
        const auto u1 = final_host.const_array(mfi);
        for (int n = 0; n < layout.ncomp(); ++n) {
            const Real low_edge = u0(0, 1, 1, n) -
                tau * dx_inv * c.carrier * (u0(0, 1, 1, n) - z_in[static_cast<std::size_t>(n)]);
            const Real high_edge = u0(nx - 1, 1, 1, n) -
                tau * dx_inv * c.carrier * (u0(nx - 1, 1, 1, n) - u0(nx - 2, 1, 1, n));
            const Real tolerance = Real(16.0) * std::numeric_limits<Real>::epsilon() * u0(0, 1, 1, n);
            EXPECT_NEAR(u1(0, 1, 1, n), low_edge, tolerance) << "component " << n;
            EXPECT_NEAR(u1(nx - 1, 1, 1, n), high_edge, tolerance) << "component " << n;
            // An interior cell keeps the high-order candidate.
            const Real interior_donor = u0(3, 1, 1, n) -
                tau * dx_inv * c.carrier * (u0(3, 1, 1, n) - u0(2, 1, 1, n));
            EXPECT_GT(std::abs(u1(3, 1, 1, n) - interior_donor), Real(1.0e3) * tolerance)
                << "component " << n;
        }
    }
}

TEST(SBMBoundary, ConfigurationMustMatchTheGeometry)
{
    BoundaryCase c;
    expect_amrex_abort(
        [&] {
            static_cast<void>(run_boundary_step(c, [] (erf_sbm::SBMTransport& transport) {
                transport.set_physical_boundaries(erf_sbm::SBMPhysicalBoundaryConfig{});
            }));
        },
        "SBM physical boundary xlo does not match the periodicity of the level geometry");
    erf_sbm::SBMPhysicalBoundaryConfig incomplete;
    incomplete.kind[xlo] = AdvectionFaceKind::PrescribedInflow;
    const auto layout = make_layout(erf_sbm::MomentMode::OneMoment);
    erf_sbm::SBMTransport transport(layout, 1);
    EXPECT_THROW(transport.set_physical_boundaries(incomplete), std::invalid_argument);
}
