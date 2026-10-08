#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>

#include "ERF_SBMAMRTransfer.H"
#include "ERF_SBMRestart.H"

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

struct ComponentOffsets {
    int one_mass;
    int one_property;
    int two_mass;
    int two_number;
    int two_property;
};

erf_sbm::SBMLayout make_transfer_layout ()
{
    erf_sbm::SBMLayoutSpec spec;
    erf_sbm::SpectralPopulationSpec one;
    one.population_id = 0;
    one.semantic_id = "one_moment";
    one.phase = erf_sbm::PopulationPhase::Liquid;
    one.moment_mode = erf_sbm::MomentMode::OneMoment;
    one.grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    one.grid.coordinate_units = "kg";
    one.grid.edges = {Real(0.0), Real(1.0), Real(2.0)};
    one.grid.pivots = {Real(0.5), Real(1.5)};
    spec.populations.push_back(one);

    erf_sbm::SpectralPopulationSpec two;
    two.population_id = 1;
    two.semantic_id = "two_moment";
    two.phase = erf_sbm::PopulationPhase::Aerosol;
    two.moment_mode = erf_sbm::MomentMode::TwoMoment;
    two.grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    two.grid.coordinate_units = "kg";
    two.grid.edges = {Real(0.0), Real(1.0), Real(2.0)};
    two.grid.pivots = {Real(0.5), Real(1.5)};
    spec.populations.push_back(two);
    spec.liquid_projection = {0, 1};

    for (const auto& property_info : std::vector<std::pair<std::string, int>>{
             {"one_moment_extensive", 0}, {"two_moment_carried", 1}}) {
        erf_sbm::AttachedPropertyDescriptor property;
        property.name = property_info.first;
        property.semantic_id = property_info.first + ".semantic";
        property.units = "kg m^-3";
        property.carrier_population = property_info.second;
        property.kind = property_info.second == 0
                            ? erf_sbm::PropertyKind::ExtensiveMass
                            : erf_sbm::PropertyKind::NumberCarried;
        property.support = erf_sbm::SupportRequirement::None;
        property.remap_policy = erf_sbm::PropertyRemapPolicy::CarrierBinConservative;
        property.transported = true;
        property.support_min = Real(0.0);
        spec.attached_properties.push_back(property);
    }
    return erf_sbm::SBMLayout(std::move(spec));
}

ComponentOffsets offsets_for (const erf_sbm::SBMLayout& layout)
{
    return {layout.populations()[0].mass_offset,
            layout.property_offset(0),
            layout.populations()[1].mass_offset,
            layout.populations()[1].number_offset,
            layout.property_offset(1)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real amplitude (const int i, const int j, const int k) noexcept
{
    return Real(1.0) + Real(0.002) * i + Real(0.001) * j + Real(0.0005) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real base_component (const int component, const int i, const int j, const int k,
                     const ComponentOffsets offsets) noexcept
{
    const Real a = amplitude(i, j, k);
    if (component >= offsets.one_mass && component < offsets.one_mass + 2) {
        return Real(0.2) * a * static_cast<Real>(component - offsets.one_mass + 1);
    }
    if (component >= offsets.one_property && component < offsets.one_property + 2) {
        return Real(0.02) * a * static_cast<Real>(component - offsets.one_property + 1);
    }
    if (component >= offsets.two_mass && component < offsets.two_mass + 2) {
        const int bin = component - offsets.two_mass;
        const Real number = Real(0.3) * a * static_cast<Real>(bin + 1);
        const Real mean = bin == 0 ? Real(0.35) : Real(1.35);
        return number * mean;
    }
    if (component >= offsets.two_number && component < offsets.two_number + 2) {
        return Real(0.3) * a * static_cast<Real>(component - offsets.two_number + 1);
    }
    if (component >= offsets.two_property && component < offsets.two_property + 2) {
        return Real(0.15) * a * static_cast<Real>(component - offsets.two_property + 1);
    }
    return Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real density_value (const int i, const int j, const int k) noexcept
{
    return Real(1.0) + Real(0.003) * i + Real(0.002) * j + Real(0.001) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real measure_value (const int i, const int j, const int k) noexcept
{
    return Real(0.7) + Real(0.004) * i + Real(0.003) * j + Real(0.002) * k;
}

Geometry make_geometry (const Box& domain)
{
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    return Geometry(domain, &physical, amrex::CoordSys::cartesian, periodic);
}

DistributionMapping shifted_mapping (const BoxArray& boxes, const int offset)
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    amrex::Vector<int> owners(static_cast<std::size_t>(boxes.size()));
    for (int box = 0; box < boxes.size(); ++box) {
        owners[box] = (box + offset) % nprocs;
    }
    return DistributionMapping(std::move(owners));
}

void fill_carriers (MultiFab& density, MultiFab& measure)
{
    for (amrex::MFIter mfi(density, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto rho = density.array(mfi);
        const auto omega = measure.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            rho(i, j, k, 0) = density_value(i, j, k);
            omega(i, j, k, 0) = measure_value(i, j, k);
        });
    }
}

void fill_spectrum (MultiFab& spectrum, const MultiFab& density,
                    const erf_sbm::SBMLayout& layout)
{
    const ComponentOffsets offsets = offsets_for(layout);
    const int ncomp = layout.ncomp();
    for (amrex::MFIter mfi(spectrum, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto rho = density.const_array(mfi);
        const auto state = spectrum.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real carrier = rho(i, j, k, 0);
            for (int component = 0; component < ncomp; ++component) {
                state(i, j, k, component) = carrier * base_component(
                    component, i, j, k, offsets);
            }
        });
    }
}

erf_sbm::SBMAMRStateView timed_view (const MultiFab& spectrum,
                                    const MultiFab& density,
                                    const MultiFab& measure,
                                    const double time)
{
    return {&spectrum, time, &density, 0, time, &measure, 0, time};
}

bool globally_admissible (const MultiFab& state,
                          const erf_sbm::SBMLayout& layout)
{
    std::string diagnostic;
    const bool local_ok = erf_sbm::authoritative_state_admissible(
        state, layout, 1, &diagnostic);
    int any_bad = local_ok ? 0 : 1;
    amrex::ParallelDescriptor::ReduceIntMax(any_bad);
    return any_bad == 0;
}

void append_fingerprint (const MultiFab& state, std::vector<Real>& result)
{
    for (int component = 0; component < state.nComp(); ++component) {
        MultiFab moments(state.boxArray(), state.DistributionMap(), 2, 0);
        for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box box = mfi.tilebox();
            const auto values = state.const_array(mfi);
            const auto out = moments.array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const Real value = values(i, j, k, component);
                const Real weight = Real(1.0) + Real(0.01) * i +
                                    Real(0.001) * j + Real(0.0001) * k;
                out(i, j, k, 0) = weight * value;
                out(i, j, k, 1) = value * value;
            });
        }
        result.push_back(state.sum(component));
        result.push_back(moments.sum(0));
        result.push_back(moments.sum(1));
    }
}

std::vector<Real> run_transfer_decomposition (const int coarse_max_size,
                                              const int fine_max_size,
                                              const int coarse_owner_shift,
                                              const int fine_owner_shift)
{
    Box coarse_domain(IntVect(0), IntVect(7));
    Box fine_domain = coarse_domain;
    fine_domain.refine(IntVect(2));
    const Geometry coarse_geometry = make_geometry(coarse_domain);
    const Geometry fine_geometry = make_geometry(fine_domain);
    BoxArray coarse_boxes(coarse_domain);
    BoxArray fine_boxes(fine_domain);
    coarse_boxes.maxSize(coarse_max_size);
    fine_boxes.maxSize(fine_max_size);
    const DistributionMapping coarse_mapping =
        shifted_mapping(coarse_boxes, coarse_owner_shift);
    const DistributionMapping fine_mapping =
        shifted_mapping(fine_boxes, fine_owner_shift);
    const auto layout = make_transfer_layout();

    MultiFab coarse_density(coarse_boxes, coarse_mapping, 1, 0);
    MultiFab coarse_measure(coarse_boxes, coarse_mapping, 1, 0);
    MultiFab fine_density(fine_boxes, fine_mapping, 1, 0);
    MultiFab fine_measure(fine_boxes, fine_mapping, 1, 0);
    fill_carriers(coarse_density, coarse_measure);
    fill_carriers(fine_density, fine_measure);
    MultiFab coarse_state(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    MultiFab fine_state(fine_boxes, fine_mapping, layout.ncomp(), 0);
    MultiFab fine_target_state(fine_boxes, fine_mapping, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_density, layout);
    fill_spectrum(fine_state, fine_density, layout);
    fine_target_state.setVal(Real(0.0));

    MultiFab coarse_candidate(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    coarse_candidate.setVal(Real(-11.0));
    const auto fine_view = timed_view(fine_state, fine_density, fine_measure, 0.5);
    const auto coarse_view = timed_view(coarse_state, coarse_density,
                                        coarse_measure, 0.5);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, IntVect(2), 0,
        coarse_candidate, diagnostic)) << diagnostic;
    EXPECT_TRUE(globally_admissible(coarse_candidate, layout));

    MultiFab coarse_mapped(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    MultiFab fine_mapped(fine_boxes, fine_mapping, layout.ncomp(), 0);
    for (amrex::MFIter mfi(coarse_candidate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto u = coarse_candidate.const_array(mfi);
        const auto omega = coarse_measure.const_array(mfi);
        const auto out = coarse_mapped.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(box, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int component) noexcept {
                out(i, j, k, component) = omega(i, j, k, 0) * u(i, j, k, component);
            });
    }
    for (amrex::MFIter mfi(fine_state, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto u = fine_state.const_array(mfi);
        const auto omega = fine_measure.const_array(mfi);
        const auto out = fine_mapped.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(box, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int component) noexcept {
                out(i, j, k, component) = omega(i, j, k, 0) * u(i, j, k, component);
            });
    }

    std::vector<Real> result;
    for (int component = 0; component < layout.ncomp(); ++component) {
        const Real coarse_inventory = coarse_mapped.sum(component);
        const Real fine_inventory = fine_mapped.sum(component) / Real(8.0);
        EXPECT_NEAR(coarse_inventory, fine_inventory,
                    Real(512.0) * std::numeric_limits<Real>::epsilon() *
                        std::max(Real(1.0), std::abs(fine_inventory)))
            << "restriction component=" << component;
        result.push_back(coarse_inventory);
    }
    append_fingerprint(coarse_candidate, result);

    MultiFab fine_candidate(fine_boxes, fine_mapping, layout.ncomp(), 0);
    fine_candidate.setVal(Real(-17.0));
    const auto fine_target_view = timed_view(
        fine_target_state, fine_density, fine_measure, 0.5);
    EXPECT_TRUE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_target_view, coarse_geometry, fine_geometry,
        IntVect(2), 1, fine_candidate, diagnostic)) << diagnostic;
    EXPECT_TRUE(globally_admissible(fine_candidate, layout));

    MultiFab prolongation_error(fine_boxes, fine_mapping, layout.ncomp(), 0);
    const ComponentOffsets offsets = offsets_for(layout);
    for (amrex::MFIter mfi(fine_candidate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto candidate = fine_candidate.const_array(mfi);
        const auto rho = fine_density.const_array(mfi);
        const auto error = prolongation_error.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(box, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int component) noexcept {
                const Real expected = rho(i, j, k, 0) * base_component(
                    component, i / 2, j / 2, k / 2, offsets);
                error(i, j, k, component) =
                    amrex::Math::abs(candidate(i, j, k, component) - expected);
            });
    }
    for (int component = 0; component < layout.ncomp(); ++component) {
        EXPECT_LE(prolongation_error.norm0(component),
                  Real(256.0) * std::numeric_limits<Real>::epsilon())
            << "prolongation component=" << component;
    }
    append_fingerprint(fine_candidate, result);
    return result;
}

} // namespace

TEST(SBMAMRTransferParallel, DecompositionInvariantRestrictionAndProlongation)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "this transfer qualification requires at least two MPI ranks";
    }

    const auto one_box = run_transfer_decomposition(8, 16, 0, 0);
    const auto split_boxes = run_transfer_decomposition(4, 8, 1, 0);
    ASSERT_EQ(one_box.size(), split_boxes.size());
    for (std::size_t index = 0; index < one_box.size(); ++index) {
        EXPECT_NEAR(one_box[index], split_boxes[index],
                    Real(1024.0) * std::numeric_limits<Real>::epsilon() *
                        std::max({Real(1.0), std::abs(one_box[index]),
                                  std::abs(split_boxes[index])}))
            << "decomposition fingerprint index=" << index;
    }
}

TEST(SBMAMRTransferParallel, RankLocalInvalidSpectrumFailsCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "the rank-local invalid-state control requires at least two MPI ranks";
    }

    const int invalid_owner = 1;
    Box coarse_domain(IntVect(0), IntVect(7));
    Box fine_domain = coarse_domain;
    fine_domain.refine(IntVect(2));
    BoxArray coarse_boxes(coarse_domain);
    BoxArray fine_boxes(fine_domain);
    coarse_boxes.maxSize(4);
    fine_boxes.maxSize(8);
    const DistributionMapping coarse_mapping = shifted_mapping(coarse_boxes, 1);
    const DistributionMapping fine_mapping = shifted_mapping(fine_boxes, 0);
    const auto layout = make_transfer_layout();

    MultiFab coarse_density(coarse_boxes, coarse_mapping, 1, 0);
    MultiFab coarse_measure(coarse_boxes, coarse_mapping, 1, 0);
    MultiFab fine_density(fine_boxes, fine_mapping, 1, 0);
    MultiFab fine_measure(fine_boxes, fine_mapping, 1, 0);
    fill_carriers(coarse_density, coarse_measure);
    fill_carriers(fine_density, fine_measure);
    MultiFab coarse_state(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    MultiFab fine_state(fine_boxes, fine_mapping, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_density, layout);
    fill_spectrum(fine_state, fine_density, layout);

    int invalid_box = -1;
    for (int box = 0; box < fine_boxes.size(); ++box) {
        if (fine_mapping[box] == invalid_owner) {
            invalid_box = box;
            break;
        }
    }
    ASSERT_GE(invalid_box, 0);
    const Box invalid_box_region = fine_boxes[invalid_box];
    const IntVect invalid_cell = invalid_box_region.smallEnd();
    int local_wrote_invalid = 0;
    const int invalid_component = layout.populations()[0].mass_offset;
    for (amrex::MFIter mfi(fine_state); mfi.isValid(); ++mfi) {
        if (mfi.index() != invalid_box) continue;
        local_wrote_invalid = 1;
        const Box cell(invalid_cell, invalid_cell);
        const auto state = fine_state.array(mfi);
        amrex::ParallelFor(cell,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                state(i, j, k, invalid_component) = Real(-1.0);
            });
    }
    amrex::ParallelDescriptor::ReduceIntSum(local_wrote_invalid);
    EXPECT_EQ(local_wrote_invalid, 1);

    MultiFab candidate(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    candidate.setVal(Real(123.0));
    const auto fine_view = timed_view(fine_state, fine_density, fine_measure, 0.5);
    const auto coarse_view = timed_view(coarse_state, coarse_density,
                                        coarse_measure, 0.5);
    std::string diagnostic;
    const bool accepted = erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, IntVect(2), 0, candidate, diagnostic);

    int minimum_accepted = accepted ? 1 : 0;
    int maximum_accepted = minimum_accepted;
    amrex::ParallelDescriptor::ReduceIntMin(minimum_accepted);
    amrex::ParallelDescriptor::ReduceIntMax(maximum_accepted);
    EXPECT_EQ(minimum_accepted, 0) << diagnostic;
    EXPECT_EQ(maximum_accepted, 0) << diagnostic;
    for (int component = 0; component < layout.ncomp(); ++component) {
        EXPECT_DOUBLE_EQ(candidate.min(component), Real(123.0))
            << "candidate component=" << component;
        EXPECT_DOUBLE_EQ(candidate.max(component), Real(123.0))
            << "candidate component=" << component;
    }
}

TEST(SBMAMRTransferParallel, RankLocalPositiveAverageUnderflowFailsCollectively)
{
    if (amrex::ParallelDescriptor::NProcs() < 2) {
        GTEST_SKIP() << "collective averaging-underflow qualification requires "
                        "at least two MPI ranks";
    }

    Box coarse_domain(IntVect(0), IntVect(3));
    Box fine_domain = coarse_domain;
    fine_domain.refine(IntVect(2));
    BoxArray coarse_boxes(coarse_domain);
    BoxArray fine_boxes(fine_domain);
    coarse_boxes.maxSize(2);
    fine_boxes.maxSize(2);
    const DistributionMapping coarse_mapping = shifted_mapping(coarse_boxes, 0);
    const DistributionMapping fine_mapping = shifted_mapping(fine_boxes, 0);
    const auto layout = make_transfer_layout();

    MultiFab coarse_density(coarse_boxes, coarse_mapping, 1, 0);
    MultiFab coarse_measure(coarse_boxes, coarse_mapping, 1, 0);
    MultiFab fine_density(fine_boxes, fine_mapping, 1, 0);
    MultiFab fine_measure(fine_boxes, fine_mapping, 1, 0);
    coarse_density.setVal(Real(1.0));
    coarse_measure.setVal(Real(1.0));
    fine_density.setVal(Real(1.0));
    fine_measure.setVal(Real(1.0));
    MultiFab coarse_state(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    MultiFab fine_state(fine_boxes, fine_mapping, layout.ncomp(), 0);
    coarse_state.setVal(Real(0.0));
    fine_state.setVal(Real(0.0));

    const Real smallest = std::numeric_limits<Real>::denorm_min();
    ASSERT_GT(smallest, Real(0.0));
    ASSERT_EQ(smallest / Real(8.0), Real(0.0));
    int invalid_box = -1;
    for (int box = 0; box < fine_boxes.size(); ++box) {
        if (fine_mapping[box] == 1) {
            invalid_box = box;
            break;
        }
    }
    ASSERT_GE(invalid_box, 0);
    const Box invalid_region = fine_boxes[invalid_box];
    const IntVect invalid_cell = invalid_region.smallEnd();
    const int component = layout.populations()[0].mass_offset;
    int local_wrote_invalid = 0;
    for (amrex::MFIter mfi(fine_state, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        if (mfi.index() != invalid_box)
            continue;
        local_wrote_invalid = 1;
        const Box cell(invalid_cell, invalid_cell);
        const auto state = fine_state.array(mfi);
        amrex::ParallelFor(cell,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                               state(i, j, k, component) = smallest;
                           });
    }
    amrex::ParallelDescriptor::ReduceIntSum(local_wrote_invalid);
    EXPECT_EQ(local_wrote_invalid, 1);
    std::string diagnostic;
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(fine_state, layout, 1,
                                                        &diagnostic))
        << diagnostic;

    MultiFab fine_before(fine_boxes, fine_mapping, layout.ncomp(), 0);
    MultiFab coarse_before(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    MultiFab::Copy(fine_before, fine_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_before, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab candidate(coarse_boxes, coarse_mapping, layout.ncomp(), 0);
    candidate.setVal(Real(31.0));
    const auto fine_view =
        timed_view(fine_state, fine_density, fine_measure, 0.25);
    const auto coarse_view =
        timed_view(coarse_state, coarse_density, coarse_measure, 0.25);
    const bool accepted = erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, IntVect(2), 0, candidate, diagnostic);
    int minimum_accepted = accepted ? 1 : 0;
    int maximum_accepted = minimum_accepted;
    amrex::ParallelDescriptor::ReduceIntMin(minimum_accepted);
    amrex::ParallelDescriptor::ReduceIntMax(maximum_accepted);
    EXPECT_EQ(minimum_accepted, 0) << diagnostic;
    EXPECT_EQ(maximum_accepted, 0) << diagnostic;
    EXPECT_NE(diagnostic.find("positive"), std::string::npos) << diagnostic;

    MultiFab fine_source_error(fine_boxes, fine_mapping, layout.ncomp(), 0);
    MultiFab coarse_source_error(coarse_boxes, coarse_mapping, layout.ncomp(),
                                 0);
    MultiFab::Copy(fine_source_error, fine_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Subtract(fine_source_error, fine_before, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_source_error, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Subtract(coarse_source_error, coarse_before, 0, 0, layout.ncomp(),
                       0);
    for (int comp = 0; comp < layout.ncomp(); ++comp) {
        EXPECT_EQ(fine_source_error.norm0(comp), Real(0.0))
            << "fine component=" << comp;
        EXPECT_EQ(coarse_source_error.norm0(comp), Real(0.0))
            << "coarse component=" << comp;
    }
}
