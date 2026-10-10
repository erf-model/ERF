#include <gtest/gtest.h>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_SBMAMRTransfer.H"
#include "ERF_SBMRestart.H"

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
    return Real(1.0) + Real(0.025) * i + Real(0.015) * j + Real(0.01) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real base_component (const int component, const int i, const int j, const int k,
                     const ComponentOffsets offsets,
                     const bool empty_upper_bin) noexcept
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
        if (empty_upper_bin && bin == 1 && i == 1 && j == 1 && k == 1) return Real(0.0);
        const Real number = Real(0.3) * a * static_cast<Real>(bin + 1);
        const Real mean = bin == 0 ? Real(0.35) : Real(1.35);
        return number * mean;
    }
    if (component >= offsets.two_number && component < offsets.two_number + 2) {
        const int bin = component - offsets.two_number;
        if (empty_upper_bin && bin == 1 && i == 1 && j == 1 && k == 1) return Real(0.0);
        return Real(0.3) * a * static_cast<Real>(bin + 1);
    }
    if (component >= offsets.two_property && component < offsets.two_property + 2) {
        const int bin = component - offsets.two_property;
        if (empty_upper_bin && bin == 1 && i == 1 && j == 1 && k == 1) return Real(0.0);
        return Real(0.15) * a * static_cast<Real>(bin + 1);
    }
    return Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real coarse_density (const int i, const int j, const int k) noexcept
{
    return Real(1.0) + Real(0.01) * i + Real(0.02) * j + Real(0.03) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real coarse_measure (const int i, const int j, const int k) noexcept
{
    return Real(0.7) + Real(0.02) * i + Real(0.01) * j + Real(0.015) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real restriction_fine_density (const int i, const int j, const int k) noexcept
{
    return Real(0.8) + Real(0.004) * i + Real(0.006) * j + Real(0.003) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real restriction_fine_measure (const int i, const int j, const int k) noexcept
{
    return Real(0.45) + Real(0.012) * i + Real(0.009) * j + Real(0.007) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real prolong_fine_measure (const int i, const int j, const int k) noexcept
{
    constexpr Real delta = Real(0.2);
    const int ci = i / 2;
    const int cj = j / 2;
    const int ck = k / 2;
    const int parity = (i % 2 + j % 2 + k % 2) % 2;
    const Real sign = parity == 0 ? Real(-1.0) : Real(1.0);
    return coarse_measure(ci, cj, ck) * (Real(1.0) + delta * sign);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real prolong_fine_density (const int i, const int j, const int k) noexcept
{
    constexpr Real delta = Real(0.2);
    const int ci = i / 2;
    const int cj = j / 2;
    const int ck = k / 2;
    const int parity = (i % 2 + j % 2 + k % 2) % 2;
    const Real sign = parity == 0 ? Real(-1.0) : Real(1.0);
    return coarse_density(ci, cj, ck) * (Real(1.0) - delta * sign) /
           (Real(1.0) - delta * delta);
}

void fill_carriers (MultiFab& density, MultiFab& measure, const bool coarse,
                    const bool prolongation_fine = false)
{
    for (amrex::MFIter mfi(density, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto rho = density.array(mfi);
        const auto omega = measure.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (coarse) {
                rho(i, j, k, 0) = coarse_density(i, j, k);
                omega(i, j, k, 0) = coarse_measure(i, j, k);
            } else if (prolongation_fine) {
                rho(i, j, k, 0) = prolong_fine_density(i, j, k);
                omega(i, j, k, 0) = prolong_fine_measure(i, j, k);
            } else {
                rho(i, j, k, 0) = restriction_fine_density(i, j, k);
                omega(i, j, k, 0) = restriction_fine_measure(i, j, k);
            }
        });
    }
}

void fill_spectrum (MultiFab& spectrum, const MultiFab& density,
                    const erf_sbm::SBMLayout& layout,
                    const bool empty_upper_bin = false)
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
                    component, i, j, k, offsets, empty_upper_bin);
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

void expect_same_values (const MultiFab& actual, const MultiFab& expected)
{
    ASSERT_EQ(actual.nComp(), expected.nComp());
    MultiFab difference(actual.boxArray(), actual.DistributionMap(),
                        actual.nComp(), 0);
    MultiFab::Copy(difference, actual, 0, 0, actual.nComp(), 0);
    MultiFab::Subtract(difference, expected, 0, 0, actual.nComp(), 0);
    for (int component = 0; component < actual.nComp(); ++component) {
        EXPECT_EQ(difference.norm0(component), Real(0.0))
            << "component=" << component;
    }
}

Geometry make_geometry (const Box& domain)
{
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    return Geometry(domain, &physical, amrex::CoordSys::cartesian, periodic);
}

Geometry
make_geometry_with_nonperiodic_direction (const Box& domain,
                                          const int direction)
{
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    int periodic[AMREX_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        periodic[dir] = dir == direction ? 0 : 1;
    }
    return Geometry(domain, &physical, amrex::CoordSys::cartesian, periodic);
}

Box coarse_domain ()
{
    return Box(IntVect(0, 0, 0), IntVect(3, 3, 3));
}

Box fine_coverage ()
{
    return Box(IntVect(2, 2, 2), IntVect(5, 5, 5));
}

void run_restriction_averages_mapped_amount_and_preserves_uncovered_state ()
{
    const auto layout = make_transfer_layout();
    const BoxArray coarse_ba(coarse_domain());
    const DistributionMapping coarse_dm(coarse_ba);
    const BoxArray fine_ba(fine_coverage());
    const DistributionMapping fine_dm(fine_ba);
    MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
    fill_carriers(coarse_rho, coarse_omega, true);
    fill_carriers(fine_rho, fine_omega, false);

    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_state(fine_ba, fine_dm, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_rho, layout);
    fill_spectrum(fine_state, fine_rho, layout);
    MultiFab coarse_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(coarse_before, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab fine_before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab::Copy(fine_before, fine_state, 0, 0, layout.ncomp(), 0);
    MultiFab candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    candidate.setVal(Real(99.0));

    constexpr int ratio_value = 2;
    const IntVect ratio(ratio_value);
    auto fine_view = timed_view(fine_state, fine_rho, fine_omega, 0.5);
    auto coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    std::string diagnostic;

    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(coarse_state, layout, 0,
                                                        &diagnostic))
        << diagnostic;
    MultiFab uncovered_roundtrip_change(coarse_ba, coarse_dm, 1, 0);
    uncovered_roundtrip_change.setVal(Real(0.0));
    for (amrex::MFIter mfi(coarse_state, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto state = coarse_state.const_array(mfi);
        const auto omega = coarse_omega.const_array(mfi);
        const auto changed = uncovered_roundtrip_change.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const bool covered =
                    i >= 1 && i <= 2 && j >= 1 && j <= 2 && k >= 1 && k <= 2;
                int found_change = 0;
                if (!covered) {
                    const Real measure = omega(i, j, k, 0);
                    for (int component = 0; component < ncomp; ++component) {
                        const Real value = state(i, j, k, component);
                        if ((measure * value) / measure != value)
                            found_change = 1;
                    }
                }
                changed(i, j, k, 0) = static_cast<Real>(found_change);
            });
    }
    ASSERT_EQ(uncovered_roundtrip_change.max(0), Real(1.0))
        << "fixture must distinguish the old uncovered H round trip";

    ASSERT_TRUE(erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, ratio, 0, candidate, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
        candidate, layout, 0, &diagnostic)) << diagnostic;

    // Independently evaluate every child amount from its analytic fixture
    // values.  The expected candidate includes the preexisting uncovered U.
    MultiFab expected(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab arithmetic_u(coarse_ba, coarse_dm, layout.ncomp(), 0);
    const ComponentOffsets offsets = offsets_for(layout);
    const int ncomp = layout.ncomp();
    for (amrex::MFIter mfi(expected, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto expected_array = expected.array(mfi);
        const auto arithmetic_array = arithmetic_u.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real omega_c = coarse_measure(i, j, k);
            const Real rho_c = coarse_density(i, j, k);
            const bool covered = i >= 1 && i <= 2 && j >= 1 && j <= 2 && k >= 1 && k <= 2;
            for (int component = 0; component < ncomp; ++component) {
                if (!covered) {
                    expected_array(i, j, k, component) = rho_c *
                        base_component(component, i, j, k, offsets, false);
                    arithmetic_array(i, j, k, component) =
                        expected_array(i, j, k, component);
                    continue;
                }
                Real mapped_sum = Real(0.0);
                Real u_sum = Real(0.0);
                for (int oz = 0; oz < ratio_value; ++oz) {
                    for (int oy = 0; oy < ratio_value; ++oy) {
                        for (int ox = 0; ox < ratio_value; ++ox) {
                            const int fi = ratio_value * i + ox;
                            const int fj = ratio_value * j + oy;
                            const int fk = ratio_value * k + oz;
                            const Real rho_f = restriction_fine_density(fi, fj, fk);
                            const Real u = rho_f * base_component(
                                component, fi, fj, fk, offsets, false);
                            const Real omega_f = restriction_fine_measure(fi, fj, fk);
                            mapped_sum += omega_f * u;
                            u_sum += u;
                        }
                    }
                }
                constexpr Real child_count = Real(8.0);
                expected_array(i, j, k, component) =
                    (mapped_sum / child_count) / omega_c;
                arithmetic_array(i, j, k, component) = u_sum / child_count;
            }
        });
    }
    MultiFab difference(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(difference, candidate, 0, 0, ncomp, 0);
    MultiFab::Subtract(difference, expected, 0, 0, ncomp, 0);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_LE(difference.norm0(component),
                  Real(64.0) * std::numeric_limits<Real>::epsilon())
            << "component=" << component;
    }

    MultiFab uncovered_error(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab coarse_source_error(coarse_ba, coarse_dm, ncomp, 0);
    for (amrex::MFIter mfi(candidate, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const Box box = mfi.tilebox();
        const auto actual = candidate.const_array(mfi);
        const auto original = coarse_before.const_array(mfi);
        const auto source = coarse_state.const_array(mfi);
        const auto uncovered = uncovered_error.array(mfi);
        const auto source_error = coarse_source_error.array(mfi);
        amrex::ParallelFor(
            box, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int component) noexcept {
                const bool covered =
                    i >= 1 && i <= 2 && j >= 1 && j <= 2 && k >= 1 && k <= 2;
                uncovered(i, j, k, component) =
                    covered ? Real(0.0)
                            : amrex::Math::abs(actual(i, j, k, component) -
                                               original(i, j, k, component));
                source_error(i, j, k, component) = amrex::Math::abs(
                    source(i, j, k, component) - original(i, j, k, component));
            });
    }
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_EQ(uncovered_error.norm0(component), Real(0.0))
            << "uncovered component=" << component;
        EXPECT_EQ(coarse_source_error.norm0(component), Real(0.0))
            << "coarse source component=" << component;
    }

    MultiFab arithmetic_difference(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(arithmetic_difference, candidate, 0, 0, ncomp, 0);
    MultiFab::Subtract(arithmetic_difference, arithmetic_u, 0, 0, ncomp, 0);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_GT(arithmetic_difference.norm0(component),
                  Real(1.0e-3) * candidate.norm0(component))
            << "component=" << component;
    }

    MultiFab::Subtract(fine_state, fine_before, 0, 0, ncomp, 0);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_EQ(fine_state.norm0(component), Real(0.0))
            << "component=" << component;
    }

    // A noncanonical fine source fails before any candidate can be committed.
    MultiFab invalid_fine(fine_ba, fine_dm, ncomp, 0);
    MultiFab::Copy(invalid_fine, fine_before, 0, 0, ncomp, 0);
    invalid_fine.setVal(Real(-1.0), layout.populations()[1].mass_offset, 1, 0);
    candidate.setVal(Real(123.0));
    auto invalid_view = timed_view(invalid_fine, fine_rho, fine_omega, 0.5);
    EXPECT_FALSE(erf_sbm::RestrictMappedSpectrum(
        layout, invalid_view, coarse_view, ratio, 0, candidate, diagnostic));
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_DOUBLE_EQ(candidate.min(component), Real(123.0))
            << "component=" << component;
        EXPECT_DOUBLE_EQ(candidate.max(component), Real(123.0))
            << "component=" << component;
    }
}

TEST(SBMAMRTransfer, RestrictionAveragesMappedAmountAndPreservesUncoveredState)
{
    run_restriction_averages_mapped_amount_and_preserves_uncovered_state();
}

void
run_restriction_late_quotient_failure_is_atomic ()
{
    const auto layout = make_transfer_layout();
    const Box c_domain = coarse_domain();
    Box f_domain = c_domain;
    f_domain.refine(IntVect(2));
    const BoxArray coarse_ba(c_domain);
    const BoxArray fine_ba(f_domain);
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);

    MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
    fill_carriers(coarse_rho, coarse_omega, true);
    fill_carriers(fine_rho, fine_omega, false);
    coarse_omega.setVal(std::numeric_limits<Real>::min());

    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_state(fine_ba, fine_dm, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_rho, layout);
    fill_spectrum(fine_state, fine_rho, layout);
    fine_state.mult(Real(1.0e30), 0, layout.ncomp(), 0);
    std::string diagnostic;
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(fine_state, layout, 1,
                                                        &diagnostic))
        << diagnostic;

    MultiFab candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    candidate.setVal(Real(41.0));
    MultiFab candidate_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab coarse_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab coarse_rho_before(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega_before(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho_before(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega_before(fine_ba, fine_dm, 1, 0);
    MultiFab::Copy(candidate_before, candidate, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_before, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(fine_before, fine_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_rho_before, coarse_rho, 0, 0, 1, 0);
    MultiFab::Copy(coarse_omega_before, coarse_omega, 0, 0, 1, 0);
    MultiFab::Copy(fine_rho_before, fine_rho, 0, 0, 1, 0);
    MultiFab::Copy(fine_omega_before, fine_omega, 0, 0, 1, 0);

    const auto fine_view = timed_view(fine_state, fine_rho, fine_omega, 0.5);
    const auto coarse_view =
        timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    EXPECT_FALSE(erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, IntVect(2), 0, candidate, diagnostic));
    EXPECT_NE(diagnostic.find("recovering coarse spectrum"), std::string::npos)
        << diagnostic;
    expect_same_values(candidate, candidate_before);
    expect_same_values(coarse_state, coarse_before);
    expect_same_values(fine_state, fine_before);
    expect_same_values(coarse_rho, coarse_rho_before);
    expect_same_values(coarse_omega, coarse_omega_before);
    expect_same_values(fine_rho, fine_rho_before);
    expect_same_values(fine_omega, fine_omega_before);
}

TEST(SBMAMRTransfer, RestrictionLateQuotientFailureIsAtomic)
{
    run_restriction_late_quotient_failure_is_atomic();
}

void run_restriction_rejects_positive_mapped_average_underflow ()
{
    const auto layout = make_transfer_layout();
    const BoxArray coarse_ba(coarse_domain());
    const DistributionMapping coarse_dm(coarse_ba);
    const Box fine_box(IntVect(0, 0, 0), IntVect(1, 1, 1));
    const BoxArray fine_ba(fine_box);
    const DistributionMapping fine_dm(fine_ba);
    MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
    coarse_rho.setVal(Real(1.0));
    coarse_omega.setVal(Real(1.0));
    fine_rho.setVal(Real(1.0));
    fine_omega.setVal(Real(1.0));

    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_state(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    coarse_state.setVal(Real(0.0));
    fine_state.setVal(Real(0.0));
    candidate.setVal(Real(29.0));
    const int component = layout.populations()[0].mass_offset;
    const Real smallest = std::numeric_limits<Real>::denorm_min();
    ASSERT_GT(smallest, Real(0.0));
    volatile Real runtime_smallest = smallest;
    const Real averaged = static_cast<Real>(runtime_smallest) / Real(8.0);
    ASSERT_EQ(averaged, Real(0.0))
        << "the selected precision must round this one-child average to zero";

    const auto set_one_child = [&fine_state, component] (const Real amount) {
        fine_state.setVal(Real(0.0));
        for (amrex::MFIter mfi(fine_state, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const Box valid_box = mfi.validbox();
            const IntVect first_cell = valid_box.smallEnd();
            const Box cell(first_cell, first_cell);
            const auto state = fine_state.array(mfi);
            amrex::ParallelFor(
                cell, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    state(i, j, k, component) = amount;
                });
        }
    };

    set_one_child(smallest);
    std::string diagnostic;
    ASSERT_TRUE(fine_state.is_finite(0, layout.ncomp(), 0));
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(fine_state, layout, 1,
                                                        &diagnostic))
        << diagnostic;
    MultiFab fine_before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab coarse_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(fine_before, fine_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_before, coarse_state, 0, 0, layout.ncomp(), 0);

    const auto fine_view = timed_view(fine_state, fine_rho, fine_omega, 0.0);
    const auto coarse_view =
        timed_view(coarse_state, coarse_rho, coarse_omega, 0.0);
    EXPECT_FALSE(erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, IntVect(2), 0, candidate, diagnostic));
    EXPECT_NE(diagnostic.find("underflow"), std::string::npos) << diagnostic;

    MultiFab fine_source_error(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab coarse_source_error(coarse_ba, coarse_dm, layout.ncomp(), 0);
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

    // Eight subnormal units average to one representable subnormal unit.
    const Real representable_amount = Real(8.0) * smallest;
    ASSERT_GT(representable_amount, smallest);
    set_one_child(representable_amount);
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(fine_state, layout, 1,
                                                        &diagnostic))
        << diagnostic;
    EXPECT_TRUE(erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, IntVect(2), 0, candidate, diagnostic))
        << diagnostic;
    MultiFab expected(coarse_ba, coarse_dm, layout.ncomp(), 0);
    expected.setVal(Real(0.0));
    for (amrex::MFIter mfi(expected, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const IntVect coarse_cell(0, 0, 0);
        const Box cell(coarse_cell, coarse_cell);
        const auto state = expected.array(mfi);
        amrex::ParallelFor(cell,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                               state(i, j, k, component) = smallest;
                           });
    }
    MultiFab positive_control_error(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(positive_control_error, candidate, 0, 0, layout.ncomp(), 0);
    MultiFab::Subtract(positive_control_error, expected, 0, 0, layout.ncomp(),
                       0);
    for (int comp = 0; comp < layout.ncomp(); ++comp) {
        EXPECT_EQ(positive_control_error.norm0(comp), Real(0.0))
            << "positive-control component=" << comp;
    }
}

TEST(SBMAMRTransfer, RestrictionRejectsPositiveMappedAverageUnderflow)
{
    run_restriction_rejects_positive_mapped_average_underflow();
}

void run_prolongation_uses_same_time_piecewise_constant_carrier_ratio ()
{
    const auto layout = make_transfer_layout();
    const Box c_domain = coarse_domain();
    Box f_domain = c_domain;
    f_domain.refine(IntVect(2));
    const Geometry cgeom = make_geometry(c_domain);
    const Geometry fgeom = make_geometry(f_domain);
    const BoxArray coarse_ba(c_domain);
    const DistributionMapping coarse_dm(coarse_ba);
    const BoxArray fine_ba(fine_coverage());
    const DistributionMapping fine_dm(fine_ba);

    MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
    fill_carriers(coarse_rho, coarse_omega, true);
    fill_carriers(fine_rho, fine_omega, false, true);
    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_storage(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab candidate(fine_ba, fine_dm, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_rho, layout, true);
    fine_storage.setVal(Real(0.0));
    candidate.setVal(Real(77.0));

    auto coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    auto fine_view = timed_view(fine_storage, fine_rho, fine_omega, 0.5);
    std::string diagnostic;
    ASSERT_TRUE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
        candidate, layout, 1, &diagnostic)) << diagnostic;

    MultiFab expected(fine_ba, fine_dm, layout.ncomp(), 0);
    const ComponentOffsets offsets = offsets_for(layout);
    const int ncomp = layout.ncomp();
    for (amrex::MFIter mfi(expected, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto rho = fine_rho.const_array(mfi);
        const auto result = expected.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const int ci = i / 2;
            const int cj = j / 2;
            const int ck = k / 2;
            const Real carrier = rho(i, j, k, 0);
            for (int component = 0; component < ncomp; ++component) {
                result(i, j, k, component) = carrier * base_component(
                    component, ci, cj, ck, offsets, true);
            }
        });
    }
    MultiFab difference(fine_ba, fine_dm, ncomp, 0);
    MultiFab::Copy(difference, candidate, 0, 0, ncomp, 0);
    MultiFab::Subtract(difference, expected, 0, 0, ncomp, 0);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_LE(difference.norm0(component),
                  Real(64.0) * std::numeric_limits<Real>::epsilon())
            << "component=" << component;
    }

    // Every child carries the parent dry-air-relative state, including the
    // 2M mass/number means and attached material ratios.
    const auto& two = layout.populations()[1];
    MultiFab ratio_error(fine_ba, fine_dm, 2, 0);
    for (amrex::MFIter mfi(ratio_error, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto state = candidate.const_array(mfi);
        const auto errors = ratio_error.array(mfi);
        const int mass0 = two.mass_offset;
        const int mass1 = two.mass_offset + 1;
        const int number0 = two.number_offset;
        const int number1 = two.number_offset + 1;
        const int property0 = layout.property_offset(1);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const int ci = i / 2;
            const int cj = j / 2;
            const int ck = k / 2;
            const bool empty_upper = ci == 1 && cj == 1 && ck == 1;
            const Real n0 = state(i, j, k, number0);
            const Real n1 = state(i, j, k, number1);
            errors(i, j, k, 0) = n0 > Real(0.0)
                ? amrex::Math::abs(state(i, j, k, mass0) / n0 - Real(0.35))
                : Real(0.0);
            errors(i, j, k, 1) = empty_upper
                ? amrex::Math::abs(n1) + amrex::Math::abs(state(i, j, k, mass1)) +
                      amrex::Math::abs(state(i, j, k, property0 + 1))
                : (n1 > Real(0.0)
                       ? amrex::Math::abs(state(i, j, k, mass1) / n1 - Real(1.35)) +
                             amrex::Math::abs(state(i, j, k, property0 + 1) / n1 - Real(0.5))
                       : Real(0.0));
        });
    }
    EXPECT_LE(ratio_error.norm0(0),
              Real(64.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_LE(ratio_error.norm0(1),
              Real(64.0) * std::numeric_limits<Real>::epsilon());

    // Compare mapped inventories over the refined patch.  Fine cell volumes
    // are 1/8 of coarse computational volumes at this refinement ratio.
    MultiFab coarse_inventory(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab fine_inventory(fine_ba, fine_dm, ncomp, 0);
    coarse_inventory.setVal(Real(0.0));
    for (amrex::MFIter mfi(coarse_inventory, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto u = coarse_state.const_array(mfi);
        const auto omega = coarse_omega.const_array(mfi);
        const auto h = coarse_inventory.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const bool covered = i >= 1 && i <= 2 && j >= 1 && j <= 2 && k >= 1 && k <= 2;
            for (int component = 0; component < ncomp; ++component) {
                h(i, j, k, component) = covered
                    ? omega(i, j, k, 0) * u(i, j, k, component)
                    : Real(0.0);
            }
        });
    }
    for (amrex::MFIter mfi(fine_inventory, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto u = candidate.const_array(mfi);
        const auto omega = fine_omega.const_array(mfi);
        const auto h = fine_inventory.array(mfi);
        amrex::ParallelFor(box, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int component) noexcept {
            h(i, j, k, component) = omega(i, j, k, 0) * u(i, j, k, component);
        });
    }
    for (int component = 0; component < ncomp; ++component) {
        const Real coarse_sum = coarse_inventory.sum(component);
        const Real fine_sum = fine_inventory.sum(component) / Real(8.0);
        EXPECT_NEAR(coarse_sum, fine_sum,
                    Real(128.0) * std::numeric_limits<Real>::epsilon() *
                        std::max(Real(1.0), amrex::Math::abs(coarse_sum)));
    }

    // A tuple with a wrong-time density is rejected before candidate writes.
    candidate.setVal(Real(77.0));
    auto stale_fine_view = fine_view;
    stale_fine_view.density_time = 0.75;
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, stale_fine_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic));
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_DOUBLE_EQ(candidate.min(component), Real(77.0))
            << "component=" << component;
        EXPECT_DOUBLE_EQ(candidate.max(component), Real(77.0))
            << "component=" << component;
    }
    EXPECT_NE(diagnostic.find("times do not match"), std::string::npos);

    candidate.setVal(Real(77.0));
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 0,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("fine level"), std::string::npos);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_DOUBLE_EQ(candidate.min(component), Real(77.0))
            << "component=" << component;
        EXPECT_DOUBLE_EQ(candidate.max(component), Real(77.0))
            << "component=" << component;
    }

    auto incomplete_view = fine_view;
    incomplete_view.dry_air_density = nullptr;
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, incomplete_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("incomplete"), std::string::npos);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_DOUBLE_EQ(candidate.min(component), Real(77.0))
            << "component=" << component;
        EXPECT_DOUBLE_EQ(candidate.max(component), Real(77.0))
            << "component=" << component;
    }

    // Exact same-time views also work at a later semantic time after the
    // carrier and conservative spectrum are both changed consistently.
    coarse_rho.mult(Real(2.0), 0, 1, 0);
    fine_rho.mult(Real(2.0), 0, 1, 0);
    coarse_state.mult(Real(2.0), 0, ncomp, 0);
    coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 1.5);
    fine_view = timed_view(fine_storage, fine_rho, fine_omega, 1.5);
    ASSERT_TRUE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic)) << diagnostic;
    candidate.mult(Real(0.5), 0, ncomp, 0);
    MultiFab::Subtract(candidate, expected, 0, 0, ncomp, 0);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_LE(candidate.norm0(component),
                  Real(64.0) * std::numeric_limits<Real>::epsilon())
            << "component=" << component;
    }
}

TEST(SBMAMRTransfer, ProlongationUsesSameTimePiecewiseConstantCarrierRatio)
{
    run_prolongation_uses_same_time_piecewise_constant_carrier_ratio();
}

TEST(SBMAMRTransfer, ProlongationRejectsQuotientAndProductUnderflow)
{
    const auto layout = make_transfer_layout();
    const Box c_domain = coarse_domain();
    Box f_domain = c_domain;
    f_domain.refine(IntVect(2));
    const Geometry cgeom = make_geometry(c_domain);
    const Geometry fgeom = make_geometry(f_domain);
    const BoxArray coarse_ba(c_domain);
    const DistributionMapping coarse_dm(coarse_ba);
    const BoxArray fine_ba(fine_coverage());
    const DistributionMapping fine_dm(fine_ba);
    MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
    fill_carriers(coarse_rho, coarse_omega, true);
    fill_carriers(fine_rho, fine_omega, false, true);
    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_storage(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab candidate(fine_ba, fine_dm, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_rho, layout);
    fine_storage.setVal(Real(0.0));
    candidate.setVal(Real(33.0));
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    const Real tiny_scale = Real(1000.0) * std::numeric_limits<Real>::min();
    ASSERT_GT(tiny, Real(0.0));

    coarse_state.mult(tiny_scale, 0, layout.ncomp(), 0);
    coarse_rho.setVal(std::numeric_limits<Real>::max());
    auto coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    auto fine_view = timed_view(fine_storage, fine_rho, fine_omega, 0.5);
    std::string diagnostic;
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("underflowed"), std::string::npos);

    coarse_rho.setVal(Real(1.0));
    fill_spectrum(coarse_state, coarse_rho, layout);
    coarse_state.mult(tiny_scale, 0, layout.ncomp(), 0);
    fine_rho.setVal(tiny);
    fine_omega.setVal(Real(1.0));
    candidate.setVal(Real(34.0));
    MultiFab candidate_before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab coarse_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_rho_before(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega_before(fine_ba, fine_dm, 1, 0);
    MultiFab coarse_state_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(candidate_before, candidate, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_before, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(fine_rho_before, fine_rho, 0, 0, 1, 0);
    MultiFab::Copy(fine_omega_before, fine_omega, 0, 0, 1, 0);
    MultiFab::Copy(coarse_state_before, coarse_state, 0, 0, layout.ncomp(), 0);
    coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    fine_view = timed_view(fine_storage, fine_rho, fine_omega, 0.5);
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("underflowed"), std::string::npos);
    expect_same_values(candidate, candidate_before);
    expect_same_values(coarse_state, coarse_state_before);
    expect_same_values(fine_rho, fine_rho_before);
    expect_same_values(fine_omega, fine_omega_before);

    // A finite but large fine carrier makes rho_f*z overflow only during the
    // last reconstruction pass, after the old implementation had written the
    // public destination.
    coarse_rho.setVal(Real(1.0));
    fill_spectrum(coarse_state, coarse_rho, layout);
    coarse_state.mult(Real(1.0e30), 0, layout.ncomp(), 0);
    fine_rho.setVal(std::numeric_limits<Real>::max());
    fine_omega.setVal(Real(1.0));
    candidate.setVal(Real(35.0));
    MultiFab::Copy(candidate_before, candidate, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_state_before, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(fine_rho_before, fine_rho, 0, 0, 1, 0);
    MultiFab::Copy(fine_omega_before, fine_omega, 0, 0, 1, 0);
    coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    fine_view = timed_view(fine_storage, fine_rho, fine_omega, 0.5);
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1, candidate,
        diagnostic));
    EXPECT_NE(diagnostic.find("overflowed or underflowed"), std::string::npos)
        << diagnostic;
    expect_same_values(candidate, candidate_before);
    expect_same_values(coarse_state, coarse_state_before);
    expect_same_values(fine_rho, fine_rho_before);
    expect_same_values(fine_omega, fine_omega_before);
}

TEST(SBMAMRTransfer, RestrictionRejectsNoncoarsenableFineBoxesBeforeAverage)
{
    const auto layout = make_transfer_layout();
    const Box coarse_domain(IntVect(0), IntVect(2));
    const Box fine_domain(IntVect(0), IntVect(7));
    const BoxArray coarse_ba(coarse_domain);
    const BoxArray fine_ba(fine_domain);
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const IntVect ratio(3);
    ASSERT_FALSE(fine_ba.coarsenable(ratio));

    MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
    fill_carriers(coarse_rho, coarse_omega, true);
    fill_carriers(fine_rho, fine_omega, false);
    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_state(fine_ba, fine_dm, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_rho, layout);
    fill_spectrum(fine_state, fine_rho, layout);
    MultiFab coarse_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab::Copy(coarse_before, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(fine_before, fine_state, 0, 0, layout.ncomp(), 0);

    MultiFab candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab candidate_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    candidate.setVal(Real(-23.0));
    MultiFab::Copy(candidate_before, candidate, 0, 0, layout.ncomp(), 0);
    const auto fine_view = timed_view(fine_state, fine_rho, fine_omega, 0.3);
    const auto coarse_view =
        timed_view(coarse_state, coarse_rho, coarse_omega, 0.3);
    std::string diagnostic;
    EXPECT_FALSE(erf_sbm::RestrictMappedSpectrum(
        layout, fine_view, coarse_view, ratio, 0, candidate, diagnostic));
    EXPECT_NE(diagnostic.find("coarsenable"), std::string::npos) << diagnostic;
    expect_same_values(candidate, candidate_before);
    expect_same_values(coarse_state, coarse_before);
    expect_same_values(fine_state, fine_before);
}

void
run_prolongation_diagnostic_precedence ()
{
    const auto layout = make_transfer_layout();
    const Box c_domain = coarse_domain();
    Box f_domain = c_domain;
    f_domain.refine(IntVect(2));
    const BoxArray coarse_ba(c_domain);
    const BoxArray fine_ba(f_domain);
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const Geometry cgeom = make_geometry(c_domain);
    const Geometry nonperiodic_fgeom =
        make_geometry_with_nonperiodic_direction(f_domain, 0);
    Box bad_f_domain = f_domain;
    IntVect bad_hi = bad_f_domain.bigEnd();
    bad_hi[0] -= 1;
    bad_f_domain = Box(bad_f_domain.smallEnd(), bad_hi);
    const Geometry bad_fgeom = make_geometry(bad_f_domain);

    MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
    fill_carriers(coarse_rho, coarse_omega, true);
    fill_carriers(fine_rho, fine_omega, false, true);
    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_target_state(fine_ba, fine_dm, layout.ncomp(), 0);
    fill_spectrum(coarse_state, coarse_rho, layout, true);
    fine_target_state.setVal(Real(0.0));
    const auto coarse_view =
        timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    const auto fine_view =
        timed_view(fine_target_state, fine_rho, fine_omega, 0.5);
    MultiFab coarse_before(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab coarse_rho_before(coarse_ba, coarse_dm, 1, 0);
    MultiFab coarse_omega_before(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_rho_before(fine_ba, fine_dm, 1, 0);
    MultiFab fine_omega_before(fine_ba, fine_dm, 1, 0);
    MultiFab::Copy(coarse_before, coarse_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(fine_before, fine_target_state, 0, 0, layout.ncomp(), 0);
    MultiFab::Copy(coarse_rho_before, coarse_rho, 0, 0, 1, 0);
    MultiFab::Copy(coarse_omega_before, coarse_omega, 0, 0, 1, 0);
    MultiFab::Copy(fine_rho_before, fine_rho, 0, 0, 1, 0);
    MultiFab::Copy(fine_omega_before, fine_omega, 0, 0, 1, 0);

    std::string diagnostic;
    // The earlier alias/layout defect must survive a later geometry failure.
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, nonperiodic_fgeom, IntVect(2), 1,
        fine_target_state, diagnostic));
    EXPECT_NE(diagnostic.find("separate from both input tuples"),
              std::string::npos)
        << diagnostic;
    expect_same_values(fine_target_state, fine_before);

    MultiFab candidate(fine_ba, fine_dm, layout.ncomp(), 0);
    candidate.setVal(Real(83.0));
    MultiFab candidate_before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab::Copy(candidate_before, candidate, 0, 0, layout.ncomp(), 0);
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, nonperiodic_fgeom, IntVect(2), 1,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("periodic geometry"), std::string::npos)
        << diagnostic;
    expect_same_values(candidate, candidate_before);

    candidate.setVal(Real(84.0));
    MultiFab::Copy(candidate_before, candidate, 0, 0, layout.ncomp(), 0);
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, bad_fgeom, IntVect(2), 1,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("does not map the coarse domain"),
              std::string::npos)
        << diagnostic;
    expect_same_values(candidate, candidate_before);
    expect_same_values(coarse_state, coarse_before);
    expect_same_values(fine_target_state, fine_before);
    expect_same_values(coarse_rho, coarse_rho_before);
    expect_same_values(coarse_omega, coarse_omega_before);
    expect_same_values(fine_rho, fine_rho_before);
    expect_same_values(fine_omega, fine_omega_before);
}

TEST(SBMAMRTransfer, ProlongationPreservesDiagnosticPrecedence)
{
    run_prolongation_diagnostic_precedence();
}

} // namespace
