#include <gtest/gtest.h>

#include <AMReX_BoxArray.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_YAFluxRegister.H>

#include "ERF_SBMAMRTransfer.H"
#include "ERF_SBMAMRFluxRegister.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMAMRTime.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace {

using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::Real;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ya_oracle_coarse_flux (int dir, int i, int j, int k, int n) noexcept
{
    return Real(1.0) + Real(3.0) * dir + Real(0.1) * i + Real(0.2) * j +
           Real(0.3) * k + Real(0.01) * n;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ya_oracle_fine_flux (int dir, int i, int j, int k, int n) noexcept
{
    return Real(2.0) + Real(5.0) * dir + Real(0.01) * i + Real(0.02) * j +
           Real(0.03) * k + Real(0.001) * n;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real accepted_chunk_rate (int dir, int i, int j, int k, int global_component,
                          int /*method*/, int stage, int substep) noexcept
{
    const Real low = Real(0.4) + Real(0.03) * dir + Real(0.001) * i +
                     Real(0.002) * j + Real(0.003) * k +
                     Real(0.004) * global_component;
    const Real high = low + Real(0.8) + Real(0.01) * stage +
                      Real(0.02) * substep;
    const Real lambda = Real(0.35) + Real(0.01) * (global_component % 4);
    return low + lambda * (high - low);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real signed_mapping_oracle_rate (int dir, int i, int j, int k,
                                 int global_component, bool fine_side) noexcept
{
    const Real sign = global_component % 2 == 0 ? Real(1.0) : Real(-1.0);
    const Real level_offset = fine_side ? Real(0.0018) : Real(0.0008);
    const Real magnitude = Real(0.001) + level_offset + Real(0.0003) * dir +
        Real(0.00002) * i + Real(0.00003) * j + Real(0.00004) * k +
        Real(0.00001) * global_component;
    return sign * magnitude;
}

erf_auxiliary::MappedFaceFluxRate make_signed_mapped_rate (
    const BoxArray& cell_ba, const DistributionMapping& mapping,
    const std::vector<int>& global_components, const bool fine_side)
{
    erf_auxiliary::MappedFaceFluxRate rate;
    rate.define(cell_ba, mapping, static_cast<int>(global_components.size()), 0);
    amrex::Gpu::DeviceVector<int> device_components(global_components.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, global_components.begin(),
                     global_components.end(), device_components.begin());
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto& field = rate.dir(dir);
        for (MFIter mfi(field, false); mfi.isValid(); ++mfi) {
            const Box box = mfi.validbox();
            const auto values = field.array(mfi);
            const int* globals = device_components.data();
            const int ncomp = static_cast<int>(global_components.size());
            ParallelFor(box, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int local) noexcept {
                values(i,j,k,local) = signed_mapping_oracle_rate(
                    dir, i, j, k, globals[local], fine_side);
            });
        }
    }
    amrex::Gpu::streamSynchronize();
    return rate;
}

void accumulate_mapped_reference_integral (
    erf_auxiliary::IntegratedMappedFaceFlux& integral,
    const erf_auxiliary::MappedFaceFluxRate& rate,
    const std::vector<int>& global_components, const double weight)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int local = 0; local < rate.nComp(); ++local) {
            MultiFab::Saxpy(integral.dir(dir), static_cast<Real>(weight),
                            rate.dir(dir), local,
                            global_components[static_cast<std::size_t>(local)],
                            1, 0);
        }
    }
}

std::uint64_t contiguous_component_runs (const std::vector<int>& components)
{
    std::uint64_t runs = components.empty() ? 0U : 1U;
    for (std::size_t index = 1; index < components.size(); ++index) {
        if (components[index] != components[index - 1U] + 1) ++runs;
    }
    return runs;
}

erf_auxiliary::MappedFaceFluxRate make_rate_chunk (
    const BoxArray& cell_ba, const DistributionMapping& mapping,
    const int first_global_component, const int ncomp,
    const int method, const int stage, const int substep)
{
    erf_auxiliary::MappedFaceFluxRate rate;
    rate.define(cell_ba, mapping, ncomp, 0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto& field = rate.dir(dir);
        for (MFIter mfi(field, false); mfi.isValid(); ++mfi) {
            const Box box = mfi.validbox();
            const auto flux = field.array(mfi);
            ParallelFor(box, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int local) noexcept {
                flux(i,j,k,local) = accepted_chunk_rate(
                    dir, i, j, k, first_global_component + local,
                    method, stage, substep);
            });
        }
    }
    amrex::Gpu::streamSynchronize();
    return rate;
}

void accumulate_reference_integral (
    erf_auxiliary::IntegratedMappedFaceFlux& integral,
    const erf_auxiliary::MappedFaceFluxRate& rate,
    const int first_global_component, const double weight)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int local = 0; local < rate.nComp(); ++local) {
            MultiFab::Saxpy(integral.dir(dir), static_cast<Real>(weight),
                            rate.dir(dir), local,
                            first_global_component + local, 1, 0);
        }
    }
}

void add_integral_to_reference_register (
    amrex::YAFluxRegister& flux_register,
    const erf_auxiliary::IntegratedMappedFaceFlux& integral,
    const BoxArray& cell_ba, const DistributionMapping& mapping,
    const amrex::Geometry& geometry, const bool coarse_side)
{
    MultiFab cell_layout(cell_ba, mapping, 1, 0,
                         amrex::MFInfo().SetAlloc(false));
    const auto dx = geometry.CellSizeArray();
    for (MFIter mfi(cell_layout, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        std::array<amrex::FArrayBox const*, AMREX_SPACEDIM> faces{};
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            faces[static_cast<std::size_t>(dir)] = &integral.dir(dir)[mfi];
        }
        if (coarse_side) {
            flux_register.CrseAdd(mfi, faces, dx.data(), Real(1.0),
                                  amrex::RunOn::Host);
        } else {
            flux_register.FineAdd(mfi, faces, dx.data(), Real(1.0),
                                  amrex::RunOn::Host);
        }
    }
}

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

erf_sbm::SBMLayout make_one_moment_layout (const int ncomp)
{
    erf_sbm::SBMLayoutSpec spec;
    erf_sbm::SpectralPopulationSpec population;
    population.population_id = 0;
    population.semantic_id = "memory_census_liquid";
    population.phase = erf_sbm::PopulationPhase::Liquid;
    population.moment_mode = erf_sbm::MomentMode::OneMoment;
    population.mass_state_units = "kg m^-3";
    population.number_state_units = "m^-3";
    population.grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    population.grid.coordinate_units = "kg";
    population.grid.edges.resize(static_cast<std::size_t>(ncomp) + 1U);
    population.grid.pivots.resize(static_cast<std::size_t>(ncomp));
    for (int bin = 0; bin <= ncomp; ++bin) {
        population.grid.edges[static_cast<std::size_t>(bin)] =
            static_cast<Real>(bin);
        if (bin < ncomp) {
            population.grid.pivots[static_cast<std::size_t>(bin)] =
                static_cast<Real>(bin) + Real(0.5);
        }
    }
    spec.populations.push_back(std::move(population));
    spec.liquid_projection = {0, ncomp / 2};
    return erf_sbm::SBMLayout(std::move(spec));
}

void fill_positive_one_moment_state (MultiFab& state)
{
    const int ncomp = state.nComp();
    for (MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto values = state.array(mfi);
        ParallelFor(box, ncomp, [=] AMREX_GPU_DEVICE
                    (int i, int j, int k, int n) noexcept {
            values(i,j,k,n) = Real(0.25) + Real(0.001) * static_cast<Real>(n);
        });
    }
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

Real value_at (const MultiFab& field, const IntVect& cell, const int component)
{
    for (MFIter mfi(field); mfi.isValid(); ++mfi) {
        const auto& fab = field[mfi];
        if (!mfi.validbox().contains(cell)) continue;
        Real value = Real(0.0);
        const auto offset = fab.box().index(cell);
        amrex::Gpu::dtoh_memcpy_async(
            &value, fab.dataPtr(component) + offset, sizeof(Real));
        amrex::Gpu::streamSynchronize();
        return value;
    }
    ADD_FAILURE() << "cell is not owned by this rank: " << cell;
    return Real(0.0);
}

void expect_roundoff_equivalent (const MultiFab& actual,
                                 const MultiFab& expected)
{
    ASSERT_EQ(actual.nComp(), expected.nComp());
    MultiFab difference(actual.boxArray(), actual.DistributionMap(),
                        actual.nComp(), 0);
    MultiFab::Copy(difference, actual, 0, 0, actual.nComp(), 0);
    MultiFab::Subtract(difference, expected, 0, 0, actual.nComp(), 0);
    constexpr Real operation_reordering_bound =
        Real(64.0) * std::numeric_limits<Real>::epsilon();
    for (int component = 0; component < actual.nComp(); ++component) {
        const Real scale = std::max(
            Real(1.0), std::max(actual.norm0(component), expected.norm0(component)));
        EXPECT_LE(difference.norm0(component), operation_reordering_bound * scale)
            << "component=" << component;
    }
}

Geometry make_geometry (const Box& domain)
{
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    return Geometry(domain, &physical, amrex::CoordSys::cartesian, periodic);
}

void fill_ya_oracle_measure (MultiFab& mapped_measure)
{
    for (MFIter mfi(mapped_measure, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto omega = mapped_measure.array(mfi);
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            omega(i,j,k,0) = Real(1.5) + Real(0.01) * i +
                             Real(0.02) * j + Real(0.03) * k;
        });
    }
}

void fill_ya_oracle_flux (MultiFab& flux, const int dir, const bool fine)
{
    const int ncomp = flux.nComp();
    for (MFIter mfi(flux, false); mfi.isValid(); ++mfi) {
        const Box box = mfi.validbox();
        const auto values = flux.array(mfi);
        ParallelFor(box, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            values(i,j,k,n) = fine ? ya_oracle_fine_flux(dir, i, j, k, n)
                                   : ya_oracle_coarse_flux(dir, i, j, k, n);
        });
    }
}

void set_cell_value (MultiFab& field, const IntVect& cell,
                     const int component, const Real value)
{
    for (MFIter mfi(field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto values = field.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (i == cell[0] && j == cell[1] && k == cell[2]) {
                values(i,j,k,component) = value;
            }
        });
    }
}

void roundtrip_through_measure (MultiFab& values, const MultiFab& measure)
{
    const int ncomp = values.nComp();
    for (MFIter mfi(values, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto state = values.array(mfi);
        const auto omega = measure.const_array(mfi);
        ParallelFor(box, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            state(i,j,k,n) = omega(i,j,k,0) * state(i,j,k,n) / omega(i,j,k,0);
        });
    }
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

    // Reproduce ERF::timeStep's last fine start (time + (i-1)*dt_fine)
    // followed by Advance's target (fine_start + dt_fine), against the
    // independently formed coarse endpoint (time + dt_coarse).
    const double coarse_begin = 1000.0;
    const double coarse_dt = 0.1;
    const double fine_dt = coarse_dt / 2.0;
    const double coarse_end = coarse_begin + coarse_dt;
    const double final_fine_begin = coarse_begin + 1.0 * fine_dt;
    const double final_fine_end = final_fine_begin + fine_dt;
    ASSERT_NE(final_fine_end, coarse_end);
    auto independently_timed_fine = fine_view;
    independently_timed_fine.spectrum_time = final_fine_end;
    independently_timed_fine.density_time = final_fine_end;
    independently_timed_fine.measure_time = final_fine_end;
    auto independently_timed_coarse = coarse_view;
    independently_timed_coarse.spectrum_time = coarse_end;
    independently_timed_coarse.density_time = coarse_end;
    independently_timed_coarse.measure_time = coarse_end;
    MultiFab exact_time_candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    EXPECT_FALSE(erf_sbm::RestrictMappedSpectrumPrivate(
        layout, independently_timed_fine, independently_timed_coarse,
        ratio, 0, exact_time_candidate, diagnostic));
    MultiFab bounded_time_candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    ASSERT_TRUE(erf_sbm::RestrictMappedSpectrumPrivate(
        layout, independently_timed_fine, independently_timed_coarse,
        ratio, 0, bounded_time_candidate, diagnostic, 16, false, nullptr,
        nullptr, fine_dt)) << diagnostic;
    expect_roundoff_equivalent(bounded_time_candidate, candidate);

    MultiFab private_candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    private_candidate.setVal(Real(99.0));
    erf_sbm::SBMTransferAllocationCensus forward_trace;
    ASSERT_TRUE(erf_sbm::RestrictMappedSpectrumPrivate(
        layout, fine_view, coarse_view, ratio, 0, private_candidate,
        diagnostic, 1, false, nullptr, &forward_trace)) << diagnostic;
    MultiFab private_difference(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(private_difference, private_candidate, 0, 0,
                   layout.ncomp(), 0);
    MultiFab::Subtract(private_difference, candidate, 0, 0,
                       layout.ncomp(), 0);
    for (int component = 0; component < layout.ncomp(); ++component) {
        EXPECT_LE(private_difference.norm0(component),
                  Real(64.0) * std::numeric_limits<Real>::epsilon())
            << "private synchronization candidate differs from the independent "
            << "transactional restriction oracle at component " << component;
    }
    erf_sbm::SBMTransferAllocationCensus reverse_trace;
    MultiFab reverse_candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    ASSERT_TRUE(erf_sbm::RestrictMappedSpectrumPrivate(
        layout, fine_view, coarse_view, ratio, 0, reverse_candidate,
        diagnostic, 1, true, nullptr, &reverse_trace)) << diagnostic;
    ASSERT_GE(forward_trace.processed_chunk_first_components.size(), 3U);
    ASSERT_EQ(reverse_trace.processed_chunk_first_components.size(),
              forward_trace.processed_chunk_first_components.size());
    for (std::size_t index = 0;
         index < forward_trace.processed_chunk_first_components.size(); ++index) {
        EXPECT_EQ(reverse_trace.processed_chunk_first_components[index],
                  forward_trace.processed_chunk_first_components[
                      forward_trace.processed_chunk_first_components.size() - 1U - index])
            << "reverse restriction must change actual complete-group visitation";
    }
    expect_roundoff_equivalent(reverse_candidate, candidate);

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

TEST(SBMAMRTime, SemanticComparisonsRejectNonfiniteAndAdjacentSubsteps)
{
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(erf_sbm::semantic_time_matches(inf, 1.0, 0.1));
    EXPECT_FALSE(erf_sbm::semantic_time_matches(-inf, 1.0, 0.1));
    EXPECT_FALSE(erf_sbm::semantic_time_matches(nan, nan, 0.1));
    EXPECT_FALSE(erf_sbm::semantic_time_matches_exact(inf, inf));

    const double large = 1.0e16;
    const double adjacent = std::nextafter(large, inf);
    EXPECT_FALSE(erf_sbm::semantic_time_matches(large, adjacent, 2.0));
    EXPECT_FALSE(erf_sbm::semantic_time_matches(large, large + 64.0, 2.0));

    // Two ordinary host addition orders can round a shared semantic timestamp
    // differently while remaining far inside the dt_min/32 cap.
    const double first_order = (1.0e12 + 0.1) + 0.2;
    const double second_order = 1.0e12 + (0.1 + 0.2);
    EXPECT_NE(first_order, second_order);
    EXPECT_TRUE(erf_sbm::semantic_time_matches(first_order, second_order, 1.0));
    EXPECT_FALSE(erf_sbm::semantic_time_matches_exact(first_order, second_order));

    // This is ERF's actual recursive fine-step expression at dt=0.1 with two
    // substeps, followed by the coarse timeStep endpoint expression.
    const double coarse_begin = 1000.0;
    const double coarse_dt = 0.1;
    const int fine_substeps = 2;
    const double fine_dt = coarse_dt / static_cast<double>(fine_substeps);
    const double coarse_end = coarse_begin + coarse_dt;
    const double last_fine_begin = coarse_begin +
        static_cast<double>(fine_substeps - 1) * fine_dt;
    const double last_fine_end = last_fine_begin + fine_dt;
    ASSERT_NE(last_fine_end, coarse_end);
    EXPECT_TRUE(erf_sbm::semantic_time_matches(last_fine_end, coarse_end, fine_dt));
    EXPECT_FALSE(erf_sbm::semantic_time_matches_exact(last_fine_end, coarse_end));
    EXPECT_FALSE(erf_sbm::semantic_time_matches(
        last_fine_end, last_fine_end + fine_dt, fine_dt));
    EXPECT_FALSE(erf_sbm::semantic_time_matches(last_fine_end, coarse_end, 0.0));
}

TEST(SBMAMRTime, FluxRegisterRejectsNonfiniteAndNeighboringCoarseTimes)
{
    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(3));
    const Box fine_domain(IntVect(0), IntVect(7));
    const BoxArray coarse_ba(coarse_domain);
    const BoxArray fine_ba(amrex::refine(Box(IntVect(0), IntVect(3)), ratio));
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const Geometry coarse_geom = make_geometry(coarse_domain);
    const Geometry fine_geom = make_geometry(fine_domain);
    erf_sbm::SBMAMRFluxRegister flux_register;
    flux_register.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                         fine_geom, coarse_geom, ratio, 1, 1);
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::string diagnostic;
    EXPECT_FALSE(flux_register.begin_interval(inf, 1.0, 1, diagnostic));
    EXPECT_FALSE(flux_register.begin_interval(0.0, nan, 1, diagnostic));

    const double large = 1.0e16;
    const double end = large + 64.0;
    ASSERT_TRUE(flux_register.begin_interval(large, end, 2, diagnostic))
        << diagnostic;
    EXPECT_FALSE(flux_register.begin_level_step(
        true, std::nextafter(large, inf), end, diagnostic));
    EXPECT_NE(diagnostic.find("does not match"), std::string::npos)
        << diagnostic;
    ASSERT_TRUE(flux_register.begin_level_step(true, large, end, diagnostic))
        << diagnostic;
    ASSERT_TRUE(flux_register.discard_interval(diagnostic)) << diagnostic;
}

TEST(SBMAMRTransfer, DeviceAdmissionUsesHostOracleOnlyAfterRejection)
{
    const auto layout = make_transfer_layout();
    const BoxArray grids(coarse_domain());
    const DistributionMapping mapping(grids);
    MultiFab density(grids, mapping, 1, 0);
    MultiFab measure(grids, mapping, 1, 0);
    MultiFab state(grids, mapping, layout.ncomp(), 0);
    fill_carriers(density, measure, true);
    fill_spectrum(state, density, layout);
    MultiFab valid_state(grids, mapping, layout.ncomp(), 0);
    MultiFab::Copy(valid_state, state, 0, 0, layout.ncomp(), 0);

    erf_sbm::SBMStateAdmission admission(layout);
    admission.define_level(0, grids, mapping);
    const int one_mass = layout.populations()[0].mass_offset;
    const int one_property = layout.property_offset(0);
    const int two_mass = layout.populations()[1].mass_offset;
    const int two_number = layout.populations()[1].number_offset;
    std::string diagnostic;
    const auto compare_case = [&] (const char* label, const bool expected,
                                   auto&& mutate) {
        MultiFab::Copy(state, valid_state, 0, 0, layout.ncomp(), 0);
        mutate(state);
        std::string host_diagnostic;
        const bool host_ok = erf_sbm::authoritative_state_admissible(
            state, layout, 0, &host_diagnostic);
        std::string device_diagnostic;
        const bool device_ok = admission.admissible(
            state, 0, device_diagnostic);
        EXPECT_EQ(host_ok, expected) << label << ": " << host_diagnostic;
        EXPECT_EQ(device_ok, host_ok)
            << label << ": device=" << device_diagnostic
            << ", host=" << host_diagnostic;
    };

    compare_case("valid 1M/2M state with attached properties", true,
                 [] (MultiFab&) {});
    compare_case("global terminal upper edge is admitted", true,
                 [=] (MultiFab& values) {
        values.setVal(Real(1.0), two_number + 1, 1, 0);
        values.setVal(Real(2.0), two_mass + 1, 1, 0);
    });
    compare_case("representable positive subnormal with no orphan material",
                 true, [=] (MultiFab& values) {
        values.setVal(std::numeric_limits<Real>::denorm_min(), one_mass, 1, 0);
        values.setVal(Real(0.0), one_property, 1, 0);
    });
    compare_case("interior 2M half-open upper edge is rejected", false,
                 [=] (MultiFab& values) {
        values.setVal(Real(1.0), two_number, 1, 0);
        values.setVal(Real(1.0), two_mass, 1, 0);
    });
    compare_case("finite negative 1M amount is rejected", false,
                 [=] (MultiFab& values) {
        values.setVal(Real(-1.0e-12), one_mass, 1, 0);
    });
    compare_case("NaN state is rejected", false, [=] (MultiFab& values) {
        values.setVal(std::numeric_limits<Real>::quiet_NaN(), one_mass, 1, 0);
    });
    compare_case("infinite state is rejected", false, [=] (MultiFab& values) {
        values.setVal(std::numeric_limits<Real>::infinity(), one_mass, 1, 0);
    });
    compare_case("1M attached material without carrier is rejected", false,
                 [=] (MultiFab& values) {
        values.setVal(Real(0.0), one_mass, 1, 0);
        values.setVal(Real(0.1), one_property, 1, 0);
    });
    compare_case("2M attached material without number carrier is rejected",
                 false, [=] (MultiFab& values) {
        values.setVal(Real(0.0), two_mass + 1, 1, 0);
        values.setVal(Real(0.0), two_number + 1, 1, 0);
        values.setVal(Real(0.1), layout.property_offset(1) + 1, 1, 0);
    });

    auto counts = admission.counters();
    EXPECT_EQ(counts.device_admissions, 9U);
    EXPECT_EQ(counts.successful_admissions, 3U);
    EXPECT_EQ(counts.host_scans, 6U);
    EXPECT_GT(counts.host_scan_fab_bytes, 0U);
#ifndef AMREX_USE_GPU
    EXPECT_EQ(counts.host_d2h_bytes, 0U);
#endif
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
    coarse_omega.setVal(Real(1.0));
    for (MFIter mfi(fine_omega, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto omega = fine_omega.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const int parity = (i % 2 + j % 2 + k % 2) % 2;
            const Real sign = parity == 0 ? Real(-1.0) : Real(1.0);
            omega(i, j, k, 0) = Real(1.0) + Real(0.2) * sign;
        });
    }
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
    MultiFab private_candidate(fine_ba, fine_dm, layout.ncomp(), 0);
    erf_sbm::SBMTransferAllocationCensus prolong_forward_trace;
    ASSERT_TRUE(erf_sbm::ProlongCarrierRelativeSpectrumPrivate(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        private_candidate, diagnostic, 1, false, nullptr,
        &prolong_forward_trace)) << diagnostic;
    expect_roundoff_equivalent(private_candidate, candidate);
    MultiFab reverse_candidate(fine_ba, fine_dm, layout.ncomp(), 0);
    erf_sbm::SBMTransferAllocationCensus prolong_reverse_trace;
    ASSERT_TRUE(erf_sbm::ProlongCarrierRelativeSpectrumPrivate(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        reverse_candidate, diagnostic, 1, true, nullptr,
        &prolong_reverse_trace)) << diagnostic;
    ASSERT_GE(prolong_forward_trace.processed_chunk_first_components.size(), 3U);
    ASSERT_EQ(prolong_reverse_trace.processed_chunk_first_components.size(),
              prolong_forward_trace.processed_chunk_first_components.size());
    for (std::size_t index = 0;
         index < prolong_forward_trace.processed_chunk_first_components.size(); ++index) {
        EXPECT_EQ(prolong_reverse_trace.processed_chunk_first_components[index],
                  prolong_forward_trace.processed_chunk_first_components[
                      prolong_forward_trace.processed_chunk_first_components.size() - 1U - index])
            << "reverse prolongation must change actual complete-group visitation";
    }
    expect_roundoff_equivalent(reverse_candidate, candidate);

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

TEST(SBMAMRTransfer, ProlongationRejectsMappedCarrierInventoryMismatchAtomically)
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
    MultiFab coarse_state(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab fine_storage(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab candidate(fine_ba, fine_dm, layout.ncomp(), 0);
    coarse_rho.setVal(Real(1.0));
    coarse_omega.setVal(Real(1.0));
    fine_rho.setVal(Real(2.0));
    fine_omega.setVal(Real(1.0));
    fine_storage.setVal(Real(0.0));
    fill_spectrum(coarse_state, coarse_rho, layout);
    candidate.setVal(Real(-17.0));
    MultiFab before(fine_ba, fine_dm, layout.ncomp(), 0);
    MultiFab::Copy(before, candidate, 0, 0, layout.ncomp(), 0);
    const auto coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    const auto fine_view = timed_view(fine_storage, fine_rho, fine_omega, 0.5);

    std::string diagnostic;
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("inventory differs"), std::string::npos)
        << diagnostic;
    expect_same_values(candidate, before);
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
    fine_rho.setVal(std::numeric_limits<Real>::max());
    coarse_omega.setVal(Real(1.0) / std::numeric_limits<Real>::max());
    fine_omega.setVal(Real(1.0) / std::numeric_limits<Real>::max());
    auto coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
    auto fine_view = timed_view(fine_storage, fine_rho, fine_omega, 0.5);
    std::string diagnostic;
    EXPECT_FALSE(erf_sbm::ProlongCarrierRelativeSpectrum(
        layout, coarse_view, fine_view, cgeom, fgeom, IntVect(2), 1,
        candidate, diagnostic));
    EXPECT_NE(diagnostic.find("underflowed"), std::string::npos);

    coarse_rho.setVal(Real(1.0));
    coarse_omega.setVal(Real(1.0));
    fill_spectrum(coarse_state, coarse_rho, layout);
    coarse_state.mult(tiny_scale, 0, layout.ncomp(), 0);
    fine_rho.setVal(tiny_scale);
    fine_omega.setVal(Real(1.0) / tiny_scale);
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
    coarse_omega.setVal(Real(1.0));
    fill_spectrum(coarse_state, coarse_rho, layout);
    coarse_state.mult(Real(1.0e30), 0, layout.ncomp(), 0);
    fine_rho.setVal(std::numeric_limits<Real>::max());
    fine_omega.setVal(Real(1.0) / std::numeric_limits<Real>::max());
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
    EXPECT_NE(diagnostic.find("matching coarse and fine periodicity"), std::string::npos)
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

TEST(SBMAMRYA, SignedMappedCorrectionMatchesIndependentFaceAmounts)
{
    ASSERT_EQ(AMREX_SPACEDIM, 3);
    const auto layout = make_transfer_layout();
    const int ncomp = layout.ncomp();
    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(7));
    const Box fine_domain(IntVect(0), IntVect(15));
    const Box covered_fine_box = amrex::refine(Box(IntVect(2), IntVect(5)), ratio);
    const BoxArray coarse_ba(coarse_domain);
    const BoxArray fine_ba(covered_fine_box);
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const Geometry coarse_geom = make_geometry(coarse_domain);
    const Geometry fine_geom = make_geometry(fine_domain);

    MultiFab coarse_cells(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab fine_cells(fine_ba, fine_dm, ncomp, 0);
    MultiFab mapped_measure(coarse_ba, coarse_dm, 1, 0);
    MultiFab covered(coarse_ba, coarse_dm, 1, 1);
    MultiFab fine_marker(fine_ba, fine_dm, 1, 0);
    MultiFab mapped_correction(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab correction_invalid(coarse_ba, coarse_dm, 1, 0);
    for (int n = 0; n < ncomp; ++n) {
        coarse_cells.setVal(Real(100.0) * (n + 1), n, 1, 0);
    }
    fine_marker.setVal(Real(1.0));
    covered.setVal(Real(0.0));
    amrex::average_down(fine_marker, covered, 0, 1, ratio);
    covered.FillBoundary(coarse_geom.periodicity());
    mapped_correction.setVal(Real(0.0));
    correction_invalid.setVal(Real(0.0));
    fill_ya_oracle_measure(mapped_measure);

    std::array<MultiFab, AMREX_SPACEDIM> coarse_flux;
    std::array<MultiFab, AMREX_SPACEDIM> fine_flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        BoxArray coarse_faces = coarse_ba;
        BoxArray fine_faces = fine_ba;
        coarse_faces.surroundingNodes(dir);
        fine_faces.surroundingNodes(dir);
        coarse_flux[dir].define(coarse_faces, coarse_dm, ncomp, 0);
        fine_flux[dir].define(fine_faces, fine_dm, ncomp, 0);
        fill_ya_oracle_flux(coarse_flux[dir], dir, false);
        fill_ya_oracle_flux(fine_flux[dir], dir, true);
    }
    amrex::Gpu::streamSynchronize();

    amrex::YAFluxRegister register_for_sbm(
        fine_ba, coarse_ba, fine_dm, coarse_dm, fine_geom, coarse_geom,
        ratio, 1, ncomp);
    register_for_sbm.reset();
    const auto coarse_dx = coarse_geom.CellSizeArray();
    const auto fine_dx = fine_geom.CellSizeArray();
    for (MFIter mfi(coarse_cells, false); mfi.isValid(); ++mfi) {
        std::array<amrex::FArrayBox const*, AMREX_SPACEDIM> fluxes{
            &coarse_flux[0][mfi], &coarse_flux[1][mfi], &coarse_flux[2][mfi]};
        register_for_sbm.CrseAdd(mfi, fluxes, coarse_dx.data(), Real(1.0),
                                 amrex::RunOn::Host);
    }
    for (MFIter mfi(fine_cells, false); mfi.isValid(); ++mfi) {
        std::array<amrex::FArrayBox const*, AMREX_SPACEDIM> fluxes{
            &fine_flux[0][mfi], &fine_flux[1][mfi], &fine_flux[2][mfi]};
        register_for_sbm.FineAdd(mfi, fluxes, fine_dx.data(), Real(1.0),
                                 amrex::RunOn::Host);
    }
    register_for_sbm.Reflux(mapped_correction, 0, 0, ncomp);
    amrex::Gpu::streamSynchronize();
    std::vector<int> all_components(static_cast<std::size_t>(ncomp));
    for (int component = 0; component < ncomp; ++component) {
        all_components[static_cast<std::size_t>(component)] = component;
    }
    std::string diagnostic;
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        coarse_cells, mapped_correction, mapped_measure, covered,
        all_components, correction_invalid, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(
        correction_invalid, diagnostic)) << diagnostic;

    const auto coarse_dx_product =
        coarse_dx[0] * coarse_dx[1] * coarse_dx[2];
    const int tangential[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    for (MFIter mfi(coarse_cells, false); mfi.isValid(); ++mfi) {
        const Box box = mfi.validbox();
        const auto state = coarse_cells.const_array(mfi);
        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    const int c[3] = {i, j, k};
                    const bool cell_covered = i >= 2 && i <= 5 &&
                                         j >= 2 && j <= 5 &&
                                         k >= 2 && k <= 5;
                    const Real omega = Real(1.5) + Real(0.01) * i +
                                       Real(0.02) * j + Real(0.03) * k;
                    for (int n = 0; n < ncomp; ++n) {
                        const Real initial = Real(100.0) * (n + 1);
                        Real expected_delta_h = Real(0.0);
                        Real expected_delta_q = Real(0.0);
                        if (!cell_covered) {
                            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                                for (int side : {-1, 1}) {
                                    const int neighbor = c[dir] + side;
                                    const bool neighbor_covered = neighbor >= 2 &&
                                        neighbor <= 5 &&
                                        c[(dir + 1) % 3] >= 2 &&
                                        c[(dir + 1) % 3] <= 5 &&
                                        c[(dir + 2) % 3] >= 2 &&
                                        c[(dir + 2) % 3] <= 5;
                                    if (!neighbor_covered) continue;

                                    const int face = c[dir] + (side > 0 ? 1 : 0);
                                    int coarse_face[3] = {i, j, k};
                                    coarse_face[dir] = face;
                                    const Real coarse_integral =
                                        ya_oracle_coarse_flux(
                                            dir, coarse_face[0], coarse_face[1],
                                            coarse_face[2], n);

                                    int fine_face[3] = {2*i, 2*j, 2*k};
                                    fine_face[dir] = 2 * face;
                                    Real fine_sum = Real(0.0);
                                    for (int a = 0; a < 2; ++a) {
                                        for (int b = 0; b < 2; ++b) {
                                            int subface[3] = {fine_face[0],
                                                              fine_face[1],
                                                              fine_face[2]};
                                            subface[tangential[dir][0]] += a;
                                            subface[tangential[dir][1]] += b;
                                            fine_sum += ya_oracle_fine_flux(
                                                dir, subface[0], subface[1],
                                                subface[2], n);
                                        }
                                    }
                                    const Real sigma = side > 0
                                        ? Real(1.0) : Real(-1.0);
                                    expected_delta_h += sigma / coarse_dx[dir] *
                                        (coarse_integral - fine_sum / Real(4.0));

                                    const Real coarse_area =
                                        coarse_dx_product / coarse_dx[dir];
                                    const Real fine_cell_volume =
                                        coarse_dx_product / Real(8.0);
                                    const Real fine_area =
                                        fine_cell_volume / fine_dx[dir];
                                    expected_delta_q += sigma *
                                        (coarse_area * coarse_integral -
                                         fine_area * fine_sum);
                                }
                            }
                        }
                        const Real expected = initial + expected_delta_h / omega;
                        EXPECT_NEAR(state(i,j,k,n), expected, Real(2.0e-11))
                            << "cell=" << i << ',' << j << ',' << k
                            << " component=" << n;
                        const Real delta_u = state(i,j,k,n) - initial;
                        const Real mapped_amount = coarse_dx_product * omega * delta_u;
                        EXPECT_NEAR(mapped_amount,
                                     expected_delta_q,
                                     Real(2.0e-11));
                        if (cell_covered) {
                            EXPECT_EQ(state(i,j,k,n), initial)
                                << "covered coarse parent changed at "
                                << i << ',' << j << ',' << k
                                << " component=" << n;
                        }
                    }
                }
            }
        }
    }
}

TEST(SBMAMRYA, LocalizedCorrectionPreservesCoveredAndFarfieldBits)
{
    const Box domain(IntVect(0), IntVect(3));
    const IntVect ratio(2, 3, 2);
    const Geometry geometry = make_geometry(domain);
    BoxArray coarse_ba(domain);
    coarse_ba.maxSize(2);
    const DistributionMapping coarse_dm(coarse_ba);
    const BoxArray fine_ba(amrex::refine(Box(IntVect(0), IntVect(0)), ratio));
    const DistributionMapping fine_dm(fine_ba);
    MultiFab marker(fine_ba, fine_dm, 1, 0);
    MultiFab covered(coarse_ba, coarse_dm, 1, 1);
    marker.setVal(Real(1.0));
    covered.setVal(Real(0.0));
    amrex::average_down(marker, covered, 0, 1, ratio);
    covered.FillBoundary(geometry.periodicity());

    MultiFab candidate(coarse_ba, coarse_dm, 2, 0);
    MultiFab before(coarse_ba, coarse_dm, 2, 0);
    MultiFab measure(coarse_ba, coarse_dm, 1, 0);
    MultiFab correction(coarse_ba, coarse_dm, 2, 0);
    MultiFab invalid(coarse_ba, coarse_dm, 1, 0);
    candidate.setVal(Real(0.1));
    measure.setVal(Real(0.1));
    correction.setVal(Real(0.0));
    invalid.setVal(Real(0.0));
    set_cell_value(candidate, IntVect(2, 2, 2), 1, -Real(0.0));
    set_cell_value(candidate, IntVect(0, 0, 0), 1, -Real(0.0));
    set_cell_value(correction, IntVect(3, 0, 0), 0, Real(0.005));
    set_cell_value(correction, IntVect(0, 3, 0), 1, -Real(0.002));
    MultiFab::Copy(before, candidate, 0, 0, 2, 0);
    std::string diagnostic;
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        candidate, correction, measure, covered, {0, 1}, invalid,
        diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(invalid, diagnostic))
        << diagnostic;
    EXPECT_EQ(value_at(candidate, IntVect(3,0,0), 0),
              (Real(0.1) * Real(0.1) + Real(0.005)) / Real(0.1));
    EXPECT_EQ(value_at(candidate, IntVect(0,3,0), 1),
              (Real(0.1) * Real(0.1) - Real(0.002)) / Real(0.1));
    EXPECT_TRUE(std::signbit(value_at(candidate, IntVect(2,2,2), 1)));
    EXPECT_TRUE(std::signbit(value_at(candidate, IntVect(0,0,0), 1)));
    for (int k = 0; k <= 3; ++k) for (int j = 0; j <= 3; ++j)
        for (int i = 0; i <= 3; ++i) for (int component = 0; component < 2; ++component) {
            const bool corrected =
                (i == 3 && j == 0 && k == 0 && component == 0) ||
                (i == 0 && j == 3 && k == 0 && component == 1);
            if (corrected) continue;
            const Real actual = value_at(candidate, IntVect(i,j,k), component);
            const Real original = value_at(before, IntVect(i,j,k), component);
            EXPECT_EQ(std::memcmp(&actual, &original, sizeof(Real)), 0)
                << "bitwise identity failed at " << i << ',' << j << ',' << k
                << " component=" << component;
        }

    MultiFab zero_delta(coarse_ba, coarse_dm, 1, 0);
    MultiFab zero_invalid(coarse_ba, coarse_dm, 1, 0);
    MultiFab zero_candidate(coarse_ba, coarse_dm, 2, 0);
    zero_delta.setVal(Real(0.0));
    zero_invalid.setVal(Real(0.0));
    MultiFab::Copy(zero_candidate, before, 0, 0, 2, 0);
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        zero_candidate, zero_delta, measure, covered, {0}, zero_invalid,
        diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(
        zero_invalid, diagnostic)) << diagnostic;
    const Real zero_actual = value_at(zero_candidate, IntVect(2,2,2), 1);
    const Real zero_original = value_at(before, IntVect(2,2,2), 1);
    EXPECT_EQ(std::memcmp(&zero_actual, &zero_original, sizeof(Real)), 0);

    MultiFab bad_delta(coarse_ba, coarse_dm, 1, 0);
    MultiFab bad_invalid(coarse_ba, coarse_dm, 1, 0);
    MultiFab bad_candidate(coarse_ba, coarse_dm, 2, 0);
    bad_delta.setVal(Real(0.0));
    bad_invalid.setVal(Real(0.0));
    MultiFab::Copy(bad_candidate, before, 0, 0, 2, 0);
    set_cell_value(bad_delta, IntVect(2, 2, 2), 0, Real(0.001));
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        bad_candidate, bad_delta, measure, covered, {0}, bad_invalid,
        diagnostic)) << diagnostic;
    EXPECT_FALSE(erf_sbm::MappedRefluxCorrectionMaskIsValid(
        bad_invalid, diagnostic));
    EXPECT_NE(diagnostic.find("outside"), std::string::npos) << diagnostic;

    MultiFab global_roundtrip(coarse_ba, coarse_dm, 2, 0);
    MultiFab::Copy(global_roundtrip, before, 0, 0, 2, 0);
    roundtrip_through_measure(global_roundtrip, measure);
    MultiFab::Subtract(global_roundtrip, before, 0, 0, 2, 0);
    EXPECT_GT(global_roundtrip.norm0(0), Real(0.0));
}

TEST(SBMAMRYA, DirectAcceptedStageRatesMatchCompletedIntegralReference)
{
    const auto layout = make_transfer_layout();
    const int ncomp = layout.ncomp();
    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(7));
    const Box fine_domain(IntVect(0), IntVect(15));
    const BoxArray coarse_ba(coarse_domain);
    const BoxArray fine_ba(amrex::refine(Box(IntVect(2), IntVect(5)), ratio));
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const Geometry coarse_geom = make_geometry(coarse_domain);
    const Geometry fine_geom = make_geometry(fine_domain);
    const double coarse_begin = 0.0;
    const double coarse_end = 0.3;
    constexpr int fine_substeps = 3;
    const double fine_dt = (coarse_end - coarse_begin) / fine_substeps;

    for (const auto method : {erf_auxiliary::HostIntegrator::CompressibleRK3,
                              erf_auxiliary::HostIntegrator::AnelasticHeun,
                              erf_auxiliary::HostIntegrator::AnelasticMidPoint}) {
        const int method_id =
            method == erf_auxiliary::HostIntegrator::CompressibleRK3 ? 0
            : method == erf_auxiliary::HostIntegrator::AnelasticHeun ? 1 : 2;
        const int stages = method_id == 0 ? 3 : 2;
        // RK3 and MidPoint restart every stage from the old state, so only
        // the final stage's rate is the step's flux; Heun averages its two.
        const auto stage_weight = [=] (int stage, double dt) {
            if (method_id == 0) return stage == 2 ? dt : 0.0;
            if (method_id == 2) return stage == 1 ? dt : 0.0;
            return 0.5 * dt;
        };

        erf_sbm::SBMAMRFluxRegister direct_register;
        direct_register.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                               fine_geom, coarse_geom, ratio, 1, ncomp);
        std::string diagnostic;
        ASSERT_TRUE(direct_register.begin_interval(
            coarse_begin, coarse_end, fine_substeps, diagnostic)) << diagnostic;

        erf_auxiliary::IntegratedMappedFaceFlux coarse_integral;
        coarse_integral.define(coarse_ba, coarse_dm, ncomp, 0);
        coarse_integral.setVal(Real(0.0));
        erf_auxiliary::IntegratedMappedFaceFlux fine_integral;
        fine_integral.define(fine_ba, fine_dm, ncomp, 0);
        fine_integral.setVal(Real(0.0));

        auto advance_reference_step = [&] (
            const bool coarse_side, const double step_begin,
            const double step_end, const double dt, const int substep) {
            const auto& ba = coarse_side ? coarse_ba : fine_ba;
            const auto& dm = coarse_side ? coarse_dm : fine_dm;
            auto& integral = coarse_side ? coarse_integral : fine_integral;
            for (int stage = 0; stage < stages; ++stage) {
                EXPECT_TRUE(direct_register.begin_stage(
                    coarse_side, method, stage, step_begin, diagnostic))
                    << diagnostic;
                const double weight = stage_weight(stage, dt);
                erf_auxiliary::AuxiliaryStageRecipe recipe;
                ASSERT_TRUE(erf_auxiliary::MakeAuxiliaryStageRecipe(
                    method, stage, dt, recipe, diagnostic)) << diagnostic;
                EXPECT_EQ(recipe.completed_ledger_time, weight)
                    << "method=" << method_id << " stage=" << stage;
                // Reverse chunk order on alternating stages/substeps.  The
                // global component identity must be independent of traversal.
                for (int chunk = 1; chunk >= 0; --chunk) {
                    const int first = chunk * 5;
                    auto accepted = make_rate_chunk(
                        ba, dm, first, 5, method_id, stage, substep);
                    const std::vector<int> component_map{
                        first, first + 1, first + 2, first + 3, first + 4};
                    ASSERT_TRUE(direct_register.add_stage_rates(
                        coarse_side, accepted, component_map, weight,
                        diagnostic)) << diagnostic;
                    accumulate_reference_integral(
                        integral, accepted, first, weight);
                }
                EXPECT_TRUE(direct_register.finish_stage(
                    coarse_side, diagnostic)) << diagnostic;
            }
            EXPECT_TRUE(direct_register.accept_level_step(
                coarse_side, step_begin, step_end, diagnostic)) << diagnostic;
        };

        ASSERT_TRUE(direct_register.begin_level_step(
            true, coarse_begin, coarse_end, diagnostic)) << diagnostic;
        advance_reference_step(true, coarse_begin, coarse_end,
                               coarse_end - coarse_begin, 0);

        MultiFab direct_state(coarse_ba, coarse_dm, ncomp, 0);
        MultiFab reference_state(coarse_ba, coarse_dm, ncomp, 0);
        MultiFab before_reflux(coarse_ba, coarse_dm, ncomp, 0);
        for (int n = 0; n < ncomp; ++n) {
            direct_state.setVal(Real(0.25) * (n + 1), n, 1, 0);
            reference_state.setVal(Real(0.25) * (n + 1), n, 1, 0);
        }

        for (int substep = 0; substep < fine_substeps; ++substep) {
            const double step_begin = coarse_begin + substep * fine_dt;
            const double step_end = step_begin + fine_dt;
            ASSERT_TRUE(direct_register.begin_level_step(
                false, step_begin, step_end, diagnostic)) << diagnostic;
            advance_reference_step(false, step_begin, step_end,
                                   fine_dt, substep + 1);
            if (substep == 0) {
                EXPECT_FALSE(direct_register.ready_to_reflux());
                MultiFab::Copy(before_reflux, direct_state, 0, 0, ncomp, 0);
                EXPECT_FALSE(direct_register.reflux(direct_state, diagnostic));
                expect_same_values(direct_state, before_reflux);
            }
        }
        ASSERT_TRUE(direct_register.ready_to_reflux());

        amrex::YAFluxRegister reference_register(
            fine_ba, coarse_ba, fine_dm, coarse_dm, fine_geom, coarse_geom,
            ratio, 1, ncomp);
        reference_register.reset();
        add_integral_to_reference_register(
            reference_register, coarse_integral, coarse_ba, coarse_dm,
            coarse_geom, true);
        add_integral_to_reference_register(
            reference_register, fine_integral, fine_ba, fine_dm,
            fine_geom, false);
        reference_register.Reflux(reference_state, 0, 0, ncomp);
        amrex::Gpu::streamSynchronize();

        MultiFab direct_correction(coarse_ba, coarse_dm, 5, 0);
        for (int first = 0; first < ncomp; first += 5) {
            std::vector<int> component_map{
                first, first + 1, first + 2, first + 3, first + 4};
            direct_correction.setVal(Real(0.0));
            ASSERT_TRUE(direct_register.reflux_components(
                direct_correction, component_map, diagnostic)) << diagnostic;
            MultiFab::Add(direct_state, direct_correction, 0, first, 5, 0);
        }
        ASSERT_TRUE(direct_register.can_accept_reflux(diagnostic)) << diagnostic;
        expect_roundoff_equivalent(direct_state, reference_state);
        // Independent local finite-volume oracle at the uncovered coarse cell
        // immediately left of the x-high refined interface.  The expected
        // correction is assembled directly from the accepted, stage-weighted
        // face integrals and an explicit outward-normal sign; it does not
        // inspect the final candidate difference or reuse YA's sign formula.
        const IntVect interface_cell(1, 3, 3);
        const Box cell_region(interface_cell, interface_cell);
        const Box coarse_face(IntVect(2, 3, 3), IntVect(2, 3, 3),
                              IntVect::TheDimensionVector(0));
        const Real inv_dx = Real(1.0) / coarse_geom.CellSize(0);
        const ComponentOffsets semantic = offsets_for(layout);
        const std::array<int, 3> sample_components{
            semantic.two_mass, semantic.two_number, semantic.two_property};
        const int operation_count = 2 * stages * (fine_substeps + 1) + 12;
        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real q = static_cast<Real>(operation_count) * eps;
        const Real gamma = q / (Real(1.0) - q);
        for (const int component : sample_components) {
            const Real coarse_amount =
                coarse_integral.dir(0).sum(coarse_face, component, false);
            Real fine_amount_sum = Real(0.0);
            for (int y_subface = 0; y_subface < 2; ++y_subface) {
                for (int z_subface = 0; z_subface < 2; ++z_subface) {
                    const IntVect fine_face_cell(4, 6 + y_subface,
                                                 6 + z_subface);
                    const Box fine_face(fine_face_cell, fine_face_cell,
                                        IntVect::TheDimensionVector(0));
                    fine_amount_sum +=
                        fine_integral.dir(0).sum(fine_face, component, false);
                }
            }
            const Real fine_face_average = fine_amount_sum / Real(4.0);
            const Real expected_delta_h =
                (coarse_amount - fine_face_average) * inv_dx;
            const Real initial_h = Real(0.25) * (component + 1);
            const Real actual_delta_h =
                direct_state.sum(cell_region, component, false) - initial_h;
            const Real magnitude =
                (std::abs(coarse_amount) +
                 (std::abs(fine_amount_sum) / Real(4.0))) * inv_dx;
            const Real oracle_tolerance =
                (gamma + Real(2.0) * eps) * magnitude;
            const Real wrong_sign_delta = -expected_delta_h;
            const Real omit_fine_delta = coarse_amount * inv_dx;
            const Real double_interval_delta =
                expected_delta_h * static_cast<Real>(coarse_end - coarse_begin);
            ASSERT_GT(std::abs(expected_delta_h), Real(100.0) * oracle_tolerance)
                << "fixture lacks a discriminating local interface amount; component="
                << component << " method=" << method_id;
            EXPECT_LE(std::abs(actual_delta_h - expected_delta_h),
                      oracle_tolerance)
                << "component=" << component << " method=" << method_id
                << " expected_delta_h=" << expected_delta_h
                << " actual_delta_h=" << actual_delta_h
                << " bound=" << oracle_tolerance;
            EXPECT_GT(std::abs(actual_delta_h - wrong_sign_delta),
                      oracle_tolerance)
                << "local FV oracle failed to distinguish reversed face sign";
            EXPECT_GT(std::abs(actual_delta_h - omit_fine_delta),
                      oracle_tolerance)
                << "local FV oracle failed to distinguish omitted fine integral";
            EXPECT_GT(std::abs(actual_delta_h - double_interval_delta),
                      oracle_tolerance)
                << "local FV oracle failed to distinguish double interval weighting";
            std::cout << "SBM_LOCAL_FV_ORACLE method=" << method_id
                      << " component=" << component
                      << " coarse_integral=" << coarse_amount
                      << " fine_subface_sum=" << fine_amount_sum
                      << " delta_h=" << expected_delta_h
                      << " actual_delta_h=" << actual_delta_h
                      << " bound=" << oracle_tolerance
                      << " operation_count=" << operation_count
                      << " wrong_sign_separation="
                      << std::abs(actual_delta_h - wrong_sign_delta)
                      << " omit_fine_separation="
                      << std::abs(actual_delta_h - omit_fine_delta)
                      << " double_interval_separation="
                      << std::abs(actual_delta_h - double_interval_delta)
                      << '\n';
        }
        EXPECT_TRUE(direct_register.has_outstanding_amount());
        ASSERT_TRUE(direct_register.accept_reflux(diagnostic)) << diagnostic;
        EXPECT_FALSE(direct_register.has_outstanding_amount());
        MultiFab consumed_snapshot(coarse_ba, coarse_dm, ncomp, 0);
        MultiFab::Copy(consumed_snapshot, direct_state, 0, 0, ncomp, 0);
        EXPECT_FALSE(direct_register.reflux(direct_state, diagnostic));
        expect_same_values(direct_state, consumed_snapshot);
    }
}

TEST(SBMAMRYA, InterleavedTwoMomentPropertiesMatchSignedOracleAndReduceRuns)
{
    ASSERT_EQ(AMREX_SPACEDIM, 3);
    const auto layout = make_transfer_layout();
    const int ncomp = layout.ncomp();
    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(7));
    const Box fine_domain(IntVect(0), IntVect(15));
    const BoxArray coarse_ba(coarse_domain);
    const BoxArray fine_ba(amrex::refine(Box(IntVect(2), IntVect(5)), ratio));
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const Geometry coarse_geom = make_geometry(coarse_domain);
    const Geometry fine_geom = make_geometry(fine_domain);
    const auto chunks = erf_sbm::make_constraint_closure_chunks(layout, 16);
    ASSERT_EQ(chunks.size(), 1U);
    const std::vector<int> interleaved = chunks.front().components;
    ASSERT_EQ(interleaved.size(), static_cast<std::size_t>(ncomp));
    std::vector<int> sorted = interleaved;
    std::sort(sorted.begin(), sorted.end());
    ASSERT_NE(interleaved, sorted);
    ASSERT_GT(contiguous_component_runs(interleaved),
              contiguous_component_runs(sorted));
    EXPECT_EQ(contiguous_component_runs(interleaved), 10U);
    EXPECT_EQ(contiguous_component_runs(sorted), 1U);
    for (int component = 0; component < ncomp; ++component) {
        EXPECT_EQ(sorted[static_cast<std::size_t>(component)], component);
    }

    auto populate_accepted_interval = [&] (
        erf_sbm::SBMAMRFluxRegister& flux_register,
        erf_auxiliary::IntegratedMappedFaceFlux& coarse_integral,
        erf_auxiliary::IntegratedMappedFaceFlux& fine_integral,
        const std::vector<int>& component_map) {
        std::string diagnostic;
        ASSERT_TRUE(flux_register.begin_interval(0.0, 1.0, 1, diagnostic))
            << diagnostic;
        coarse_integral.define(coarse_ba, coarse_dm, ncomp, 0);
        coarse_integral.setVal(Real(0.0));
        fine_integral.define(fine_ba, fine_dm, ncomp, 0);
        fine_integral.setVal(Real(0.0));
        for (const bool coarse_side : {true, false}) {
            const auto& ba = coarse_side ? coarse_ba : fine_ba;
            const auto& dm = coarse_side ? coarse_dm : fine_dm;
            auto& integral = coarse_side ? coarse_integral : fine_integral;
            ASSERT_TRUE(flux_register.begin_level_step(
                coarse_side, 0.0, 1.0, diagnostic)) << diagnostic;
            for (int stage = 0; stage < 3; ++stage) {
                ASSERT_TRUE(flux_register.begin_stage(
                    coarse_side, erf_auxiliary::HostIntegrator::CompressibleRK3,
                    stage, 0.0, diagnostic)) << diagnostic;
                auto rate = make_signed_mapped_rate(
                    ba, dm, component_map, !coarse_side);
                const double completed_weight = stage == 2 ? 1.0 : 0.0;
                ASSERT_TRUE(flux_register.add_stage_rates(
                    coarse_side, rate, component_map, completed_weight,
                    diagnostic)) << diagnostic;
                accumulate_mapped_reference_integral(
                    integral, rate, component_map, completed_weight);
                ASSERT_TRUE(flux_register.finish_stage(coarse_side, diagnostic))
                    << diagnostic;
            }
            ASSERT_TRUE(flux_register.accept_level_step(
                coarse_side, 0.0, 1.0, diagnostic)) << diagnostic;
        }
        ASSERT_TRUE(flux_register.ready_to_reflux());
    };

    erf_sbm::SBMAMRFluxRegister interleaved_register;
    interleaved_register.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                                fine_geom, coarse_geom, ratio, 1, ncomp);
    erf_auxiliary::IntegratedMappedFaceFlux coarse_integral;
    erf_auxiliary::IntegratedMappedFaceFlux fine_integral;
    populate_accepted_interval(interleaved_register, coarse_integral,
                               fine_integral, interleaved);

    // Duplicate mapping is rejected before YA reads or writes the compact
    // correction. Omitting a component leaves the interval unconsumable.
    MultiFab duplicate_correction(coarse_ba, coarse_dm, ncomp, 0);
    duplicate_correction.setVal(Real(77.0));
    std::vector<int> duplicate_map = interleaved;
    duplicate_map[1] = duplicate_map[0];
    std::string diagnostic;
    EXPECT_FALSE(interleaved_register.reflux_components(
        duplicate_correction, duplicate_map, diagnostic));
    EXPECT_DOUBLE_EQ(duplicate_correction.min(0), Real(77.0));
    EXPECT_DOUBLE_EQ(duplicate_correction.max(ncomp - 1), Real(77.0));

    MultiFab marker(fine_ba, fine_dm, 1, 0);
    MultiFab covered(coarse_ba, coarse_dm, 1, 1);
    marker.setVal(Real(1.0));
    covered.setVal(Real(0.0));
    amrex::average_down(marker, covered, 0, 1, ratio);
    covered.FillBoundary(coarse_geom.periodicity());
    MultiFab measure(coarse_ba, coarse_dm, 1, 0);
    MultiFab initial(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab density(coarse_ba, coarse_dm, 1, 0);
    fill_carriers(density, measure, true);
    fill_spectrum(initial, density, layout);
    MultiFab interleaved_state(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab sorted_state(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab::Copy(interleaved_state, initial, 0, 0, ncomp, 0);
    MultiFab::Copy(sorted_state, initial, 0, 0, ncomp, 0);
    MultiFab correction(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab invalid(coarse_ba, coarse_dm, 1, 0);
    correction.setVal(Real(0.0));
    invalid.setVal(Real(0.0));
    ASSERT_TRUE(interleaved_register.reflux_components(
        correction, interleaved, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        interleaved_state, correction, measure, covered, interleaved,
        invalid, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(invalid, diagnostic))
        << diagnostic;
    ASSERT_TRUE(interleaved_register.can_accept_reflux(diagnostic)) << diagnostic;
    EXPECT_EQ(interleaved_register.reflux_call_count(),
              contiguous_component_runs(interleaved));
    EXPECT_EQ(interleaved_register.reflux_sync_count(),
              interleaved_register.reflux_call_count());

    // Fault injection: swap one property/global slot with another member while
    // reusing the same compact YA output. The independent reference must
    // distinguish this map error.
    std::vector<int> swapped_map = interleaved;
    std::swap(swapped_map[1], swapped_map[2]);
    MultiFab swapped_state(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab swapped_invalid(coarse_ba, coarse_dm, 1, 0);
    MultiFab::Copy(swapped_state, initial, 0, 0, ncomp, 0);
    swapped_invalid.setVal(Real(0.0));
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        swapped_state, correction, measure, covered, swapped_map,
        swapped_invalid, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(
        swapped_invalid, diagnostic)) << diagnostic;
    MultiFab swapped_difference(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab::Copy(swapped_difference, swapped_state, 0, 0, ncomp, 0);
    MultiFab::Subtract(swapped_difference, interleaved_state, 0, 0, ncomp, 0);
    Real swapped_mapping_error = Real(0.0);
    for (int component = 0; component < ncomp; ++component) {
        swapped_mapping_error = std::max(swapped_mapping_error,
                                         swapped_difference.norm0(component));
    }
    EXPECT_GT(swapped_mapping_error, Real(0.0));

    erf_sbm::SBMAMRFluxRegister sorted_register;
    sorted_register.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                           fine_geom, coarse_geom, ratio, 1, ncomp);
    erf_auxiliary::IntegratedMappedFaceFlux sorted_coarse_integral;
    erf_auxiliary::IntegratedMappedFaceFlux sorted_fine_integral;
    populate_accepted_interval(sorted_register, sorted_coarse_integral,
                               sorted_fine_integral, sorted);
    correction.setVal(Real(0.0));
    invalid.setVal(Real(0.0));
    ASSERT_TRUE(sorted_register.reflux_components(
        correction, sorted, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        sorted_state, correction, measure, covered, sorted, invalid,
        diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(invalid, diagnostic))
        << diagnostic;
    ASSERT_TRUE(sorted_register.can_accept_reflux(diagnostic)) << diagnostic;
    EXPECT_EQ(sorted_register.reflux_call_count(),
              contiguous_component_runs(sorted));
    EXPECT_EQ(sorted_register.reflux_sync_count(),
              sorted_register.reflux_call_count());
    EXPECT_LT(sorted_register.reflux_call_count(),
              interleaved_register.reflux_call_count());

    // Pinned YA consumes the same completed global integrated face fields.
    MultiFab reference_state(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab::Copy(reference_state, initial, 0, 0, ncomp, 0);
    amrex::YAFluxRegister reference_register(
        fine_ba, coarse_ba, fine_dm, coarse_dm, fine_geom, coarse_geom,
        ratio, 1, ncomp);
    reference_register.reset();
    add_integral_to_reference_register(reference_register, coarse_integral,
                                       coarse_ba, coarse_dm, coarse_geom, true);
    add_integral_to_reference_register(reference_register, fine_integral,
                                       fine_ba, fine_dm, fine_geom, false);
    MultiFab reference_correction(coarse_ba, coarse_dm, ncomp, 0);
    reference_correction.setVal(Real(0.0));
    reference_register.Reflux(reference_correction, 0, 0, ncomp);
    amrex::Gpu::streamSynchronize();
    MultiFab reference_invalid(coarse_ba, coarse_dm, 1, 0);
    reference_invalid.setVal(Real(0.0));
    std::vector<int> all_components(static_cast<std::size_t>(ncomp));
    for (int component = 0; component < ncomp; ++component) {
        all_components[static_cast<std::size_t>(component)] = component;
    }
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        reference_state, reference_correction, measure, covered,
        all_components, reference_invalid, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(
        reference_invalid, diagnostic)) << diagnostic;
    expect_roundoff_equivalent(interleaved_state, reference_state);
    expect_roundoff_equivalent(sorted_state, reference_state);

    // Independently sum signed coarse-minus-fine face amounts over each
    // uncovered coarse cell adjacent to the covered fine patch.
    const auto coarse_dx = coarse_geom.CellSizeArray();
    const int tangential[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    MultiFab fv_expected(coarse_ba, coarse_dm, ncomp, 0);
    MultiFab::Copy(fv_expected, initial, 0, 0, ncomp, 0);
    for (MFIter mfi(fv_expected, false); mfi.isValid(); ++mfi) {
        const Box box = mfi.validbox();
        const auto expected = fv_expected.array(mfi);
        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    const int cell[3] = {i, j, k};
                    const bool is_covered = i >= 2 && i <= 5 &&
                                            j >= 2 && j <= 5 &&
                                            k >= 2 && k <= 5;
                    const Real omega = coarse_measure(i, j, k);
                    for (int component = 0; component < ncomp; ++component) {
                        Real delta_h = Real(0.0);
                        if (!is_covered) {
                            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                                for (const int side : {-1, 1}) {
                                    const int neighbor = cell[dir] + side;
                                    const bool neighbor_covered = neighbor >= 2 &&
                                        neighbor <= 5 &&
                                        cell[(dir + 1) % 3] >= 2 &&
                                        cell[(dir + 1) % 3] <= 5 &&
                                        cell[(dir + 2) % 3] >= 2 &&
                                        cell[(dir + 2) % 3] <= 5;
                                    if (!neighbor_covered) continue;
                                    const int face = cell[dir] + (side > 0 ? 1 : 0);
                                    int coarse_face[3] = {i, j, k};
                                    coarse_face[dir] = face;
                                    const Real coarse_flux = signed_mapping_oracle_rate(
                                        dir, coarse_face[0], coarse_face[1],
                                        coarse_face[2], component, false);
                                    int fine_face[3] = {2*i, 2*j, 2*k};
                                    fine_face[dir] = 2 * face;
                                    Real fine_sum = Real(0.0);
                                    for (int a = 0; a < 2; ++a) {
                                        for (int b = 0; b < 2; ++b) {
                                            int subface[3] = {fine_face[0],
                                                              fine_face[1],
                                                              fine_face[2]};
                                            subface[tangential[dir][0]] += a;
                                            subface[tangential[dir][1]] += b;
                                            fine_sum += signed_mapping_oracle_rate(
                                                dir, subface[0], subface[1],
                                                subface[2], component, true);
                                        }
                                    }
                                    const Real sigma = side > 0
                                        ? Real(1.0) : Real(-1.0);
                                    delta_h += sigma / coarse_dx[dir] *
                                        (coarse_flux - fine_sum / Real(4.0));
                                }
                            }
                        }
                        expected(i,j,k,component) += delta_h / omega;
                    }
                }
            }
        }
    }
    expect_roundoff_equivalent(interleaved_state, fv_expected);
    expect_roundoff_equivalent(sorted_state, fv_expected);

    // A missing global component cannot complete or consume the pending YA
    // interval, even though its already formed compact correction is private.
    erf_sbm::SBMAMRFluxRegister omitted_register;
    omitted_register.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                            fine_geom, coarse_geom, ratio, 1, ncomp);
    erf_auxiliary::IntegratedMappedFaceFlux omitted_coarse;
    erf_auxiliary::IntegratedMappedFaceFlux omitted_fine;
    populate_accepted_interval(omitted_register, omitted_coarse,
                              omitted_fine, interleaved);
    std::vector<int> omitted_map = interleaved;
    omitted_map.pop_back();
    MultiFab omitted_correction(coarse_ba, coarse_dm, ncomp, 0);
    omitted_correction.setVal(Real(0.0));
    EXPECT_TRUE(omitted_register.reflux_components(
        omitted_correction, omitted_map, diagnostic)) << diagnostic;
    EXPECT_FALSE(omitted_register.can_accept_reflux(diagnostic));
    EXPECT_FALSE(omitted_register.accept_reflux(diagnostic));
    EXPECT_TRUE(omitted_register.has_outstanding_amount());
    ASSERT_TRUE(omitted_register.discard_interval(diagnostic)) << diagnostic;
    EXPECT_FALSE(omitted_register.has_outstanding_amount());
}

TEST(SBMAMRTransfer, MeasuresRealTransferOperatorAndDenseYAAllocations)
{
    ASSERT_EQ(AMREX_SPACEDIM, 3);
    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(15));
    const Box fine_domain(IntVect(0), IntVect(31));
    const Geometry coarse_geom = make_geometry(coarse_domain);
    const Geometry fine_geom = make_geometry(fine_domain);
    const BoxArray coarse_ba(coarse_domain);
    const BoxArray fine_ba(amrex::refine(Box(IntVect(4), IntVect(11)), ratio));
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    auto local_cells = [] (const BoxArray& boxes,
                           const DistributionMapping& mapping) {
        std::uint64_t count = 0;
        for (int box = 0; box < boxes.size(); ++box) {
            if (mapping[box] == amrex::ParallelDescriptor::MyProc()) {
                count += static_cast<std::uint64_t>(boxes[box].numPts());
            }
        }
        return count;
    };
    const std::uint64_t coarse_cells = local_cells(coarse_ba, coarse_dm);
    const std::uint64_t fine_cells = local_cells(fine_ba, fine_dm);
    constexpr std::uint64_t bytes_per_real = sizeof(Real);

    for (const int ncomp : {32, 128, 256}) {
        const auto layout = make_one_moment_layout(ncomp);
        MultiFab coarse_rho(coarse_ba, coarse_dm, 1, 0);
        MultiFab coarse_omega(coarse_ba, coarse_dm, 1, 0);
        MultiFab fine_rho(fine_ba, fine_dm, 1, 0);
        MultiFab fine_omega(fine_ba, fine_dm, 1, 0);
        MultiFab coarse_state(coarse_ba, coarse_dm, ncomp, 0);
        MultiFab fine_state(fine_ba, fine_dm, ncomp, 0);
        MultiFab coarse_candidate(coarse_ba, coarse_dm, ncomp, 0);
        MultiFab fine_candidate(fine_ba, fine_dm, ncomp, 0);
        coarse_rho.setVal(Real(1.0));
        coarse_omega.setVal(Real(1.0));
        fine_rho.setVal(Real(1.0));
        fine_omega.setVal(Real(1.0));
        fill_positive_one_moment_state(coarse_state);
        fill_positive_one_moment_state(fine_state);

        erf_sbm::SBMStateAdmission admission(layout);
        admission.define_level(0, coarse_ba, coarse_dm);
        admission.define_level(1, fine_ba, fine_dm);
        const auto chunks = erf_sbm::make_constraint_closure_chunks(layout, 16);
        int max_chunk_components = 0;
        for (const auto& chunk : chunks) {
            max_chunk_components = std::max(
                max_chunk_components, static_cast<int>(chunk.components.size()));
        }
        ASSERT_GT(max_chunk_components, 0);

        const std::uint64_t ya_before = amrex::TotalBytesAllocatedInFabs();
        erf_sbm::SBMAMRFluxRegister owner;
        owner.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                     fine_geom, coarse_geom, ratio, 1, ncomp);
        const std::uint64_t ya_total =
            static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs()) - ya_before;
        const std::uint64_t ya_dense_coarse =
            coarse_cells * static_cast<std::uint64_t>(ncomp) * bytes_per_real;
        const std::uint64_t caller_coarse_candidate =
            coarse_cells * static_cast<std::uint64_t>(ncomp) * bytes_per_real;
        BoxArray coarse_on_fine_ba = fine_ba;
        coarse_on_fine_ba.coarsen(ratio);
        const std::uint64_t fill_scratch_before =
            static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs());
        std::uint64_t fill_stage_scratch_bytes = 0;
        {
            // Match fill_coarse_parent_ghosts: two coarse spectral endpoint
            // chunks plus old/new density, with two coarse ghost cells.
            MultiFab fill_scratch(coarse_on_fine_ba, fine_dm,
                                  2 * max_chunk_components + 2,
                                  IntVect(2));
            amrex::Gpu::streamSynchronize();
            fill_stage_scratch_bytes =
                static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs()) -
                fill_scratch_before;
        }

        const auto fine_view = timed_view(fine_state, fine_rho, fine_omega, 0.5);
        const auto coarse_view = timed_view(coarse_state, coarse_rho, coarse_omega, 0.5);
        std::string diagnostic;
        erf_sbm::SBMTransferAllocationCensus restriction_census;
        ASSERT_TRUE(erf_sbm::RestrictMappedSpectrumPrivate(
            layout, fine_view, coarse_view, ratio, 0, coarse_candidate,
            diagnostic, 16, false, &admission, &restriction_census)) << diagnostic;
        const std::uint64_t restriction_peak =
            restriction_census.peak_extra_bytes();
        EXPECT_EQ(restriction_census.support_loss_collective_votes, 1U)
            << "restriction must make one global support-loss decision, independent of ncomp";
        const std::uint64_t restriction_bound =
            (2ULL * static_cast<std::uint64_t>(max_chunk_components) *
                 (fine_cells + coarse_cells) +
             2ULL * fine_cells + 3ULL * coarse_cells) * bytes_per_real;
        EXPECT_LE(restriction_peak, restriction_bound + 64ULL * 1024ULL)
            << "ncomp=" << ncomp << " K=" << max_chunk_components;

        erf_sbm::SBMTransferAllocationCensus prolongation_census;
        ASSERT_TRUE(erf_sbm::ProlongCarrierRelativeSpectrumPrivate(
            layout, coarse_view, timed_view(fine_state, fine_rho, fine_omega, 0.5),
            coarse_geom, fine_geom, ratio, 1, fine_candidate, diagnostic,
            16, false, &admission, &prolongation_census)) << diagnostic;
        const std::uint64_t prolongation_peak =
            prolongation_census.peak_extra_bytes();
        const std::uint64_t prolongation_bound = std::max<std::uint64_t>(
            3ULL * fine_cells + 5ULL * coarse_cells,
            static_cast<std::uint64_t>(max_chunk_components) *
                    (fine_cells + coarse_cells) + fine_cells + coarse_cells) *
            bytes_per_real;
        EXPECT_LE(prolongation_peak, prolongation_bound + 64ULL * 1024ULL)
            << "ncomp=" << ncomp << " K=" << max_chunk_components;
        EXPECT_EQ(admission.counters().host_scans, 0U);
        EXPECT_EQ(admission.counters().successful_admissions, 5U);

        const std::uint64_t full_fine_payload =
            fine_cells * static_cast<std::uint64_t>(ncomp) * bytes_per_real;
        const std::uint64_t before_mutant =
            static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs());
        std::uint64_t observed_full_fine_payload = 0;
        {
            MultiFab artificial_full_fine(fine_ba, fine_dm, ncomp, 0);
            amrex::Gpu::streamSynchronize();
            observed_full_fine_payload =
                static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs()) -
                before_mutant;
            EXPECT_GE(observed_full_fine_payload, full_fine_payload);
        }

        // Recreate the three full-component prolongation temporaries in the
        // audited baseline source (coarse_z, fine_z, and its private candidate)
        // to measure their actual FAB payload on this same layout.  This is a
        // source-shaped allocation control, not a run of the old algorithm.
        const auto before_legacy_prolongation =
            static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs());
        std::uint64_t legacy_prolongation_staging_bytes = 0;
        {
            MultiFab legacy_coarse_z(coarse_ba, coarse_dm, ncomp, 0);
            MultiFab legacy_fine_z(fine_ba, fine_dm, ncomp, 0);
            MultiFab legacy_fine_candidate(fine_ba, fine_dm, ncomp, 0);
            amrex::Gpu::streamSynchronize();
            legacy_prolongation_staging_bytes =
                static_cast<std::uint64_t>(amrex::TotalBytesAllocatedInFabs()) -
                before_legacy_prolongation;
        }
        const std::uint64_t expected_legacy_prolongation_staging =
            (coarse_cells + 2ULL * fine_cells) *
            static_cast<std::uint64_t>(ncomp) * bytes_per_real;
        EXPECT_GE(legacy_prolongation_staging_bytes,
                  expected_legacy_prolongation_staging);
        std::cout << "SBM_TRANSFER_CENSUS ncomp=" << ncomp
                  << " K=" << max_chunk_components
                  << " local_coarse_cells=" << coarse_cells
                  << " local_fine_cells=" << fine_cells
                  << " restriction_peak_scratch_bytes=" << restriction_peak
                  << " prolongation_peak_scratch_bytes=" << prolongation_peak
                  << " support_loss_collective_votes="
                  << restriction_census.support_loss_collective_votes
                  << " legacy_restriction_full_fine_staging_bytes="
                  << observed_full_fine_payload
                  << " legacy_prolongation_source_shaped_staging_bytes="
                  << legacy_prolongation_staging_bytes
                  << " caller_coarse_candidate_bytes=" << caller_coarse_candidate
                  << " temporal_fill_stage_scratch_bytes="
                  << fill_stage_scratch_bytes
                  << " temporal_fill_stage_scratch_components="
                  << 2 * max_chunk_components + 2
                  << " ya_dense_coarse_bytes=" << ya_dense_coarse
                  << " ya_total_interface_payload_bytes=" << ya_total
                  << " full_fine_spectrum_mutant_bytes=" << full_fine_payload
                  << " host_scans=" << admission.counters().host_scans
                  << " success_host_d2h_bytes=" << admission.counters().host_d2h_bytes
                  << " admission_mpi_reductions=" << admission.counters().mpi_reductions
                  << '\n';
    }
}

TEST(SBMAMRTransfer, MeasuresRealProductionTransportAndLedgerAllocations)
{
    ASSERT_EQ(AMREX_SPACEDIM, 3);
    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(15));
    const Box fine_domain = amrex::refine(coarse_domain, ratio);
    const Geometry coarse_geom = make_geometry(coarse_domain);
    const Geometry fine_geom = make_geometry(fine_domain);
    BoxArray coarse_ba(coarse_domain);
    BoxArray fine_ba(amrex::refine(Box(IntVect(4), IntVect(11)), ratio));
    coarse_ba.maxSize(8);
    fine_ba.maxSize(16);
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const int local_rank = amrex::ParallelDescriptor::MyProc();
    auto local_cells = [local_rank] (const BoxArray& boxes,
                                     const DistributionMapping& dm) {
        std::uint64_t cells = 0;
        for (int box = 0; box < boxes.size(); ++box) {
            if (dm[box] == local_rank) {
                cells += static_cast<std::uint64_t>(boxes[box].numPts());
            }
        }
        return cells;
    };
    const std::uint64_t coarse_cells = local_cells(coarse_ba, coarse_dm);
    const std::uint64_t fine_cells = local_cells(fine_ba, fine_dm);
    constexpr std::uint64_t bytes_per_real = sizeof(Real);

    for (const int ncomp : {32, 128, 256}) {
        const auto layout = make_one_moment_layout(ncomp);
        const auto groups = erf_sbm::make_constraint_groups(layout);
        const auto chunks = erf_sbm::make_constraint_closure_chunks(layout, 16);
        int max_chunk_components = 0;
        int cell_ratio_components = 0;
        for (const auto& chunk : chunks) {
            max_chunk_components = std::max(
                max_chunk_components, static_cast<int>(chunk.components.size()));
            int chunk_ratios = 0;
            for (const int group_index : chunk.group_indices) {
                const auto& group = groups[static_cast<std::size_t>(group_index)];
                // Matches SBMTransport's ratio_count: linear constraints,
                // two sides for every property support, plus any 2M reserve.
                chunk_ratios += static_cast<int>(group.constraints.size()) +
                    2 * static_cast<int>(group.attached_property_indices.size()) +
                    (group.moment_mode == erf_sbm::MomentMode::TwoMoment &&
                             group.bin < layout.populations()[static_cast<std::size_t>(
                                 group.population_id)].grid.nbins() - 1
                         ? 1 : 0);
            }
            cell_ratio_components = std::max(cell_ratio_components,
                                             chunk_ratios);
        }
        ASSERT_EQ(max_chunk_components, 16);
        ASSERT_GT(cell_ratio_components, 0);

        std::uint64_t state_bytes = 0;
        std::uint64_t production_transport_bytes = 0;
        {
            const auto state_before = amrex::TotalBytesAllocatedInFabs();
            erf_sbm::SBMStateManager states(layout, 2);
            states.define(0, coarse_ba, coarse_dm, 0.0);
            states.define(1, fine_ba, fine_dm, 0.0);
            state_bytes = static_cast<std::uint64_t>(
                amrex::TotalBytesAllocatedInFabs() - state_before);
            EXPECT_EQ(state_bytes,
                2ULL * (coarse_cells + fine_cells) *
                    static_cast<std::uint64_t>(ncomp) * bytes_per_real);

            const auto transport_before = amrex::TotalBytesAllocatedInFabs();
            erf_sbm::SBMTransport transport(layout, 2, 16, false);
            transport.define(0, coarse_ba, coarse_dm);
            transport.define(1, fine_ba, fine_dm);
            amrex::Gpu::streamSynchronize();
            production_transport_bytes = static_cast<std::uint64_t>(
                amrex::TotalBytesAllocatedInFabs() - transport_before);
        }
        const std::uint64_t production_state_transport_bytes =
            state_bytes + production_transport_bytes;

        std::uint64_t reference_state_transport_bytes = 0;
        {
            const auto reference_before = amrex::TotalBytesAllocatedInFabs();
            erf_sbm::SBMStateManager states(layout, 2);
            states.define(0, coarse_ba, coarse_dm, 0.0);
            states.define(1, fine_ba, fine_dm, 0.0);
            erf_sbm::SBMTransport transport(layout, 2, 16, true);
            transport.define(0, coarse_ba, coarse_dm);
            transport.define(1, fine_ba, fine_dm);
            amrex::Gpu::streamSynchronize();
            reference_state_transport_bytes = static_cast<std::uint64_t>(
                amrex::TotalBytesAllocatedInFabs() - reference_before);
        }
        ASSERT_GE(reference_state_transport_bytes,
                  production_state_transport_bytes);
        const std::uint64_t completed_ledger_extra_bytes =
            reference_state_transport_bytes - production_state_transport_bytes;
        EXPECT_GT(completed_ledger_extra_bytes, 0U);

        std::uint64_t ghost_bytes = 0;
        const auto accumulate_transport_ghosts = [&] (
            const BoxArray& boxes, const DistributionMapping& dm) {
            for (int box = 0; box < boxes.size(); ++box) {
                if (dm[box] != local_rank) continue;
                Box grow2 = boxes[box];
                grow2.grow(2);
                Box grow1 = boxes[box];
                grow1.grow(1);
                const auto valid = static_cast<std::uint64_t>(
                    boxes[box].numPts());
                const auto two_ghost = static_cast<std::uint64_t>(
                    grow2.numPts()) - valid;
                const auto one_ghost = static_cast<std::uint64_t>(
                    grow1.numPts()) - valid;
                ghost_bytes += 4ULL * two_ghost *
                    static_cast<std::uint64_t>(ncomp) * bytes_per_real;
                ghost_bytes += one_ghost *
                    static_cast<std::uint64_t>(cell_ratio_components) *
                    bytes_per_real;
            }
        };
        accumulate_transport_ghosts(coarse_ba, coarse_dm);
        accumulate_transport_ghosts(fine_ba, fine_dm);

        const auto ya_before = amrex::TotalBytesAllocatedInFabs();
        std::uint64_t ya_bytes = 0;
        {
            erf_sbm::SBMAMRFluxRegister owner;
            owner.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                         fine_geom, coarse_geom, ratio, 1, ncomp);
            amrex::Gpu::streamSynchronize();
            ya_bytes = static_cast<std::uint64_t>(
                amrex::TotalBytesAllocatedInFabs() - ya_before);
        }
        const auto ya_dense_coarse = coarse_cells *
            static_cast<std::uint64_t>(ncomp) * bytes_per_real;

        for (int rank = 0;
             rank < amrex::ParallelDescriptor::NProcs(); ++rank) {
            if (local_rank == rank) {
                std::cout << "SBM_REAL_OBJECT_CENSUS ranks="
                          << amrex::ParallelDescriptor::NProcs()
                          << " rank=" << local_rank
                          << " ncomp=" << ncomp
                          << " K=" << max_chunk_components
                          << " cell_ratio_components=" << cell_ratio_components
                          << " coarse_cells_local=" << coarse_cells
                          << " fine_cells_local=" << fine_cells
                          << " old_new_state_bytes=" << state_bytes
                          << " production_transport_bytes="
                          << production_transport_bytes
                          << " transport_ghost_region_bytes=" << ghost_bytes
                          << " YA_dense_coarse_bytes=" << ya_dense_coarse
                          << " YA_total_bytes=" << ya_bytes
                          << " completed_integral_ledger_extra_bytes="
                          << completed_ledger_extra_bytes
                          << std::endl;
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
}

TEST(SBMAMRYA, ReportsDirectAndReferenceAllocationCensus)
{
    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(31));
    const Box fine_domain(IntVect(0), IntVect(63));
    BoxArray coarse_ba(coarse_domain);
    coarse_ba.maxSize(16);
    BoxArray fine_ba(amrex::refine(Box(IntVect(8), IntVect(23)), ratio));
    fine_ba.maxSize(16);
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const Geometry coarse_geom = make_geometry(coarse_domain);
    const Geometry fine_geom = make_geometry(fine_domain);
    constexpr int chunk_components = 16;

    auto payload_bytes = [] (const MultiFab& field) {
        std::uint64_t bytes = 0;
        for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
            const auto& fab = field[mfi];
            bytes += static_cast<std::uint64_t>(fab.box().numPts()) *
                     static_cast<std::uint64_t>(fab.nComp()) * sizeof(Real);
        }
        return bytes;
    };
    auto valid_payload_bytes = [] (const MultiFab& field) {
        std::uint64_t bytes = 0;
        for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
            bytes += static_cast<std::uint64_t>(mfi.validbox().numPts()) *
                     static_cast<std::uint64_t>(field.nComp()) * sizeof(Real);
        }
        return bytes;
    };

    for (const int ncomp : {32, 128, 256}) {
        const auto layout = make_one_moment_layout(ncomp);
        const auto groups = erf_sbm::make_constraint_groups(layout);
        const auto chunks = erf_sbm::make_constraint_closure_chunks(layout, 16);
        int cell_ratio_components = 0;
        int max_chunk_components = 0;
        for (const auto& chunk : chunks) {
            max_chunk_components = std::max(
                max_chunk_components, static_cast<int>(chunk.components.size()));
            int chunk_ratios = 0;
            for (const int group_index : chunk.group_indices) {
                const auto& group = groups[static_cast<std::size_t>(group_index)];
                chunk_ratios += static_cast<int>(group.constraints.size()) +
                    2 * static_cast<int>(group.attached_property_indices.size()) +
                    (group.moment_mode == erf_sbm::MomentMode::TwoMoment &&
                             group.bin < layout.populations()[static_cast<std::size_t>(
                                 group.population_id)].grid.nbins() - 1
                         ? 1 : 0);
            }
            cell_ratio_components = std::max(cell_ratio_components,
                                             chunk_ratios);
        }
        ASSERT_EQ(max_chunk_components, chunk_components);
        ASSERT_GT(cell_ratio_components, 0);

        std::uint64_t local_coarse_cells = 0;
        std::uint64_t local_fine_cells = 0;
        for (int box = 0; box < coarse_ba.size(); ++box) {
            if (coarse_dm[box] == amrex::ParallelDescriptor::MyProc()) {
                local_coarse_cells +=
                    static_cast<std::uint64_t>(coarse_ba[box].numPts());
            }
        }
        for (int box = 0; box < fine_ba.size(); ++box) {
            if (fine_dm[box] == amrex::ParallelDescriptor::MyProc()) {
                local_fine_cells +=
                    static_cast<std::uint64_t>(fine_ba[box].numPts());
            }
        }
        const std::uint64_t dense_coarse_bytes = local_coarse_cells *
            static_cast<std::uint64_t>(ncomp) * sizeof(Real);
        const auto before = amrex::TotalBytesAllocatedInFabs();
        std::uint64_t direct_allocated = 0;
        std::uint64_t reference_extra = 0;
        {
            erf_sbm::SBMAMRFluxRegister owner;
            owner.define(fine_ba, coarse_ba, fine_dm, coarse_dm,
                         fine_geom, coarse_geom, ratio, 1, ncomp);
            erf_auxiliary::MappedFaceFluxRate low, high, accepted, lambda,
                projected;
            low.define(fine_ba, fine_dm, chunk_components, 0);
            high.define(fine_ba, fine_dm, chunk_components, 0);
            accepted.define(fine_ba, fine_dm, chunk_components, 0);
            lambda.define(fine_ba, fine_dm, 1, 0);
            projected.define(fine_ba, fine_dm, 2, 0);
            low.setVal(Real(0.0));
            high.setVal(Real(0.0));
            accepted.setVal(Real(0.0));
            lambda.setVal(Real(0.0));
            projected.setVal(Real(0.0));
            erf_auxiliary::IntegratedMappedFaceFlux projected_ledger;
            projected_ledger.define(fine_ba, fine_dm, 2, 0);
            projected_ledger.setVal(Real(0.0));
            erf_auxiliary::MappedFaceFluxRate coarse_low, coarse_high,
                coarse_accepted, coarse_lambda, coarse_projected;
            coarse_low.define(coarse_ba, coarse_dm, chunk_components, 0);
            coarse_high.define(coarse_ba, coarse_dm, chunk_components, 0);
            coarse_accepted.define(coarse_ba, coarse_dm, chunk_components, 0);
            coarse_lambda.define(coarse_ba, coarse_dm, 1, 0);
            coarse_projected.define(coarse_ba, coarse_dm, 2, 0);
            coarse_low.setVal(Real(0.0));
            coarse_high.setVal(Real(0.0));
            coarse_accepted.setVal(Real(0.0));
            coarse_lambda.setVal(Real(0.0));
            coarse_projected.setVal(Real(0.0));
            erf_auxiliary::IntegratedMappedFaceFlux coarse_projected_ledger;
            coarse_projected_ledger.define(coarse_ba, coarse_dm, 2, 0);
            coarse_projected_ledger.setVal(Real(0.0));
            amrex::Gpu::streamSynchronize();
            direct_allocated = static_cast<std::uint64_t>(
                amrex::TotalBytesAllocatedInFabs() - before);

            // Count the existing M3 transport storage on both active levels:
            // anchor/target, four full-component WENO input views with their
            // two-cell ghost regions, and the old/new authoritative views.
            const auto persistent_before = amrex::TotalBytesAllocatedInFabs();
            std::vector<std::unique_ptr<MultiFab>> persistent_fields;
            auto add_persistent_field = [&] (const BoxArray& ba,
                                             const DistributionMapping& dm,
                                             const int components,
                                             const int ghosts) {
                auto field = std::make_unique<MultiFab>(
                    ba, dm, components, ghosts);
                field->setVal(Real(0.0));
                persistent_fields.push_back(std::move(field));
            };
            for (const auto& level : std::array{
                     std::pair<const BoxArray*, const DistributionMapping*>{
                         &coarse_ba, &coarse_dm},
                     std::pair<const BoxArray*, const DistributionMapping*>{
                         &fine_ba, &fine_dm}}) {
                for (int state_view = 0; state_view < 2; ++state_view) {
                    add_persistent_field(*level.first, *level.second,
                                         ncomp, 0);
                }
                for (int transport_view = 0; transport_view < 2;
                     ++transport_view) {
                    add_persistent_field(*level.first, *level.second,
                                         ncomp, 0);
                }
                for (int ghost_view = 0; ghost_view < 4; ++ghost_view) {
                    add_persistent_field(*level.first, *level.second,
                                         ncomp, 2);
                }
                add_persistent_field(*level.first, *level.second,
                                     chunk_components, 0); // low_trial_h
                add_persistent_field(*level.first, *level.second,
                                     cell_ratio_components, 1);
                for (int scalar_view = 0; scalar_view < 3; ++scalar_view) {
                    add_persistent_field(*level.first, *level.second, 1, 0);
                }
            }
            amrex::Gpu::streamSynchronize();
            const auto state_and_transport_bytes = static_cast<std::uint64_t>(
                amrex::TotalBytesAllocatedInFabs() - persistent_before);
            const auto persistent_payload = [&] () {
                std::uint64_t bytes = 0;
                for (const auto& field : persistent_fields) {
                    bytes += payload_bytes(*field);
                }
                return bytes;
            }();
            std::uint64_t ghost_region_bytes = 0;
            for (std::size_t field = 0; field < persistent_fields.size(); ++field) {
                if (field % 13 >= 4 && field % 13 < 8) {
                    ghost_region_bytes += payload_bytes(*persistent_fields[field]) -
                                          valid_payload_bytes(*persistent_fields[field]);
                }
            }
            EXPECT_EQ(state_and_transport_bytes, persistent_payload);

            BoxArray coarse_on_fine_ba = fine_ba;
            coarse_on_fine_ba.coarsen(ratio);
            const auto fill_scratch_before = amrex::TotalBytesAllocatedInFabs();
            std::uint64_t temporal_fill_stage_scratch_bytes = 0;
            std::uint64_t parent_fill_source_pack_bytes = 0;
            {
                MultiFab source_pack(coarse_ba, coarse_dm,
                                     2 * chunk_components + 2, 0);
                MultiFab fill_scratch(coarse_on_fine_ba, fine_dm,
                                      2 * chunk_components + 2,
                                      amrex::IntVect(2));
                amrex::Gpu::streamSynchronize();
                for (amrex::MFIter mfi(source_pack, false);
                     mfi.isValid(); ++mfi) {
                    parent_fill_source_pack_bytes +=
                        static_cast<std::uint64_t>(
                            source_pack[mfi].nBytesOwned());
                }
                temporal_fill_stage_scratch_bytes = payload_bytes(fill_scratch);
                EXPECT_EQ(static_cast<std::uint64_t>(
                              amrex::TotalBytesAllocatedInFabs() -
                              fill_scratch_before),
                          temporal_fill_stage_scratch_bytes +
                              parent_fill_source_pack_bytes);
            }
            const auto direct_componentwise_stage_peak_bytes =
                direct_allocated + state_and_transport_bytes +
                temporal_fill_stage_scratch_bytes;
            const auto packed_comparison_stage_peak_bytes =
                direct_componentwise_stage_peak_bytes +
                parent_fill_source_pack_bytes;

            // Keep the one private full-spectrum coarse candidate live while
            // measuring the two mutually exclusive sync phases. Restriction
            // uses only complete-group K-sized H/support chunks; reflux uses
            // one K-sized correction. Neither phase allocates a second full
            // fine-spectrum staging field.
            const auto sync_before = amrex::TotalBytesAllocatedInFabs();
            std::uint64_t sync_scratch_bytes = 0;
            {
                MultiFab coarse_candidate(coarse_ba, coarse_dm, ncomp, 0);
                coarse_candidate.setVal(Real(0.0));
                std::uint64_t restriction_phase_bytes = 0;
                {
                    MultiFab fine_h_chunk(fine_ba, fine_dm,
                                          chunk_components, 0);
                    MultiFab fine_support_chunk(fine_ba, fine_dm,
                                                chunk_components, 0);
                    MultiFab coarse_h_chunk(coarse_ba, coarse_dm,
                                            chunk_components, 0);
                    MultiFab coarse_support_chunk(coarse_ba, coarse_dm,
                                                 chunk_components, 0);
                    MultiFab fine_coverage(fine_ba, fine_dm, 1, 0);
                    MultiFab coarse_coverage(coarse_ba, coarse_dm, 1, 0);
                    MultiFab lost_support(coarse_ba, coarse_dm, 1, 0);
                    MultiFab fine_invalid(fine_ba, fine_dm, 1, 0);
                    MultiFab coarse_invalid(coarse_ba, coarse_dm, 1, 0);
                    fine_h_chunk.setVal(Real(0.0));
                    fine_support_chunk.setVal(Real(0.0));
                    coarse_h_chunk.setVal(Real(0.0));
                    coarse_support_chunk.setVal(Real(0.0));
                    fine_coverage.setVal(Real(0.0));
                    coarse_coverage.setVal(Real(0.0));
                    lost_support.setVal(Real(0.0));
                    fine_invalid.setVal(Real(0.0));
                    coarse_invalid.setVal(Real(0.0));
                    amrex::Gpu::streamSynchronize();
                    restriction_phase_bytes = static_cast<std::uint64_t>(
                        amrex::TotalBytesAllocatedInFabs() - sync_before);
                }

                std::uint64_t reflux_phase_bytes = 0;
                {
                    MultiFab fine_coverage(fine_ba, fine_dm, 1, 0);
                    MultiFab covered(coarse_ba, coarse_dm, 1, 1);
                    MultiFab correction(coarse_ba, coarse_dm,
                                        chunk_components, 0);
                    MultiFab correction_invalid(coarse_ba, coarse_dm, 1, 0);
                    fine_coverage.setVal(Real(1.0));
                    covered.setVal(Real(0.0));
                    correction.setVal(Real(0.0));
                    correction_invalid.setVal(Real(0.0));
                    amrex::Gpu::streamSynchronize();
                    reflux_phase_bytes = static_cast<std::uint64_t>(
                        amrex::TotalBytesAllocatedInFabs() - sync_before);
                }
                sync_scratch_bytes = std::max(restriction_phase_bytes,
                                              reflux_phase_bytes);
            }

            {
                const auto ledger_before = amrex::TotalBytesAllocatedInFabs();
                erf_auxiliary::IntegratedMappedFaceFlux reference_ledger;
                reference_ledger.define(fine_ba, fine_dm, ncomp, 0);
                reference_ledger.setVal(Real(0.0));
                erf_auxiliary::IntegratedMappedFaceFlux coarse_reference_ledger;
                coarse_reference_ledger.define(coarse_ba, coarse_dm, ncomp, 0);
                coarse_reference_ledger.setVal(Real(0.0));
                amrex::Gpu::streamSynchronize();
                reference_extra = static_cast<std::uint64_t>(
                    amrex::TotalBytesAllocatedInFabs() - ledger_before);
                EXPECT_EQ(payload_bytes(reference_ledger.dir(0)) +
                              payload_bytes(reference_ledger.dir(1)) +
                              payload_bytes(reference_ledger.dir(2)) +
                              payload_bytes(coarse_reference_ledger.dir(0)) +
                              payload_bytes(coarse_reference_ledger.dir(1)) +
                              payload_bytes(coarse_reference_ledger.dir(2)),
                          reference_extra);
            }

            const auto direct_peak_bytes = direct_allocated +
                state_and_transport_bytes + sync_scratch_bytes;
            for (int rank = 0;
                 rank < amrex::ParallelDescriptor::NProcs(); ++rank) {
                if (amrex::ParallelDescriptor::MyProc() == rank) {
                    std::cout << "SBM_ALLOC_CENSUS ranks="
                              << amrex::ParallelDescriptor::NProcs()
                              << " ncomp=" << ncomp
                              << " rank=" << amrex::ParallelDescriptor::MyProc()
                              << " coarse_cells_local="
                              << local_coarse_cells
                              << " fine_cells_local=" << local_fine_cells
                              << " YA_dense_coarse_bytes=" << dense_coarse_bytes
                              << " YA_plus_face_chunks_bytes=" << direct_allocated
                              << " old_new_transport_and_ghost_bytes="
                              << state_and_transport_bytes
                              << " temporal_fill_stage_scratch_bytes="
                              << temporal_fill_stage_scratch_bytes
                              << " parent_fill_source_pack_bytes="
                              << parent_fill_source_pack_bytes
                              << " packed_comparison_stage_peak_total_bytes="
                              << packed_comparison_stage_peak_bytes
                              << " ghost_region_bytes=" << ghost_region_bytes
                              << " sync_candidate_peak_bytes=" << sync_scratch_bytes
                              << " direct_sync_peak_total_bytes=" << direct_peak_bytes
                              << " direct_componentwise_stage_peak_total_bytes="
                              << direct_componentwise_stage_peak_bytes
                              << " completed_integral_ledger_bytes=" << reference_extra
                              << " rate_chunk_components=" << chunk_components
                              << " cell_ratio_components="
                              << cell_ratio_components
                              << std::endl;
                }
                amrex::ParallelDescriptor::Barrier();
            }
        }
    }
}

// Gate 4B regrid fixture: the fine level moves one parent to the right.
// Parents i=1 are released, i=2 retained, and i=3 newly refined.
constexpr double regrid_time = 0.5;

Box regrid_old_fine_box () { return Box(IntVect(2, 2, 2), IntVect(5, 5, 5)); }
Box regrid_new_fine_box () { return Box(IntVect(4, 2, 2), IntVect(7, 5, 5)); }
Box regrid_retained_box () { return Box(IntVect(4, 2, 2), IntVect(5, 5, 5)); }
Box regrid_prolonged_box () { return Box(IntVect(6, 2, 2), IntVect(7, 5, 5)); }

// Every child's mapped carrier equals its parent's: omega_c rho_c = rho_P.
void fill_regrid_fine_carriers (MultiFab& density, MultiFab& measure)
{
    for (MFIter mfi(density, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto rho = density.array(mfi);
        const auto omega = measure.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const int parity = (i % 2 + j % 2 + k % 2) % 2;
            const Real sign = parity == 0 ? Real(-1.0) : Real(1.0);
            rho(i, j, k, 0) = prolong_fine_density(i, j, k);
            omega(i, j, k, 0) = Real(1.0) + Real(0.2) * sign;
        });
    }
}

struct RegridFixture {
    erf_sbm::SBMLayout layout = make_transfer_layout();
    Geometry coarse_geom = make_geometry(coarse_domain());
    Geometry fine_geom = make_geometry(amrex::refine(coarse_domain(), IntVect(2)));
    MultiFab coarse_rho, coarse_omega, coarse_state, unsynchronized_coarse_state;
    MultiFab old_rho, old_omega, old_state;
    MultiFab new_rho, new_omega;

    RegridFixture (const BoxArray& old_ba, const DistributionMapping& old_dm,
                   const BoxArray& new_ba, const DistributionMapping& new_dm)
    {
        const BoxArray coarse_ba(coarse_domain());
        const DistributionMapping coarse_dm(coarse_ba);
        const int ncomp = layout.ncomp();
        coarse_rho.define(coarse_ba, coarse_dm, 1, 0);
        coarse_omega.define(coarse_ba, coarse_dm, 1, 0);
        coarse_state.define(coarse_ba, coarse_dm, ncomp, 0);
        unsynchronized_coarse_state.define(coarse_ba, coarse_dm, ncomp, 0);
        fill_carriers(coarse_rho, coarse_omega, true);
        coarse_omega.setVal(Real(1.0));
        fill_spectrum(unsynchronized_coarse_state, coarse_rho, layout);

        old_rho.define(old_ba, old_dm, 1, 0);
        old_omega.define(old_ba, old_dm, 1, 0);
        old_state.define(old_ba, old_dm, ncomp, 0);
        fill_regrid_fine_carriers(old_rho, old_omega);
        // Fine values use fine indices, so they differ from any parent prolongation.
        fill_spectrum(old_state, old_rho, layout);

        new_rho.define(new_ba, new_dm, 1, 0);
        new_omega.define(new_ba, new_dm, 1, 0);
        fill_regrid_fine_carriers(new_rho, new_omega);

        // The last synchronization restricted the old fine level onto its parents.
        std::string diagnostic;
        if (!erf_sbm::RestrictMappedSpectrum(
                layout, timed_view(old_state, old_rho, old_omega, regrid_time),
                timed_view(unsynchronized_coarse_state, coarse_rho, coarse_omega,
                           regrid_time),
                IntVect(2), 0, coarse_state, diagnostic)) {
            ADD_FAILURE() << diagnostic;
        }
    }

    [[nodiscard]] erf_sbm::SBMAMRStateView coarse_view () const
    {
        return timed_view(coarse_state, coarse_rho, coarse_omega, regrid_time);
    }

    bool remap (const MultiFab& retained, MultiFab& candidate,
                std::string& diagnostic, const double fine_time = regrid_time,
                erf_sbm::SBMTransferAllocationCensus* census = nullptr) const
    {
        return erf_sbm::RemapRegriddedSpectrumPrivate(
            layout, coarse_view(), retained, new_rho, 0, new_omega, 0, fine_time,
            coarse_geom, fine_geom, IntVect(2), 1, candidate, diagnostic, 16,
            nullptr, nullptr, census);
    }
};

// Copy the valid cells of field inside region into a single-box MultiFab.
MultiFab restrict_to_box (const MultiFab& field, const Box& region)
{
    const BoxArray region_ba(region);
    MultiFab result(region_ba, DistributionMapping(region_ba), field.nComp(), 0);
    result.setVal(Real(0.0));
    result.ParallelCopy(field, 0, 0, field.nComp());
    return result;
}

// Composite mapped inventory; a fine child carries 1/8 of a parent volume.
std::vector<Real> composite_mapped_inventory (const MultiFab& coarse_state,
                                              const MultiFab& coarse_omega,
                                              const MultiFab& fine_state,
                                              const MultiFab& fine_omega)
{
    MultiFab fine_ones(fine_state.boxArray(), fine_state.DistributionMap(), 1, 0);
    MultiFab covered(coarse_state.boxArray(), coarse_state.DistributionMap(), 1, 0);
    fine_ones.setVal(Real(1.0));
    covered.setVal(Real(0.0));
    amrex::average_down(fine_ones, covered, 0, 1, IntVect(2));
    const int ncomp = coarse_state.nComp();
    MultiFab coarse_amount(coarse_state.boxArray(), coarse_state.DistributionMap(),
                           ncomp, 0);
    MultiFab fine_amount(fine_state.boxArray(), fine_state.DistributionMap(),
                         ncomp, 0);
    for (MFIter mfi(coarse_amount, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto u = coarse_state.const_array(mfi);
        const auto omega = coarse_omega.const_array(mfi);
        const auto mask = covered.const_array(mfi);
        const auto h = coarse_amount.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                h(i, j, k, n) = (Real(1.0) - mask(i, j, k, 0)) *
                                omega(i, j, k, 0) * u(i, j, k, n);
            });
    }
    for (MFIter mfi(fine_amount, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto u = fine_state.const_array(mfi);
        const auto omega = fine_omega.const_array(mfi);
        const auto h = fine_amount.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                h(i, j, k, n) = omega(i, j, k, 0) * u(i, j, k, n) / Real(8.0);
            });
    }
    std::vector<Real> inventory(static_cast<std::size_t>(ncomp));
    for (int n = 0; n < ncomp; ++n) {
        inventory[static_cast<std::size_t>(n)] =
            coarse_amount.sum(n) + fine_amount.sum(n);
    }
    return inventory;
}

// Each new fine cell's expected value from its parent's dry-air-relative
// state, and the parent amount U copied unchanged (the negative control).
void fill_regrid_parent_expectations (MultiFab& ratio_expected,
                                      MultiFab& copied_parent_u,
                                      const MultiFab& fine_rho,
                                      const ComponentOffsets& offsets)
{
    const int ncomp = ratio_expected.nComp();
    for (MFIter mfi(ratio_expected, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto rho = fine_rho.const_array(mfi);
        const auto z_parent = ratio_expected.array(mfi);
        const auto u_parent = copied_parent_u.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            for (int n = 0; n < ncomp; ++n) {
                const Real z = base_component(n, i / 2, j / 2, k / 2, offsets, false);
                z_parent(i, j, k, n) = rho(i, j, k, 0) * z;
                u_parent(i, j, k, n) = coarse_density(i / 2, j / 2, k / 2) * z;
            }
        });
    }
}

TEST(SBMAMRRegrid, CellCountsClassifyRetainedProlongedAndReleasedCells)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    const auto counts = erf_sbm::CountRegridCells(old_ba, new_ba);
    EXPECT_EQ(counts.cells, 64);
    EXPECT_EQ(counts.retained, regrid_retained_box().numPts());
    EXPECT_EQ(counts.prolonged, regrid_prolonged_box().numPts());
    EXPECT_EQ(counts.released, 32);
    const auto created = erf_sbm::CountRegridCells(BoxArray(), new_ba);
    EXPECT_EQ(created.retained, 0);
    EXPECT_EQ(created.prolonged, 64);
}

TEST(SBMAMRRegrid, RetainedCellsAreBitwiseAndNewCellsCarryTheParentRatio)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    const RegridFixture f(old_ba, DistributionMapping(old_ba),
                          new_ba, DistributionMapping(new_ba));
    const int ncomp = f.layout.ncomp();
    MultiFab candidate(new_ba, DistributionMapping(new_ba), ncomp, 0);
    candidate.setVal(Real(77.0));
    std::string diagnostic;
    ASSERT_TRUE(f.remap(f.old_state, candidate, diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
        candidate, f.layout, 1, &diagnostic)) << diagnostic;

    // Retained cells keep the accepted fine values bit for bit.
    expect_same_values(restrict_to_box(candidate, regrid_retained_box()),
                       restrict_to_box(f.old_state, regrid_retained_box()));

    // New cells carry their never-refined parent's dry-air-relative state.
    const ComponentOffsets offsets = offsets_for(f.layout);
    MultiFab expected(new_ba, DistributionMapping(new_ba), ncomp, 0);
    MultiFab copied_parent_u(new_ba, DistributionMapping(new_ba), ncomp, 0);
    fill_regrid_parent_expectations(expected, copied_parent_u, f.new_rho, offsets);
    const MultiFab prolonged = restrict_to_box(candidate, regrid_prolonged_box());
    const MultiFab ratio_expected = restrict_to_box(expected, regrid_prolonged_box());
    const MultiFab copy_expected = restrict_to_box(copied_parent_u, regrid_prolonged_box());
    for (int n = 0; n < ncomp; ++n) {
        MultiFab error(prolonged.boxArray(), prolonged.DistributionMap(), 1, 0);
        MultiFab::Copy(error, prolonged, n, 0, 1, 0);
        MultiFab::Subtract(error, ratio_expected, n, 0, 1, 0);
        const Real scale = std::max(Real(1.0), ratio_expected.norm0(n));
        EXPECT_LE(error.norm0(0), Real(64.0) * std::numeric_limits<Real>::epsilon() * scale)
            << "component=" << n;
        // Negative control: copying the parent amount U, not its ratio,
        // misses the children's 20% density variation.
        MultiFab::Copy(error, prolonged, n, 0, 1, 0);
        MultiFab::Subtract(error, copy_expected, n, 0, 1, 0);
        if (ratio_expected.norm0(n) > Real(0.0)) {
            EXPECT_GT(error.norm0(0), Real(0.1) * ratio_expected.norm0(n))
                << "component=" << n;
        }
    }
}

TEST(SBMAMRRegrid, CompositeInventoryIsConservedAcrossRetainedNewAndReleasedParents)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    const RegridFixture f(old_ba, DistributionMapping(old_ba),
                          new_ba, DistributionMapping(new_ba));
    MultiFab candidate(new_ba, DistributionMapping(new_ba), f.layout.ncomp(), 0);
    std::string diagnostic;
    ASSERT_TRUE(f.remap(f.old_state, candidate, diagnostic)) << diagnostic;

    const auto before = composite_mapped_inventory(f.coarse_state, f.coarse_omega,
                                                   f.old_state, f.old_omega);
    const auto after = composite_mapped_inventory(f.coarse_state, f.coarse_omega,
                                                  candidate, f.new_omega);
    // Negative control: without the synchronization that restricted the
    // released parents, the same regrid does not conserve the composite.
    const auto unsynchronized = composite_mapped_inventory(
        f.unsynchronized_coarse_state, f.coarse_omega, candidate, f.new_omega);
    for (std::size_t n = 0; n < before.size(); ++n) {
        const Real scale = std::max(Real(1.0), std::abs(before[n]));
        EXPECT_NEAR(after[n], before[n],
                    Real(256.0) * std::numeric_limits<Real>::epsilon() * scale)
            << "component=" << n;
        if (std::abs(before[n]) > Real(0.0)) {
            EXPECT_GT(std::abs(unsynchronized[n] - before[n]),
                      Real(1.0e-6) * std::abs(before[n]))
                << "component=" << n;
        }
    }
}

TEST(SBMAMRRegrid, RejectionLeavesEveryInputUnchanged)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    RegridFixture f(old_ba, DistributionMapping(old_ba),
                    new_ba, DistributionMapping(new_ba));
    const int ncomp = f.layout.ncomp();
    const auto snapshot = [] (const MultiFab& field) {
        MultiFab copy(field.boxArray(), field.DistributionMap(), field.nComp(), 0);
        MultiFab::Copy(copy, field, 0, 0, field.nComp(), 0);
        return copy;
    };
    const MultiFab coarse_before = snapshot(f.coarse_state);
    const MultiFab old_before = snapshot(f.old_state);
    MultiFab candidate(new_ba, DistributionMapping(new_ba), ncomp, 0);
    std::string diagnostic;

    EXPECT_FALSE(f.remap(f.old_state, candidate, diagnostic, 0.75));
    EXPECT_NE(diagnostic.find("same-time"), std::string::npos) << diagnostic;

    // A retained cell that is inadmissible must reject the whole candidate.
    MultiFab corrupt_retained = snapshot(f.old_state);
    set_cell_value(corrupt_retained, IntVect(4, 2, 2), offsets_for(f.layout).one_mass,
                   Real(-1.0));
    EXPECT_FALSE(f.remap(corrupt_retained, candidate, diagnostic));
    EXPECT_NE(diagnostic.find("inadmissible"), std::string::npos) << diagnostic;

    // A new-grid carrier that breaks the parent identity is rejected.
    f.new_rho.mult(Real(2.0), 0, 1, 0);
    EXPECT_FALSE(f.remap(f.old_state, candidate, diagnostic));
    EXPECT_NE(diagnostic.find("inventory differs"), std::string::npos) << diagnostic;
    f.new_rho.mult(Real(0.5), 0, 1, 0);

    expect_same_values(f.coarse_state, coarse_before);
    expect_same_values(f.old_state, old_before);
}

TEST(SBMAMRRegrid, RetainedCellsDoNotConstrainTheProlongation)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    const Box retained_box = regrid_retained_box();
    const Box prolonged_box = regrid_prolonged_box();
    // Breaks its parent's carrier identity and underflows the prolonged ratio.
    const Real tiny_carrier = std::numeric_limits<Real>::denorm_min();
    std::string diagnostic;

    RegridFixture retained_fault(old_ba, DistributionMapping(old_ba),
                                 new_ba, DistributionMapping(new_ba));
    set_cell_value(retained_fault.new_rho, retained_box.smallEnd(), 0,
                   tiny_carrier);
    MultiFab candidate(new_ba, DistributionMapping(new_ba),
                       retained_fault.layout.ncomp(), 0);
    candidate.setVal(Real(0.0));
    ASSERT_TRUE(retained_fault.remap(retained_fault.old_state, candidate, diagnostic))
        << diagnostic;
    expect_same_values(restrict_to_box(candidate, retained_box),
                       restrict_to_box(retained_fault.old_state, retained_box));

    RegridFixture prolonged_fault(old_ba, DistributionMapping(old_ba),
                                  new_ba, DistributionMapping(new_ba));
    set_cell_value(prolonged_fault.new_rho, prolonged_box.smallEnd(), 0,
                   tiny_carrier);
    EXPECT_FALSE(prolonged_fault.remap(prolonged_fault.old_state, candidate, diagnostic));
    EXPECT_FALSE(diagnostic.empty());
}

TEST(SBMAMRRegrid, CandidateIsDecompositionInvariant)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    const RegridFixture whole(old_ba, DistributionMapping(old_ba),
                              new_ba, DistributionMapping(new_ba));
    MultiFab reference(new_ba, DistributionMapping(new_ba), whole.layout.ncomp(), 0);
    std::string diagnostic;
    ASSERT_TRUE(whole.remap(whole.old_state, reference, diagnostic)) << diagnostic;

    BoxArray split_old(regrid_old_fine_box());
    split_old.maxSize(2);
    BoxArray split_new(regrid_new_fine_box());
    split_new.maxSize(IntVect(2, 4, 2));
    const RegridFixture split(split_old, DistributionMapping(split_old),
                              split_new, DistributionMapping(split_new));
    MultiFab split_candidate(split_new, DistributionMapping(split_new),
                             split.layout.ncomp(), 0);
    ASSERT_TRUE(split.remap(split.old_state, split_candidate, diagnostic)) << diagnostic;
    expect_same_values(restrict_to_box(split_candidate, regrid_new_fine_box()),
                       restrict_to_box(reference, regrid_new_fine_box()));
}

TEST(SBMAMRRegrid, ScratchStaysWithinTheChunkedProlongationModel)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    const RegridFixture f(old_ba, DistributionMapping(old_ba),
                          new_ba, DistributionMapping(new_ba));
    MultiFab candidate(new_ba, DistributionMapping(new_ba), f.layout.ncomp(), 0);
    erf_sbm::SBMTransferAllocationCensus census;
    std::string diagnostic;
    ASSERT_TRUE(f.remap(f.old_state, candidate, diagnostic, regrid_time, &census))
        << diagnostic;
    int max_chunk_components = 0;
    for (const auto& chunk : erf_sbm::make_constraint_closure_chunks(f.layout, 16)) {
        max_chunk_components = std::max(max_chunk_components,
                                        static_cast<int>(chunk.components.size()));
    }
    const auto fine_cells = static_cast<std::uint64_t>(new_ba.numPts());
    const auto coarse_cells = static_cast<std::uint64_t>(coarse_domain().numPts());
    // Never a second full-spectrum fine buffer: the candidate is the only
    // one.  The prolonged-cell mask and its parent mask add one fine and one
    // coarse component.
    const std::uint64_t bound = (std::max<std::uint64_t>(
        3ULL * fine_cells + 5ULL * coarse_cells,
        static_cast<std::uint64_t>(max_chunk_components) *
            (fine_cells + coarse_cells) + fine_cells + coarse_cells) +
        fine_cells + coarse_cells) * sizeof(Real);
    EXPECT_LE(census.peak_extra_bytes(), bound + 64ULL * 1024ULL);
    EXPECT_FALSE(census.processed_chunk_first_components.empty());
}

TEST(SBMAMRRegrid, ReplacedLevelInvalidatesOldViewsAndMeasure)
{
    const BoxArray old_ba(regrid_old_fine_box());
    const BoxArray new_ba(regrid_new_fine_box());
    const DistributionMapping old_dm(old_ba);
    const DistributionMapping new_dm(new_ba);
    const RegridFixture f(old_ba, old_dm, new_ba, new_dm);
    const int ncomp = f.layout.ncomp();

    erf_sbm::SBMStateManager manager(f.layout, 2);
    manager.define(1, old_ba, old_dm, regrid_time);
    MultiFab::Copy(manager.new_state_for_initialization(1), f.old_state, 0, 0, ncomp, 0);
    std::string diagnostic;
    ASSERT_TRUE(manager.begin_step(1, regrid_time, diagnostic)) << diagnostic;
    MultiFab candidate(new_ba, new_dm, ncomp, 0);
    ASSERT_TRUE(f.remap(f.old_state, candidate, diagnostic)) << diagnostic;
    // A level inside a physical step cannot be replaced.
    EXPECT_FALSE(manager.can_replace_level(1, candidate, regrid_time, diagnostic));
    ASSERT_TRUE(manager.accept_stage_target(1, regrid_time, true, diagnostic)) << diagnostic;
    ASSERT_TRUE(manager.old_valid(1));

    EXPECT_FALSE(manager.can_replace_level(1, candidate, 0.75, diagnostic));
    MultiFab wrong_components(new_ba, new_dm, ncomp - 1, 0);
    EXPECT_FALSE(manager.can_replace_level(1, wrong_components, regrid_time, diagnostic));
    ASSERT_TRUE(manager.can_replace_level(1, candidate, regrid_time, diagnostic))
        << diagnostic;
    MultiFab expected(new_ba, new_dm, ncomp, 0);
    MultiFab::Copy(expected, candidate, 0, 0, ncomp, 0);
    manager.replace_level(1, std::move(candidate), regrid_time);
    EXPECT_EQ(manager.new_state(1).boxArray(), new_ba);
    EXPECT_FALSE(manager.old_valid(1));
    EXPECT_THROW(static_cast<void>(manager.old_state(1)), std::logic_error);
    EXPECT_EQ(manager.new_time(1), regrid_time);
    expect_same_values(manager.new_state(1), expected);
    ASSERT_TRUE(manager.begin_step(1, regrid_time, diagnostic)) << diagnostic;
    EXPECT_EQ(manager.old_state(1).boxArray(), new_ba);

    // A recreated transport level has no measure until it is rebuilt on the new grids.
    erf_sbm::SBMTransport transport(f.layout, 2);
    const auto planar = [] (const BoxArray& cells) {
        amrex::BoxList boxes = cells.boxList();
        for (auto& box : boxes) { box.setRange(2, 0); }
        return BoxArray(std::move(boxes));
    };
    const auto build_measure = [&] (const BoxArray& cells, const DistributionMapping& dm) {
        MultiFab detj(cells, dm, 1, 0);
        MultiFab mx(planar(cells), dm, 1, 0);
        MultiFab my(planar(cells), dm, 1, 0);
        detj.setVal(Real(1.0));
        mx.setVal(Real(1.0));
        my.setVal(Real(1.0));
        std::string measure_diagnostic;
        EXPECT_TRUE(transport.rebuild_static_measure(1, detj, mx, my,
                                                     measure_diagnostic))
            << measure_diagnostic;
    };
    transport.define(1, old_ba, old_dm);
    build_measure(old_ba, old_dm);
    ASSERT_TRUE(transport.measure_is_ready(1));
    transport.destroy(1);
    EXPECT_FALSE(transport.is_defined(1));
    transport.define(1, new_ba, new_dm);
    EXPECT_FALSE(transport.measure_is_ready(1));
    EXPECT_THROW(static_cast<void>(transport.static_measure(1)), std::logic_error);
    build_measure(new_ba, new_dm);
    ASSERT_TRUE(transport.measure_is_ready(1));
    EXPECT_EQ(transport.static_measure(1).boxArray(), new_ba);
}

} // namespace
