#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_YAFluxRegister.H>

#include "ERF_AdvectionSrcForScalars.H"
#include "ERF_SBMAdvectionBoundary.H"
#include "ERF_SBMAMRFluxRegister.H"
#include "ERF_SBMAMRTransfer.H"
#include "ERF_IndexDefines.H"
#include "ERF_SBMRemapping.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMStageOwnership.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
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

constexpr int transport_nx = 8;

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

#if !defined(AMREX_USE_FLOAT)
struct A4FixturePoint
{
    Real rho;
    Real mass0;
    Real number0;
    Real mass1;
    Real number1;
    Real property0;
    Real property1;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE A4FixturePoint
a4_fixture_point (const Real x, const Real y, const Real z, const Real time)
{
    A4FixturePoint point{};
    point.rho = Real(1.0) + Real(0.07) * x + Real(0.025) * y +
        Real(0.015) * z + Real(0.12) * time + Real(0.018) * x * time;
    point.number0 = Real(0.012) + Real(0.0011) * x + Real(0.0003) * y +
        Real(0.0002) * z + Real(0.0007) * time + Real(0.00015) * x * time;
    const Real mean0 = Real(0.25) + Real(0.012) * x + Real(0.004) * y +
        Real(0.002) * z + Real(0.014) * time + Real(0.003) * x * time;
    point.mass0 = point.number0 * mean0;
    point.number1 = Real(0.008) + Real(0.0004) * x + Real(0.0002) * y +
        Real(0.0001) * z + Real(0.0003) * time + Real(0.00008) * x * time;
    const Real mean1 = Real(0.70) + Real(0.009) * x + Real(0.003) * y +
        Real(0.001) * z + Real(0.008) * time + Real(0.002) * x * time;
    point.mass1 = point.number1 * mean1;
    point.property0 = point.number0 *
        (Real(0.18) + Real(0.006) * x + Real(0.002) * y + Real(0.005) * time);
    point.property1 = point.number1 *
        (Real(0.32) + Real(0.004) * x + Real(0.001) * y + Real(0.003) * time);
    return point;
}

Real a4_weno_z3_left (const Real sm2, const Real sm1, const Real s)
{
    const Real b1 = (sm1 - sm2) * (sm1 - sm2);
    const Real b2 = (s - sm1) * (s - sm1);
    const Real tau = std::abs(b2 - b1);
    const Real eps = sizeof(Real) == sizeof(float)
        ? Real(1.0e-12) : Real(1.0e-40);
    const Real w1 = Real(1.0 / 3.0) *
        (Real(1.0) + tau * tau / ((eps + b1) * (eps + b1)));
    const Real w2 = Real(2.0 / 3.0) *
        (Real(1.0) + tau * tau / ((eps + b2) * (eps + b2)));
    return (w1 * (-sm2 + Real(3.0) * sm1) +
            w2 * (sm1 + s)) / (Real(2.0) * (w1 + w2));
}

void fill_a4_fixture_state (MultiFab& state, MultiFab& density,
                            const bool fine, const Real time,
                            const int mass0, const int number0,
                            const int property0)
{
    const int scale = fine ? 2 : 1;
    for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto values = state.array(mfi);
        const auto rho = density.array(mfi);
        const int m0 = mass0, n0 = number0;
        const int m1 = mass0 + 1, n1 = number0 + 1, p0 = property0;
        amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const Real x = static_cast<Real>(i) / scale;
                const Real y = static_cast<Real>(j) / scale;
                const Real z = static_cast<Real>(k) / scale;
                const A4FixturePoint q = a4_fixture_point(x, y, z, time);
                const Real state_rho = fine ? Real(1.0) : q.rho;
                rho(i, j, k, Rho_comp) = state_rho;
                values(i, j, k, m0) = state_rho * q.mass0;
                values(i, j, k, n0) = state_rho * q.number0;
                values(i, j, k, m1) = state_rho * q.mass1;
                values(i, j, k, n1) = state_rho * q.number1;
                values(i, j, k, p0) = state_rho * q.property0;
                values(i, j, k, p0 + 1) = state_rho * q.property1;
            });
    }
}

void set_a4_x_face_velocity (MultiFab& avg_xmom, const IntVect& face,
                             const Real velocity)
{
    for (amrex::MFIter mfi(avg_xmom, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto u = avg_xmom.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j,
                                                      int k) noexcept {
            if ((i == face[0] || i == face[0] + 1) &&
                j == face[1] && k == face[2])
                u(i, j, k, 0) = velocity;
        });
    }
}
#endif

BoxArray
project_to_xy (const BoxArray& cell_ba)
{
    amrex::BoxList boxes = cell_ba.boxList();
    for (auto& box : boxes)
        box.setRange(2, 0);
    return BoxArray(std::move(boxes));
}

struct RunOptions
{
    erf_sbm::MomentMode mode{erf_sbm::MomentMode::OneMoment};
    bool varying_density{false};
    bool discontinuity{false};
    bool attached_property{false};
    Real carrier_x{Real(0.0)};
    bool noncanonical_shared_edge{false};
    int max_groups_per_chunk{16};
    double dt{0.01};
    bool smooth_profile{false};
    bool complex_profile{false};
    bool adversarial_one_moment_profile{false};
    bool canonical_upper_edge_profile{false};
    int canonical_upper_edge_seed{0};
    bool multidirectional_carrier{false};
    bool compare_native_candidate{false};
    bool compare_direct_moment_candidate{false};
    bool compare_historical_p2{false};
    bool expect_limiter_active{false};
    bool anelastic_heun{false};
    bool mapped_geometry{false};
    bool advance_all_rk3_stages{false};
    bool verify_spectral_ledger_identity{false};
    bool verify_direct_register_reference{false};
    int completed_steps{1};
    bool amplify_corrector_input{false};
    bool distinct_density_roles{false};
    Real rho_anchor_slope{Real(0.0)};
    Real rho_input_slope{Real(0.0)};
    Real rho_target_slope{Real(0.0)};
};

struct RunSummary
{
    std::vector<Real> inventory;
    std::vector<Real> weighted_inventory;
    std::vector<Real> squared_inventory;
    Real minimum_interior_upper_gap_fraction{std::numeric_limits<Real>::max()};
    Real minimum_interior_number{std::numeric_limits<Real>::max()};
    Real native_candidate_interior_upper_gap_fraction{
        std::numeric_limits<Real>::max()};
    std::vector<Real> spectral_ledger_relative_errors;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
canonical_upper_edge_probe_low_endpoint (const int index,
                                         const int seed) noexcept
{
    const unsigned int hash =
        (static_cast<unsigned int>(index + 1) * 1664525u) +
        (static_cast<unsigned int>(seed + 1) * 1013904223u);
    const int level = static_cast<int>((hash >> 16) & 7u);
    int multiplier = 1;
    switch (level) {
    case 0: multiplier = 1; break;
    case 1: multiplier = 2; break;
    case 2: multiplier = 4; break;
    case 3: multiplier = 8; break;
    case 4: multiplier = 16; break;
    case 5: multiplier = 32; break;
    case 6: multiplier = 64; break;
    default: multiplier = 128; break;
    }
    const Real eta =
        Real(128.0) * std::numeric_limits<Real>::epsilon();
    return Real(3.0) * eta * Real(0.01) * static_cast<Real>(multiplier);
}

std::unique_ptr<MultiFab>
make_native_candidate (const erf_sbm::SBMLayout& layout,
                       const MultiFab& spectrum,
                       const MultiFab& conserved,
                       const MultiFab& avg_xmom,
                       const MultiFab& avg_ymom,
                       const MultiFab& avg_zmom,
                       const MultiFab& measure,
                       const Geometry& geom,
                       const double interval,
                       const bool direct_physical_moment_weno = false,
                       const MultiFab* anchor_state = nullptr)
{
    MultiFab intensive(spectrum.boxArray(), spectrum.DistributionMap(),
                       layout.ncomp(), 2);
    for (amrex::MFIter mfi(intensive); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = spectrum.const_array(mfi);
        const auto rho = conserved.const_array(mfi);
        const auto z = intensive.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                z(i, j, k, n) = state(i, j, k, n) / rho(i, j, k, Rho_comp);
            });
    }
    intensive.FillBoundary(geom.periodicity());

    erf_auxiliary::MappedFaceFluxRate high_rate;
    high_rate.define(spectrum.boxArray(), spectrum.DistributionMap(),
                     layout.ncomp(), 0);
    for (const auto& population : layout.populations()) {
        if (population.moment_mode != erf_sbm::MomentMode::TwoMoment)
            continue;
        if (direct_physical_moment_weno)
            continue;
        for (int bin = 0; bin < population.grid.nbins(); ++bin) {
            const int mass = population.mass_offset + bin;
            const int number = population.number_offset + bin;
            const Real lower =
                population.grid.edges()[static_cast<std::size_t>(bin)];
            const Real upper =
                population.grid.edges()[static_cast<std::size_t>(bin + 1)];
            for (amrex::MFIter mfi(intensive); mfi.isValid(); ++mfi) {
                const Box bx = mfi.validbox();
                const auto z = intensive.array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const Real mass_value = z(i, j, k, mass);
                        const Real number_value = z(i, j, k, number);
                        const Real width = upper - lower;
                        z(i, j, k, mass) =
                            (upper * number_value - mass_value) / width;
                        z(i, j, k, number) =
                            (mass_value - lower * number_value) / width;
                    });
            }
        }
    }
    intensive.FillBoundary(geom.periodicity());

    for (amrex::MFIter mfi(intensive); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const amrex::GpuArray<const amrex::Array4<Real>, AMREX_SPACEDIM>
            flux_views{{high_rate.dir(0).array(mfi),
                        high_rate.dir(1).array(mfi),
                        high_rate.dir(2).array(mfi)}};
        for (int component = 0; component < layout.ncomp(); ++component) {
            BuildScalarAdvectionFluxes(
                bx, intensive.const_array(mfi), component, flux_views,
                component, avg_xmom.const_array(mfi), avg_ymom.const_array(mfi),
                avg_zmom.const_array(mfi), AdvType::Weno_3Z, AdvType::Weno_3Z,
                Real(0.0), Real(0.0));
        }
    }

    for (const auto& population : layout.populations()) {
        if (population.moment_mode != erf_sbm::MomentMode::TwoMoment)
            continue;
        if (direct_physical_moment_weno)
            continue;
        for (int bin = 0; bin < population.grid.nbins(); ++bin) {
            const int mass = population.mass_offset + bin;
            const int number = population.number_offset + bin;
            const Real lower =
                population.grid.edges()[static_cast<std::size_t>(bin)];
            const Real upper =
                population.grid.edges()[static_cast<std::size_t>(bin + 1)];
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                auto& face = high_rate.dir(dir);
                for (amrex::MFIter mfi(face); mfi.isValid(); ++mfi) {
                    const Box bx = mfi.validbox();
                    const auto flux = face.array(mfi);
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real blo = flux(i, j, k, mass);
                            const Real bhi = flux(i, j, k, number);
                            flux(i, j, k, mass) = lower * blo + upper * bhi;
                            flux(i, j, k, number) = blo + bhi;
                        });
                }
            }
        }
    }

    auto candidate = std::make_unique<MultiFab>(
        spectrum.boxArray(), spectrum.DistributionMap(), layout.ncomp(), 0);
    const auto inv_dx = geom.InvCellSizeArray();
    const Real tau = static_cast<Real>(interval);
    const MultiFab& anchor = anchor_state != nullptr ? *anchor_state : spectrum;
    for (amrex::MFIter mfi(*candidate); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto initial = anchor.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto fx = high_rate.dir(0).const_array(mfi);
        const auto fy = high_rate.dir(1).const_array(mfi);
        const auto fz = high_rate.dir(2).const_array(mfi);
        const auto out = candidate->array(mfi);
        const Real dx = inv_dx[0], dy = inv_dx[1], dz = inv_dx[2];
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                const Real divergence =
                    (fx(i + 1, j, k, n) - fx(i, j, k, n)) * dx +
                    (fy(i, j + 1, k, n) - fy(i, j, k, n)) * dy +
                    (fz(i, j, k + 1, n) - fz(i, j, k, n)) * dz;
                out(i, j, k, n) = (omega(i, j, k, 0) * initial(i, j, k, n) -
                                   tau * divergence) /
                                  omega(i, j, k, 0);
            });
    }
    return candidate;
}

// Frozen test-only reconstruction copied from the P2 Cartesian implementation
// at 7f5742d9ac1edb3e6eca1be4f8c0f2d357b0781a. This deliberately does not call
// ERF's production WENO helper: it is an independent oracle for retained
// periodic Cartesian endpoint-number behavior.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
historical_p2_weno_z3_face (const Real qm2, const Real qm1, const Real q,
                            const Real qp1, const Real carrier) noexcept
{
    const bool positive = carrier >= Real(0.0);
    const Real q0 = positive ? Real(0.5) * (-qm2 + Real(3.0) * qm1)
                             : Real(0.5) * (Real(3.0) * q - qp1);
    const Real q1 = Real(0.5) * (qm1 + q);
    const Real beta0 = positive ? (qm1 - qm2) * (qm1 - qm2)
                               : (qp1 - q) * (qp1 - q);
    const Real beta1 = (q - qm1) * (q - qm1);
#ifdef AMREX_USE_FLOAT
    constexpr Real epsilon = Real(1.0e-12);
#else
    constexpr Real epsilon = Real(1.0e-40);
#endif
    const Real tau = amrex::Math::abs(beta1 - beta0);
    const Real w0 = (Real(1.0) / Real(3.0)) *
                    (Real(1.0) + tau * tau /
                                     ((epsilon + beta0) * (epsilon + beta0)));
    const Real w1 = (Real(2.0) / Real(3.0)) *
                    (Real(1.0) + tau * tau /
                                     ((epsilon + beta1) * (epsilon + beta1)));
    return (w0 * q0 + w1 * q1) / (w0 + w1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
historical_p2_weno_face (const amrex::Array4<const Real>& values,
                         const int i, const int j, const int k,
                         const int component, const int dir,
                         const Real carrier) noexcept
{
    if (dir == 0) {
        return historical_p2_weno_z3_face(values(i - 2, j, k, component),
                                           values(i - 1, j, k, component),
                                           values(i, j, k, component),
                                           values(i + 1, j, k, component),
                                           carrier);
    }
    if (dir == 1) {
        return historical_p2_weno_z3_face(values(i, j - 2, k, component),
                                           values(i, j - 1, k, component),
                                           values(i, j, k, component),
                                           values(i, j + 1, k, component),
                                           carrier);
    }
    return historical_p2_weno_z3_face(values(i, j, k - 2, component),
                                       values(i, j, k - 1, component),
                                       values(i, j, k, component),
                                       values(i, j, k + 1, component),
                                       carrier);
}

std::unique_ptr<MultiFab>
make_historical_p2_cartesian_candidate (
    const erf_sbm::SBMLayout& layout, const MultiFab& spectrum,
    const MultiFab& conserved, const MultiFab& avg_xmom,
    const MultiFab& avg_ymom, const MultiFab& avg_zmom,
    const MultiFab& measure, const Geometry& geom, const double interval)
{
    MultiFab intensive(spectrum.boxArray(), spectrum.DistributionMap(),
                       layout.ncomp(), 2);
    for (amrex::MFIter mfi(intensive); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = spectrum.const_array(mfi);
        const auto rho = conserved.const_array(mfi);
        const auto z = intensive.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                z(i, j, k, n) = state(i, j, k, n) / rho(i, j, k, Rho_comp);
            });
    }
    intensive.FillBoundary(geom.periodicity());

    for (const auto& population : layout.populations()) {
        if (population.moment_mode != erf_sbm::MomentMode::TwoMoment) continue;
        for (int bin = 0; bin < population.grid.nbins(); ++bin) {
            const int mass = population.mass_offset + bin;
            const int number = population.number_offset + bin;
            const Real lower = population.grid.edges()[static_cast<std::size_t>(bin)];
            const Real upper = population.grid.edges()[static_cast<std::size_t>(bin + 1)];
            for (amrex::MFIter mfi(intensive); mfi.isValid(); ++mfi) {
                const Box bx = mfi.validbox();
                const auto z = intensive.array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const Real mass_value = z(i, j, k, mass);
                        const Real number_value = z(i, j, k, number);
                        const Real width = upper - lower;
                        z(i, j, k, mass) =
                            (upper * number_value - mass_value) / width;
                        z(i, j, k, number) =
                            (mass_value - lower * number_value) / width;
                    });
            }
        }
    }
    intensive.FillBoundary(geom.periodicity());

    erf_auxiliary::MappedFaceFluxRate high_rate;
    high_rate.define(spectrum.boxArray(), spectrum.DistributionMap(),
                     layout.ncomp(), 0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto& face = high_rate.dir(dir);
        const MultiFab& carrier = dir == 0 ? avg_xmom :
                                  (dir == 1 ? avg_ymom : avg_zmom);
        for (amrex::MFIter mfi(face, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            const Box bx = mfi.tilebox();
            const auto q = intensive.const_array(mfi);
            const auto mass_rate = carrier.const_array(mfi);
            const auto output = face.array(mfi);
            const int ncomp = layout.ncomp();
            amrex::ParallelFor(
                bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    output(i, j, k, n) =
                        mass_rate(i, j, k) * historical_p2_weno_face(
                                                 q, i, j, k, n, dir,
                                                 mass_rate(i, j, k));
                });
        }
    }
    for (const auto& population : layout.populations()) {
        if (population.moment_mode != erf_sbm::MomentMode::TwoMoment) continue;
        for (int bin = 0; bin < population.grid.nbins(); ++bin) {
            const int mass = population.mass_offset + bin;
            const int number = population.number_offset + bin;
            const Real lower = population.grid.edges()[static_cast<std::size_t>(bin)];
            const Real upper = population.grid.edges()[static_cast<std::size_t>(bin + 1)];
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                auto& face = high_rate.dir(dir);
                for (amrex::MFIter mfi(face); mfi.isValid(); ++mfi) {
                    const Box bx = mfi.validbox();
                    const auto flux = face.array(mfi);
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real blo = flux(i, j, k, mass);
                            const Real bhi = flux(i, j, k, number);
                            flux(i, j, k, mass) = lower * blo + upper * bhi;
                            flux(i, j, k, number) = blo + bhi;
                        });
                }
            }
        }
    }

    auto candidate = std::make_unique<MultiFab>(
        spectrum.boxArray(), spectrum.DistributionMap(), layout.ncomp(), 0);
    const auto inv_dx = geom.InvCellSizeArray();
    const Real tau = static_cast<Real>(interval);
    for (amrex::MFIter mfi(*candidate); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto initial = spectrum.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto fx = high_rate.dir(0).const_array(mfi);
        const auto fy = high_rate.dir(1).const_array(mfi);
        const auto fz = high_rate.dir(2).const_array(mfi);
        const auto out = candidate->array(mfi);
        const Real dx = inv_dx[0], dy = inv_dx[1], dz = inv_dx[2];
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                const Real divergence =
                    (fx(i + 1, j, k, n) - fx(i, j, k, n)) * dx +
                    (fy(i, j + 1, k, n) - fy(i, j, k, n)) * dy +
                    (fz(i, j, k + 1, n) - fz(i, j, k, n)) * dz;
                out(i, j, k, n) = (omega(i, j, k, 0) * initial(i, j, k, n) -
                                   tau * divergence) / omega(i, j, k, 0);
            });
    }
    return candidate;
}

Real
mapped_inventory (const MultiFab& state,
                  const MultiFab& measure,
                  const int component)
{
    MultiFab product(state.boxArray(), state.DistributionMap(), 1, 0);
    for (amrex::MFIter mfi(product); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto value = state.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto out = product.array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                out(i, j, k, 0) = omega(i, j, k, 0) * value(i, j, k, component);
            });
    }
    return product.sum(0);
}

Real
minimum_normalized_upper_gap (const MultiFab& state,
                              const int mass_component,
                              const int number_component,
                              const Real upper_edge)
{
    MultiFab gap(state.boxArray(), state.DistributionMap(), 1, 0);
    for (amrex::MFIter mfi(gap); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto values = state.const_array(mfi);
        const auto out = gap.array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const Real mass = values(i, j, k, mass_component);
                const Real number = values(i, j, k, number_component);
                out(i, j, k, 0) =
                    (upper_edge * number - mass) / (upper_edge * number + mass);
            });
    }
    return gap.min(0);
}

RunSummary
run_transport (const RunOptions& options)
{
    constexpr int nx = transport_nx;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const Box domain(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    const Geometry geom(domain, &physical, amrex::CoordSys::cartesian,
                        periodic);
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const auto layout = make_layout(options.mode, options.attached_property);

    erf_sbm::SBMStateManager state_manager(layout, 1);
    state_manager.define(0, ba, dm, 0.0);
    erf_sbm::SBMTransport transport(state_manager.layout(),
                                    (options.verify_spectral_ledger_identity ||
                                     options.verify_direct_register_reference) ? 2 : 1,
                                    options.max_groups_per_chunk);
    transport.define(0, ba, dm);

    std::unique_ptr<erf_sbm::SBMAMRFluxRegister> direct_register;
    BoxArray fine_interface_ba;
    DistributionMapping fine_interface_dm;
    std::unique_ptr<Geometry> fine_interface_geom;
    const IntVect interface_ratio(2);
    if (options.verify_direct_register_reference) {
        const Box fine_domain = amrex::refine(domain, interface_ratio);
        const Box covered_coarse(
            IntVect(nx / 4, 0, 0), IntVect(3 * nx / 4 - 1, ny - 1, nz - 1));
        fine_interface_ba = BoxArray(amrex::refine(covered_coarse,
                                                   interface_ratio));
        fine_interface_dm = DistributionMapping(fine_interface_ba);
        const amrex::RealBox fine_physical({0.0, 0.0, 0.0},
                                           {1.0, 1.0, 1.0});
        const int fine_periodic[AMREX_SPACEDIM] = {1, 1, 1};
        fine_interface_geom = std::make_unique<Geometry>(
            fine_domain, &fine_physical, amrex::CoordSys::cartesian,
            fine_periodic);
        direct_register =
            std::make_unique<erf_sbm::SBMAMRFluxRegister>();
        direct_register->define(fine_interface_ba, ba, fine_interface_dm, dm,
                                *fine_interface_geom, geom, interface_ratio,
                                1, layout.ncomp());
        transport.set_amr_flux_registers(0, direct_register.get(), nullptr);
    }

    MultiFab detj(ba, dm, 1, 0);
    const auto map_ba = project_to_xy(ba);
    MultiFab mx(map_ba, dm, 1, 0);
    MultiFab my(map_ba, dm, 1, 0);
    detj.setVal(Real(1.0));
    mx.setVal(Real(1.0));
    my.setVal(Real(1.0));
    if (options.mapped_geometry) {
        for (amrex::MFIter mfi(detj); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto det = detj.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                        int k) noexcept {
                // Synthetic static mapped measure; this is operator evidence,
                // not a claim that native terrain runtime was exercised.
                det(i, j, k, 0) = Real(1.0) + Real(0.02) * k + Real(0.005) * i;
            });
        }
        for (amrex::MFIter mfi(mx); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto x = mx.array(mfi);
            const auto y = my.array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
                    x(i, j, 0, 0) = Real(1.0) + Real(0.03) * i;
                    y(i, j, 0, 0) = Real(1.0) + Real(0.02) * j;
                });
        }
    }
    std::string diagnostic;
    if (!transport.rebuild_static_measure(0, detj, mx, my, diagnostic)) {
        ADD_FAILURE() << diagnostic;
        return {};
    }

    MultiFab conserved_anchor(ba, dm, 3, 0);
    MultiFab conserved_input(ba, dm, 3, 0);
    MultiFab conserved_target(ba, dm, 3, 0);
    conserved_anchor.setVal(Real(0.0));
    const Real default_density_slope =
        options.varying_density ? Real(0.1) : Real(0.0);
    const Real anchor_slope = options.distinct_density_roles
                                  ? options.rho_anchor_slope
                                  : default_density_slope;
    const Real input_slope = options.distinct_density_roles
                                 ? options.rho_input_slope
                                 : default_density_slope;
    const Real target_slope = options.distinct_density_roles
                                  ? options.rho_target_slope
                                  : default_density_slope;
    for (amrex::MFIter mfi(conserved_anchor); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto old = conserved_anchor.array(mfi);
        const bool separate_roles = options.distinct_density_roles;
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const Real coordinate = separate_roles
                                            ? static_cast<Real>(i) / Real(nx)
                                            : static_cast<Real>(j);
                old(i, j, k, Rho_comp) = Real(1.0) + anchor_slope * coordinate;
            });
    }
    MultiFab::Copy(conserved_input, conserved_anchor, 0, 0, 3, 0);
    MultiFab::Copy(conserved_target, conserved_anchor, 0, 0, 3, 0);
    if (options.distinct_density_roles) {
        for (amrex::MFIter mfi(conserved_input); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto input = conserved_input.array(mfi);
            const auto target = conserved_target.array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Real coordinate = static_cast<Real>(i) / Real(nx);
                    input(i, j, k, Rho_comp) =
                        Real(1.0) + input_slope * coordinate;
                    target(i, j, k, Rho_comp) =
                        Real(1.0) + target_slope * coordinate;
                });
        }
    }

    auto& spectrum = state_manager.new_state_for_initialization(0);
    spectrum.setVal(Real(0.0));
    const auto& population = state_manager.layout().populations().front();
    const int mass0 = population.mass_offset;
    const int mass1 = mass0 + 1;
    const int number0 = population.number_offset;
    const int number1 = number0 < 0 ? -1 : number0 + 1;
    const Real lower_edge0 = population.grid.edges()[0];
    const Real upper_edge0 = population.grid.edges()[1];
    const int property0 = options.attached_property
                              ? state_manager.layout().property_offset(0)
                              : -1;
    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto rho = conserved_input.const_array(mfi);
        const auto state = spectrum.array(mfi);
        const bool two_moment = options.mode == erf_sbm::MomentMode::TwoMoment;
        const bool discontinuity = options.discontinuity;
        const bool has_property = options.attached_property;
        const bool smooth_profile = options.smooth_profile;
        const bool complex_profile = options.complex_profile;
        const bool adversarial_one_moment_profile =
            options.adversarial_one_moment_profile;
        const bool canonical_upper_edge_profile =
            options.canonical_upper_edge_profile;
        const int canonical_upper_edge_seed =
            options.canonical_upper_edge_seed;
        const bool noncanonical_shared_edge = options.noncanonical_shared_edge;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                    int k) noexcept {
            const Real density = rho(i, j, k, Rho_comp);
            const Real phase = static_cast<Real>(i) / Real(nx);
            const bool high_cell = discontinuity && (i % 2 == 0);
            const int pattern = (17 * i + 31 * j + 13 * k) % 7;
            const Real scale =
                discontinuity && !high_cell ? Real(0.02) : Real(1.0);
            if (two_moment) {
                if (noncanonical_shared_edge) {
                    const Real number = Real(0.01);
                    state(i, j, k, mass0) = density * Real(0.5) * number;
                    state(i, j, k, number0) = density * number;
                    return;
                }
                if (canonical_upper_edge_profile) {
                    const Real low_endpoint =
                        canonical_upper_edge_probe_low_endpoint(
                            i, canonical_upper_edge_seed);
                    const Real high_endpoint = Real(0.01);
                    const Real number_bin0 = low_endpoint + high_endpoint;
                    const Real mass_bin0 = lower_edge0 * low_endpoint +
                                           upper_edge0 * high_endpoint;
                    const Real number_bin1 = Real(0.02);
                    const Real mass_bin1 =
                        Real(0.75) * number_bin1;
                    state(i, j, k, mass0) = density * mass_bin0;
                    state(i, j, k, number0) = density * number_bin0;
                    state(i, j, k, mass1) = density * mass_bin1;
                    state(i, j, k, number1) = density * number_bin1;
                    if (has_property) {
                        state(i, j, k, property0) =
                            density * number_bin0 * Real(0.25);
                        state(i, j, k, property0 + 1) =
                            density * number_bin1 * Real(0.25);
                    }
                    return;
                }
                const Real n0 =
                    complex_profile
                        ? (pattern == 0 ? Real(0.01) : Real(0.00001))
                        : scale *
                              (smooth_profile
                                   ? Real(0.01) *
                                         (Real(1.0) +
                                          Real(0.1) * amrex::Math::sinpi(
                                                          Real(2.0) * phase))
                                   : Real(0.01));
                const Real mean0 =
                    complex_profile
                        ? (pattern % 2 == 0 ? Real(0.11) : Real(0.45))
                    : smooth_profile
                        ? Real(0.30) +
                              Real(0.04) * amrex::Math::cospi(Real(2.0) * phase)
                        : (discontinuity ? (high_cell ? Real(0.15) : Real(0.40))
                                         : Real(0.30));
                const Real mean1 =
                    complex_profile
                        ? (pattern % 3 == 0 ? Real(0.51) : Real(0.99))
                    : smooth_profile
                        ? Real(0.75) +
                              Real(0.05) * amrex::Math::sinpi(Real(2.0) * phase)
                        : (discontinuity ? (high_cell ? Real(0.75) : Real(0.70))
                                         : Real(0.75));
                const Real m0 = n0 * mean0;
                const Real n1 =
                    complex_profile
                        ? (pattern == 1 ? Real(0.02) : Real(0.00002))
                        : scale *
                              (smooth_profile
                                   ? Real(0.02) *
                                         (Real(1.0) +
                                          Real(0.08) * amrex::Math::cospi(
                                                           Real(2.0) * phase))
                                   : Real(0.02));
                const Real m1 = n1 * mean1;
                state(i, j, k, mass0) = density * m0;
                state(i, j, k, mass1) = density * m1;
                state(i, j, k, number0) = density * n0;
                state(i, j, k, number1) = density * n1;
                if (has_property) {
                    const Real ratio =
                        complex_profile
                            ? ((i + j + k) % 2 == 0 ? Real(0.0) : Real(0.5))
                            : (high_cell ? Real(0.0) : Real(0.5));
                    state(i, j, k, property0) = density * n0 * ratio;
                    state(i, j, k, property0 + 1) = density * n1 * ratio;
                }
            } else {
                const int profile_cell = i % 8;
                // At CFL 0.88, this asymmetric positive stencil triggers
                // WENO overshoot without a roundoff-sensitive FCT result.
                Real adversarial_m0 = Real(1.0e-5);
                if (profile_cell == 1 || profile_cell == 4 ||
                    profile_cell == 7) {
                    adversarial_m0 = Real(5.0e-5);
                } else if (profile_cell == 2 || profile_cell == 6) {
                    adversarial_m0 = Real(5.0e-4);
                } else if (profile_cell == 3) {
                    adversarial_m0 = Real(1.0e-3);
                } else if (profile_cell == 5) {
                    adversarial_m0 = Real(1.0e-4);
                }
                const Real m0 =
                    adversarial_one_moment_profile
                        ? adversarial_m0
                        : scale * (discontinuity && high_cell ? Real(0.04)
                                                              : Real(0.001));
                const Real smooth_m0 =
                    Real(0.02) +
                    Real(0.005) * amrex::Math::sinpi(Real(2.0) * phase);
                const Real smooth_m1 =
                    Real(0.03) +
                    Real(0.002) * amrex::Math::cospi(Real(2.0) * phase);
                const Real m1 =
                    smooth_profile
                        ? smooth_m1
                        : scale * (discontinuity && high_cell ? Real(0.02)
                                                              : Real(0.003));
                state(i, j, k, mass0) =
                    density * (smooth_profile ? smooth_m0 : m0);
                state(i, j, k, mass1) = density * m1;
                if (has_property) {
                    const Real ratio = high_cell ? Real(0.0) : Real(0.5);
                    state(i, j, k, property0) =
                        density * (m0 / Real(0.3)) * ratio;
                    state(i, j, k, property0 + 1) =
                        density * (m1 / Real(0.75)) * ratio;
                }
            }
        });
    }
    MultiFab initial(ba, dm, layout.ncomp(), 0);
    MultiFab::Copy(initial, spectrum, 0, 0, layout.ncomp(), 0);
    if (options.canonical_upper_edge_profile &&
        !erf_sbm::authoritative_state_admissible(initial, layout, 0,
                                                  &diagnostic)) {
        ADD_FAILURE() << "canonical edge probe starts inadmissible: "
                      << diagnostic;
        return {};
    }

    MultiFab avg_xmom(amrex::convert(ba, IntVect::TheDimensionVector(0)), dm, 1,
                      0);
    MultiFab avg_ymom(amrex::convert(ba, IntVect::TheDimensionVector(1)), dm, 1,
                      0);
    MultiFab avg_zmom(amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1,
                      0);
    avg_xmom.setVal(options.carrier_x);
    avg_ymom.setVal(options.multidirectional_carrier ? Real(0.11) : Real(0.0));
    avg_zmom.setVal(options.multidirectional_carrier ? Real(0.09) : Real(0.0));

    // In the full-RK3 path, dt is the completed physical-step interval;
    // each stage below supplies ERF's own stage interval and times.
    const double dt = options.dt;
    const double stage0_interval =
        options.advance_all_rk3_stages ? dt / 3.0 : dt;
    const double stage0_target_time = stage0_interval;
    std::unique_ptr<MultiFab> native_candidate;
    std::unique_ptr<MultiFab> direct_moment_candidate;
    std::unique_ptr<MultiFab> historical_candidate;
    if (options.compare_native_candidate) {
        native_candidate = make_native_candidate(
            layout, initial, conserved_input, avg_xmom, avg_ymom, avg_zmom,
            transport.static_measure(0), geom, stage0_interval);
    }
    Real native_candidate_interior_upper_gap_fraction =
        std::numeric_limits<Real>::max();
    if (options.canonical_upper_edge_profile && native_candidate) {
        native_candidate_interior_upper_gap_fraction =
            minimum_normalized_upper_gap(*native_candidate, mass0, number0,
                                         upper_edge0);
    }
    if (options.adversarial_one_moment_profile) {
        if (options.mode != erf_sbm::MomentMode::OneMoment ||
            !native_candidate) {
            ADD_FAILURE() << "the adversarial 1M fixture requires an "
                             "independent native candidate";
            return {};
        }
        if (!erf_sbm::authoritative_state_admissible(initial, layout, 0,
                                                      &diagnostic)) {
            ADD_FAILURE() << "adversarial initial state is inadmissible: "
                          << diagnostic;
            return {};
        }
        const Real native_min = native_candidate->min(mass0);
        const Real native_scale = std::max(
            native_candidate->norm0(), std::numeric_limits<Real>::min());
        // Require a physical positivity violation with a large roundoff margin.
        EXPECT_LT(native_min, -Real(1.0e-6) * native_scale)
            << "unrestricted native WENO candidate did not materially violate "
               "1M positivity";
    }
    if (options.compare_direct_moment_candidate) {
        direct_moment_candidate = make_native_candidate(
            layout, initial, conserved_input, avg_xmom, avg_ymom, avg_zmom,
            transport.static_measure(0), geom, dt, true);
    }
    if (options.compare_historical_p2) {
        historical_candidate = make_historical_p2_cartesian_candidate(
            layout, initial, conserved_input, avg_xmom, avg_ymom, avg_zmom,
            transport.static_measure(0), geom, dt);
    }
    if (direct_register) {
        if (options.completed_steps != 1 ||
            (!options.advance_all_rk3_stages && !options.anelastic_heun)) {
            ADD_FAILURE() << "direct YA reference check requires one complete RK3 or Heun step";
            return {};
        }
        if (!direct_register->begin_interval(0.0, dt, 2, diagnostic) ||
            !transport.begin_physical_step(0, 0.0, dt, diagnostic)) {
            ADD_FAILURE() << diagnostic;
            return {};
        }
    }
    if (!state_manager.begin_step(0, 0.0, diagnostic)) {
        ADD_FAILURE() << diagnostic;
        return {};
    }
    transport.advance_stage_from_host(
        0,
        options.anelastic_heun ? erf_auxiliary::HostIntegrator::AnelasticHeun
                               : erf_auxiliary::HostIntegrator::CompressibleRK3,
        0, 0.0, 0.0, stage0_target_time, stage0_interval, state_manager,
        conserved_anchor, conserved_input, conserved_target, avg_xmom,
        avg_ymom, avg_zmom, geom, 1, 2);

    const auto& final_state = state_manager.new_state(0);
    if (!erf_sbm::authoritative_state_admissible(final_state, layout, 0,
                                                 &diagnostic)) {
        ADD_FAILURE() << diagnostic;
        return {};
    }
    Real minimum_interior_upper_gap_fraction =
        std::numeric_limits<Real>::max();
    Real minimum_interior_number = std::numeric_limits<Real>::max();
    if (options.canonical_upper_edge_profile) {
        MultiFab canonical_metrics(ba, dm, 2, 0);
        for (amrex::MFIter mfi(canonical_metrics); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto state = final_state.const_array(mfi);
            const auto metrics = canonical_metrics.array(mfi);
            const Real upper = upper_edge0;
            const int mass_component = mass0;
            const int number_component = number0;
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Real mass = state(i, j, k, mass_component);
                    const Real number = state(i, j, k, number_component);
                    const Real scale = upper * number + mass;
                    metrics(i, j, k, 0) = (upper * number - mass) / scale;
                    metrics(i, j, k, 1) = number;
                });
        }
        minimum_interior_upper_gap_fraction = canonical_metrics.min(0);
        minimum_interior_number = canonical_metrics.min(1);
    }
    EXPECT_TRUE(transport.projected_ledger(0).step_active());
    EXPECT_EQ(transport.projected_ledger(0).next_stage(), 1);

    for (int component = 0; component < layout.ncomp(); ++component) {
        const Real final_inventory =
            options.mapped_geometry
                ? mapped_inventory(final_state, transport.static_measure(0),
                                   component)
                : final_state.sum(component);
        const Real initial_inventory =
            options.mapped_geometry
                ? mapped_inventory(initial, transport.static_measure(0),
                                   component)
                : initial.sum(component);
        EXPECT_NEAR(final_inventory, initial_inventory,
                    Real(128.0) * std::numeric_limits<Real>::epsilon() *
                        std::max(std::abs(initial_inventory),
                                 std::numeric_limits<Real>::min()));
    }

    MultiFab difference(ba, dm, layout.ncomp(), 0);
    MultiFab::Copy(difference, final_state, 0, 0, layout.ncomp(), 0);
    difference.minus(initial, 0, layout.ncomp(), 0);
    if (options.carrier_x == Real(0.0) ||
        (options.varying_density && !options.smooth_profile)) {
        // Zero carrier is the identity. A constant dry-air ratio under a
        // constant divergence-free carrier remains unchanged even when rho
        // varies.
        EXPECT_LE(difference.norm0(),
                  Real(256.0) * std::numeric_limits<Real>::epsilon());
    } else if (options.discontinuity) {
        EXPECT_GT(
            difference.norm0(),
            Real(2.0) * std::numeric_limits<Real>::epsilon() *
                std::max(initial.norm0(), std::numeric_limits<Real>::min()));
    } else if (options.smooth_profile || options.complex_profile) {
        EXPECT_GT(
            difference.norm0(),
            Real(2.0) * std::numeric_limits<Real>::epsilon() *
                std::max(initial.norm0(), std::numeric_limits<Real>::min()));
    }

    if (native_candidate) {
        MultiFab candidate_error(ba, dm, layout.ncomp(), 0);
        MultiFab::Copy(candidate_error, final_state, 0, 0, layout.ncomp(), 0);
        candidate_error.minus(*native_candidate, 0, layout.ncomp(), 0);
        if (options.adversarial_one_moment_profile) {
            EXPECT_GT(candidate_error.norm0(),
                      Real(1.0e-6) *
                          std::max(native_candidate->norm0(),
                                   std::numeric_limits<Real>::min()))
                << "accepted transport should materially correct the "
                   "inadmissible native candidate";
        } else if (options.expect_limiter_active) {
            EXPECT_GT(candidate_error.norm0(),
                      Real(0.25) * std::numeric_limits<Real>::epsilon() *
                          std::max(native_candidate->norm0(),
                                   std::numeric_limits<Real>::min()));
        } else {
            EXPECT_LE(candidate_error.norm0(),
                      Real(2.0) * std::numeric_limits<Real>::epsilon() *
                          std::max(native_candidate->norm0(),
                                   std::numeric_limits<Real>::min()));
        }
        if (options.distinct_density_roles) {
            auto wrong_target_candidate = make_native_candidate(
                layout, initial, conserved_target, avg_xmom, avg_ymom,
                avg_zmom, transport.static_measure(0), geom,
                stage0_interval);
            MultiFab wrong_density_difference(ba, dm, layout.ncomp(), 0);
            MultiFab::Copy(wrong_density_difference, *native_candidate, 0, 0,
                           layout.ncomp(), 0);
            wrong_density_difference.minus(*wrong_target_candidate, 0,
                                          layout.ncomp(), 0);
            EXPECT_GT(wrong_density_difference.norm0(),
                      Real(64.0) * std::numeric_limits<Real>::epsilon() *
                          std::max(native_candidate->norm0(),
                                   std::numeric_limits<Real>::min()));
        }
    }

    if (direct_moment_candidate) {
        if (native_candidate) {
            MultiFab noncommutation(ba, dm, layout.ncomp(), 0);
            MultiFab::Copy(noncommutation, *native_candidate, 0, 0,
                           layout.ncomp(), 0);
            noncommutation.minus(*direct_moment_candidate, 0, layout.ncomp(),
                                 0);
            EXPECT_GT(noncommutation.norm0(),
                      Real(256.0) * std::numeric_limits<Real>::epsilon() *
                          std::max(native_candidate->norm0(),
                                   std::numeric_limits<Real>::min()));
        } else {
            ADD_FAILURE() << "direct-moment comparison needs the native "
                             "endpoint candidate";
        }
    }

    if (historical_candidate) {
        MultiFab historical_error(ba, dm, layout.ncomp(), 0);
        MultiFab::Copy(historical_error, final_state, 0, 0, layout.ncomp(), 0);
        historical_error.minus(*historical_candidate, 0, layout.ncomp(), 0);
        EXPECT_LE(historical_error.norm0(),
                  Real(2048.0) * std::numeric_limits<Real>::epsilon() *
                      std::max(historical_candidate->norm0(),
                               std::numeric_limits<Real>::min()));
    }

    if (options.anelastic_heun) {
        auto heun_forward_euler = make_native_candidate(
            layout, final_state, conserved_input, avg_xmom, avg_ymom, avg_zmom,
            transport.static_measure(0), geom, dt);
        MultiFab heun_expected(ba, dm, layout.ncomp(), 0);
        for (amrex::MFIter mfi(heun_expected); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto old = initial.const_array(mfi);
            const auto fe = heun_forward_euler->const_array(mfi);
            const auto out = heun_expected.array(mfi);
            const int ncomp = layout.ncomp();
            amrex::ParallelFor(
                bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    out(i, j, k, n) =
                        Real(0.5) * (old(i, j, k, n) + fe(i, j, k, n));
                });
        }
        transport.advance_stage_from_host(
            0, erf_auxiliary::HostIntegrator::AnelasticHeun, 1, 0.0, dt, dt, dt,
            state_manager, conserved_anchor, conserved_input, conserved_target,
            avg_xmom, avg_ymom, avg_zmom, geom, 1, 2);
        MultiFab heun_error(ba, dm, layout.ncomp(), 0);
        MultiFab::Copy(heun_error, final_state, 0, 0, layout.ncomp(), 0);
        heun_error.minus(heun_expected, 0, layout.ncomp(), 0);
        if (options.expect_limiter_active) {
            EXPECT_GT(heun_error.norm0(),
                      Real(0.25) * std::numeric_limits<Real>::epsilon() *
                          std::max(heun_expected.norm0(),
                                   std::numeric_limits<Real>::min()));
        } else {
            EXPECT_LE(heun_error.norm0(),
                      Real(8.0) * std::numeric_limits<Real>::epsilon() *
                          std::max(heun_expected.norm0(),
                                   std::numeric_limits<Real>::min()));
        }
        EXPECT_TRUE(transport.projected_ledger(0).step_complete());
        EXPECT_EQ(transport.projected_ledger(0).next_stage(), 2);
    }

    if (options.advance_all_rk3_stages) {
        if (options.completed_steps < 1) {
            ADD_FAILURE() << "RK3 transport test requires at least one step";
            return {};
        }
        for (int step = 0; step < options.completed_steps; ++step) {
            const double step_old_time = static_cast<double>(step) * dt;
            const double step_new_time = step_old_time + dt;
            const double stage0_time = step_old_time + dt / 3.0;
            const double stage1_time = step_old_time + dt / 2.0;
            if (step > 0) {
                if (!state_manager.begin_step(0, step_old_time, diagnostic)) {
                    ADD_FAILURE() << diagnostic;
                    return {};
                }
                transport.advance_stage_from_host(
                    0, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
                    step_old_time, step_old_time, stage0_time, dt / 3.0,
                    state_manager, conserved_anchor, conserved_input,
                    conserved_target, avg_xmom, avg_ymom, avg_zmom, geom, 1, 2);
            }
            if (step == 0 && options.amplify_corrector_input) {
                auto& stage_input = state_manager.new_target_storage(0);
                for (amrex::MFIter mfi(stage_input); mfi.isValid(); ++mfi) {
                    const Box bx = mfi.validbox();
                    const auto state = stage_input.array(mfi);
                    const int ncomp = layout.ncomp();
                    amrex::ParallelFor(
                        bx, ncomp,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k,
                                             int n) noexcept {
                            state(i, j, k, n) *= Real(10.0);
                        });
                }
            }
            transport.advance_stage_from_host(
                0, erf_auxiliary::HostIntegrator::CompressibleRK3, 1,
                step_old_time, stage0_time, stage1_time, dt / 2.0,
                state_manager, conserved_anchor, conserved_input,
                conserved_target, avg_xmom, avg_ymom, avg_zmom, geom, 1, 2);
            transport.advance_stage_from_host(
                0, erf_auxiliary::HostIntegrator::CompressibleRK3, 2,
                step_old_time, stage1_time, step_new_time, dt, state_manager,
                conserved_anchor, conserved_input, conserved_target, avg_xmom,
                avg_ymom, avg_zmom, geom, 1, 2);
            EXPECT_TRUE(transport.projected_ledger(0).step_complete());
            EXPECT_EQ(transport.projected_ledger(0).next_stage(), 3);
            const auto& accepted_spectrum = state_manager.new_state(0);
            EXPECT_TRUE(accepted_spectrum.is_finite(0, accepted_spectrum.nComp(), 0));
            diagnostic.clear();
            EXPECT_TRUE(erf_sbm::authoritative_state_admissible(
                accepted_spectrum, layout, 0, &diagnostic))
                << "after completed step " << step << ": " << diagnostic;
            for (int component = 0; component < layout.ncomp(); ++component) {
                EXPECT_NEAR(accepted_spectrum.sum(component), initial.sum(component),
                            Real(512.0) *
                                std::numeric_limits<Real>::epsilon() *
                                std::max(std::abs(initial.sum(component)),
                                         std::numeric_limits<Real>::min()))
                    << "component=" << component << " step=" << step;
            }
        }
    }

    if (direct_register) {
        if (!transport.accept_physical_step(0, 0.0, dt, diagnostic)) {
            ADD_FAILURE() << "coarse transport step rejected by direct YA owner: "
                          << diagnostic;
            return {};
        }

        // This fixture supplies real coarse-side rates from the FCT
        // transport. Complete the interface with two zero fine substeps so
        // the owner can form a reflux candidate without inventing any fine
        // transfer. The reference below uses the full accepted integral
        // ledger produced by the same transport stages.
        erf_auxiliary::MappedFaceFluxRate zero_fine_rate;
        zero_fine_rate.define(fine_interface_ba, fine_interface_dm,
                              layout.ncomp(), 0);
        zero_fine_rate.setVal(Real(0.0));
        std::vector<int> global_components(
            static_cast<std::size_t>(layout.ncomp()));
        for (int component = 0; component < layout.ncomp(); ++component) {
            global_components[static_cast<std::size_t>(component)] = component;
        }
        const auto method = options.anelastic_heun
            ? erf_auxiliary::HostIntegrator::AnelasticHeun
            : erf_auxiliary::HostIntegrator::CompressibleRK3;
        const int stages = options.anelastic_heun ? 2 : 3;
        for (int substep = 0; substep < 2; ++substep) {
            const double fine_begin = dt * static_cast<double>(substep) / 2.0;
            const double fine_end = dt * static_cast<double>(substep + 1) / 2.0;
            const double fine_dt = fine_end - fine_begin;
            if (!direct_register->begin_level_step(false, fine_begin, fine_end,
                                                   diagnostic)) {
                ADD_FAILURE() << diagnostic;
                return {};
            }
            for (int stage = 0; stage < stages; ++stage) {
                if (!direct_register->begin_stage(false, method, stage,
                                                  fine_begin, diagnostic)) {
                    ADD_FAILURE() << diagnostic;
                    return {};
                }
                const double ledger_time = options.anelastic_heun
                    ? 0.5 * fine_dt : (stage == 2 ? fine_dt : 0.0);
                if (!direct_register->add_stage_rates(
                        false, zero_fine_rate, global_components, ledger_time,
                        diagnostic) ||
                    !direct_register->finish_stage(false, diagnostic)) {
                    ADD_FAILURE() << diagnostic;
                    return {};
                }
            }
            if (!direct_register->accept_level_step(false, fine_begin, fine_end,
                                                    diagnostic)) {
                ADD_FAILURE() << diagnostic;
                return {};
            }
        }

        amrex::YAFluxRegister reference_register(
            fine_interface_ba, ba, fine_interface_dm, dm,
            *fine_interface_geom, geom, interface_ratio, 1, layout.ncomp());
        reference_register.reset();
        const auto& integrated =
            transport.completed_spectral_ledger(0).integrated_flux();
        MultiFab no_allocation_cells(
            ba, dm, 1, 0, amrex::MFInfo().SetAlloc(false));
        const auto dx = geom.CellSizeArray();
        const amrex::RunOn run_on = amrex::Gpu::inLaunchRegion()
            ? amrex::RunOn::Gpu : amrex::RunOn::Host;
        for (amrex::MFIter mfi(no_allocation_cells,
                               amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            std::array<amrex::FArrayBox const*, AMREX_SPACEDIM> faces{};
            for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
                faces[static_cast<std::size_t>(direction)] =
                    &integrated.dir(direction)[mfi];
            }
            reference_register.CrseAdd(
                mfi, faces, dx.data(), Real(1.0), 0, 0, layout.ncomp(),
                run_on);
        }
        MultiFab direct_candidate(ba, dm, layout.ncomp(), 0);
        MultiFab reference_candidate(ba, dm, layout.ncomp(), 0);
        direct_candidate.setVal(Real(0.0));
        reference_candidate.setVal(Real(0.0));
        if (!direct_register->reflux(direct_candidate, diagnostic)) {
            ADD_FAILURE() << diagnostic;
            return {};
        }
        reference_register.Reflux(reference_candidate, 0, 0, layout.ncomp());
        amrex::Gpu::streamSynchronize();
        // norm0 skips NaN, so check finiteness before comparing.
        EXPECT_TRUE(direct_candidate.is_finite(0, layout.ncomp(), 0));
        EXPECT_TRUE(reference_candidate.is_finite(0, layout.ncomp(), 0));
        MultiFab direct_reference_error(ba, dm, layout.ncomp(), 0);
        MultiFab::Copy(direct_reference_error, direct_candidate, 0, 0,
                       layout.ncomp(), 0);
        direct_reference_error.minus(reference_candidate, 0,
                                     layout.ncomp(), 0);
        const Real comparison_scale = std::max(
            {direct_candidate.norm0(), reference_candidate.norm0(), Real(1.0)});
        EXPECT_LE(direct_reference_error.norm0(),
                  Real(4096.0) * std::numeric_limits<Real>::epsilon() *
                      comparison_scale)
            << "limiter-accepted direct YA contributions differ from the "
               "completed-integral reference";
        EXPECT_TRUE(direct_register->accept_reflux(diagnostic)) << diagnostic;
    }

    const auto& current_spectrum = state_manager.new_state(0);
    std::vector<Real> spectral_ledger_relative_errors;
    if (options.verify_spectral_ledger_identity) {
        if (!options.advance_all_rk3_stages || options.completed_steps != 1) {
            ADD_FAILURE() << "completed spectral-ledger identity requires one full RK3 step";
            return {};
        }
        const auto& ledger = transport.completed_spectral_ledger(0);
        if (!ledger.step_complete()) {
            ADD_FAILURE() << "completed spectral ledger did not finish the RK3 step";
            return {};
        }
        MultiFab identity_error(ba, dm, layout.ncomp(), 0);
        const auto& integrated_flux = ledger.integrated_flux();
        const auto& omega_field = transport.static_measure(0);
        const auto inv_dx = geom.InvCellSizeArray();
        for (amrex::MFIter mfi(identity_error, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const Box bx = mfi.tilebox();
            const auto old_state = initial.const_array(mfi);
            const auto new_state = current_spectrum.const_array(mfi);
            const auto omega = omega_field.const_array(mfi);
            const auto fx = integrated_flux.dir(0).const_array(mfi);
            const auto fy = integrated_flux.dir(1).const_array(mfi);
            const auto fz = integrated_flux.dir(2).const_array(mfi);
            const auto error = identity_error.array(mfi);
            const int ncomp = layout.ncomp();
            amrex::ParallelFor(
                bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k,
                                     int component) noexcept {
                    const Real divergence =
                        erf_auxiliary::ComputationalMappedDivergence(
                            fx(i + 1, j, k, component), fx(i, j, k, component),
                            fy(i, j + 1, k, component), fy(i, j, k, component),
                            fz(i, j, k + 1, component), fz(i, j, k, component),
                            inv_dx[0], inv_dx[1], inv_dx[2]);
                    const Real predicted =
                        (omega(i, j, k, 0) * old_state(i, j, k, component) -
                         divergence) / omega(i, j, k, 0);
                    error(i, j, k, component) = amrex::Math::abs(
                        predicted - new_state(i, j, k, component));
                });
        }
        for (int component = 0; component < layout.ncomp(); ++component) {
            const Real scale = std::max(
                {initial.norm0(component), current_spectrum.norm0(component),
                 std::numeric_limits<Real>::min()});
            spectral_ledger_relative_errors.push_back(
                identity_error.norm0(component) / scale);
        }
    }

    MultiFab expected_core(ba, dm, 3, 0);
    MultiFab::Copy(expected_core, conserved_target, 0, 0, 3, 0);
    for (amrex::MFIter mfi(expected_core); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = current_spectrum.const_array(mfi);
        const auto core = expected_core.array(mfi);
        const int m0_comp = mass0;
        const int m1_comp = mass1;
        amrex::ParallelFor(bx,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                               core(i, j, k, 1) = state(i, j, k, m0_comp);
                               core(i, j, k, 2) = state(i, j, k, m1_comp);
                           });
    }
    for (int comp = 1; comp <= 2; ++comp) {
        MultiFab projection_error(ba, dm, 1, 0);
        for (amrex::MFIter mfi(expected_core); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto expected = expected_core.const_array(mfi);
            const auto actual = conserved_target.const_array(mfi);
            const auto error = projection_error.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                        int k) noexcept {
                error(i, j, k, 0) = amrex::Math::abs(expected(i, j, k, comp) -
                                                     actual(i, j, k, comp));
            });
        }
        EXPECT_LE(projection_error.norm0(),
                  Real(128.0) * std::numeric_limits<Real>::epsilon());
    }

    if (options.attached_property) {
        MultiFab invalid(ba, dm, 1, 0);
        invalid.setVal(Real(0.0));
        const int property_offset = layout.property_offset(0);
        const auto& final_population = layout.populations().front();
        for (amrex::MFIter mfi(invalid); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto state = current_spectrum.const_array(mfi);
            const auto bad = invalid.array(mfi);
            const int carrier0 = final_population.number_offset >= 0
                                     ? final_population.number_offset
                                     : final_population.mass_offset;
            const int carrier1 = carrier0 + 1;
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (state(i, j, k, property_offset) < Real(0.0) ||
                        state(i, j, k, property_offset) >
                            Real(0.5) * state(i, j, k, carrier0) ||
                        state(i, j, k, property_offset + 1) < Real(0.0) ||
                        state(i, j, k, property_offset + 1) >
                            Real(0.5) * state(i, j, k, carrier1)) {
                        bad(i, j, k, 0) = Real(1.0);
                    }
                });
        }
        EXPECT_EQ(invalid.norm0(), Real(0.0));
    }

    MultiFab fingerprints(ba, dm, 2 * layout.ncomp(), 0);
    for (amrex::MFIter mfi(current_spectrum); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = current_spectrum.const_array(mfi);
        const auto out = fingerprints.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                const Real value = state(i, j, k, n);
                const Real weight = Real(1.0) + Real(0.01) * i +
                                    Real(0.001) * j + Real(0.0001) * k;
                out(i, j, k, n) = weight * value;
                out(i, j, k, ncomp + n) = value * value;
            });
    }
    RunSummary summary;
    summary.minimum_interior_upper_gap_fraction =
        minimum_interior_upper_gap_fraction;
    summary.minimum_interior_number = minimum_interior_number;
    summary.native_candidate_interior_upper_gap_fraction =
        native_candidate_interior_upper_gap_fraction;
    summary.spectral_ledger_relative_errors =
        std::move(spectral_ledger_relative_errors);
    for (int component = 0; component < layout.ncomp(); ++component) {
        summary.inventory.push_back(current_spectrum.sum(component));
        summary.weighted_inventory.push_back(fingerprints.sum(component));
        summary.squared_inventory.push_back(
            fingerprints.sum(layout.ncomp() + component));
    }
    return summary;
}

#if !defined(AMREX_USE_FLOAT)
struct SmoothTranslationErrors
{
    int nz;
    Real mass_l1;
    Real number_l1;
    Real native_mass_l1;
    Real native_number_l1;
    Real max_mass_difference;
    Real max_number_difference;
};

void
run_smooth_vertical_mapped_density_weighted_translation_test (
    std::vector<SmoothTranslationErrors>& measurements)
{
    constexpr int nx = 4;
    constexpr int ny = 4;
    constexpr Real carrier = Real(0.2);
    constexpr double final_time = 0.25;
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    const auto& population = layout.populations().front();
    const int mass_component = population.mass_offset;
    const int number_component = population.number_offset;
    const Real lower_edge = population.grid.edges().front();
    const Real upper_edge = population.grid.edges()[1];

    for (const int nz : {16, 32, 64}) {
        const Box domain(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
        const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
        const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
        const Geometry geom(domain, &physical, amrex::CoordSys::cartesian,
                            periodic);
        const BoxArray ba(domain);
        const DistributionMapping dm(ba);

        erf_sbm::SBMStateManager state_manager(layout, 1);
        state_manager.define(0, ba, dm, 0.0);
        erf_sbm::SBMTransport transport(layout, 1);
        transport.define(0, ba, dm);

        const auto map_ba = project_to_xy(ba);
        MultiFab detj(ba, dm, 1, 0);
        MultiFab mx(map_ba, dm, 1, 0);
        MultiFab my(map_ba, dm, 1, 0);
        for (amrex::MFIter mfi(detj); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto det = detj.array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Real zeta = (static_cast<Real>(k) + Real(0.5)) /
                                     static_cast<Real>(nz);
                    det(i, j, k, 0) =
                        Real(1.0) + Real(0.15) *
                                       amrex::Math::sinpi(Real(2.0) * zeta);
                });
        }
        mx.setVal(Real(1.0));
        my.setVal(Real(1.0));
        std::string diagnostic;
        ASSERT_TRUE(transport.rebuild_static_measure(0, detj, mx, my,
                                                      diagnostic))
            << diagnostic;
        const MultiFab& measure = transport.static_measure(0);
        EXPECT_LT(measure.min(0), Real(1.0));
        EXPECT_GT(measure.max(0), Real(1.0));

        MultiFab conserved_anchor(ba, dm, 3, 0);
        MultiFab conserved_input(ba, dm, 3, 0);
        MultiFab conserved_target(ba, dm, 3, 0);
        conserved_anchor.setVal(Real(0.0));
        for (amrex::MFIter mfi(conserved_anchor); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto omega = measure.const_array(mfi);
            const auto conserved = conserved_anchor.array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    conserved(i, j, k, Rho_comp) =
                        Real(1.0) / omega(i, j, k, 0);
                });
        }
        MultiFab::Copy(conserved_input, conserved_anchor, 0, 0, 3, 0);
        MultiFab::Copy(conserved_target, conserved_anchor, 0, 0, 3, 0);

        auto& spectrum = state_manager.new_state_for_initialization(0);
        spectrum.setVal(Real(0.0));
        const Real cell_width = Real(1.0) / static_cast<Real>(nz);
        const Real pi = amrex::Math::pi<Real>();
        // Initialize the analytic profile as finite-volume cell averages.
        for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto rho = conserved_anchor.const_array(mfi);
            const auto state = spectrum.array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Real left = static_cast<Real>(k) * cell_width;
                    const Real right = static_cast<Real>(k + 1) * cell_width;
                    const Real sin_average =
                        (amrex::Math::cospi(Real(2.0) * left) -
                         amrex::Math::cospi(Real(2.0) * right)) /
                        (Real(2.0) * pi * cell_width);
                    const Real cos_average =
                        (amrex::Math::sinpi(Real(2.0) * right) -
                         amrex::Math::sinpi(Real(2.0) * left)) /
                        (Real(2.0) * pi * cell_width);
                    const Real blo = Real(1.0) + Real(0.10) * sin_average;
                    const Real bhi = Real(0.8) + Real(0.08) * cos_average;
                    const Real density = rho(i, j, k, Rho_comp);
                    const Real number = blo + bhi;
                    const Real mass = lower_edge * blo + upper_edge * bhi;
                    state(i, j, k, mass_component) = density * mass;
                    state(i, j, k, number_component) = density * number;
                });
        }
        ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
            spectrum, layout, 0, &diagnostic)) << diagnostic;

        MultiFab avg_xmom(amrex::convert(
                              ba, IntVect::TheDimensionVector(0)),
                          dm, 1, 0);
        MultiFab avg_ymom(amrex::convert(
                              ba, IntVect::TheDimensionVector(1)),
                          dm, 1, 0);
        MultiFab avg_zmom(amrex::convert(
                              ba, IntVect::TheDimensionVector(2)),
                          dm, 1, 0);
        avg_xmom.setVal(Real(0.0));
        avg_ymom.setVal(Real(0.0));
        avg_zmom.setVal(carrier);

        MultiFab initial(ba, dm, layout.ncomp(), 0);
        MultiFab::Copy(initial, spectrum, 0, 0, layout.ncomp(), 0);
        MultiFab native_reference(ba, dm, layout.ncomp(), 0);
        MultiFab::Copy(native_reference, initial, 0, 0, layout.ncomp(), 0);
        const Real initial_mass_inventory =
            mapped_inventory(initial, measure, mass_component);
        const Real initial_number_inventory =
            mapped_inventory(initial, measure, number_component);

        const double dt = 1.0 / static_cast<double>(nz);
        const int step_count = nz / 4;
        for (int step = 0; step < step_count; ++step) {
            const double step_old_time = static_cast<double>(step) * dt;
            const double stage1_time = step_old_time + dt / 3.0;
            const double stage2_time = step_old_time + dt / 2.0;
            const double step_new_time = step_old_time + dt;
            if (!state_manager.begin_step(0, step_old_time, diagnostic)) {
                ADD_FAILURE() << diagnostic;
                return;
            }
            transport.advance_stage_from_host(
                0, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
                step_old_time, step_old_time, stage1_time,
                stage1_time - step_old_time, state_manager, conserved_anchor,
                conserved_input, conserved_target, avg_xmom, avg_ymom,
                avg_zmom, geom, 1, 2);
            transport.advance_stage_from_host(
                0, erf_auxiliary::HostIntegrator::CompressibleRK3, 1,
                step_old_time, stage1_time, stage2_time,
                stage2_time - step_old_time, state_manager, conserved_anchor,
                conserved_input, conserved_target, avg_xmom, avg_ymom,
                avg_zmom, geom, 1, 2);
            transport.advance_stage_from_host(
                0, erf_auxiliary::HostIntegrator::CompressibleRK3, 2,
                step_old_time, stage2_time, step_new_time,
                step_new_time - step_old_time, state_manager, conserved_anchor,
                conserved_input, conserved_target, avg_xmom, avg_ymom,
                avg_zmom, geom, 1, 2);

            // Independent unlimited ERF-native WENO-Z3 reference for the full
            // three-stage compressible RK recurrence. Each stage reconstructs
            // from its input state while retaining the H^n anchor.
            auto native_stage1 = make_native_candidate(
                layout, native_reference, conserved_anchor, avg_xmom,
                avg_ymom, avg_zmom, measure, geom, dt / 3.0);
            auto native_stage2 = make_native_candidate(
                layout, *native_stage1, conserved_anchor, avg_xmom,
                avg_ymom, avg_zmom, measure, geom, dt / 2.0, false,
                &native_reference);
            auto native_stage3 = make_native_candidate(
                layout, *native_stage2, conserved_anchor, avg_xmom,
                avg_ymom, avg_zmom, measure, geom, dt, false,
                &native_reference);
            MultiFab::Copy(native_reference, *native_stage3, 0, 0,
                           layout.ncomp(), 0);
        }

        const auto& final_state = state_manager.new_state(0);
        ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
            final_state, layout, 0, &diagnostic)) << diagnostic;
        EXPECT_TRUE(transport.projected_ledger(0).step_complete());

        const Real final_mass_inventory =
            mapped_inventory(final_state, measure, mass_component);
        const Real final_number_inventory =
            mapped_inventory(final_state, measure, number_component);
        const Real epsilon = std::numeric_limits<Real>::epsilon();
        EXPECT_NEAR(final_mass_inventory, initial_mass_inventory,
                    Real(128.0) * epsilon *
                        std::max(std::abs(initial_mass_inventory),
                                 std::numeric_limits<Real>::min()));
        EXPECT_NEAR(final_number_inventory, initial_number_inventory,
                    Real(128.0) * epsilon *
                        std::max(std::abs(initial_number_inventory),
                                 std::numeric_limits<Real>::min()));

        MultiFab error(ba, dm, 6, 0);
        for (amrex::MFIter mfi(error); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto omega = measure.const_array(mfi);
            const auto state = final_state.const_array(mfi);
            const auto native = native_reference.const_array(mfi);
            const auto out = error.array(mfi);
            const Real translated_distance =
                carrier * static_cast<Real>(final_time);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Real left = static_cast<Real>(k) * cell_width -
                                      translated_distance;
                    const Real right = static_cast<Real>(k + 1) * cell_width -
                                       translated_distance;
                    const Real sin_average =
                        (amrex::Math::cospi(Real(2.0) * left) -
                         amrex::Math::cospi(Real(2.0) * right)) /
                        (Real(2.0) * pi * cell_width);
                    const Real cos_average =
                        (amrex::Math::sinpi(Real(2.0) * right) -
                         amrex::Math::sinpi(Real(2.0) * left)) /
                        (Real(2.0) * pi * cell_width);
                    const Real blo = Real(1.0) + Real(0.10) * sin_average;
                    const Real bhi = Real(0.8) + Real(0.08) * cos_average;
                    const Real exact_number = blo + bhi;
                    const Real exact_mass =
                        lower_edge * blo + upper_edge * bhi;
                    out(i, j, k, 0) = amrex::Math::abs(
                        omega(i, j, k, 0) * state(i, j, k, mass_component) -
                        exact_mass);
                    out(i, j, k, 1) = amrex::Math::abs(
                        omega(i, j, k, 0) * state(i, j, k, number_component) -
                        exact_number);
                    out(i, j, k, 2) = amrex::Math::abs(
                        omega(i, j, k, 0) *
                            native(i, j, k, mass_component) -
                        exact_mass);
                    out(i, j, k, 3) = amrex::Math::abs(
                        omega(i, j, k, 0) *
                            native(i, j, k, number_component) -
                        exact_number);
                    out(i, j, k, 4) = amrex::Math::abs(
                        state(i, j, k, mass_component) -
                        native(i, j, k, mass_component));
                    out(i, j, k, 5) = amrex::Math::abs(
                        state(i, j, k, number_component) -
                        native(i, j, k, number_component));
                });
        }
        const Real cell_count = static_cast<Real>(domain.numPts());
        const Real mass_l1 = error.sum(0) / cell_count;
        const Real number_l1 = error.sum(1) / cell_count;
        const Real native_mass_l1 = error.sum(2) / cell_count;
        const Real native_number_l1 = error.sum(3) / cell_count;
        const Real max_mass_difference = error.norm0(4);
        const Real max_number_difference = error.norm0(5);
        // Roundoff parity over every stored spectral component shows that the
        // group limiter made no material correction during the full evolution.
        for (int component = 0; component < layout.ncomp(); ++component) {
            MultiFab component_difference(ba, dm, 1, 0);
            for (amrex::MFIter mfi(component_difference); mfi.isValid(); ++mfi) {
                const Box bx = mfi.validbox();
                const auto m3 = final_state.const_array(mfi);
                const auto native = native_reference.const_array(mfi);
                const auto diff = component_difference.array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        diff(i, j, k, 0) = amrex::Math::abs(
                            m3(i, j, k, component) -
                            native(i, j, k, component));
                    });
            }
            const Real parity_tolerance =
                Real(512.0) * epsilon *
                std::max(native_reference.norm0(component),
                         std::numeric_limits<Real>::min());
            EXPECT_LE(component_difference.norm0(), parity_tolerance)
                << "complete-step M3/native WENO-Z3 parity failed for N="
                << nz << " component=" << component;
        }
        const Real error_norm_tolerance =
            measure.max(0) * Real(512.0) * epsilon *
                std::max(native_reference.norm0(mass_component),
                         native_reference.norm0(number_component)) +
            Real(16.0) * epsilon;
        EXPECT_NEAR(mass_l1, native_mass_l1, error_norm_tolerance)
            << "analytic mass error differs from native WENO-Z3 at N=" << nz;
        EXPECT_NEAR(number_l1, native_number_l1, error_norm_tolerance)
            << "analytic number error differs from native WENO-Z3 at N="
            << nz;
        measurements.push_back({nz, mass_l1, number_l1, native_mass_l1,
                                native_number_l1, max_mass_difference,
                                max_number_difference});
    }
}

#endif

TEST(SBMTransport, SmoothVerticalMappedDensityWeightedTranslationConverges)
{
#if defined(AMREX_USE_FLOAT)
    GTEST_SKIP() << "Observed smooth-transport order is qualified in DOUBLE; "
                    "SINGLE remains covered by the functional, realizability, "
                    "and parity tests.";
#else
    std::vector<SmoothTranslationErrors> measurements;
    run_smooth_vertical_mapped_density_weighted_translation_test(
        measurements);
    ASSERT_EQ(measurements.size(), 3U);
    ASSERT_EQ(measurements[0].nz, 16);
    ASSERT_EQ(measurements[1].nz, 32);
    ASSERT_EQ(measurements[2].nz, 64);

    const double mass16 = static_cast<double>(measurements[0].mass_l1);
    const double mass32 = static_cast<double>(measurements[1].mass_l1);
    const double mass64 = static_cast<double>(measurements[2].mass_l1);
    const double number16 = static_cast<double>(measurements[0].number_l1);
    const double number32 = static_cast<double>(measurements[1].number_l1);
    const double number64 = static_cast<double>(measurements[2].number_l1);
    const double native_mass16 =
        static_cast<double>(measurements[0].native_mass_l1);
    const double native_mass32 =
        static_cast<double>(measurements[1].native_mass_l1);
    const double native_mass64 =
        static_cast<double>(measurements[2].native_mass_l1);
    const double native_number16 =
        static_cast<double>(measurements[0].native_number_l1);
    const double native_number32 =
        static_cast<double>(measurements[1].native_number_l1);
    const double native_number64 =
        static_cast<double>(measurements[2].native_number_l1);
    const auto observed_order = [](const double coarse, const double fine) {
        return coarse > 0.0 && fine > 0.0
                   ? std::log(coarse / fine) / std::log(2.0)
                   : -std::numeric_limits<double>::infinity();
    };
    const double mass_order_16_32 = observed_order(mass16, mass32);
    const double mass_order_32_64 = observed_order(mass32, mass64);
    const double number_order_16_32 = observed_order(number16, number32);
    const double number_order_32_64 = observed_order(number32, number64);
    const double native_mass_order_16_32 =
        observed_order(native_mass16, native_mass32);
    const double native_mass_order_32_64 =
        observed_order(native_mass32, native_mass64);
    const double native_number_order_16_32 =
        observed_order(native_number16, native_number32);
    const double native_number_order_32_64 =
        observed_order(native_number32, native_number64);

    std::cout << std::setprecision(12)
              << "SBM/native smooth mapped/density-weighted errors (N, M3 mass "
                 "L1, native mass L1, M3 number L1, native number L1, max "
                 "mass U difference, max number U difference):\n";
    for (const auto& result : measurements) {
        std::cout << result.nz << ", " << result.mass_l1 << ", "
                  << result.native_mass_l1 << ", " << result.number_l1
                  << ", " << result.native_number_l1 << ", "
                  << result.max_mass_difference << ", "
                  << result.max_number_difference << '\n';
    }
    std::cout << "observed orders (mass 16-32, 32-64; number 16-32, 32-64): "
              << mass_order_16_32 << ", " << mass_order_32_64 << "; "
              << number_order_16_32 << ", " << number_order_32_64 << '\n';
    std::cout << "native observed orders (mass 16-32, 32-64; number 16-32, "
                 "32-64): "
              << native_mass_order_16_32 << ", "
              << native_mass_order_32_64 << "; "
              << native_number_order_16_32 << ", "
              << native_number_order_32_64 << '\n';

    std::ostringstream mass_report;
    mass_report << std::setprecision(12)
                << "mass L1 errors [N=16,32,64] = [" << mass16 << ", "
                << mass32 << ", " << mass64 << "], observed orders = ["
                << mass_order_16_32 << ", " << mass_order_32_64 << ']';
    std::ostringstream number_report;
    number_report << std::setprecision(12)
                  << "number L1 errors [N=16,32,64] = [" << number16 << ", "
                  << number32 << ", " << number64 << "], observed orders = ["
                  << number_order_16_32 << ", " << number_order_32_64 << ']';

    EXPECT_GT(mass16, mass32) << mass_report.str();
    EXPECT_GT(mass32, mass64) << mass_report.str();
    EXPECT_GT(number16, number32) << number_report.str();
    EXPECT_GT(number32, number64) << number_report.str();
    EXPECT_TRUE(std::isfinite(mass16) && std::isfinite(mass32) &&
                std::isfinite(mass64) && std::isfinite(number16) &&
                std::isfinite(number32) && std::isfinite(number64));
    EXPECT_TRUE(std::isfinite(native_mass16) &&
                std::isfinite(native_mass32) &&
                std::isfinite(native_mass64) &&
                std::isfinite(native_number16) &&
                std::isfinite(native_number32) &&
                std::isfinite(native_number64));
    EXPECT_TRUE(std::isfinite(mass_order_16_32) &&
                std::isfinite(mass_order_32_64) &&
                std::isfinite(number_order_16_32) &&
                std::isfinite(number_order_32_64));
    EXPECT_TRUE(std::isfinite(native_mass_order_16_32) &&
                std::isfinite(native_mass_order_32_64) &&
                std::isfinite(native_number_order_16_32) &&
                std::isfinite(native_number_order_32_64));
#endif
}

TEST(SBMTransport, ZeroCarrierIdentityOneAndTwoMoment)
{
    run_transport({erf_sbm::MomentMode::OneMoment});
    run_transport({erf_sbm::MomentMode::TwoMoment});
}

TEST(SBMTransport, PeriodicAdvectionConservesAndProjectsBothRepresentations)
{
    run_transport(
        {erf_sbm::MomentMode::OneMoment, false, true, false, Real(0.2)});
    run_transport(
        {erf_sbm::MomentMode::TwoMoment, false, true, false, Real(0.2)});
}

TEST(SBMTransport, ConstantDryAirRatioSurvivesVariableDensity)
{
    run_transport(
        {erf_sbm::MomentMode::OneMoment, true, false, false, Real(0.2)});
    run_transport(
        {erf_sbm::MomentMode::TwoMoment, true, false, false, Real(0.2)});
}

TEST(SBMTransport, ERFHostSeamUsesPredictorDensityForCompressibleStage)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.smooth_profile = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.08;
    options.compare_native_candidate = true;
    options.distinct_density_roles = true;
    options.rho_anchor_slope = Real(0.1);
    options.rho_input_slope = Real(0.1);
    options.rho_target_slope = Real(0.55);
    run_transport(options);
}

// Keep device lambdas out of the private TestBody generated by TEST for NVCC.
void run_host_slow_copy_preserves_accepted_liquid_projection_test ()
{
    const Box domain(IntVect(0, 0, 0), IntVect(3, 3, 3));
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    const Geometry geom(domain, &physical, amrex::CoordSys::cartesian,
                        periodic);
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const auto layout = make_layout(erf_sbm::MomentMode::OneMoment);
    erf_sbm::SBMStateManager state_manager(layout, 1);
    state_manager.define(0, ba, dm, 0.0);
    erf_sbm::SBMTransport transport(layout, 1);
    transport.define(0, ba, dm);

    MultiFab detj(ba, dm, 1, 0);
    const auto map_ba = project_to_xy(ba);
    MultiFab mx(map_ba, dm, 1, 0);
    MultiFab my(map_ba, dm, 1, 0);
    detj.setVal(Real(1.0));
    mx.setVal(Real(1.0));
    my.setVal(Real(1.0));
    std::string diagnostic;
    ASSERT_TRUE(transport.rebuild_static_measure(0, detj, mx, my, diagnostic))
        << diagnostic;

    // Place qc, qr and an unrelated component right after RhoQ1_comp so they
    // never overlap it, whatever scalars the build adds before the moisture
    // components (5, 6, 7 and 8 components when RhoQ1_comp is 4).
    constexpr int qc = RhoQ1_comp + 1;
    constexpr int qr = RhoQ1_comp + 2;
    constexpr int unrelated = RhoQ1_comp + 3;
    constexpr int ncomp = RhoQ1_comp + 4;
    MultiFab anchor(ba, dm, ncomp, 0);
    MultiFab predictor(ba, dm, ncomp, 0);
    MultiFab target(ba, dm, ncomp, 0);
    anchor.setVal(Real(0.0));
    anchor.setVal(Real(1.0), Rho_comp, 1, 0);
    anchor.setVal(Real(300.0), RhoTheta_comp, 1, 0);
    anchor.setVal(Real(0.01), RhoQ1_comp, 1, 0);
    MultiFab::Copy(predictor, anchor, 0, 0, ncomp, 0);
    MultiFab::Copy(target, anchor, 0, 0, ncomp, 0);
    predictor.setVal(Real(0.03), RhoQ1_comp, 1, 0);
    predictor.setVal(Real(11.0), qc, 1, 0);
    predictor.setVal(Real(13.0), qr, 1, 0);
    predictor.setVal(Real(17.0), unrelated, 1, 0);

    auto& spectrum = state_manager.new_state_for_initialization(0);
    spectrum.setVal(Real(0.0));
    const auto& population = layout.populations().front();
    spectrum.setVal(Real(0.12), population.mass_offset, 1, 0);
    spectrum.setVal(Real(0.34), population.mass_offset + 1, 1, 0);

    MultiFab avg_xmom(convert(ba, IntVect::TheDimensionVector(0)), dm, 1, 0);
    MultiFab avg_ymom(convert(ba, IntVect::TheDimensionVector(1)), dm, 1, 0);
    MultiFab avg_zmom(convert(ba, IntVect::TheDimensionVector(2)), dm, 1, 0);
    avg_xmom.setVal(Real(0.0));
    avg_ymom.setVal(Real(0.0));
    avg_zmom.setVal(Real(0.0));
    ASSERT_TRUE(state_manager.begin_step(0, 0.0, diagnostic)) << diagnostic;
    transport.advance_stage_from_host(
        0, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
        0.0, 0.0, 0.1, 0.1, state_manager, anchor, predictor, target,
        avg_xmom, avg_ymom, avg_zmom, geom, qc, qr);

    // Exercise the production-used component copy operation immediately after
    // the host stage transaction, before any later end-of-step projection.
    for (amrex::MFIter mfi(target, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto current = target.array(mfi);
        const auto source = predictor.const_array(mfi);
        amrex::ParallelFor(bx, target.nComp() - 2,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int nn) noexcept {
                const int component = 2 + nn;
                erf_sbm::copy_host_slow_component(
                    current, source, i, j, k, component, true, qc, qr);
            });
    }
    EXPECT_EQ(target.min(RhoQ1_comp), Real(0.03));
    EXPECT_EQ(target.min(unrelated), Real(17.0));
    EXPECT_EQ(target.min(qc), Real(0.12));
    EXPECT_EQ(target.min(qr), Real(0.34));

    // A non-SBM host copy continues to own every slow component, including
    // the compact liquid indices.
    MultiFab non_sbm_target(ba, dm, ncomp, 0);
    MultiFab::Copy(non_sbm_target, anchor, 0, 0, ncomp, 0);
    for (amrex::MFIter mfi(non_sbm_target, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto current = non_sbm_target.array(mfi);
        const auto source = predictor.const_array(mfi);
        amrex::ParallelFor(bx, non_sbm_target.nComp() - 2,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int nn) noexcept {
                const int component = 2 + nn;
                erf_sbm::copy_host_slow_component(
                    current, source, i, j, k, component, false, qc, qr);
            });
    }
    EXPECT_EQ(non_sbm_target.min(qc), Real(11.0));
    EXPECT_EQ(non_sbm_target.min(qr), Real(13.0));
    EXPECT_EQ(non_sbm_target.min(RhoQ1_comp), Real(0.03));

    MultiFab projection_error(ba, dm, 2, 0);
    for (amrex::MFIter mfi(projection_error); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto accepted = spectrum.const_array(mfi);
        const auto core = target.const_array(mfi);
        const auto error = projection_error.array(mfi);
        const int mass_offset = population.mass_offset;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            error(i, j, k, 0) = core(i, j, k, qc) - accepted(i, j, k, mass_offset);
            error(i, j, k, 1) = core(i, j, k, qr) - accepted(i, j, k, mass_offset + 1);
        });
    }
    EXPECT_EQ(projection_error.norm0(0), Real(0.0));
    EXPECT_EQ(projection_error.norm0(1), Real(0.0));
}

struct G1GhostSample
{
    Real physical{Real(-1.0)};
    Real endpoint{Real(-1.0)};
    int physical_hit{0};
    int endpoint_hit{0};
};

void run_coarse_stage_g1_ghost_test (const bool inject_stale_write_fault = false,
                                    const bool rounded_three_substeps = false,
                                    const int quotient_fault = 0,
                                    const bool two_moment_attached = false,
                                    const IntVect probe_cell = IntVect(1, 3, 3),
                                    const int probe_component = 0,
                                    G1GhostSample* sample_result = nullptr,
                                    const int probe_stage = -1,
                                    const int probe_role = 0,
                                    const bool full_rk3_step = false,
                                    const bool componentwise_parent_fill = true,
                                    const bool omit_write_tally = false)
{
    const erf_sbm::SBMLayout layout =
        make_layout(two_moment_attached ? erf_sbm::MomentMode::TwoMoment
                                        : erf_sbm::MomentMode::OneMoment,
                    two_moment_attached);
    const Box coarse_domain(IntVect(0, 0, 0), IntVect(3, 3, 3));
    const Box fine_domain = amrex::refine(coarse_domain, IntVect(2));
    const BoxArray coarse_ba(coarse_domain);
    const DistributionMapping coarse_dm(coarse_ba);
    const BoxArray fine_ba(amrex::refine(
        Box(IntVect(1, 1, 1), IntVect(2, 2, 2)), IntVect(2)));
    const DistributionMapping fine_dm(fine_ba);
    const amrex::RealBox fine_physical({0.0, 0.0, 0.0},
                                       {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    const Geometry fine_geom(fine_domain, &fine_physical,
                             amrex::CoordSys::cartesian, periodic);
    const double coarse_begin =
        (rounded_three_substeps || two_moment_attached) ? 0.1 : 0.0;
    const double coarse_end =
        (rounded_three_substeps || two_moment_attached) ? 0.2 : 1.0;
    const double fine_dt = rounded_three_substeps
        ? (coarse_end - coarse_begin) / 3.0 : (coarse_end - coarse_begin);

    erf_sbm::SBMStateManager state_manager(layout, 2);
    std::string diagnostic;
    state_manager.define(0, coarse_ba, coarse_dm, coarse_begin);
    state_manager.define(1, fine_ba, fine_dm, coarse_begin);
    const auto set_two_moment_state = [&layout] (MultiFab& state,
                                                  const bool future) {
        state.setVal(Real(0.0));
        const auto& population = layout.populations().front();
        const int property = layout.property_offset(0);
        const std::array<Real, 6> values = future
            ? std::array<Real, 6>{Real(1.2), Real(1.4), Real(4.0),
                                  Real(2.0), Real(0.8), Real(0.4)}
            : std::array<Real, 6>{Real(0.3), Real(0.7), Real(2.0),
                                  Real(1.0), Real(0.2), Real(0.1)};
        state.setVal(values[0], population.mass_offset, 1, 0);
        state.setVal(values[1], population.mass_offset + 1, 1, 0);
        state.setVal(values[2], population.number_offset, 1, 0);
        state.setVal(values[3], population.number_offset + 1, 1, 0);
        state.setVal(values[4], property, 1, 0);
        state.setVal(values[5], property + 1, 1, 0);
    };
    if (two_moment_attached) {
        set_two_moment_state(state_manager.new_state_for_initialization(0), false);
        set_two_moment_state(state_manager.new_state_for_initialization(1), false);
        ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
            state_manager.new_state(0), layout, 0, &diagnostic)) << diagnostic;
        ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
            state_manager.new_state(1), layout, 1, &diagnostic)) << diagnostic;
    } else {
        state_manager.new_state_for_initialization(0).setVal(Real(1.0));
        state_manager.new_state_for_initialization(1).setVal(Real(1.0));
    }
    if (quotient_fault == 1) {
        // Make the first requested stage itself unrepresentable.  A huge
        // future endpoint alone can drive WENO's smoothness arithmetic out
        // of range before the next stage reaches the intended quotient check.
        state_manager.new_state_for_initialization(0).setVal(
            std::numeric_limits<Real>::max(), 0, 1, 0);
    }
    MultiFab coarse_rho_old(coarse_ba, coarse_dm, 3, 0);
    MultiFab coarse_rho_new(coarse_ba, coarse_dm, 3, 0);
    coarse_rho_old.setVal(Real(0.0));
    coarse_rho_new.setVal(Real(0.0));
    Real rho_old_value = Real(1.0);
    Real rho_new_value = Real(3.0);
    if (quotient_fault == 1) {
        rho_old_value = rho_new_value = std::numeric_limits<Real>::min();
    } else if (quotient_fault == 2) {
        rho_old_value = rho_new_value = std::numeric_limits<Real>::max();
    } else if (quotient_fault == 3) {
        rho_old_value = Real(0.0);
    }
    coarse_rho_old.setVal(rho_old_value, Rho_comp, 1, 0);
    coarse_rho_new.setVal(rho_new_value, Rho_comp, 1, 0);
    if (quotient_fault == 2) {
        state_manager.new_state_for_initialization(0).setVal(
            std::numeric_limits<Real>::denorm_min(), 0, 1, 0);
    }
    ASSERT_TRUE(state_manager.begin_step(0, coarse_begin, diagnostic)) << diagnostic;
    if (two_moment_attached) {
        set_two_moment_state(state_manager.new_target_storage(0), true);
        ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
            state_manager.new_target_storage(0), layout, 0, &diagnostic))
            << diagnostic;
    } else {
        state_manager.new_target_storage(0).setVal(Real(9.0));
    }
    if (quotient_fault == 1) {
        state_manager.new_target_storage(0).setVal(
            std::numeric_limits<Real>::max(), 0, 1, 0);
    } else if (quotient_fault == 2) {
        state_manager.new_target_storage(0).setVal(
            std::numeric_limits<Real>::denorm_min(), 0, 1, 0);
    } else if (quotient_fault == 4) {
        state_manager.new_target_storage(0).setVal(
            std::numeric_limits<Real>::quiet_NaN(), 0, 1, 0);
    }
    ASSERT_TRUE(state_manager.accept_stage_target(0, coarse_end, true, diagnostic))
        << diagnostic;

    erf_sbm::SBMTransport transport(layout, 2, 16, false);
    transport.define(0, coarse_ba, coarse_dm);
    transport.define(1, fine_ba, fine_dm);
    MultiFab coarse_detj(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_detj(fine_ba, fine_dm, 1, 0);
    const auto coarse_map_ba = project_to_xy(coarse_ba);
    const auto fine_map_ba = project_to_xy(fine_ba);
    MultiFab coarse_mx(coarse_map_ba, coarse_dm, 1, 0);
    MultiFab coarse_my(coarse_map_ba, coarse_dm, 1, 0);
    MultiFab fine_mx(fine_map_ba, fine_dm, 1, 0);
    MultiFab fine_my(fine_map_ba, fine_dm, 1, 0);
    coarse_detj.setVal(Real(1.0));
    fine_detj.setVal(Real(1.0));
    coarse_mx.setVal(Real(1.0));
    coarse_my.setVal(Real(1.0));
    fine_mx.setVal(Real(1.0));
    fine_my.setVal(Real(1.0));
    ASSERT_TRUE(transport.rebuild_static_measure(
        1, fine_detj, fine_mx, fine_my, diagnostic)) << diagnostic;

    MultiFab fine_rho_anchor(fine_ba, fine_dm, 3, 0);
    MultiFab fine_rho_input(fine_ba, fine_dm, 3, 0);
    MultiFab fine_rho_target(fine_ba, fine_dm, 3, 0);
    fine_rho_anchor.setVal(Real(0.0));
    fine_rho_input.setVal(Real(0.0));
    fine_rho_target.setVal(Real(0.0));
    fine_rho_anchor.setVal(Real(1.0), Rho_comp, 1, 0);
    fine_rho_input.setVal(Real(1.0), Rho_comp, 1, 0);
    fine_rho_target.setVal(Real(1.0), Rho_comp, 1, 0);
    MultiFab avg_xmom(amrex::convert(fine_ba,
                                    IntVect::TheDimensionVector(0)),
                      fine_dm, 1, 0);
    MultiFab avg_ymom(amrex::convert(fine_ba,
                                    IntVect::TheDimensionVector(1)),
                      fine_dm, 1, 0);
    MultiFab avg_zmom(amrex::convert(fine_ba,
                                    IntVect::TheDimensionVector(2)),
                      fine_dm, 1, 0);
    avg_xmom.setVal(Real(0.0));
    avg_ymom.setVal(Real(0.0));
    avg_zmom.setVal(Real(0.0));

    amrex::Gpu::DeviceVector<Real> sampled_intensive(1, Real(-1.0));
    amrex::Gpu::DeviceVector<Real> sampled_endpoint(1, Real(-1.0));
    amrex::Gpu::DeviceVector<int> sample_hit(1, 0);
    amrex::Gpu::DeviceVector<int> endpoint_hit(1, 0);
    erf_sbm::SBMAMRStageGhostProbe probe;
    probe.fine_cell = probe_cell;
    probe.global_component = probe_component;
    probe.requested_stage = probe_stage >= 0 ? probe_stage
        : (rounded_three_substeps ? 0 : -1);
    probe.requested_role = probe_role;
    probe.sampled_intensive_device = sampled_intensive.data();
    probe.sampled_endpoint_device = two_moment_attached && probe_role == 0
        ? sampled_endpoint.data() : nullptr;
    probe.hit_device = sample_hit.data();
    probe.endpoint_hit_device = two_moment_attached && probe_role == 0
        ? endpoint_hit.data() : nullptr;
    probe.seed_finite_stale_value = inject_stale_write_fault;
    probe.skipped_global_component =
        (inject_stale_write_fault || omit_write_tally) ? probe_component : -1;
    probe.skip_write_tally_only = omit_write_tally;
    erf_sbm::SBMAMRStageParent parent;
    parent.conserved_old = &coarse_rho_old;
    parent.conserved_new = &coarse_rho_new;
    parent.old_time = coarse_begin;
    parent.new_time = coarse_end;
    parent.fine_step_dt = fine_dt;
    parent.refinement_ratio = IntVect(2);
    parent.ghost_probe = &probe;
    parent.componentwise_parent_fill = componentwise_parent_fill;

    if (rounded_three_substeps) {
        for (int substep = 0; substep < 3; ++substep) {
            const double step_old_time =
                coarse_begin + static_cast<double>(substep) * fine_dt;
            if (substep == 2) {
                ASSERT_NE(state_manager.new_time(1), step_old_time)
                    << "the fixture must reproduce ERF's non-bitwise third fine start";
            }
            ASSERT_TRUE(state_manager.begin_step(
                1, step_old_time, diagnostic, fine_dt)) << diagnostic;
            EXPECT_DOUBLE_EQ(state_manager.old_time(1), step_old_time);
            const double stage0_time = step_old_time + fine_dt / 3.0;
            const double stage1_time = step_old_time + fine_dt / 2.0;
            const double step_end_time = step_old_time + fine_dt;
            transport.advance_stage_from_host(
                1, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
                step_old_time, step_old_time, stage0_time, fine_dt / 3.0,
                state_manager, fine_rho_anchor, fine_rho_input,
                fine_rho_target, avg_xmom, avg_ymom, avg_zmom, fine_geom,
                1, 2, &parent);
            transport.advance_stage_from_host(
                1, erf_auxiliary::HostIntegrator::CompressibleRK3, 1,
                step_old_time, stage0_time, stage1_time, fine_dt / 2.0,
                state_manager, fine_rho_anchor, fine_rho_input,
                fine_rho_target, avg_xmom, avg_ymom, avg_zmom, fine_geom,
                1, 2, &parent);
            transport.advance_stage_from_host(
                1, erf_auxiliary::HostIntegrator::CompressibleRK3, 2,
                step_old_time, stage1_time, step_end_time, fine_dt,
                state_manager, fine_rho_anchor, fine_rho_input,
                fine_rho_target, avg_xmom, avg_ymom, avg_zmom, fine_geom,
                1, 2, &parent);
            EXPECT_TRUE(transport.projected_ledger(1).step_complete());
        }
    } else if (full_rk3_step) {
        ASSERT_TRUE(state_manager.begin_step(1, coarse_begin, diagnostic,
                                             fine_dt)) << diagnostic;
        const double stage0_time = coarse_begin + fine_dt / 3.0;
        const double stage1_time = coarse_begin + fine_dt / 2.0;
        const double step_end_time = coarse_begin + fine_dt;
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
            coarse_begin, coarse_begin, stage0_time, fine_dt / 3.0,
            state_manager, fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 1,
            coarse_begin, stage0_time, stage1_time, fine_dt / 2.0,
            state_manager, fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 2,
            coarse_begin, stage1_time, step_end_time, fine_dt,
            state_manager, fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);
    } else {
        ASSERT_TRUE(state_manager.begin_step(1, coarse_begin, diagnostic,
                                             fine_dt)) << diagnostic;
        const double stage_mid_time =
            coarse_begin + 0.5 * (coarse_end - coarse_begin);
        const double half_dt = 0.5 * (coarse_end - coarse_begin);
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
            coarse_begin, coarse_begin, stage_mid_time, half_dt,
            state_manager, fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 1,
            coarse_begin, stage_mid_time, coarse_end, half_dt,
            state_manager, fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);
    }

    Real sampled_value = Real(-1.0);
    int sample_was_hit = 0;
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     sampled_intensive.begin(), sampled_intensive.end(),
                     &sampled_value);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, sample_hit.begin(),
                     sample_hit.end(), &sample_was_hit);
    Real sampled_endpoint_value = Real(-1.0);
    int endpoint_was_hit = 0;
    if (two_moment_attached) {
        amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                         sampled_endpoint.begin(), sampled_endpoint.end(),
                         &sampled_endpoint_value);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, endpoint_hit.begin(),
                         endpoint_hit.end(), &endpoint_was_hit);
    }
    amrex::Gpu::streamSynchronize();
    EXPECT_EQ(sample_was_hit, 1)
        << "the probe cell must be a two-ghost coarse/fine donor cell";
    if (two_moment_attached && probe_role == 0) {
        EXPECT_EQ(endpoint_was_hit, 1)
            << "the probe cell must be transformed after physical G1 fill";
    }
    if (sample_result != nullptr) {
        sample_result->physical = sampled_value;
        sample_result->endpoint = sampled_endpoint_value;
        sample_result->physical_hit = sample_was_hit;
        sample_result->endpoint_hit = endpoint_was_hit;
    }
    if (two_moment_attached) {
        double requested = coarse_begin;
        if (probe_role == 0 && probe_stage == 1) {
            requested = coarse_begin + fine_dt / 3.0;
        } else if (probe_role == 0 && probe_stage == 2) {
            requested = coarse_begin + fine_dt / 2.0;
        }
        const Real theta = static_cast<Real>((requested - coarse_begin) /
                                             (coarse_end - coarse_begin));
        Real old_value = Real(0.0);
        Real new_value = Real(0.0);
        switch (probe_component) {
        case 0: old_value = Real(0.3); new_value = Real(1.2); break;
        case 1: old_value = Real(0.7); new_value = Real(1.4); break;
        case 2: old_value = Real(2.0); new_value = Real(4.0); break;
        case 3: old_value = Real(1.0); new_value = Real(2.0); break;
        case 4: old_value = Real(0.2); new_value = Real(0.8); break;
        default: old_value = Real(0.1); new_value = Real(0.4); break;
        }
        const Real rho_t = (Real(1.0) - theta) * Real(1.0) +
                           theta * Real(3.0);
        const Real expected_physical =
            ((Real(1.0) - theta) * old_value + theta * new_value) / rho_t;
        EXPECT_NEAR(sampled_value, expected_physical, Real(64.0) *
                    std::numeric_limits<Real>::epsilon());
        Real expected_endpoint = Real(-1.0);
        if (probe_role == 0) {
            expected_endpoint = expected_physical;
        }
        if (probe_role == 0 && (probe_component == 0 || probe_component == 2)) {
            const Real mass = probe_component == 0
                ? expected_physical
                : ((Real(1.0) - theta) * Real(0.3) +
                   theta * Real(1.2)) / rho_t;
            const Real number = probe_component == 2
                ? expected_physical
                : ((Real(1.0) - theta) * Real(2.0) +
                   theta * Real(4.0)) / rho_t;
            expected_endpoint = probe_component == 0
                ? (Real(0.5) * number - mass) / Real(0.4)
                : (mass - Real(0.1) * number) / Real(0.4);
        }
        if (probe_role == 0) {
            EXPECT_NEAR(sampled_endpoint_value, expected_endpoint, Real(96.0) *
                        std::numeric_limits<Real>::epsilon());
        }
        if (probe_component == 0 && probe_role == 0 && theta > Real(0.0)) {
            const Real endpoint_ratio_interpolation =
                (Real(1.0) - theta) * (Real(0.3) / Real(1.0)) +
                theta * (Real(1.2) / Real(3.0));
            const Real future_density_substitution =
                ((Real(1.0) - theta) * Real(0.3) +
                 theta * Real(1.2)) / Real(3.0);
            EXPECT_GT(std::abs(sampled_value - endpoint_ratio_interpolation),
                      Real(100.0) * std::numeric_limits<Real>::epsilon());
            EXPECT_GT(std::abs(sampled_value - future_density_substitution),
                      Real(100.0) * std::numeric_limits<Real>::epsilon());
            EXPECT_GT(std::abs(sampled_endpoint_value - sampled_value),
                      Real(100.0) * std::numeric_limits<Real>::epsilon())
                << "the physical M,N ghost must be endpoint-transformed before WENO";
        }
    } else if (rounded_three_substeps) {
        const double requested = coarse_begin + 2.0 * fine_dt;
        const Real theta = static_cast<Real>((requested - coarse_begin) /
                                             (coarse_end - coarse_begin));
        const Real expected = (Real(1.0) + Real(8.0) * theta) /
                              (Real(1.0) + Real(2.0) * theta);
        EXPECT_NEAR(sampled_value, expected, Real(64.0) *
                    std::numeric_limits<Real>::epsilon());
    } else {
        EXPECT_NEAR(sampled_value, Real(2.5), Real(32.0) *
                    std::numeric_limits<Real>::epsilon());
        const Real endpoint_ratio_interpolation =
            Real(0.5) * (Real(1.0) / Real(1.0) + Real(9.0) / Real(3.0));
        const Real future_density_substitution =
            (Real(0.5) * Real(1.0) + Real(0.5) * Real(9.0)) / Real(3.0);
        EXPECT_NE(sampled_value, endpoint_ratio_interpolation);
        EXPECT_NE(sampled_value, future_density_substitution);
    }

    const auto& fine_state = state_manager.new_state(1);
    for (int component = 0; component < layout.ncomp(); ++component) {
        if (two_moment_attached) {
            const std::array<Real, 6> initial{
                Real(0.3), Real(0.7), Real(2.0), Real(1.0),
                Real(0.2), Real(0.1)};
            const Real expected = initial[static_cast<std::size_t>(component)];
            EXPECT_EQ(fine_state.min(component), expected);
            EXPECT_EQ(fine_state.max(component), expected);
        } else {
            EXPECT_EQ(fine_state.min(component), Real(1.0));
            EXPECT_EQ(fine_state.max(component), Real(1.0));
        }
    }
}

TEST(SBMTransport, FineStageCoarseGhostUsesG1PhysicalEndpointInterpolation)
{
    run_coarse_stage_g1_ghost_test();
}

TEST(SBMTransport, RoundedThirdFineSubstepCompletesRealStageZeroAndRK3)
{
    run_coarse_stage_g1_ghost_test(false, true);
}

TEST(SBMTransport, TwoMomentAttachedG1AndEndpointGhostsCoverFaceEdgeCornerLayers)
{
    const std::array<IntVect, 6> ghost_cells{
        IntVect(1, 3, 3), IntVect(0, 3, 3),
        IntVect(1, 1, 3), IntVect(0, 0, 3),
        IntVect(1, 1, 1), IntVect(0, 0, 0)};
    for (const IntVect& cell : ghost_cells) {
        for (const int component : {0, 2, 4}) {
            SCOPED_TRACE(::testing::Message()
                         << "ghost=" << cell << " component=" << component);
            for (int stage = 0; stage < 3; ++stage) {
                for (int role = 0; role < 3; ++role) {
                    SCOPED_TRACE(::testing::Message()
                                 << "stage=" << stage << " role=" << role);
                    G1GhostSample componentwise_reference;
                    G1GhostSample packed_comparison;
                    run_coarse_stage_g1_ghost_test(
                        false, false, 0, true, cell, component,
                        &componentwise_reference, stage, role, true, true);
                    run_coarse_stage_g1_ghost_test(
                        false, false, 0, true, cell, component,
                        &packed_comparison, stage, role, true, false);
                    EXPECT_EQ(componentwise_reference.physical_hit, 1);
                    EXPECT_EQ(packed_comparison.physical_hit, 1);
                    EXPECT_EQ(componentwise_reference.endpoint_hit,
                              role == 0 ? 1 : 0);
                    EXPECT_EQ(packed_comparison.endpoint_hit,
                              role == 0 ? 1 : 0);
                    EXPECT_EQ(packed_comparison.physical,
                              componentwise_reference.physical);
                    EXPECT_EQ(packed_comparison.endpoint,
                              componentwise_reference.endpoint);
                }
            }
        }
    }
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

TEST(SBMTransport, MissingParentGhostWriteRejectsFiniteStaleSentinel)
{
    expect_amrex_abort([] { run_coarse_stage_g1_ghost_test(true, false); },
                       "left a coarse-only parent ghost component unwritten or invalid");
}

TEST(SBMTransport, MissingParentGhostWriteTallyRejectsCompleteStores)
{
    expect_amrex_abort(
        [] {
            run_coarse_stage_g1_ghost_test(
                false, false, 0, false, IntVect(1, 3, 3), 0, nullptr, -1,
                0, false, true, true);
        },
        "omitted a required parent ghost component write");
}

TEST(SBMTransport, AssignmentOnlyMissingParentStoreFailsAcrossModesRolesAndHaloDepths)
{
    const std::array<IntVect, 2> ghost_cells{
        IntVect(1, 3, 3), IntVect(0, 3, 3)};
    for (const bool two_moment_attached : {false, true}) {
        for (const int role : {0, 1}) {
            for (const IntVect& ghost : ghost_cells) {
                SCOPED_TRACE(::testing::Message()
                             << "two_moment_attached=" << two_moment_attached
                             << " role=" << role << " ghost=" << ghost);
                expect_amrex_abort(
                    [&] {
                        run_coarse_stage_g1_ghost_test(
                            true, false, 0, two_moment_attached, ghost, 0,
                            nullptr, 0, role);
                    },
                    "left a coarse-only parent ghost component unwritten or invalid");
            }
        }
    }
}

TEST(SBMTransport, ParentGhostRejectsOverflowingUOverRho)
{
    expect_amrex_abort([] { run_coarse_stage_g1_ghost_test(false, false, 1); },
                       "unrepresentable U(t)/rho_d(t) quotient");
}

TEST(SBMTransport, ParentGhostRejectsPositiveSupportUnderflow)
{
    expect_amrex_abort([] { run_coarse_stage_g1_ghost_test(false, false, 2); },
                       "unrepresentable U(t)/rho_d(t) quotient");
}

TEST(SBMTransport, ParentGhostRejectsZeroDensityEndpoint)
{
    expect_amrex_abort([] { run_coarse_stage_g1_ghost_test(false, false, 3); },
                       "endpoint data are not finite with positive density");
}

TEST(SBMTransport, ParentGhostRejectsNonfiniteSpectrumEndpoint)
{
    expect_amrex_abort([] { run_coarse_stage_g1_ghost_test(false, false, 4); },
                       "endpoint data are not finite with positive density");
}

struct SpectralLedgerResult {
    std::vector<Real> spectrum;
    std::vector<Real> projected;
    std::array<long, AMREX_SPACEDIM> face_counts{};
};

SpectralLedgerResult run_full_spectral_ledger (
    const int max_groups_per_chunk,
    const erf_sbm::MomentMode mode = erf_sbm::MomentMode::TwoMoment,
    const bool with_property = true)
{
    const Box domain(IntVect(0, 0, 0), IntVect(3, 3, 3));
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    const Geometry geom(domain, &physical, amrex::CoordSys::cartesian, periodic);
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const auto layout = make_layout(mode, with_property);
    erf_sbm::SBMStateManager manager(layout, 1);
    manager.define(0, ba, dm, 0.0);
    erf_sbm::SBMTransport transport(layout, 2, max_groups_per_chunk);
    transport.define(0, ba, dm);

    MultiFab detj(ba, dm, 1, 0);
    const auto map_ba = project_to_xy(ba);
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
    const auto& population = layout.populations().front();
    auto& initial = manager.new_state_for_initialization(0);
    initial.setVal(Real(0.0));
    initial.setVal(Real(0.18), population.mass_offset, 1, 0);
    initial.setVal(Real(0.45), population.mass_offset + 1, 1, 0);
    if (mode == erf_sbm::MomentMode::TwoMoment) {
        initial.setVal(Real(0.6), population.number_offset, 2, 0);
        if (with_property) {
            initial.setVal(Real(0.15), layout.property_offset(0), 2, 0);
        }
    }

    MultiFab avg_xmom(amrex::convert(ba, IntVect::TheDimensionVector(0)), dm, 1, 0);
    MultiFab avg_ymom(amrex::convert(ba, IntVect::TheDimensionVector(1)), dm, 1, 0);
    MultiFab avg_zmom(amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1, 0);
    avg_xmom.setVal(Real(0.1));
    avg_ymom.setVal(Real(0.0));
    avg_zmom.setVal(Real(0.0));

    constexpr double dt = 0.12;
    if (!manager.begin_step(0, 0.0, diagnostic)) {
        ADD_FAILURE() << diagnostic;
        return {};
    }
    transport.advance_stage_from_host(
        0, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
        0.0, 0.0, dt / 3.0, dt / 3.0, manager, conserved_anchor,
        conserved_input, conserved_target, avg_xmom, avg_ymom, avg_zmom,
        geom, 1, 2);
    transport.advance_stage_from_host(
        0, erf_auxiliary::HostIntegrator::CompressibleRK3, 1,
        0.0, dt / 3.0, dt / 2.0, dt / 2.0, manager, conserved_anchor,
        conserved_input, conserved_target, avg_xmom, avg_ymom, avg_zmom,
        geom, 1, 2);
    transport.advance_stage_from_host(
        0, erf_auxiliary::HostIntegrator::CompressibleRK3, 2,
        0.0, dt / 2.0, dt, dt, manager, conserved_anchor,
        conserved_input, conserved_target, avg_xmom, avg_ymom, avg_zmom,
        geom, 1, 2);

    const auto& spectral_ledger = transport.completed_spectral_ledger(0);
    const auto& projected_ledger = transport.projected_ledger(0);
    EXPECT_TRUE(spectral_ledger.step_complete());
    EXPECT_EQ(spectral_ledger.integrated_flux().nComp(), layout.ncomp());
    EXPECT_TRUE(projected_ledger.step_complete());
    EXPECT_EQ(projected_ledger.integrated_flux().nComp(), 2);

    SpectralLedgerResult result;
    for (int component = 0; component < layout.ncomp(); ++component) {
        const Real actual = spectral_ledger.integrated_flux().dir(0).norm0(component);
        const Real expected = dt * Real(0.1) * manager.new_state(0).norm0(component);
        EXPECT_NEAR(actual, expected,
                    Real(128.0) * std::numeric_limits<Real>::epsilon());
        EXPECT_GT(actual, Real(0.0));
        result.spectrum.push_back(actual);
    }
    const Real projected_expected0 = dt * Real(0.1) * Real(0.18);
    const Real projected_expected1 = dt * Real(0.1) * Real(0.45);
    EXPECT_NEAR(projected_ledger.integrated_flux().dir(0).norm0(0),
                projected_expected0,
                Real(128.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_NEAR(projected_ledger.integrated_flux().dir(0).norm0(1),
                projected_expected1,
                Real(128.0) * std::numeric_limits<Real>::epsilon());
    for (int component = 0; component < 2; ++component) {
        result.projected.push_back(
            projected_ledger.integrated_flux().dir(0).norm0(component));
    }
    for (int direction = 0; direction < AMREX_SPACEDIM; ++direction) {
        const BoxArray face_ba = amrex::convert(
            ba, IntVect::TheDimensionVector(direction));
        long count = 0;
        for (int box_index = 0; box_index < face_ba.size(); ++box_index) {
            count += face_ba[box_index].numPts();
        }
        result.face_counts[static_cast<std::size_t>(direction)] = count;
    }
    return result;
}

TEST(SBMTransport, AMRCapableCompletedLedgerRetainsEveryAcceptedComponent)
{
    const auto single_level_layout = make_layout(erf_sbm::MomentMode::TwoMoment, true);
    erf_sbm::SBMTransport single_level(single_level_layout, 1);
    const Box domain(IntVect(0, 0, 0), IntVect(3, 3, 3));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    single_level.define(0, ba, dm);
    EXPECT_THROW(static_cast<void>(single_level.completed_spectral_ledger(0)),
                 std::logic_error);

    const auto small_chunks = run_full_spectral_ledger(1);
    const auto large_chunks = run_full_spectral_ledger(16);
    ASSERT_EQ(small_chunks.spectrum.size(), large_chunks.spectrum.size());
    for (std::size_t component = 0; component < small_chunks.spectrum.size(); ++component) {
        EXPECT_NEAR(small_chunks.spectrum[component], large_chunks.spectrum[component],
                    Real(128.0) * std::numeric_limits<Real>::epsilon());
    }
    ASSERT_EQ(small_chunks.projected.size(), large_chunks.projected.size());
    for (std::size_t component = 0; component < small_chunks.projected.size(); ++component) {
        EXPECT_NEAR(small_chunks.projected[component], large_chunks.projected[component],
                    Real(128.0) * std::numeric_limits<Real>::epsilon());
    }
    EXPECT_EQ(small_chunks.face_counts[0], 80);
    EXPECT_EQ(small_chunks.face_counts[1], 80);
    EXPECT_EQ(small_chunks.face_counts[2], 80);

    const auto one_moment_small = run_full_spectral_ledger(
        1, erf_sbm::MomentMode::OneMoment, false);
    const auto one_moment_large = run_full_spectral_ledger(
        16, erf_sbm::MomentMode::OneMoment, false);
    ASSERT_EQ(one_moment_small.spectrum.size(), one_moment_large.spectrum.size());
    for (std::size_t component = 0;
         component < one_moment_small.spectrum.size(); ++component) {
        EXPECT_NEAR(one_moment_small.spectrum[component],
                    one_moment_large.spectrum[component],
                    Real(128.0) * std::numeric_limits<Real>::epsilon());
    }
    ASSERT_EQ(one_moment_small.projected.size(), one_moment_large.projected.size());
    for (std::size_t component = 0;
         component < one_moment_small.projected.size(); ++component) {
        EXPECT_NEAR(one_moment_small.projected[component],
                    one_moment_large.projected[component],
                    Real(128.0) * std::numeric_limits<Real>::epsilon());
    }
}

TEST(SBMTransport, LimiterActiveCompletedSpectralLedgerReconstructsAcceptedUpdate)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.canonical_upper_edge_profile = true;
    options.canonical_upper_edge_seed = 2;
    options.attached_property = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.55;
    options.advance_all_rk3_stages = true;
    options.expect_limiter_active = true;
    options.compare_native_candidate = true;
    options.verify_spectral_ledger_identity = true;
    options.verify_direct_register_reference = true;
    options.max_groups_per_chunk = 1;

    const auto summary = run_transport(options);
    const auto layout = make_layout(options.mode, options.attached_property);
    ASSERT_EQ(summary.spectral_ledger_relative_errors.size(),
              static_cast<std::size_t>(layout.ncomp()));
    for (std::size_t component = 0;
         component < summary.spectral_ledger_relative_errors.size(); ++component) {
        EXPECT_LE(summary.spectral_ledger_relative_errors[component],
                  Real(512.0) * std::numeric_limits<Real>::epsilon())
            << "accepted-ledger component=" << component;
    }
}

TEST(SBMTransport, HostSlowCopyPreservesAcceptedLiquidProjection)
{
    run_host_slow_copy_preserves_accepted_liquid_projection_test();
}

TEST(SBMTransport, AttachedPropertiesStayInTheirAtomicSupportGroup)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.discontinuity = true;
    options.attached_property = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.6;
    run_transport(options);
}

TEST(SBMTransport, CompleteGroupChunkPoliciesAreInvariant)
{
    RunOptions one_group_options;
    one_group_options.mode = erf_sbm::MomentMode::TwoMoment;
    one_group_options.canonical_upper_edge_profile = true;
    one_group_options.canonical_upper_edge_seed = 2;
    one_group_options.attached_property = true;
    one_group_options.carrier_x = Real(0.2);
    one_group_options.multidirectional_carrier = true;
    one_group_options.max_groups_per_chunk = 1;
    one_group_options.dt = 0.3;
    one_group_options.compare_native_candidate = true;
    one_group_options.expect_limiter_active = true;
    RunOptions many_groups_options = one_group_options;
    many_groups_options.max_groups_per_chunk = 16;
    const auto one_group = run_transport(one_group_options);
    const auto many_groups = run_transport(many_groups_options);
    ASSERT_EQ(one_group.inventory.size(), many_groups.inventory.size());
    for (std::size_t component = 0; component < one_group.inventory.size();
         ++component) {
        EXPECT_NEAR(one_group.inventory[component],
                    many_groups.inventory[component],
                    Real(256.0) * std::numeric_limits<Real>::epsilon());
        EXPECT_NEAR(one_group.weighted_inventory[component],
                    many_groups.weighted_inventory[component],
                    Real(256.0) * std::numeric_limits<Real>::epsilon());
        EXPECT_NEAR(one_group.squared_inventory[component],
                    many_groups.squared_inventory[component],
                    Real(256.0) * std::numeric_limits<Real>::epsilon());
    }
}

TEST(SBMTransport, InactiveLimiterMatchesNativeWENOZ3)
{
    RunOptions one_moment;
    one_moment.mode = erf_sbm::MomentMode::OneMoment;
    one_moment.carrier_x = Real(0.2);
    one_moment.smooth_profile = true;
    one_moment.compare_native_candidate = true;
    run_transport(one_moment);

    RunOptions two_moment = one_moment;
    two_moment.mode = erf_sbm::MomentMode::TwoMoment;
    run_transport(two_moment);
}

TEST(SBMTransport, HistoricalP2CartesianEndpointWENOParity)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.carrier_x = Real(0.2);
    options.dt = 0.01;
    options.smooth_profile = true;
    options.compare_historical_p2 = true;
    run_transport(options);
}

TEST(SBMTransport, DirectMomentWENOIsNotEndpointWENO)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.varying_density = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.3;
    options.smooth_profile = true;
    options.compare_native_candidate = true;
    options.compare_direct_moment_candidate = true;
    run_transport(options);
}

TEST(SBMTransport, ActiveLimiterRestrictsNativeWENOZ3Proposal)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::OneMoment;
    options.adversarial_one_moment_profile = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.55;
    options.compare_native_candidate = true;
    options.expect_limiter_active = true;
    options.anelastic_heun = true;
    options.verify_direct_register_reference = true;
    run_transport(options);
}

TEST(SBMTransport, TwoMomentAtomicGroupLimiterRestrictsNativeProposal)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    options.complex_profile = true;
    options.multidirectional_carrier = true;
    options.carrier_x = Real(0.13);
    options.dt = 0.3;
    options.compare_native_candidate = true;
    options.expect_limiter_active = true;
    run_transport(options);
}

TEST(SBMTransport, InteriorTwoMomentCanonicalReserveLimitsClosedEdgeSaturation)
{
    constexpr Real cells_per_unit_length = static_cast<Real>(transport_nx);
    constexpr int seed = 2;
    constexpr Real carrier = Real(0.2);
    constexpr Real timestep = Real(0.55);
    const Real courant = carrier * timestep * cells_per_unit_length;
    ASSERT_LT(courant, Real(1.0));
    // With constant density and unit mapped measure, donor low order is a
    // convex combination of positive endpoint counts and remains canonical.
    for (int i = 0; i < transport_nx; ++i) {
        const Real current = canonical_upper_edge_probe_low_endpoint(i, seed);
        // The +x upwind donor wraps from cell zero to the final periodic cell.
        const int donor_index = (i + transport_nx - 1) % transport_nx;
        const Real donor = canonical_upper_edge_probe_low_endpoint(donor_index,
                                                                   seed);
        const Real donor_low = (Real(1.0) - courant) * current +
                               courant * donor;
        EXPECT_GT(donor_low, Real(0.0)) << "cell=" << i;
    }

    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.canonical_upper_edge_profile = true;
    options.canonical_upper_edge_seed = seed;
    options.attached_property = true;
    options.carrier_x = carrier;
    options.dt = timestep;
    options.compare_native_candidate = true;
    options.expect_limiter_active = true;
    const auto summary = run_transport(options);
    const Real eta =
        Real(128.0) * std::numeric_limits<Real>::epsilon();
    EXPECT_LT(summary.native_candidate_interior_upper_gap_fraction, Real(0.0));
    EXPECT_GT(summary.minimum_interior_number, Real(0.0));
    EXPECT_GT(summary.minimum_interior_upper_gap_fraction, Real(0.0));
    EXPECT_GE(summary.minimum_interior_upper_gap_fraction, eta * Real(0.5));
}

TEST(SBMTransport, RepeatedNonuniformTwoMomentTransportStaysCanonical)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.canonical_upper_edge_profile = true;
    options.canonical_upper_edge_seed = 2;
    options.attached_property = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.15;
    options.advance_all_rk3_stages = true;
    options.completed_steps = 3;
    options.max_groups_per_chunk = 1;
    options.expect_limiter_active = true;
    options.compare_native_candidate = true;
    run_transport(options);
}

TEST(SBMTransport, HeunLimiterUsesFullDtTrialWhenActive)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    options.complex_profile = true;
    options.multidirectional_carrier = true;
    options.carrier_x = Real(0.13);
    options.dt = 0.45;
    options.compare_native_candidate = true;
    options.expect_limiter_active = true;
    options.anelastic_heun = true;
    options.verify_direct_register_reference = true;
    run_transport(options);
}

TEST(SBMTransport, CanonicalGateRejectsClosedLinearSharedEdge)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    const Box domain(IntVect(0, 0, 0), IntVect(0, 0, 0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab candidate(ba, dm, layout.ncomp(), 0);
    candidate.setVal(Real(0.0));
    const auto& population = layout.populations().front();
    candidate.setVal(Real(0.01), population.number_offset, 1, 0);
    candidate.setVal(Real(0.005), population.mass_offset, 1, 0);
    std::string diagnostic;
    EXPECT_FALSE(erf_sbm::authoritative_state_admissible(
        candidate, layout, 0, &diagnostic));
    EXPECT_NE(diagnostic.find("constraint=canonical-two-moment-bin-state"),
              std::string::npos);
}

TEST(SBMTransport, CanonicalAdmissionAllowsFinalGlobalUpperEdge)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    const Box domain(IntVect(0, 0, 0), IntVect(0, 0, 0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab candidate(ba, dm, layout.ncomp(), 0);
    candidate.setVal(Real(0.0));
    const auto& population = layout.populations().front();
    const Real number = Real(0.01);
    candidate.setVal(number, population.number_offset + 1, 1, 0);
    candidate.setVal(population.grid.edges().back() * number,
                     population.mass_offset + 1, 1, 0);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::authoritative_state_admissible(
        candidate, layout, 0, &diagnostic))
        << diagnostic;
}

TEST(SBMTransport, CanonicalAdmissionAllowsInteriorLowerEdge)
{
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment);
    const Box domain(IntVect(0, 0, 0), IntVect(0, 0, 0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab candidate(ba, dm, layout.ncomp(), 0);
    candidate.setVal(Real(0.0));
    const auto& population = layout.populations().front();
    const Real number = Real(0.01);
    candidate.setVal(number, population.number_offset, 1, 0);
    candidate.setVal(population.grid.edges().front() * number,
                     population.mass_offset, 1, 0);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::authoritative_state_admissible(
        candidate, layout, 0, &diagnostic))
        << diagnostic;
}

TEST(SBMTransport, PhysicalAdvectionFacePolicyIsFailClosedAndAtomic)
{
    using erf_sbm::AdvectionFaceAction;
    using erf_sbm::AdvectionFaceKind;
    using erf_sbm::AdvectionFacePolicyInput;
    using erf_sbm::make_advection_face_policy;

    const auto periodic = make_advection_face_policy(
        {AdvectionFaceKind::Periodic, Real(0.3), false, false, true});
    EXPECT_EQ(periodic.action, AdvectionFaceAction::PeriodicSharedFace);
    EXPECT_TRUE(periodic.use_high_order_candidate);

    const auto wall = make_advection_face_policy(
        {AdvectionFaceKind::ImpermeableWall, Real(-0.4), false, false, true});
    EXPECT_EQ(wall.action, AdvectionFaceAction::ZeroNormalFlux);
    EXPECT_FALSE(wall.use_interior_donor);
    EXPECT_FALSE(wall.use_explicit_spectral_state);
    EXPECT_FALSE(wall.use_high_order_candidate);

    const auto zero_carrier = make_advection_face_policy(
        {AdvectionFaceKind::AdvectiveOutflow, Real(0.0), true, true, true});
    EXPECT_EQ(zero_carrier.action, AdvectionFaceAction::ZeroNormalFlux);
    EXPECT_TRUE(zero_carrier.accepted());
    EXPECT_FALSE(zero_carrier.use_interior_donor);
    EXPECT_FALSE(zero_carrier.use_explicit_spectral_state);
    EXPECT_FALSE(zero_carrier.use_high_order_candidate);

    const auto outward = make_advection_face_policy(
        {AdvectionFaceKind::AdvectiveOutflow, Real(0.4), false, false, true});
    EXPECT_EQ(outward.action, AdvectionFaceAction::InteriorDonor);
    EXPECT_TRUE(outward.use_interior_donor);
    EXPECT_TRUE(outward.use_high_order_candidate);

    const auto inward_outflow = make_advection_face_policy(
        {AdvectionFaceKind::AdvectiveOutflow, Real(-0.2), false, false, false});
    EXPECT_EQ(inward_outflow.action,
              AdvectionFaceAction::RejectInwardOutflow);
    EXPECT_FALSE(inward_outflow.accepted());
    EXPECT_FALSE(inward_outflow.use_interior_donor);
    EXPECT_FALSE(inward_outflow.use_explicit_spectral_state);
    EXPECT_FALSE(inward_outflow.use_high_order_candidate);

    const auto inward_outflow_with_complete_state = make_advection_face_policy(
        {AdvectionFaceKind::AdvectiveOutflow, Real(-0.2), true, true, true});
    EXPECT_EQ(inward_outflow_with_complete_state.action,
              AdvectionFaceAction::RejectInwardOutflow);
    EXPECT_FALSE(inward_outflow_with_complete_state.accepted());
    EXPECT_FALSE(inward_outflow_with_complete_state.use_interior_donor);
    EXPECT_FALSE(inward_outflow_with_complete_state.use_explicit_spectral_state);
    EXPECT_FALSE(inward_outflow_with_complete_state.use_high_order_candidate);

    const auto prescribed = make_advection_face_policy(
        {AdvectionFaceKind::PrescribedInflow, Real(-0.2), true, true, true});
    EXPECT_EQ(prescribed.action, AdvectionFaceAction::ExplicitSpectralInflow);
    EXPECT_TRUE(prescribed.use_explicit_spectral_state);
    EXPECT_TRUE(prescribed.use_high_order_candidate);

    const auto missing_prescribed_state = make_advection_face_policy(
        {AdvectionFaceKind::PrescribedInflow, Real(-0.2), false, false, true});
    EXPECT_EQ(missing_prescribed_state.action,
              AdvectionFaceAction::RejectMissingInflowState);
    EXPECT_FALSE(missing_prescribed_state.accepted());
    EXPECT_FALSE(missing_prescribed_state.use_interior_donor);
    EXPECT_FALSE(missing_prescribed_state.use_explicit_spectral_state);
    EXPECT_FALSE(missing_prescribed_state.use_high_order_candidate);

    const auto incomplete_group = make_advection_face_policy(
        {AdvectionFaceKind::PrescribedInflow, Real(-0.2), true, false, true});
    EXPECT_EQ(incomplete_group.action,
              AdvectionFaceAction::RejectIncompleteInflowState);
    EXPECT_FALSE(incomplete_group.accepted());

    const auto donor_fallback = make_advection_face_policy(
        {AdvectionFaceKind::AdvectiveOutflow, Real(0.4), false, false, false});
    EXPECT_EQ(donor_fallback.action, AdvectionFaceAction::InteriorDonor);
    EXPECT_TRUE(donor_fallback.use_interior_donor);
    EXPECT_FALSE(donor_fallback.use_high_order_candidate);
}

TEST(SBMTransport, HeunCorrectorUsesFullTrialIntervalAndHalfRecurrenceWeight)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::OneMoment;
    options.discontinuity = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.6;
    options.compare_native_candidate = true;
    options.anelastic_heun = true;
    run_transport(options);
}

TEST(SBMTransport, HeunInputTrialUsesItsMatchingPredictorDensity)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.smooth_profile = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.04;
    options.compare_native_candidate = true;
    options.anelastic_heun = true;
    options.distinct_density_roles = true;
    options.rho_anchor_slope = Real(0.15);
    options.rho_input_slope = Real(0.15);
    options.rho_target_slope = Real(0.5);
    run_transport(options);
}

TEST(SBMTransport, StaticMappedMeasureConservesHAcrossJacobianAndMapFactors)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::OneMoment;
    options.discontinuity = true;
    options.carrier_x = Real(0.2);
    options.mapped_geometry = true;
    run_transport(options);
}

TEST(SBMTransport, CompressibleCorrectorsKeepTheOldStepLimiterBase)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.discontinuity = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.3;
    options.advance_all_rk3_stages = true;
    run_transport(options);
}

TEST(SBMTransport, CompressiblePredictorCannotReplaceAnchorLimiterBase)
{
    RunOptions options;
    options.mode = erf_sbm::MomentMode::OneMoment;
    options.discontinuity = true;
    options.carrier_x = Real(0.2);
    options.dt = 0.3;
    options.advance_all_rk3_stages = true;
    options.amplify_corrector_input = true;
    run_transport(options);
}


TEST(SBMTransport, FineStageG1WenoAcceptedRateReachesMappedYAOracle)
{
#if defined(AMREX_USE_FLOAT)
    GTEST_SKIP() << "The coupled WENO/FCT/YA separations are qualified in "
                    "DOUBLE; SINGLE cannot resolve them.";
#else
    ASSERT_EQ(AMREX_SPACEDIM, 3);
    const auto layout = make_layout(erf_sbm::MomentMode::TwoMoment, true);
    const auto& population = layout.populations().front();
    const int mass0 = population.mass_offset;
    const int number0 = population.number_offset;
    const int property0 = layout.property_offset(0);
    ASSERT_GE(property0, 0);

    const IntVect ratio(2);
    const Box coarse_domain(IntVect(0), IntVect(3));
    const Box fine_domain = amrex::refine(coarse_domain, ratio);
    const BoxArray coarse_ba(coarse_domain);
    const BoxArray fine_ba(amrex::refine(Box(IntVect(1), IntVect(2)), ratio));
    const DistributionMapping coarse_dm(coarse_ba);
    const DistributionMapping fine_dm(fine_ba);
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    const Geometry coarse_geom(coarse_domain, &physical,
                               amrex::CoordSys::cartesian, periodic);
    const Geometry fine_geom(fine_domain, &physical,
                             amrex::CoordSys::cartesian, periodic);
    constexpr double t0 = 0.1;
    constexpr double t1 = 0.13;
    constexpr int fine_steps = 3;
    const double fine_dt = (t1 - t0) / fine_steps;

    erf_sbm::SBMStateManager state_manager(layout, 2);
    std::string diagnostic;
    state_manager.define(0, coarse_ba, coarse_dm, t0);
    state_manager.define(1, fine_ba, fine_dm, t0);
    MultiFab rho_old(coarse_ba, coarse_dm, 3, 0);
    MultiFab rho_new(coarse_ba, coarse_dm, 3, 0);
    MultiFab unused_density(fine_ba, fine_dm, 3, 0);
    rho_old.setVal(Real(0.0));
    rho_new.setVal(Real(0.0));
    unused_density.setVal(Real(0.0));

    auto fill_state = [&] (MultiFab& state, MultiFab& density,
                           const bool fine, const Real time) {
        fill_a4_fixture_state(state, density, fine, time,
                              mass0, number0, property0);
    };

    fill_state(state_manager.new_state_for_initialization(0), rho_old,
               false, Real(0.0));
    ASSERT_TRUE(state_manager.begin_step(0, t0, diagnostic)) << diagnostic;
    fill_state(state_manager.new_target_storage(0), rho_new,
               false, Real(1.0));
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
        state_manager.new_target_storage(0), layout, 0, &diagnostic))
        << diagnostic;
    ASSERT_TRUE(state_manager.accept_stage_target(0, t1, true, diagnostic))
        << diagnostic;
    fill_state(state_manager.new_state_for_initialization(1), unused_density,
               true, Real(0.0));

    MultiFab fine_rho_anchor(fine_ba, fine_dm, 3, 0);
    MultiFab fine_rho_input(fine_ba, fine_dm, 3, 0);
    MultiFab fine_rho_target(fine_ba, fine_dm, 3, 0);
    for (MultiFab* rho : {&fine_rho_anchor, &fine_rho_input,
                          &fine_rho_target}) {
        rho->setVal(Real(0.0));
        rho->setVal(Real(1.0), Rho_comp, 1, 0);
    }

    erf_sbm::SBMTransport transport(layout, 2, 16, false);
    transport.define(0, coarse_ba, coarse_dm);
    transport.define(1, fine_ba, fine_dm);
    const auto fine_map_ba = project_to_xy(fine_ba);
    MultiFab detj(fine_ba, fine_dm, 1, 0);
    MultiFab mx(fine_map_ba, fine_dm, 1, 0);
    MultiFab my(fine_map_ba, fine_dm, 1, 0);
    detj.setVal(Real(1.0));
    mx.setVal(Real(1.0));
    my.setVal(Real(1.0));
    ASSERT_TRUE(transport.rebuild_static_measure(
        1, detj, mx, my, diagnostic)) << diagnostic;

    const IntVect face(2, 3, 3);
    MultiFab avg_xmom(amrex::convert(
                          fine_ba, IntVect::TheDimensionVector(0)),
                      fine_dm, 1, 0);
    MultiFab avg_ymom(amrex::convert(
                          fine_ba, IntVect::TheDimensionVector(1)),
                      fine_dm, 1, 0);
    MultiFab avg_zmom(amrex::convert(
                          fine_ba, IntVect::TheDimensionVector(2)),
                      fine_dm, 1, 0);
    avg_xmom.setVal(Real(0.0));
    avg_ymom.setVal(Real(0.0));
    avg_zmom.setVal(Real(0.0));
    constexpr Real velocity = Real(0.001);
    set_a4_x_face_velocity(avg_xmom, face, velocity);

    amrex::Gpu::DeviceVector<Real> physical_stencil(9, Real(-1.0));
    amrex::Gpu::DeviceVector<Real> endpoint_stencil(9, Real(-1.0));
    amrex::Gpu::DeviceVector<Real> secondary_physical_stencil(9, Real(-1.0));
    amrex::Gpu::DeviceVector<Real> secondary_endpoint_stencil(9, Real(-1.0));
    amrex::Gpu::DeviceVector<Real> accepted_rate(9, Real(-1.0));
    amrex::Gpu::DeviceVector<Real> accepted_lambda(1, Real(-1.0));
    amrex::Gpu::DeviceVector<int> accepted_hit(1, 0);
    amrex::Gpu::DeviceVector<Real> secondary_rate(9, Real(-1.0));
    amrex::Gpu::DeviceVector<Real> secondary_lambda(1, Real(-1.0));
    amrex::Gpu::DeviceVector<int> secondary_hit(1, 0);
    erf_sbm::SBMAMRStageGhostProbe probe;
    probe.sampled_face = face;
    probe.sampled_face_direction = 0;
    probe.sampled_stage = 2;
    probe.sampled_mass_component = mass0;
    probe.sampled_number_component = number0;
    probe.sampled_property_component = property0;
    probe.sampled_physical_stencil_device = physical_stencil.data();
    probe.sampled_endpoint_stencil_device = endpoint_stencil.data();
    probe.sampled_mass_rate_device = accepted_rate.data();
    probe.sampled_number_rate_device = accepted_rate.data() + 1;
    probe.sampled_property_rate_device = accepted_rate.data() + 2;
    probe.sampled_high_mass_rate_device = accepted_rate.data() + 3;
    probe.sampled_high_number_rate_device = accepted_rate.data() + 4;
    probe.sampled_high_property_rate_device = accepted_rate.data() + 5;
    probe.sampled_low_mass_rate_device = accepted_rate.data() + 6;
    probe.sampled_low_number_rate_device = accepted_rate.data() + 7;
    probe.sampled_low_property_rate_device = accepted_rate.data() + 8;
    probe.sampled_lambda_device = accepted_lambda.data();
    probe.sampled_rate_hit_device = accepted_hit.data();
    erf_sbm::SBMAMRStageGhostProbe secondary_probe;
    secondary_probe.sampled_face = IntVect(3, 3, 3);
    secondary_probe.sampled_face_direction = 0;
    secondary_probe.sampled_stage = 2;
    secondary_probe.sampled_mass_component = mass0;
    secondary_probe.sampled_number_component = number0;
    secondary_probe.sampled_property_component = property0;
    secondary_probe.sampled_physical_stencil_device =
        secondary_physical_stencil.data();
    secondary_probe.sampled_endpoint_stencil_device =
        secondary_endpoint_stencil.data();
    secondary_probe.sampled_mass_rate_device = secondary_rate.data();
    secondary_probe.sampled_number_rate_device = secondary_rate.data() + 1;
    secondary_probe.sampled_property_rate_device = secondary_rate.data() + 2;
    secondary_probe.sampled_high_mass_rate_device = secondary_rate.data() + 3;
    secondary_probe.sampled_high_number_rate_device = secondary_rate.data() + 4;
    secondary_probe.sampled_high_property_rate_device = secondary_rate.data() + 5;
    secondary_probe.sampled_low_mass_rate_device = secondary_rate.data() + 6;
    secondary_probe.sampled_low_number_rate_device = secondary_rate.data() + 7;
    secondary_probe.sampled_low_property_rate_device = secondary_rate.data() + 8;
    secondary_probe.sampled_lambda_device = secondary_lambda.data();
    secondary_probe.sampled_rate_hit_device = secondary_hit.data();
    erf_sbm::SBMAMRStageParent parent;
    parent.conserved_old = &rho_old;
    parent.conserved_new = &rho_new;
    parent.old_time = t0;
    parent.new_time = t1;
    parent.fine_step_dt = fine_dt;
    parent.refinement_ratio = ratio;
    parent.ghost_probe = &probe;
    parent.secondary_ghost_probe = &secondary_probe;

    erf_sbm::SBMAMRFluxRegister ya;
    ya.define(fine_ba, coarse_ba, fine_dm, coarse_dm, fine_geom, coarse_geom,
              ratio, 1, layout.ncomp());
    ASSERT_TRUE(ya.begin_interval(t0, t1, fine_steps, diagnostic)) << diagnostic;
    transport.set_amr_flux_registers(1, nullptr, &ya);

    std::vector<int> components(static_cast<std::size_t>(layout.ncomp()));
    for (int n = 0; n < layout.ncomp(); ++n)
        components[static_cast<std::size_t>(n)] = n;
    erf_auxiliary::MappedFaceFluxRate zero_coarse;
    zero_coarse.define(coarse_ba, coarse_dm, layout.ncomp(), 0);
    zero_coarse.setVal(Real(0.0));
    ASSERT_TRUE(ya.begin_level_step(true, t0, t1, diagnostic)) << diagnostic;
    for (int stage = 0; stage < 3; ++stage) {
        ASSERT_TRUE(ya.begin_stage(
            true, erf_auxiliary::HostIntegrator::CompressibleRK3,
            stage, t0, diagnostic)) << diagnostic;
        ASSERT_TRUE(ya.add_stage_rates(
            true, zero_coarse, components, stage == 2 ? t1 - t0 : 0.0,
            diagnostic)) << diagnostic;
        ASSERT_TRUE(ya.finish_stage(true, diagnostic)) << diagnostic;
    }
    ASSERT_TRUE(ya.accept_level_step(true, t0, t1, diagnostic)) << diagnostic;

    std::array<Real, 3> fine_integral{Real(0.0), Real(0.0), Real(0.0)};
    std::array<Real, 3> ratio_wrong_integral{
        Real(0.0), Real(0.0), Real(0.0)};
    Real interface_lambda = Real(1.0);
    const Real lower = population.grid.edges()[0];
    const Real upper = population.grid.edges()[1];
    const Real width = upper - lower;
    for (int substep = 0; substep < fine_steps; ++substep) {
        const double begin = t0 + substep * fine_dt;
        const double end = begin + fine_dt;
        ASSERT_TRUE(state_manager.begin_step(1, begin, diagnostic, fine_dt))
            << diagnostic;
        ASSERT_TRUE(transport.begin_physical_step(1, begin, end, diagnostic))
            << diagnostic;
        const double stage0 = begin + fine_dt / 3.0;
        const double stage1 = begin + fine_dt / 2.0;
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
            begin, begin, stage0, fine_dt / 3.0, state_manager,
            fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 1,
            begin, stage0, stage1, fine_dt / 2.0, state_manager,
            fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);
        const int clear_hit = 0;
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, &clear_hit,
                         &clear_hit + 1, accepted_hit.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, &clear_hit,
                         &clear_hit + 1, secondary_hit.begin());
        transport.advance_stage_from_host(
            1, erf_auxiliary::HostIntegrator::CompressibleRK3, 2,
            begin, stage1, end, fine_dt, state_manager,
            fine_rho_anchor, fine_rho_input, fine_rho_target,
            avg_xmom, avg_ymom, avg_zmom, fine_geom, 1, 2, &parent);

        std::array<Real, 9> physical_values{};
        std::array<Real, 9> endpoint_values{};
        std::array<Real, 9> secondary_physical_values{};
        std::array<Real, 9> secondary_endpoint_values{};
        std::array<Real, 9> rate_values{};
        std::array<Real, 9> secondary_rate_values{};
        Real lambda = Real(-1.0);
        Real secondary_limiter = Real(-1.0);
        int hit = 0;
        int secondary_hit_value = 0;
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, physical_stencil.begin(),
                         physical_stencil.end(), physical_values.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, endpoint_stencil.begin(),
                         endpoint_stencil.end(), endpoint_values.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, accepted_rate.begin(),
                         accepted_rate.end(), rate_values.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                         secondary_physical_stencil.begin(),
                         secondary_physical_stencil.end(),
                         secondary_physical_values.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                         secondary_endpoint_stencil.begin(),
                         secondary_endpoint_stencil.end(),
                         secondary_endpoint_values.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, secondary_rate.begin(),
                         secondary_rate.end(), secondary_rate_values.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, accepted_lambda.begin(),
                         accepted_lambda.end(), &lambda);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, secondary_lambda.begin(),
                         secondary_lambda.end(), &secondary_limiter);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, accepted_hit.begin(),
                         accepted_hit.end(), &hit);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, secondary_hit.begin(),
                         secondary_hit.end(), &secondary_hit_value);
        amrex::Gpu::streamSynchronize();
        ASSERT_EQ(hit, 1) << "selected accepted interface rate was not captured";
        ASSERT_EQ(secondary_hit_value, 1)
            << "interior WENO discriminator was not captured";
        EXPECT_LT(lambda, Real(1.0))
            << "the selected interface face must exercise group FCT limiting";
        interface_lambda = lambda;

        const Real theta = static_cast<Real>((stage1 - t0) / (t1 - t0));
        const A4FixturePoint oldp =
            a4_fixture_point(Real(0.0), Real(1.0), Real(1.0), Real(0.0));
        const A4FixturePoint newp =
            a4_fixture_point(Real(0.0), Real(1.0), Real(1.0), Real(1.0));
        const Real rho_t = (Real(1.0) - theta) * oldp.rho +
                           theta * newp.rho;
        const auto interpolate_physical = [=] (Real a, Real b) {
            return ((Real(1.0) - theta) * oldp.rho * a +
                    theta * newp.rho * b) / rho_t;
        };
        const Real donor_m = interpolate_physical(oldp.mass0, newp.mass0);
        const Real donor_n = interpolate_physical(oldp.number0, newp.number0);
        const Real donor_p = interpolate_physical(oldp.property0, newp.property0);
        const Real donor_l = (upper * donor_n - donor_m) / width;
        const Real donor_h = (donor_m - lower * donor_n) / width;
        const Real endpoint_tol = Real(512.0) *
                                  std::numeric_limits<Real>::epsilon();
        EXPECT_NEAR(physical_values[0], donor_m, endpoint_tol);
        EXPECT_NEAR(physical_values[1], donor_n, endpoint_tol);
        EXPECT_NEAR(physical_values[2], donor_p, endpoint_tol);
        EXPECT_NEAR(physical_values[3], donor_m, endpoint_tol);
        EXPECT_NEAR(endpoint_values[0], donor_l, endpoint_tol);
        EXPECT_NEAR(endpoint_values[1], donor_h, endpoint_tol);
        EXPECT_NEAR(endpoint_values[3], donor_l, endpoint_tol);
        EXPECT_NEAR(endpoint_values[4], donor_h, endpoint_tol);

        const Real flux_l = a4_weno_z3_left(
            endpoint_values[0], endpoint_values[3], endpoint_values[6]);
        const Real flux_h = a4_weno_z3_left(
            endpoint_values[1], endpoint_values[4], endpoint_values[7]);
        const Real expected_high_m =
            velocity * (lower * flux_l + upper * flux_h);
        const Real expected_high_n = velocity * (flux_l + flux_h);
        const Real expected_high_p = velocity * a4_weno_z3_left(
            physical_values[2], physical_values[5], physical_values[8]);
        const Real base_theta = static_cast<Real>((begin - t0) / (t1 - t0));
        const Real base_rho = (Real(1.0) - base_theta) * oldp.rho +
                              base_theta * newp.rho;
        const auto base_physical = [=] (Real a, Real b) {
            return ((Real(1.0) - base_theta) * oldp.rho * a +
                    base_theta * newp.rho * b) / base_rho;
        };
        const std::array<Real, 3> expected_low{
            velocity * base_physical(oldp.mass0, newp.mass0),
            velocity * base_physical(oldp.number0, newp.number0),
            velocity * base_physical(oldp.property0, newp.property0)};
        const Real rate_tol = Real(2048.0) *
                              std::numeric_limits<Real>::epsilon();
        EXPECT_NEAR(rate_values[3], expected_high_m, rate_tol);
        EXPECT_NEAR(rate_values[4], expected_high_n, rate_tol);
        EXPECT_NEAR(rate_values[5], expected_high_p, rate_tol);
        for (int n = 0; n < 3; ++n) {
            EXPECT_NEAR(rate_values[6 + n], expected_low[
                static_cast<std::size_t>(n)], rate_tol);
            EXPECT_NEAR(rate_values[n], rate_values[6 + n] + lambda *
                        (rate_values[3 + n] - rate_values[6 + n]), rate_tol)
                << "accepted flux must use one shared complete-group lambda";
        }

        const Real wrong_m = (Real(1.0) - theta) * oldp.mass0 +
                             theta * newp.mass0;
        const Real wrong_n = (Real(1.0) - theta) * oldp.number0 +
                             theta * newp.number0;
        const Real wrong_p = (Real(1.0) - theta) * oldp.property0 +
                             theta * newp.property0;
        const Real wrong_l = (upper * wrong_n - wrong_m) / width;
        const Real wrong_h = (wrong_m - lower * wrong_n) / width;
        const Real wrong_ratio_m = velocity *
            (lower * a4_weno_z3_left(wrong_l, wrong_l, endpoint_values[6]) +
             upper * a4_weno_z3_left(wrong_h, wrong_h, endpoint_values[7]));
        const Real wrong_ratio_n = velocity *
            (a4_weno_z3_left(wrong_l, wrong_l, endpoint_values[6]) +
             a4_weno_z3_left(wrong_h, wrong_h, endpoint_values[7]));
        const Real wrong_ratio_p = velocity *
            a4_weno_z3_left(wrong_p, wrong_p, physical_values[8]);
        EXPECT_GT(std::abs(expected_high_m - wrong_ratio_m),
                  Real(100.0) * rate_tol);
        EXPECT_GT(std::abs(expected_high_n - wrong_ratio_n),
                  Real(100.0) * rate_tol);
        EXPECT_GT(std::abs(expected_high_p - wrong_ratio_p),
                  Real(100.0) * rate_tol);

        // One fine cell inside the interface has a non-degenerate native
        // stencil containing a coarse-only parent ghost and two distinct fine
        // cells. It separates physical-M/N-before-endpoint WENO from the
        // required endpoint-coordinate WENO independently of the interface
        // limiter, which is deliberately active at the boundary face.
        const Real secondary_flux_l = a4_weno_z3_left(
            secondary_endpoint_values[0], secondary_endpoint_values[3],
            secondary_endpoint_values[6]);
        const Real secondary_flux_h = a4_weno_z3_left(
            secondary_endpoint_values[1], secondary_endpoint_values[4],
            secondary_endpoint_values[7]);
        const Real secondary_high_m = velocity *
            (lower * secondary_flux_l + upper * secondary_flux_h);
        const Real secondary_high_n = velocity *
            (secondary_flux_l + secondary_flux_h);
        const Real secondary_high_p = velocity * a4_weno_z3_left(
            secondary_physical_values[2], secondary_physical_values[5],
            secondary_physical_values[8]);
        const Real secondary_wrong_ratio_m = velocity *
            (lower * a4_weno_z3_left(wrong_l,
                                     secondary_endpoint_values[3],
                                     secondary_endpoint_values[6]) +
             upper * a4_weno_z3_left(wrong_h,
                                     secondary_endpoint_values[4],
                                     secondary_endpoint_values[7]));
        const Real secondary_wrong_ratio_n = velocity *
            (a4_weno_z3_left(wrong_l, secondary_endpoint_values[3],
                             secondary_endpoint_values[6]) +
             a4_weno_z3_left(wrong_h, secondary_endpoint_values[4],
                             secondary_endpoint_values[7]));
        const Real secondary_wrong_ratio_p = velocity * a4_weno_z3_left(
            wrong_p, secondary_physical_values[5], secondary_physical_values[8]);
        EXPECT_NEAR(secondary_rate_values[3], secondary_high_m, rate_tol);
        EXPECT_NEAR(secondary_rate_values[4], secondary_high_n, rate_tol);
        EXPECT_NEAR(secondary_rate_values[5], secondary_high_p, rate_tol);
        for (int n = 0; n < 3; ++n) {
            EXPECT_NEAR(secondary_rate_values[n],
                secondary_rate_values[6 + n] + secondary_limiter *
                    (secondary_rate_values[3 + n] -
                     secondary_rate_values[6 + n]), rate_tol);
        }
        const Real physical_mass_weno = velocity * a4_weno_z3_left(
            secondary_physical_values[0], secondary_physical_values[3],
            secondary_physical_values[6]);
        const Real physical_number_weno = velocity * a4_weno_z3_left(
            secondary_physical_values[1], secondary_physical_values[4],
            secondary_physical_values[7]);
        EXPECT_GT(std::abs(secondary_high_m - physical_mass_weno),
                  Real(100.0) * rate_tol);
        EXPECT_GT(std::abs(secondary_high_n - physical_number_weno),
                  Real(100.0) * rate_tol);
        EXPECT_GT(std::abs(secondary_high_m - secondary_wrong_ratio_m),
                  Real(100.0) * rate_tol);
        EXPECT_GT(std::abs(secondary_high_n - secondary_wrong_ratio_n),
                  Real(100.0) * rate_tol);
        EXPECT_GT(std::abs(secondary_high_p - secondary_wrong_ratio_p),
                  Real(100.0) * rate_tol);
        std::cout << std::setprecision(17)
                  << "SBM_A4_WENO substep=" << substep
                  << " interface_lambda=" << lambda
                  << " interface_high_mass=" << expected_high_m
                  << " interface_high_ratio_wrong_separation="
                  << std::abs(expected_high_m - wrong_ratio_m)
                  << " interior_high_mass=" << secondary_high_m
                  << " interior_physical_M_weno_separation="
                  << std::abs(secondary_high_m - physical_mass_weno)
                  << " interior_physical_N_weno_separation="
                  << std::abs(secondary_high_n - physical_number_weno)
                  << " interior_wrong_ratio_separation="
                  << std::abs(secondary_high_m - secondary_wrong_ratio_m)
                  << " interior_lambda=" << secondary_limiter << '\n';

        const Real ratio_base_m = (Real(1.0) - base_theta) * oldp.mass0 +
                                  base_theta * newp.mass0;
        const Real ratio_base_n = (Real(1.0) - base_theta) * oldp.number0 +
                                  base_theta * newp.number0;
        const Real ratio_base_p = (Real(1.0) - base_theta) * oldp.property0 +
                                  base_theta * newp.property0;

        fine_integral[0] += static_cast<Real>(fine_dt) * rate_values[0];
        fine_integral[1] += static_cast<Real>(fine_dt) * rate_values[1];
        fine_integral[2] += static_cast<Real>(fine_dt) * rate_values[2];
        ratio_wrong_integral[0] += static_cast<Real>(fine_dt) * velocity * ratio_base_m;
        ratio_wrong_integral[1] += static_cast<Real>(fine_dt) * velocity * ratio_base_n;
        ratio_wrong_integral[2] += static_cast<Real>(fine_dt) * velocity * ratio_base_p;
        ASSERT_TRUE(transport.accept_physical_step(
            1, begin, end, diagnostic)) << diagnostic;
    }
    ASSERT_TRUE(ya.ready_to_reflux());

    MultiFab coarse_measure(coarse_ba, coarse_dm, 1, 0);
    MultiFab fine_measure(fine_ba, fine_dm, 1, 0);
    coarse_measure.setVal(Real(1.0));
    fine_measure.setVal(Real(1.0));
    erf_sbm::SBMAMRStateView coarse_view{
        &state_manager.new_state(0), t1, &rho_new, Rho_comp, t1,
        &coarse_measure, 0, t1};
    erf_sbm::SBMAMRStateView fine_view{
        &state_manager.new_state(1), t1, &fine_rho_target, Rho_comp, t1,
        &fine_measure, 0, t1};
    MultiFab candidate(coarse_ba, coarse_dm, layout.ncomp(), 0);
    MultiFab::Copy(candidate, state_manager.new_state(0), 0, 0,
                   layout.ncomp(), 0);
    ASSERT_TRUE(erf_sbm::RestrictMappedSpectrumPrivate(
        layout, fine_view, coarse_view, ratio, 0, candidate, diagnostic,
        16, false, nullptr, nullptr, fine_dt)) << diagnostic;

    MultiFab correction(coarse_ba, coarse_dm, layout.ncomp(), 0);
    correction.setVal(Real(0.0));
    ASSERT_TRUE(ya.reflux(correction, diagnostic)) << diagnostic;
    MultiFab marker(fine_ba, fine_dm, 1, 0);
    MultiFab covered(coarse_ba, coarse_dm, 1, 1);
    marker.setVal(Real(1.0));
    covered.setVal(Real(0.0));
    amrex::average_down(marker, covered, 0, 1, ratio);
    covered.FillBoundary(coarse_geom.periodicity());
    MultiFab invalid(coarse_ba, coarse_dm, 1, 0);
    invalid.setVal(Real(0.0));
    ASSERT_TRUE(erf_sbm::ApplyMappedRefluxCorrectionChunk(
        candidate, correction, coarse_measure, covered, components, invalid,
        diagnostic)) << diagnostic;
    ASSERT_TRUE(erf_sbm::MappedRefluxCorrectionMaskIsValid(
        invalid, diagnostic)) << diagnostic;

    const IntVect coarse_cell(0, 1, 1);
    const Box cell(coarse_cell, coarse_cell);
    const Real dx = coarse_geom.CellSize(0);
    const Real tol = Real(512.0) *
                     std::numeric_limits<Real>::epsilon();
    const std::array<int, 3> tested{mass0, number0, property0};
    for (int slot = 0; slot < 3; ++slot) {
        const int component = tested[static_cast<std::size_t>(slot)];
        const Real expected = -fine_integral[static_cast<std::size_t>(slot)] /
                              (Real(4.0) * dx);
        const Real base = state_manager.new_state(0).sum(cell, component, false);
        const Real actual = candidate.sum(cell, component, false) - base;
        EXPECT_NEAR(actual, expected, tol) << "component=" << component;
        EXPECT_GT(std::abs(actual + expected), Real(100.0) * tol)
            << "signed oracle must reject reversed FV sign";
        EXPECT_GT(std::abs(actual), Real(100.0) * tol)
            << "signed oracle must reject omitted fine amount";
        EXPECT_GT(std::abs(actual - expected * (t1 - t0)), Real(100.0) * tol)
            << "signed oracle must reject repeated interval weighting";
        EXPECT_GT(std::abs(actual +
                           ratio_wrong_integral[static_cast<std::size_t>(slot)] /
                               (Real(4.0) * dx)),
                  Real(100.0) * tol)
            << "FV result must distinguish endpoint-ratio interpolation";
    }
    const Real mass_delta = candidate.sum(cell, mass0, false) -
        state_manager.new_state(0).sum(cell, mass0, false);
    const Real wrong_component = -fine_integral[1] / (Real(4.0) * dx);
    EXPECT_GT(std::abs(mass_delta - wrong_component), Real(100.0) * tol)
        << "signed oracle must reject mass/number permutation";
    ASSERT_TRUE(erf_sbm::authoritative_state_admissible(
        candidate, layout, 0, &diagnostic)) << diagnostic;
    ASSERT_TRUE(state_manager.can_accept_synchronized_state(
        0, candidate, t1, diagnostic)) << diagnostic;
    const Real old_coarse_mass = state_manager.new_state(0).sum(
        cell, mass0, false);
    ASSERT_TRUE(ya.accept_reflux(diagnostic)) << diagnostic;
    state_manager.publish_synchronized_state(0, candidate);
    EXPECT_NEAR(state_manager.new_state(0).sum(cell, mass0, false) -
                    old_coarse_mass,
                -fine_integral[0] / (Real(4.0) * dx), tol);

    std::cout << std::setprecision(17)
              << "SBM_A4_COUPLED_ORACLE mass_integral=" << fine_integral[0]
              << " number_integral=" << fine_integral[1]
              << " property_integral=" << fine_integral[2]
              << " mass_delta_h=" << -fine_integral[0] / (Real(4.0) * dx)
              << " wrong_ratio_mass_integral=" << ratio_wrong_integral[0]
              << " interface_lambda=" << interface_lambda
              << " stages=RK3x3 direct_YA=true\n";
#endif
}

} // namespace
