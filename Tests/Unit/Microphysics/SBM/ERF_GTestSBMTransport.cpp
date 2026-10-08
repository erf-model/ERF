#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_AdvectionSrcForScalars.H"
#include "ERF_SBMAdvectionBoundary.H"
#include "ERF_IndexDefines.H"
#include "ERF_SBMRemapping.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMStageOwnership.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
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
    state_manager.define(0, ba, dm);
    erf_sbm::SBMTransport transport(state_manager.layout(), 1,
                                    options.max_groups_per_chunk);
    transport.define(0, ba, dm);

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

    auto& spectrum = state_manager.state(0);
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
    transport.advance_stage_from_host(
        0,
        options.anelastic_heun ? erf_auxiliary::HostIntegrator::AnelasticHeun
                               : erf_auxiliary::HostIntegrator::CompressibleRK3,
        0, 0.0, 0.0, stage0_target_time, stage0_interval, state_manager,
        conserved_anchor, conserved_input, conserved_target, avg_xmom,
        avg_ymom, avg_zmom, geom, 1, 2);

    const auto& final_state = state_manager.state(0);
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
                transport.advance_stage_from_host(
                    0, erf_auxiliary::HostIntegrator::CompressibleRK3, 0,
                    step_old_time, step_old_time, stage0_time, dt / 3.0,
                    state_manager, conserved_anchor, conserved_input,
                    conserved_target, avg_xmom, avg_ymom, avg_zmom, geom, 1, 2);
            }
            if (step == 0 && options.amplify_corrector_input) {
                for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
                    const Box bx = mfi.validbox();
                    const auto state = spectrum.array(mfi);
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
            EXPECT_TRUE(spectrum.is_finite(0, spectrum.nComp(), 0));
            diagnostic.clear();
            EXPECT_TRUE(erf_sbm::authoritative_state_admissible(
                spectrum, layout, 0, &diagnostic))
                << "after completed step " << step << ": " << diagnostic;
            for (int component = 0; component < layout.ncomp(); ++component) {
                EXPECT_NEAR(spectrum.sum(component), initial.sum(component),
                            Real(512.0) *
                                std::numeric_limits<Real>::epsilon() *
                                std::max(std::abs(initial.sum(component)),
                                         std::numeric_limits<Real>::min()))
                    << "component=" << component << " step=" << step;
            }
        }
    }

    MultiFab expected_core(ba, dm, 3, 0);
    MultiFab::Copy(expected_core, conserved_target, 0, 0, 3, 0);
    for (amrex::MFIter mfi(expected_core); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = final_state.const_array(mfi);
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
            const auto state = final_state.const_array(mfi);
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
    for (amrex::MFIter mfi(final_state); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = final_state.const_array(mfi);
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
    for (int component = 0; component < layout.ncomp(); ++component) {
        summary.inventory.push_back(final_state.sum(component));
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
        state_manager.define(0, ba, dm);
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

        auto& spectrum = state_manager.state(0);
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

        const auto& final_state = state_manager.state(0);
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
    state_manager.define(0, ba, dm);
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

    constexpr int qc = 5;
    constexpr int qr = 6;
    constexpr int unrelated = 7;
    MultiFab anchor(ba, dm, 8, 0);
    MultiFab predictor(ba, dm, 8, 0);
    MultiFab target(ba, dm, 8, 0);
    anchor.setVal(Real(0.0));
    anchor.setVal(Real(1.0), Rho_comp, 1, 0);
    anchor.setVal(Real(300.0), RhoTheta_comp, 1, 0);
    anchor.setVal(Real(0.01), RhoQ1_comp, 1, 0);
    MultiFab::Copy(predictor, anchor, 0, 0, 8, 0);
    MultiFab::Copy(target, anchor, 0, 0, 8, 0);
    predictor.setVal(Real(0.03), RhoQ1_comp, 1, 0);
    predictor.setVal(Real(11.0), qc, 1, 0);
    predictor.setVal(Real(13.0), qr, 1, 0);
    predictor.setVal(Real(17.0), unrelated, 1, 0);

    auto& spectrum = state_manager.state(0);
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
    MultiFab non_sbm_target(ba, dm, 8, 0);
    MultiFab::Copy(non_sbm_target, anchor, 0, 0, 8, 0);
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

} // namespace
