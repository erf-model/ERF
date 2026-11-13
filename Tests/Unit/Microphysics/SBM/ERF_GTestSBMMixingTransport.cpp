#include <gtest/gtest.h>

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include "ERF_AuxiliaryDiffusion.H"
#include "ERF_AuxiliaryMappedTransport.H"
#include "ERF_IndexDefines.H"
#include "ERF_SBMParticleMixing.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <algorithm>
#include <array>
#include <cmath>
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
using erf_auxiliary::HostIntegrator;
using erf_sbm::ParticleMixingPolicy;

constexpr Real two_pi = Real(6.283185307179586476925286766559);
constexpr int nx = 8;
constexpr int ny = 4;
constexpr int nz = 4;

Real eps () { return std::numeric_limits<Real>::epsilon(); }

erf_sbm::SBMLayout
make_layout (const erf_sbm::MomentMode mode, const bool with_property)
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

enum class Profile { Smooth, Step, Uniform };

struct Options
{
    erf_sbm::MomentMode mode{erf_sbm::MomentMode::OneMoment};
    bool attached_property{false};
    ParticleMixingPolicy policy{ParticleMixingPolicy::HostPassiveScalar};
    Real kh{0};
    Real kv{0};
    Real eddy_scale{1};
    // Walls bound z; a stretched mesh also has walls.
    bool walls{false};
    bool stretched{false};
    bool variable_maps{false};
    bool varying_density{false};
    Real carrier_x{0};
    Profile profile{Profile::Smooth};
};

struct Run
{
    Options options;
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    std::unique_ptr<erf_sbm::SBMStateManager> state;
    std::unique_ptr<erf_sbm::SBMTransport> transport;
    MultiFab rho;
    MultiFab target;
    MultiFab eddy;
    amrex::Vector<std::unique_ptr<MultiFab>> mapfac;
    MultiFab ax, ay, detj;
    amrex::Gpu::DeviceVector<Real> dz;
    MultiFab xmom, ymom, zmom;
    erf_sbm::SBMStageMixingInputs inputs;

    [[nodiscard]] const erf_sbm::SBMStageMixingInputs* mixing () const
    {
        return transport->mixes_particles() ? &inputs : nullptr;
    }
    [[nodiscard]] const MultiFab& spectrum () const { return state->new_state(0); }
};

// Physical mesh spacing of the stretched column; the computational spacing is 1/nz.
constexpr std::array<Real, nz> stretched_dz{Real(0.15), Real(0.2), Real(0.3), Real(0.35)};

void fill_map_factor (MultiFab& mf, const Real base, const Real amplitude, const int phase)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
            a(i, j, 0, 0) = base + amplitude * std::sin(two_pi * Real(i + phase) / Real(nx)) +
                            Real(0.5) * amplitude * std::cos(two_pi * Real(j) / Real(ny));
        });
    }
}

void fill_column (MultiFab& mf, const amrex::Gpu::DeviceVector<Real>& dz)
{
    const Real* dzp = dz.data();
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, 0) = dzp[k] * Real(nz);
        });
    }
}

void fill_density (MultiFab& rho, const bool varying)
{
    rho.setVal(Real(0.0));
    for (amrex::MFIter mfi(rho); mfi.isValid(); ++mfi) {
        const auto a = rho.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, Rho_comp) = varying
                ? Real(1.0) + Real(0.15) * std::sin(two_pi * Real(i) / Real(nx)) +
                      Real(0.05) * std::cos(two_pi * Real(j) / Real(ny)) + Real(0.02) * Real(k)
                : Real(1.0);
        });
    }
}

// Spatially varying, anisotropic eddy diffusivities; every other channel holds
// a distinct value that a wrong channel would pick up.
void fill_eddy (MultiFab& eddy, const Real scale)
{
    for (amrex::MFIter mfi(eddy); mfi.isValid(); ++mfi) {
        const auto a = eddy.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            for (int n = 0; n < EddyDiff::NumDiffs; ++n) {
                a(i, j, k, n) = Real(0.5) + Real(0.1) * Real(n);
            }
            a(i, j, k, EddyDiff::Scalar_h) =
                scale * (Real(0.012) + Real(0.006) * std::sin(two_pi * Real(i) / Real(nx)) +
                         Real(0.003) * std::cos(two_pi * Real(j) / Real(ny)));
            a(i, j, k, EddyDiff::Scalar_v) =
                scale * (Real(0.004) + Real(0.002) * std::cos(two_pi * Real(i) / Real(nx)) +
                         Real(0.001) * Real(k));
        });
    }
}

// U = rho z for the selected profile; 2M bins keep one mean mass per half
// domain and attached properties a quarter (or a tenth) of the number.
void fill_spectrum (MultiFab& spectrum, const MultiFab& rho, const erf_sbm::SBMLayout& layout,
                    const Profile profile, const bool walls)
{
    const auto& population = layout.populations().front();
    const bool two_moment = population.moment_mode == erf_sbm::MomentMode::TwoMoment;
    const int mass0 = population.mass_offset;
    const int number0 = two_moment ? population.number_offset : -1;
    const int property0 = layout.ncomp() > (two_moment ? 4 : 2) ? layout.property_offset(0) : -1;
    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const auto u = spectrum.array(mfi);
        const auto r = rho.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real x = (Real(i) + Real(0.5)) / Real(nx);
            const Real y = (Real(j) + Real(0.5)) / Real(ny);
            const Real z = (Real(k) + Real(0.5)) / Real(nz);
            Real factor = Real(1.0);
            bool upper = false;
            if (profile == Profile::Smooth) {
                factor = Real(1.0) + Real(0.3) * std::sin(two_pi * x) * std::cos(two_pi * y) +
                         Real(0.2) * (walls ? std::cos(Real(0.5) * two_pi * z)
                                            : std::cos(two_pi * z));
            } else if (profile == Profile::Step) {
                upper = walls ? z > Real(0.5) : x > Real(0.5);
                factor = upper ? Real(0.1) : Real(1.0);
            }
            const Real density = r(i, j, k, Rho_comp);
            for (int bin = 0; bin < 2; ++bin) {
                const Real number = factor * Real(0.01) * (Real(1.0) + Real(0.5) * Real(bin));
                const Real mean = bin == 0 ? (upper ? Real(0.4) : Real(0.2)) : Real(0.75);
                u(i, j, k, mass0 + bin) = density * number * mean;
                if (number0 >= 0) { u(i, j, k, number0 + bin) = density * number; }
                if (property0 >= 0) {
                    u(i, j, k, property0 + bin) =
                        density * number * (upper ? Real(0.1) : Real(0.4));
                }
            }
        });
    }
}

std::unique_ptr<Run>
make_run (const Options& options)
{
    auto run = std::make_unique<Run>();
    run->options = options;
    const bool walls = options.walls || options.stretched;
    const Box domain(IntVect(0), IntVect(nx - 1, ny - 1, nz - 1));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, walls ? 0 : 1};
    run->geom = Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodic);
    run->ba = BoxArray(domain);
    run->ba.maxSize(4);
    run->dm = DistributionMapping(run->ba);
    const auto& ba = run->ba;
    const auto& dm = run->dm;

    run->state = std::make_unique<erf_sbm::SBMStateManager>(
        make_layout(options.mode, options.attached_property), 1);
    run->state->define(0, ba, dm, 0.0);
    run->transport = std::make_unique<erf_sbm::SBMTransport>(run->state->layout(), 1, 16);
    run->transport->define(0, ba, dm);
    if (options.policy != ParticleMixingPolicy::None) {
        erf_sbm::ParticleMixingConfig config;
        config.policy = options.policy;
        config.reference_kh = options.kh;
        config.reference_kv = options.kv;
        run->transport->set_particle_mixing(
            config, options.stretched ? erf_auxiliary::DiffusionMesh::StretchedDz
                                      : erf_auxiliary::DiffusionMesh::ConstantDz);
    }
    if (walls) {
        erf_sbm::SBMPhysicalBoundaryConfig boundaries;
        boundaries.kind[2] = erf_sbm::AdvectionFaceKind::ImpermeableWall;
        boundaries.kind[5] = erf_sbm::AdvectionFaceKind::ImpermeableWall;
        run->transport->set_physical_boundaries(boundaries);
    }

    amrex::BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    const BoxArray ba2d(std::move(boxes));
    run->mapfac.resize(MapFacType::num);
    run->mapfac[MapFacType::m_x] = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    run->mapfac[MapFacType::u_x] = std::make_unique<MultiFab>(
        amrex::convert(ba2d, IntVect(1, 0, 0)), dm, 1, 0);
    run->mapfac[MapFacType::v_x] = std::make_unique<MultiFab>(
        amrex::convert(ba2d, IntVect(0, 1, 0)), dm, 1, 0);
    const Real amplitude = options.variable_maps ? Real(0.05) : Real(0.0);
    fill_map_factor(*run->mapfac[MapFacType::m_x], Real(1.0), amplitude, 0);
    fill_map_factor(*run->mapfac[MapFacType::u_x], Real(1.0), amplitude, 1);
    fill_map_factor(*run->mapfac[MapFacType::v_x], Real(1.0), amplitude, 2);

    std::vector<Real> dz_host(static_cast<std::size_t>(nz), Real(1.0) / Real(nz));
    if (options.stretched) {
        std::copy(stretched_dz.begin(), stretched_dz.end(), dz_host.begin());
    }
    run->dz.resize(dz_host.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, dz_host.begin(), dz_host.end(), run->dz.begin());
    run->detj.define(ba, dm, 1, 0);
    run->ax.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    run->ay.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    for (MultiFab* mf : {&run->detj, &run->ax, &run->ay}) { fill_column(*mf, run->dz); }
    std::string diagnostic;
    EXPECT_TRUE(run->transport->rebuild_static_measure(0, run->detj, *run->mapfac[MapFacType::m_x],
                                                       *run->mapfac[MapFacType::m_x], diagnostic))
        << diagnostic;

    run->rho.define(ba, dm, 3, 1);
    fill_density(run->rho, options.varying_density);
    run->rho.FillBoundary(run->geom.periodicity());
    run->target.define(ba, dm, 3, 1);
    MultiFab::Copy(run->target, run->rho, 0, 0, 3, 1);
    run->eddy.define(ba, dm, EddyDiff::NumDiffs, 1);
    fill_eddy(run->eddy, options.eddy_scale);
    run->eddy.FillBoundary(run->geom.periodicity());

    run->xmom.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    run->ymom.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    run->zmom.define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
    run->xmom.setVal(options.carrier_x);
    run->ymom.setVal(Real(0.0));
    run->zmom.setVal(Real(0.0));

    fill_spectrum(run->state->new_state_for_initialization(0), run->rho, run->state->layout(),
                  options.profile, walls);
    run->inputs = erf_sbm::make_stage_mixing_inputs(*run->transport, &run->eddy, run->mapfac,
                                                    &run->ax, &run->ay, run->dz.data());
    return run;
}

void stage (Run& run, const HostIntegrator method, const int index, const double step_old,
            const double input, const double target, const double interval)
{
    run.transport->advance_stage_from_host(
        0, method, index, step_old, input, target, interval, *run.state, run.rho, run.rho,
        run.target, run.xmom, run.ymom, run.zmom, run.geom, 1, 2, nullptr, run.mixing());
}

void begin_step (Run& run, const double time)
{
    std::string diagnostic;
    EXPECT_TRUE(run.state->begin_step(0, time, diagnostic)) << diagnostic;
}

// One step from time t; returns the step's end time, which the next step starts from.
double step_once (Run& run, const HostIntegrator method, const double t, const double dt)
{
    const double t_new = t + dt;
    begin_step(run, t);
    if (method == HostIntegrator::CompressibleRK3) {
        stage(run, method, 0, t, t, t + dt / 3.0, dt / 3.0);
        stage(run, method, 1, t, t + dt / 3.0, t + dt / 2.0, dt / 2.0);
        stage(run, method, 2, t, t + dt / 2.0, t_new, dt);
    } else if (method == HostIntegrator::AnelasticHeun) {
        stage(run, method, 0, t, t, t_new, dt);
        stage(run, method, 1, t, t_new, t_new, dt);
    } else {
        stage(run, method, 0, t, t, t + dt / 2.0, dt / 2.0);
        stage(run, method, 1, t, t + dt / 2.0, t_new, dt);
    }
    return t_new;
}

void advance (Run& run, const HostIntegrator method, const double dt, const int steps)
{
    double t = 0.0;
    for (int step = 0; step < steps; ++step) { t = step_once(run, method, t, dt); }
}

std::unique_ptr<MultiFab>
copy_of (const MultiFab& source)
{
    auto copy = std::make_unique<MultiFab>(source.boxArray(), source.DistributionMap(),
                                           source.nComp(), 0);
    MultiFab::Copy(*copy, source, 0, 0, source.nComp(), 0);
    return copy;
}

Real difference (const MultiFab& a, const MultiFab& b)
{
    // Separate runs build their own distribution maps, which need not agree.
    MultiFab d(a.boxArray(), a.DistributionMap(), a.nComp(), 0);
    d.ParallelCopy(b, 0, 0, a.nComp());
    d.minus(a, 0, a.nComp(), 0);
    return d.norm0(0, a.nComp(), IntVect(0));
}

Real mapped_inventory (const MultiFab& state, const MultiFab& measure, const int component)
{
    MultiFab product(state.boxArray(), state.DistributionMap(), 1, 0);
    MultiFab::Copy(product, state, component, 0, 1, 0);
    MultiFab::Multiply(product, measure, 0, 0, 1, 0);
    return product.sum(0);
}

// dU/dt of the M5a two-point operator on z = U / rho: an oracle independent
// of the native diffusion builders that form the transport's high-order rate.
void two_point_tendency (const Run& run, const MultiFab& state, MultiFab& tendency)
{
    const int ncomp = state.nComp();
    MultiFab z(run.ba, run.dm, ncomp, 1);
    z.setVal(Real(0.0));
    for (amrex::MFIter mfi(z); mfi.isValid(); ++mfi) {
        const auto out = z.array(mfi);
        const auto u = state.const_array(mfi);
        const auto r = run.rho.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                out(i, j, k, n) = u(i, j, k, n) / r(i, j, k, Rho_comp);
            });
    }
    z.FillBoundary(run.geom.periodicity());
    erf_auxiliary::MappedFaceTransmissibility t;
    t.define(run.ba, run.dm, 1, 0);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::build_particle_mixing_transmissibility(
        run.transport->particle_mixing(), run.geom, run.inputs.metrics, run.rho, Rho_comp,
        run.inputs.eddy_diffusivity, t, diagnostic)) << diagnostic;
    erf_auxiliary::MappedFaceFluxRate rate;
    rate.define(run.ba, run.dm, ncomp, 0);
    std::vector<int> components(static_cast<std::size_t>(ncomp));
    for (int n = 0; n < ncomp; ++n) { components[static_cast<std::size_t>(n)] = n; }
    EXPECT_TRUE(erf_auxiliary::BuildLowOrderDiffusionRate(rate, t, run.geom, z, components,
                                                          diagnostic)) << diagnostic;
    for (int n = 0; n < ncomp; ++n) {
        erf_auxiliary::ApplyMappedFluxTendency(rate, n, run.transport->static_measure(0), 0,
                                               tendency, n, run.geom.InvCellSizeArray());
    }
}

// a + tau * L(b)
std::unique_ptr<MultiFab>
oracle_update (const Run& run, const MultiFab& a, const MultiFab& b, const Real tau)
{
    MultiFab tendency(run.ba, run.dm, a.nComp(), 0);
    two_point_tendency(run, b, tendency);
    auto out = copy_of(a);
    MultiFab::Saxpy(*out, tau, tendency, 0, 0, a.nComp(), 0);
    return out;
}

// One step of the host integrator with the oracle operator at fixed density.
std::unique_ptr<MultiFab>
oracle_step (const Run& run, const MultiFab& u0, const HostIntegrator method, const Real dt)
{
    if (method == HostIntegrator::CompressibleRK3) {
        auto u1 = oracle_update(run, u0, u0, dt / Real(3.0));
        auto u2 = oracle_update(run, u0, *u1, dt / Real(2.0));
        return oracle_update(run, u0, *u2, dt);
    }
    if (method == HostIntegrator::AnelasticHeun) {
        auto u1 = oracle_update(run, u0, u0, dt);
        auto u2 = oracle_update(run, *u1, *u1, dt);
        MultiFab::LinComb(*u2, Real(0.5), u0, 0, Real(0.5), *u2, 0, 0, u0.nComp(), 0);
        return u2;
    }
    auto half = oracle_update(run, u0, u0, dt / Real(2.0));
    return oracle_update(run, u0, *half, dt);
}

// With the limiter inactive every accepted rate is ERF's diffusion of the
// stage input, so one step reproduces the integrator applied to that operator.
void expect_recovers_host_diffusion (const Options& options, const HostIntegrator method)
{
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    const Real dt = Real(0.05);
    advance(*run, method, dt, 1);
    const auto expected = oracle_step(*run, *initial, method, dt);
    const Real change = difference(*expected, *initial);
    const Real scale = initial->norm0(0, initial->nComp(), IntVect(0));
    EXPECT_GT(change, Real(1.0e-4) * scale);
    EXPECT_LE(difference(run->spectrum(), *expected), Real(1024.0) * eps() * scale);
}

void run_recovery ()
{
    Options n_grid;
    n_grid.variable_maps = true;
    n_grid.varying_density = true;
    expect_recovers_host_diffusion(n_grid, HostIntegrator::CompressibleRK3);
    n_grid.mode = erf_sbm::MomentMode::TwoMoment;
    expect_recovers_host_diffusion(n_grid, HostIntegrator::AnelasticMidPoint);

    // Distinct horizontal and vertical reference diffusivities take one
    // native pass each on a stretched column.
    Options s_grid;
    s_grid.policy = ParticleMixingPolicy::ConstantReference;
    s_grid.kh = Real(0.02);
    s_grid.kv = Real(0.004);
    s_grid.stretched = true;
    s_grid.varying_density = true;
    expect_recovers_host_diffusion(s_grid, HostIntegrator::AnelasticHeun);
}

Real final_error (const HostIntegrator method, const int steps, const MultiFab& reference)
{
    Options options;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.02);
    options.kv = Real(0.02);
    auto run = make_run(options);
    advance(*run, method, Real(0.4) / Real(steps), steps);
    return difference(run->spectrum(), reference);
}

// Diffusion evaluated on the actual stage input keeps each integrator's order;
// an old-state surrogate would leave only first order.
void run_temporal_order (const HostIntegrator method, const Real expected_order)
{
    Options options;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.02);
    options.kv = Real(0.02);
    auto reference = make_run(options);
    advance(*reference, method, Real(0.4) / Real(64), 64);
    const Real e4 = final_error(method, 4, reference->spectrum());
    const Real e8 = final_error(method, 8, reference->spectrum());
    const Real e16 = final_error(method, 16, reference->spectrum());
    const Real coarse = std::log2(e4 / e8);
    const Real fine = std::log2(e8 / e16);
    amrex::Print() << "SBM_M6A_TEMPORAL integrator=" << erf_auxiliary::HostIntegratorName(method)
                   << " e4=" << e4 << " e8=" << e8 << " e16=" << e16
                   << " order=" << coarse << "," << fine << '\n';
    EXPECT_GE(coarse, expected_order);
    EXPECT_GE(fine, expected_order);
}

// A stage input far from the anchor gives a high-order correction that the
// unlimited update cannot absorb: the limiter must bound it and leave every
// cell it drains admissible.
void run_limited_correction_with (Options options)
{
    options.profile = Profile::Uniform;
    auto run = make_run(options);
    const Real dt = Real(0.08);
    const auto anchor = copy_of(run->spectrum());
    begin_step(*run, 0.0);
    stage(*run, HostIntegrator::CompressibleRK3, 0, 0.0, 0.0, dt / 3.0, dt / 3.0);
    // Move nearly all material of the stage input into one cell.
    auto& input = run->state->new_target_storage(0);
    for (amrex::MFIter mfi(input); mfi.isValid(); ++mfi) {
        const auto u = input.array(mfi);
        amrex::ParallelFor(mfi.validbox(), input.nComp(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                u(i, j, k, n) *= (i == 5 && j == 2 && k == 1) ? Real(50.0) : Real(1.0e-6);
            });
    }
    const auto stage_input = copy_of(input);
    stage(*run, HostIntegrator::CompressibleRK3, 1, 0.0, dt / 3.0, dt / 2.0, dt / 2.0);

    if (options.policy != ParticleMixingPolicy::None && options.carrier_x == Real(0.0)) {
        const auto unlimited = oracle_update(*run, *anchor, *stage_input, dt / Real(2.0));
        EXPECT_LT(unlimited->min(0), -Real(0.1) * anchor->max(0));
        EXPECT_GT(difference(run->spectrum(), *unlimited), Real(0.1) * anchor->max(0));
    }
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::authoritative_state_admissible(run->spectrum(), run->state->layout(), 0,
                                                        &diagnostic)) << diagnostic;
    for (int n = 0; n < anchor->nComp(); ++n) {
        const Real before = mapped_inventory(*anchor, run->transport->static_measure(0), n);
        const Real after = mapped_inventory(run->spectrum(), run->transport->static_measure(0), n);
        EXPECT_NEAR(after, before, Real(256.0) * eps() * before) << "component " << n;
    }
}

void run_limited_correction (const Real carrier)
{
    Options options;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.02);
    options.kv = Real(0.02);
    options.carrier_x = carrier;
    run_limited_correction_with(options);
}

// U = density * 0.01 (n + 1) * v[i]: columns uniform in y and z.
void fill_columns (MultiFab& u, const amrex::GpuArray<Real, nx>& v, const Real density)
{
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const auto a = u.array(mfi);
        amrex::ParallelFor(mfi.validbox(), u.nComp(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                a(i, j, k, n) = density * Real(0.01) * Real(n + 1) * v[i];
            });
    }
}

// Where the limiter binds, the accepted rates fall back to the low-order
// exchange of the trial base. On x columns the anchor holds an empty cell h
// between neighbours p and q in a field of ones, and the stage input, at twice
// the density, a spike B in h with the same neighbours. With D the low-order
// exchange weight (anchor density) and Dh = 2 D the high-order one (input
// density), the limiter drains h exactly with r = D (p + q) / (Dh (2B - p - q)
// + D (p + q)), and the neighbours receive closed-form amounts that change if
// the exchange were left out or taken at the input density.
void run_limited_cells_keep_low_order_exchange ()
{
    Options options;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.01);
    options.kv = Real(0.01);
    auto run = make_run(options);
    constexpr int h = 3;
    const Real p = Real(1.0);
    const Real q = Real(0.25);
    const Real B = Real(50.0);
    const Real rho_input = Real(2.0);
    const Real dt = Real(0.08);
    amrex::GpuArray<Real, nx> v;
    for (int i = 0; i < nx; ++i) { v[i] = Real(1.0); }
    v[h - 1] = q;
    v[h] = Real(0.0);
    v[h + 1] = p;
    fill_columns(run->state->new_state_for_initialization(0), v, Real(1.0));
    const Real scale = run->spectrum().max(0);

    // D = tau t_x from the oracle: a unit column at h gives t_x to each neighbour.
    amrex::GpuArray<Real, nx> unit;
    for (int i = 0; i < nx; ++i) { unit[i] = i == h ? Real(1.0) : Real(0.0); }
    MultiFab delta(run->ba, run->dm, run->spectrum().nComp(), 0);
    fill_columns(delta, unit, Real(1.0));
    MultiFab tendency(run->ba, run->dm, delta.nComp(), 0);
    two_point_tendency(*run, delta, tendency);
    const Real D = (dt / Real(2.0)) * tendency.max(0) / Real(0.01);
    const Real Dh = rho_input * D;
    const Real r = D * (p + q) / (Dh * (Real(2.0) * B - p - q) + D * (p + q));
    const Real right = p - D * p + r * (Dh * (B - p) + D * p);
    const Real left = q - D * q + r * (Dh * (B - q) + D * q) + Dh * (Real(1.0) - q);
    const Real far = Real(1.0) - Dh * (Real(1.0) - q);
    // The answers of a limiter falling back to advection alone (r = 0) or to
    // an exchange at the input density must be far outside the tolerance, which
    // allows for the limiter's roundoff margin (128 eps of the drained cell's
    // gross budget, about 1e3 eps of the result here).
    const Real tolerance = std::max(Real(1.0e-10), Real(2048.0) * eps()) * scale;
    const Real r_input = (p + q) / (Real(2.0) * B);
    const Real right_input = p - Dh * p + r_input * (Dh * (B - p) + Dh * p);
    ASSERT_GT(Real(0.01) * std::abs(right - p), Real(10.0) * tolerance);
    ASSERT_GT(Real(0.01) * std::abs(right - right_input), Real(10.0) * tolerance);

    begin_step(*run, 0.0);
    stage(*run, HostIntegrator::CompressibleRK3, 0, 0.0, 0.0, dt / 3.0, dt / 3.0);
    MultiFab input_density(run->ba, run->dm, run->rho.nComp(), run->rho.nGrow());
    MultiFab::Copy(input_density, run->rho, 0, 0, run->rho.nComp(), run->rho.nGrow());
    input_density.setVal(rho_input, Rho_comp, 1, input_density.nGrow());
    amrex::GpuArray<Real, nx> spike = v;
    spike[h] = B;
    fill_columns(run->state->new_target_storage(0), spike, rho_input);
    run->transport->advance_stage_from_host(
        0, HostIntegrator::CompressibleRK3, 1, 0.0, dt / 3.0, dt / 2.0, dt / 2.0, *run->state,
        run->rho, input_density, run->target, run->xmom, run->ymom, run->zmom, run->geom, 1, 2,
        nullptr, run->mixing());

    amrex::GpuArray<Real, nx> accepted = v;
    accepted[h - 2] = far;
    accepted[h - 1] = left;
    accepted[h] = Real(0.0);
    accepted[h + 1] = right;
    MultiFab expected(run->ba, run->dm, delta.nComp(), 0);
    fill_columns(expected, accepted, Real(1.0));
    EXPECT_LE(difference(run->spectrum(), expected), tolerance);
}

// Forked death tests hang under MPI, so amrex::Abort throws for the check.
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
        << "expected amrex::Abort containing \"" << expected << "\", got \"" << message << '"';
}

// One Heun predictor stage with advective and diffusive demand fractions of
// the trial interval: uniform mass flux c, reference diffusivity K and unit
// density give c dx_inv dt and 2 K (dx_inv^2 + dy_inv^2 + dz_inv^2) dt.
void run_combined_stage (const Real advective, const Real diffusive)
{
    const Real dt = Real(0.01);
    Options options;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = diffusive / (dt * Real(2.0) * Real(nx * nx + ny * ny + nz * nz));
    options.kv = options.kh;
    options.carrier_x = advective / (dt * Real(nx));
    options.profile = Profile::Uniform;
    auto run = make_run(options);
    begin_step(*run, 0.0);
    stage(*run, HostIntegrator::AnelasticHeun, 0, 0.0, 0.0, dt, dt);
}

void run_combined_bound ()
{
    run_combined_stage(Real(0.4), Real(0.4));
    expect_amrex_abort([] { run_combined_stage(Real(0.6), Real(0.6)); },
                       "SBM combined advection-diffusion outgoing-demand condition failed");
}

void run_zero_coefficient ()
{
    Options mixing;
    mixing.eddy_scale = Real(0.0);
    mixing.carrier_x = Real(2.0);
    mixing.profile = Profile::Step;
    mixing.mode = erf_sbm::MomentMode::TwoMoment;
    mixing.attached_property = true;
    Options none = mixing;
    none.policy = ParticleMixingPolicy::None;
    auto with_zero_mixing = make_run(mixing);
    auto without_mixing = make_run(none);
    advance(*with_zero_mixing, HostIntegrator::CompressibleRK3, Real(0.01), 2);
    advance(*without_mixing, HostIntegrator::CompressibleRK3, Real(0.01), 2);
    EXPECT_EQ(difference(with_zero_mixing->spectrum(), without_mixing->spectrum()), Real(0.0));
}

void run_constant_ratio ()
{
    Options options;
    options.profile = Profile::Uniform;
    options.varying_density = true;
    options.variable_maps = true;
    options.carrier_x = Real(1.5);
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    advance(*run, HostIntegrator::CompressibleRK3, Real(0.02), 2);
    EXPECT_LE(difference(run->spectrum(), *initial),
              Real(256.0) * eps() * initial->norm0(0, initial->nComp(), IntVect(0)));
}

void run_walls ()
{
    Options options;
    options.stretched = true;
    options.profile = Profile::Step;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    advance(*run, HostIntegrator::CompressibleRK3, Real(0.05), 3);
    EXPECT_GT(difference(run->spectrum(), *initial), Real(1.0e-3) * initial->max(0));
    for (int n = 0; n < initial->nComp(); ++n) {
        const Real before = mapped_inventory(*initial, run->transport->static_measure(0), n);
        const Real after = mapped_inventory(run->spectrum(), run->transport->static_measure(0), n);
        EXPECT_NEAR(after, before, Real(256.0) * eps() * before) << "component " << n;
    }
    const auto& ledger = run->transport->boundary_ledger();
    for (std::size_t slot = 0; slot < ledger.outward.size(); ++slot) {
        EXPECT_EQ(ledger.outward[slot], 0.0);
        EXPECT_EQ(ledger.gross[slot], 0.0);
    }
}

void run_grouped_properties ()
{
    Options options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    options.profile = Profile::Step;
    options.carrier_x = Real(1.0);
    options.varying_density = true;
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    double t = 0.0;
    for (int step = 0; step < 3; ++step) {
        t = step_once(*run, HostIntegrator::CompressibleRK3, t, 0.03);
        std::string diagnostic;
        EXPECT_TRUE(erf_sbm::authoritative_state_admissible(run->spectrum(),
                                                            run->state->layout(), 0, &diagnostic))
            << "step " << step << ": " << diagnostic;
    }
    EXPECT_GT(difference(run->spectrum(), *initial), Real(1.0e-3) * initial->max(0));
    for (int n = 0; n < initial->nComp(); ++n) {
        const Real before = mapped_inventory(*initial, run->transport->static_measure(0), n);
        const Real after = mapped_inventory(run->spectrum(), run->transport->static_measure(0), n);
        EXPECT_NEAR(after, before, Real(256.0) * eps() * before) << "component " << n;
    }
}

void run_missing_inputs ()
{
    auto run = make_run(Options{});
    begin_step(*run, 0.0);
    run->transport->advance_stage_from_host(
        0, HostIntegrator::CompressibleRK3, 0, 0.0, 0.0, 0.01, 0.01, *run->state, run->rho,
        run->rho, run->target, run->xmom, run->ymom, run->zmom, run->geom, 1, 2);
}

} // namespace

TEST(SBMMixingTransport, InactiveLimiterRecoversHostDiffusion) { run_recovery(); }
TEST(SBMMixingTransport, DiffusionKeepsRK3Order)
{
#ifdef AMREX_USE_FLOAT
    GTEST_SKIP() << "the finest third-order error is at single-precision roundoff";
#else
    run_temporal_order(HostIntegrator::CompressibleRK3, Real(2.7));
#endif
}
TEST(SBMMixingTransport, DiffusionKeepsHeunAndMidPointOrder)
{
    run_temporal_order(HostIntegrator::AnelasticHeun, Real(1.8));
    run_temporal_order(HostIntegrator::AnelasticMidPoint, Real(1.8));
}
TEST(SBMMixingTransport, HighOrderDiffusionCorrectionIsLimited)
{
    run_limited_correction(Real(0.0));
    run_limited_correction(Real(2.0));
}
// Without mixing too, a cell the limiter drains to its bound must not round
// past it and fail admission.
TEST(SBMMixingTransport, DrainedCellStaysAdmissibleWithoutMixing)
{
    Options options;
    options.policy = ParticleMixingPolicy::None;
    options.carrier_x = Real(2.0);
    run_limited_correction_with(options);
}
TEST(SBMMixingTransport, CombinedBoundRejectsSeparatelyAdmissibleTerms) { run_combined_bound(); }
TEST(SBMMixingTransport, LimitedCellsKeepTheTrialBaseExchange)
{
    run_limited_cells_keep_low_order_exchange();
}
TEST(SBMMixingTransport, ZeroCoefficientIsBitwiseAdvection) { run_zero_coefficient(); }
TEST(SBMMixingTransport, ConstantRatioSurvivesVaryingDensityAndCoefficients)
{
    run_constant_ratio();
}
TEST(SBMMixingTransport, WallsCloseTheInventoryWithoutBoundaryTransfer) { run_walls(); }
TEST(SBMMixingTransport, GroupedPropertiesStayAdmissibleAndConserved)
{
    run_grouped_properties();
}
TEST(SBMMixingTransport, MixingWithoutStageInputsAborts)
{
    expect_amrex_abort([] { run_missing_inputs(); },
                       "SBM particle mixing HostPassiveScalar received no stage mixing inputs");
}
