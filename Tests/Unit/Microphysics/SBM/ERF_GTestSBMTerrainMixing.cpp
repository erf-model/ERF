#include <gtest/gtest.h>

#include <AMReX.H>
#include <AMReX_BCRec.H>
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
#include "ERF_ScalarDiffusion.H"
#include "ERF_Utils.H"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Terrain-fitted (VariableDz) particle mixing through the SBM stage: the
// native terrain diffusion as the high-order rate, the matched diagonal
// exchange as the low-order one, the wall and lateral-face rules, and the
// stage's conservation, admissibility and timestep checks.

namespace {

using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;
using amrex::MultiFab;
using amrex::Real;
using erf_auxiliary::DiffusionFaceKind;
using erf_auxiliary::HostIntegrator;
using erf_sbm::ParticleMixingPolicy;

constexpr Real two_pi = Real(6.283185307179586476925286766559);

Real eps () { return std::numeric_limits<Real>::epsilon(); }

// Bin edge b of nbins geometric bins from 0.1 to 1; two bins have the edges
// 0.1, 0.5, 1 instead.
Real bin_edge (const int b, const int nbins)
{
    if (nbins == 2) { return b == 0 ? Real(0.1) : (b == 1 ? Real(0.5) : Real(1.0)); }
    return Real(0.1) * std::pow(Real(10.0), Real(b) / Real(nbins));
}

erf_sbm::SBMLayout
make_layout (const erf_sbm::MomentMode mode, const bool with_property, const int nbins = 2)
{
    erf_sbm::SpectralGridSpec grid;
    grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    grid.coordinate_units = "kg particle^-1";
    if (nbins == 2) {
        grid.edges = {Real(0.1), Real(0.5), Real(1.0)};
        grid.pivots = {Real(0.3), Real(0.75)};
    } else {
        for (int b = 0; b <= nbins; ++b) { grid.edges.push_back(bin_edge(b, nbins)); }
        for (int b = 0; b < nbins; ++b) {
            grid.pivots.push_back(std::sqrt(bin_edge(b, nbins) * bin_edge(b + 1, nbins)));
        }
    }

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

// Smooth: positive and smooth in physical space. Step: a jump across
// physical height one half. DryEndpoint: in the last 2M bin, endpoint
// numbers of physical height alone, the low one empty over half of each
// wavelength and the high one wet everywhere, under one group limiter.
// Stripe: in the last 2M bin, a dry low endpoint max(cos 2 pi z, 0)^8 and a
// wet high endpoint 1 + cos(2 pi z) / 5, both of physical height alone with
// zero gradient at flat walls z = 0, 1. EdgeBin: Smooth, except that the
// interior 2M bin 0 is nearly empty and its mean mass lies within 1e-10 to
// 2e-9 (relative) of its strict upper edge.
enum class Profile { Smooth, Step, Uniform, Empty, DryEndpoint, Stripe, EdgeBin };

// Hill: a terrain-following hill decaying to a flat top. Wave: flat bottom
// and top walls with interior coordinate surfaces displaced by
// hill sin(2 pi x / lx) sin(pi zeta), so a field of height alone satisfies
// the walls' zero normal flux. Witness: the numerical analysis's curved grid,
// zeta + hill sin(2 pi xi) (1 - zeta) + plane sin(2 pi zeta).
enum class Terrain { Hill, Wave, Witness };

struct Options
{
    int nx{8};
    int ny{4};
    // Physical extents in x and y; the computational zeta runs over [0, 1].
    Real lx{1};
    Real ly{1};
    Terrain terrain{Terrain::Hill};
    erf_sbm::MomentMode mode{erf_sbm::MomentMode::OneMoment};
    // Bins of the single population; two have the fixed edges 0.1, 0.5, 1.
    int nbins{2};
    bool attached_property{false};
    // Positive: every property is this fraction of its carrier's number.
    Real property_ratio{0};
    ParticleMixingPolicy policy{ParticleMixingPolicy::HostPassiveScalar};
    Real kh{0};
    Real kv{0};
    Real eddy_scale{1};
    // HostPassiveScalar with Scalar_h = eddy_scale (1 + 0.3 sin cos) varying in
    // x and y, and Scalar_v = 0.
    bool horizontal_eddy_only{false};
    int nz{8};
    // Boxes are at most max_grid_xy x max_grid_xy x max_grid_z cells.
    int max_grid_xy{4};
    int max_grid_z{8};
    // Hill terrain z = zeta + (1 - zeta) hill (sin 2 pi xi + cos(2 pi eta) / 2) + plane xi,
    // with xi = x / lx and eta = y / ly. Wave terrain adds hill_y sin(2 pi eta) sin(pi zeta).
    Real hill{Real(0.1)};
    Real hill_y{Real(0.0)};
    Real plane{Real(0.0)};
    bool varying_density{false};
    // Map factors varying in x and y about 1, isotropic as ERF's are.
    bool variable_maps{false};
    Real carrier_x{0};
    // The x carrier falls linearly from carrier_x at xhi to -carrier_x at xlo.
    bool carrier_ramp{false};
    Profile profile{Profile::Smooth};
    // x is bounded by an inflow face at xlo, with this diffusive kind, and an
    // outflow face at xhi; otherwise x is periodic. z is bounded by walls.
    bool inflow_x{false};
    DiffusionFaceKind xlo_diffusion{DiffusionFaceKind::Dirichlet};
    bool match{true};
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
    MultiFab z_nd, z_cc, ax, ay, detj;
    MultiFab xmom, ymom, zmom;
    erf_sbm::SBMStageMixingInputs inputs;
    std::vector<Real> inflow_z;

    [[nodiscard]] const erf_sbm::SBMStageMixingInputs* mixing () const
    {
        return transport->mixes_particles() ? &inputs : nullptr;
    }
    [[nodiscard]] const MultiFab& spectrum () const { return state->new_state(0); }
    [[nodiscard]] int nz () const { return options.nz; }
};

void fill_heights (MultiFab& z_nd, const Options& options)
{
    const Real hill = options.hill;
    const Real hill_y = options.hill_y;
    const Real plane = options.plane;
    const bool wave = options.terrain == Terrain::Wave;
    const bool witness = options.terrain == Terrain::Witness;
    const Real nx = Real(options.nx);
    const Real ny = Real(options.ny);
    const Real nz = Real(options.nz);
    // The periodic parts take the node index modulo the period, so both copies
    // of a periodic node hold one height, as ERF's OverrideSync of z_phys_nd
    // gives them; the plane is translation invariant and takes the true index.
    const int px = options.nx;
    const int py = options.ny;
    for (amrex::MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
        const auto a = z_nd.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real xi = Real(((i % px) + px) % px) / nx;
            const Real eta = Real(((j % py) + py) % py) / ny;
            const Real zeta = Real(k) / nz;
            if (wave) {
                a(i, j, k, 0) = zeta + (hill * std::sin(two_pi * xi) +
                                        hill_y * std::sin(two_pi * eta)) *
                                           std::sin(Real(0.5) * two_pi * zeta);
                return;
            }
            if (witness) {
                a(i, j, k, 0) = zeta + hill * std::sin(two_pi * xi) * (Real(1.0) - zeta) +
                                plane * std::sin(two_pi * zeta);
                return;
            }
            const Real surface = hill * (std::sin(two_pi * xi) + Real(0.5) * std::cos(two_pi * eta));
            a(i, j, k, 0) = zeta + (Real(1.0) - zeta) * surface + plane * Real(i) / nx;
        });
    }
}

void fill_density (MultiFab& rho, const bool varying, const int nx, const int ny)
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

// Anisotropic eddy diffusivities; every other channel holds a distinct value
// that a wrong channel would pick up.
void fill_eddy (MultiFab& eddy, const Real scale, const int nx, const int ny,
                const bool horizontal_only)
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
                scale * (Real(0.002) + Real(0.001) * std::cos(two_pi * Real(i) / Real(nx)) +
                         Real(0.0005) * Real(k));
            if (horizontal_only) {
                a(i, j, k, EddyDiff::Scalar_h) =
                    scale * (Real(1.0) + Real(0.3) * std::sin(two_pi * Real(i) / Real(nx)) *
                                             std::cos(two_pi * Real(j) / Real(ny)));
                a(i, j, k, EddyDiff::Scalar_v) = Real(0.0);
            }
        });
    }
}

// U = rho z for the selected profile, in physical height; 2M bins keep one
// mean mass per region, 1M bins their pivot, and attached properties a
// fraction of the number. With a positive property_ratio every attached
// property is that fraction of its carrier's number in every cell, an exactly
// proportional state.
void fill_spectrum (MultiFab& spectrum, const MultiFab& rho, const MultiFab& z_cc,
                    const erf_sbm::SBMLayout& layout, const Profile profile,
                    const int nx, const int ny, const Real property_ratio)
{
    const auto& population = layout.populations().front();
    const bool two_moment = population.moment_mode == erf_sbm::MomentMode::TwoMoment;
    const int nbins = population.grid.nbins();
    const int mass0 = population.mass_offset;
    const int number0 = two_moment ? population.number_offset : -1;
    const int property0 = layout.attached_properties().empty() ? -1 : layout.property_offset(0);
    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const auto u = spectrum.array(mfi);
        const auto r = rho.const_array(mfi);
        const auto height = z_cc.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real x = (Real(i) + Real(0.5)) / Real(nx);
            const Real y = (Real(j) + Real(0.5)) / Real(ny);
            const Real z = height(i, j, k);
            const Real density = r(i, j, k, Rho_comp);
            for (int bin = 0; bin < nbins; ++bin) {
                Real factor = Real(1.0);
                bool upper = false;
                if (profile == Profile::Smooth || profile == Profile::EdgeBin) {
                    factor = Real(1.0) + Real(0.3) * std::sin(two_pi * x) * std::cos(two_pi * y) +
                             Real(0.2) * std::cos(Real(0.5) * two_pi * z);
                    if (profile == Profile::EdgeBin && two_moment && bin == 0) {
                        // Bin 0 spans [0.1, 0.5]: a mean of 0.5 (1 - gap).
                        const Real number = Real(1.0e-12) * factor;
                        const Real gap = Real(1.0e-9) * (Real(1.0) + Real(0.9) * std::cos(two_pi * x));
                        u(i, j, k, mass0) = density * number * Real(0.5) * (Real(1.0) - gap);
                        u(i, j, k, number0) = density * number;
                        if (property0 >= 0) {
                            u(i, j, k, property0) = density * number * Real(0.4);
                        }
                        continue;
                    }
                } else if (profile == Profile::Step) {
                    upper = z > Real(0.5);
                    factor = upper ? Real(0.1) : Real(1.0);
                } else if (profile == Profile::Empty) {
                    factor = Real(0.0);
                } else if ((profile == Profile::DryEndpoint || profile == Profile::Stripe) &&
                           two_moment && bin == 1) {
                    // The endpoint numbers at the bin edges 0.5 and 1.
                    const Real c = std::cos(two_pi * z);
                    const Real c2 = c * c;
                    const Real dry = profile == Profile::Stripe ? c2 * c2 * c2 * c2 : c2 * c2;
                    const Real low = Real(0.01) * (c > Real(0.0) ? dry : Real(0.0));
                    const Real high = Real(0.01) * (Real(1.0) + Real(0.2) *
                        (profile == Profile::Stripe ? c : std::sin(two_pi * z + Real(0.7))));
                    u(i, j, k, mass0 + bin) = density * (Real(0.5) * low + high);
                    u(i, j, k, number0 + bin) = density * (low + high);
                    if (property0 >= 0) {
                        u(i, j, k, property0 + bin) = density * property_ratio * (low + high);
                    }
                    continue;
                }
                const Real number = factor * Real(0.01) * (Real(1.0) + Real(0.5) * Real(bin));
                // A 1M bin's particles sit at its pivot, as do those of more
                // than two bins.
                const Real pivot = nbins == 2 ? (bin == 0 ? Real(0.3) : Real(0.75))
                    : std::sqrt(Real(0.01) * std::pow(Real(10.0), Real(2 * bin + 1) / Real(nbins)));
                const Real mean = !two_moment || nbins != 2 ? pivot
                                : bin == 0 ? (upper ? Real(0.4) : Real(0.2)) : Real(0.75);
                u(i, j, k, mass0 + bin) = density * number * mean;
                if (number0 >= 0) { u(i, j, k, number0 + bin) = density * number; }
                if (property0 >= 0) {
                    u(i, j, k, property0 + bin) = density * number *
                        (property_ratio > Real(0.0) ? property_ratio
                                                    : (upper ? Real(0.1) : Real(0.4)));
                }
            }
        });
    }
}

// A canonical, admissible inflow spectrum z = U / rho_d.
std::vector<Real> inflow_state (const erf_sbm::SBMLayout& layout)
{
    const auto& population = layout.populations().front();
    const bool two_moment = population.moment_mode == erf_sbm::MomentMode::TwoMoment;
    const int property0 = layout.ncomp() > (two_moment ? 4 : 2) ? layout.property_offset(0) : -1;
    std::vector<Real> z(static_cast<std::size_t>(layout.ncomp()), Real(0.0));
    for (int bin = 0; bin < 2; ++bin) {
        const Real number = Real(1.2) * Real(0.01) * (Real(1.0) + Real(0.5) * Real(bin));
        z[static_cast<std::size_t>(population.mass_offset + bin)] =
            number * (bin == 0 ? Real(0.3) : Real(0.75));
        if (two_moment) { z[static_cast<std::size_t>(population.number_offset + bin)] = number; }
        if (property0 >= 0) { z[static_cast<std::size_t>(property0 + bin)] = Real(0.25) * number; }
    }
    return z;
}

std::unique_ptr<Run>
make_run (const Options& options)
{
    auto run = std::make_unique<Run>();
    run->options = options;
    const int nz = options.nz;
    const int nx = options.nx;
    const int ny = options.ny;
    const Box domain(IntVect(0), IntVect(nx - 1, ny - 1, nz - 1));
    const amrex::RealBox real_box({0.0, 0.0, 0.0},
                                  {double(options.lx), double(options.ly), 1.0});
    const int periodic[AMREX_SPACEDIM] = {options.inflow_x ? 0 : 1, 1, 0};
    run->geom = Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodic);
    run->ba = BoxArray(domain);
    run->ba.maxSize(IntVect(options.max_grid_xy, options.max_grid_xy, options.max_grid_z));
    run->dm = DistributionMapping(run->ba);
    const auto& ba = run->ba;
    const auto& dm = run->dm;

    run->state = std::make_unique<erf_sbm::SBMStateManager>(
        make_layout(options.mode, options.attached_property, options.nbins), 1);
    run->state->define(0, ba, dm, 0.0);
    run->transport = std::make_unique<erf_sbm::SBMTransport>(run->state->layout(), 1, 16);
    run->transport->define(0, ba, dm);
    if (options.policy != ParticleMixingPolicy::None) {
        erf_sbm::ParticleMixingConfig config;
        config.policy = options.policy;
        config.reference_kh = options.kh;
        config.reference_kv = options.kv;
        run->transport->set_particle_mixing(config, erf_auxiliary::DiffusionMesh::VariableDz,
                                            {true, options.match});
    }
    erf_sbm::SBMPhysicalBoundaryConfig boundaries;
    boundaries.kind[2] = erf_sbm::AdvectionFaceKind::ImpermeableWall;
    boundaries.kind[5] = erf_sbm::AdvectionFaceKind::ImpermeableWall;
    if (options.inflow_x) {
        run->inflow_z = inflow_state(run->state->layout());
        boundaries.kind[0] = erf_sbm::AdvectionFaceKind::PrescribedInflow;
        boundaries.inflow_z[0] = run->inflow_z;
        boundaries.diffusion[0] = options.xlo_diffusion;
        boundaries.kind[3] = erf_sbm::AdvectionFaceKind::AdvectiveOutflow;
    }
    run->transport->set_physical_boundaries(boundaries);

    amrex::BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    const BoxArray ba2d(std::move(boxes));
    run->mapfac.resize(MapFacType::num);
    run->mapfac[MapFacType::m_x] = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    run->mapfac[MapFacType::u_x] = std::make_unique<MultiFab>(
        amrex::convert(ba2d, IntVect(1, 0, 0)), dm, 1, 0);
    run->mapfac[MapFacType::v_x] = std::make_unique<MultiFab>(
        amrex::convert(ba2d, IntVect(0, 1, 0)), dm, 1, 0);
    for (int m = 0; m < MapFacType::num; ++m) {
        const Real amplitude = options.variable_maps ? Real(0.05) : Real(0.0);
        const Real nxr = Real(nx);
        const Real nyr = Real(ny);
        for (amrex::MFIter mfi(*run->mapfac[m]); mfi.isValid(); ++mfi) {
            const auto a = run->mapfac[m]->array(mfi);
            // Periodic in the index, so both copies of a periodic face agree.
            amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
                const int ip = ((i + m) % nx + nx) % nx;
                const int jp = (j % ny + ny) % ny;
                a(i, j, 0, 0) = Real(1.0) + amplitude * std::sin(two_pi * Real(ip) / nxr) +
                                Real(0.5) * amplitude * std::cos(two_pi * Real(jp) / nyr);
            });
        }
    }

    run->z_nd.define(amrex::convert(ba, IntVect(1)), dm, 1, 2);
    fill_heights(run->z_nd, options);
    // ERF's cell-centre heights, averaged from the analytic nodes in every
    // ghost cell too: the tilted slab is translation invariant but not
    // periodic, so its periodic ghost columns must not take the images'.
    run->z_cc.define(ba, dm, 1, 2);
    for (amrex::MFIter mfi(run->z_cc); mfi.isValid(); ++mfi) {
        const auto zc = run->z_cc.array(mfi);
        const auto zn = run->z_nd.const_array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            zc(i, j, k) = Compute_Z_AtCellCenter(i, j, k, zn);
        });
    }
    run->detj.define(ba, dm, 1, 0);
    make_J(run->geom, run->z_nd, run->detj);
    run->ax.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    run->ay.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    MultiFab az(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
    make_areas(run->geom, run->z_nd, run->ax, run->ay, az);
    std::string diagnostic;
    EXPECT_TRUE(run->transport->rebuild_static_measure(0, run->detj, *run->mapfac[MapFacType::m_x],
                                                       *run->mapfac[MapFacType::m_x], diagnostic))
        << diagnostic;

    run->rho.define(ba, dm, 3, 2);
    fill_density(run->rho, options.varying_density, nx, ny);
    run->rho.FillBoundary(run->geom.periodicity());
    run->target.define(ba, dm, 3, 2);
    MultiFab::Copy(run->target, run->rho, 0, 0, 3, 2);
    run->eddy.define(ba, dm, EddyDiff::NumDiffs, 2);
    fill_eddy(run->eddy, options.eddy_scale, nx, ny, options.horizontal_eddy_only);
    run->eddy.FillBoundary(run->geom.periodicity());

    run->xmom.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    run->ymom.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    run->zmom.define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
    run->xmom.setVal(options.carrier_x);
    if (options.carrier_ramp) {
        const Real c = options.carrier_x;
        for (amrex::MFIter mfi(run->xmom); mfi.isValid(); ++mfi) {
            const auto a = run->xmom.array(mfi);
            amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                a(i, j, k, 0) = c * (Real(2 * i) / Real(nx) - Real(1.0));
            });
        }
    }
    run->ymom.setVal(Real(0.0));
    run->zmom.setVal(Real(0.0));

    fill_spectrum(run->state->new_state_for_initialization(0), run->rho, run->z_cc,
                  run->state->layout(), options.profile, nx, ny, options.property_ratio);
    run->inputs = erf_sbm::make_stage_mixing_inputs(*run->transport, &run->eddy, run->mapfac,
                                                    &run->ax, &run->ay, nullptr, &run->z_nd,
                                                    &run->z_cc);
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

double step_once (Run& run, const HostIntegrator method, const double t, const double dt)
{
    const double t_new = t + dt;
    begin_step(run, t);
    if (method == HostIntegrator::CompressibleRK3) {
        stage(run, method, 0, t, t, t + dt / 3.0, dt / 3.0);
        stage(run, method, 1, t, t + dt / 3.0, t + dt / 2.0, dt / 2.0);
        stage(run, method, 2, t, t + dt / 2.0, t_new, dt);
    } else {
        stage(run, method, 0, t, t, t_new, dt);
        stage(run, method, 1, t, t_new, t_new, dt);
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
    // Separate runs can have their own box arrays and distribution maps.
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

// dU/dt of ERF's native terrain diffusion of z = U / rho under SBM's wall
// rule, written from the library builders with this file's own quadratic
// extension past the walls and wall suppression, for a domain periodic in x
// and y: an oracle independent of the transport's terrain wrapper.
void native_tendency (const Run& run, const MultiFab& state, MultiFab& tendency)
{
    AMREX_ALWAYS_ASSERT(!run.options.inflow_x);
    const int ncomp = state.nComp();
    const int nz = run.nz();
    MultiFab z(run.ba, run.dm, ncomp, 2);
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
    for (amrex::MFIter mfi(z); mfi.isValid(); ++mfi) {
        const auto q = z.array(mfi);
        const auto h = run.z_cc.const_array(mfi);
        const Box fab = mfi.fabbox();
        for (const int ghost : {-1, nz}) {
            if (ghost < fab.smallEnd(2) || ghost > fab.bigEnd(2)) continue;
            Box layer = fab;
            layer.setRange(2, ghost, 1);
            const int a = ghost == -1 ? 0 : nz - 1;
            const int b = ghost == -1 ? 1 : nz - 2;
            const int c = ghost == -1 ? 2 : nz - 3;
            amrex::ParallelFor(layer, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int, int n) noexcept {
                const Real x = h(i, j, ghost);
                const Real xa = h(i, j, a), xb = h(i, j, b), xc = h(i, j, c);
                q(i, j, ghost, n) =
                    q(i, j, a, n) * (x - xb) * (x - xc) / ((xa - xb) * (xa - xc)) +
                    q(i, j, b, n) * (x - xa) * (x - xc) / ((xb - xa) * (xb - xc)) +
                    q(i, j, c, n) * (x - xa) * (x - xb) / ((xc - xa) * (xc - xb));
            });
        }
    }

    erf_auxiliary::MappedDiffusionCoefficients c;
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::particle_mixing_coefficients(run.transport->particle_mixing(), c,
                                                       diagnostic)) << diagnostic;
    ScalarDiffusionFluxPolicy policy;
    policy.coefficients = c.horizontal;
    policy.coefficient_mode = c.mode;
    amrex::BCRec bc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int type = run.geom.isPeriodic(d) ? ERFBCType::int_dir : ERFBCType::foextrap;
        bc.setLo(d, type);
        bc.setHi(d, type);
    }
    amrex::Gpu::DeviceVector<amrex::BCRec> bcs(1);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, &bc, &bc + 1, bcs.begin());
    MultiFab raw_x(amrex::convert(run.ba, IntVect(1, 0, 0)), run.dm, 1, IntVect(0, 0, 1));
    MultiFab raw_y(amrex::convert(run.ba, IntVect(0, 1, 0)), run.dm, 1, IntVect(0, 0, 1));
    MultiFab raw_z(amrex::convert(run.ba, IntVect(0, 0, 1)), run.dm, 1, 0);
    erf_auxiliary::MappedFaceFluxRate mapped;
    mapped.define(run.ba, run.dm, 1, 0);
    const auto inv = run.geom.InvCellSizeArray();
    const auto& m = *run.mapfac[MapFacType::m_x];
    const auto& u = *run.mapfac[MapFacType::u_x];
    const auto& v = *run.mapfac[MapFacType::v_x];
    for (int n = 0; n < ncomp; ++n) {
        for (amrex::MFIter mfi(z); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            ScalarDiffusionFieldViews field;
            field.scalar = z.const_array(mfi);
            field.scalar_comp = n;
            field.density = run.rho.const_array(mfi);
            field.rho_comp = Rho_comp;
            field.mu_turb = run.eddy.const_array(mfi);
            field.xflux = raw_x.array(mfi);
            field.yflux = raw_y.array(mfi);
            field.zflux = raw_z.array(mfi);
            BuildScalarDiffusionFluxes_T(bx, run.geom.Domain(), field, policy, c.vertical, false,
                                         run.z_nd.const_array(mfi), run.z_cc.const_array(mfi),
                                         inv, u.const_array(mfi), v.const_array(mfi),
                                         bcs.data(), 0);
            BuildScalarDiffusionMappedTransfers_T(
                bx, run.geom.Domain(), raw_x.const_array(mfi), raw_y.const_array(mfi),
                raw_z.const_array(mfi), 0, mapped.dir(0).array(mfi), 0,
                mapped.dir(1).array(mfi), 0, mapped.dir(2).array(mfi), 0,
                run.z_nd.const_array(mfi), run.ax.const_array(mfi), run.ay.const_array(mfi),
                inv, m.const_array(mfi), u.const_array(mfi), m.const_array(mfi),
                v.const_array(mfi), true, true);
        }
        erf_auxiliary::ApplyMappedFluxTendency(mapped, 0, run.transport->static_measure(0), 0,
                                               tendency, n, inv);
    }
    amrex::Gpu::streamSynchronize();
}

// a + tau * L(b)
std::unique_ptr<MultiFab>
oracle_update (const Run& run, const MultiFab& a, const MultiFab& b, const Real tau)
{
    MultiFab tendency(run.ba, run.dm, a.nComp(), 0);
    native_tendency(run, b, tendency);
    auto out = copy_of(a);
    MultiFab::Saxpy(*out, tau, tendency, 0, 0, a.nComp(), 0);
    return out;
}

std::unique_ptr<MultiFab>
oracle_rk3_step (const Run& run, const MultiFab& u0, const Real dt)
{
    auto u1 = oracle_update(run, u0, u0, dt / Real(3.0));
    auto u2 = oracle_update(run, u0, *u1, dt / Real(2.0));
    return oracle_update(run, u0, *u2, dt);
}

// Forked death tests hang under MPI, so amrex::Abort throws for the check.
template <typename Action>
std::string abort_message (Action&& action)
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
    return message;
}

void expect_conserved (const MultiFab& initial, const Run& run)
{
    for (int n = 0; n < initial.nComp(); ++n) {
        const Real before = mapped_inventory(initial, run.transport->static_measure(0), n);
        const Real after = mapped_inventory(run.spectrum(), run.transport->static_measure(0), n);
        EXPECT_NEAR(after, before, Real(256.0) * eps() * before) << "component " << n;
    }
}

// ---------------------------------------------------------------------------

// With the limiter inactive every accepted rate is the native terrain
// diffusion of the stage input, so one RK3 step reproduces RK3 applied to
// that operator, matched or not. This includes K_v = 0, where the only
// vertical transfer is horizontal diffusion along the slopes.
void expect_recovers_native_diffusion (Options options)
{
    const Real dt = Real(0.05);
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    advance(*run, HostIntegrator::CompressibleRK3, dt, 1);
    const auto expected = oracle_rk3_step(*run, *initial, dt);
    const Real change = difference(*expected, *initial);
    const Real scale = initial->norm0(0, initial->nComp(), IntVect(0));
    EXPECT_GT(change, Real(1.0e-4) * scale);
    EXPECT_LE(difference(run->spectrum(), *expected), Real(1024.0) * eps() * scale);

    options.match = false;
    auto unmatched = make_run(options);
    advance(*unmatched, HostIntegrator::CompressibleRK3, dt, 1);
    EXPECT_LE(difference(unmatched->spectrum(), *expected), Real(1024.0) * eps() * scale);
}

void run_recovery ()
{
    Options eddy;
    eddy.varying_density = true;
    eddy.variable_maps = true;
    eddy.mode = erf_sbm::MomentMode::TwoMoment;
    expect_recovers_native_diffusion(eddy);

    // No attached property: at a uniform property ratio its donor envelope has
    // no width, and roundoff in that support can keep the limiter active.
    Options horizontal_only;
    horizontal_only.policy = ParticleMixingPolicy::ConstantReference;
    horizontal_only.kh = Real(0.02);
    horizontal_only.kv = Real(0.0);
    expect_recovers_native_diffusion(horizontal_only);
}

// A step profile with a carrier, a 2M population and an attached property:
// every group stays admissible and every inventory is conserved, so the walls
// carry nothing, metric terms included.
void run_conservation ()
{
    Options options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    options.profile = Profile::Step;
    options.carrier_x = Real(1.0);
    options.varying_density = true;
    options.hill = Real(0.15);
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    double t = 0.0;
    for (int step = 0; step < 4; ++step) {
        t = step_once(*run, HostIntegrator::CompressibleRK3, t, 0.03);
        std::string diagnostic;
        EXPECT_TRUE(erf_sbm::authoritative_state_admissible(run->spectrum(),
                                                            run->state->layout(), 0, &diagnostic))
            << "step " << step << ": " << diagnostic;
    }
    EXPECT_GT(difference(run->spectrum(), *initial), Real(1.0e-3) * initial->max(0));
    expect_conserved(*initial, *run);
    const auto& ledger = run->transport->boundary_ledger();
    for (std::size_t slot = 0; slot < ledger.outward.size(); ++slot) {
        EXPECT_EQ(ledger.outward[slot], 0.0);
        EXPECT_EQ(ledger.gross[slot], 0.0);
    }
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
    options.carrier_x = Real(1.5);
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    advance(*run, HostIntegrator::CompressibleRK3, Real(0.02), 2);
    EXPECT_LE(difference(run->spectrum(), *initial),
              Real(256.0) * eps() * initial->norm0(0, initial->nComp(), IntVect(0)));
}

// Boxes split vertically read their raw side fluxes one layer beyond their
// own cells, and the lateral face states with them: the same run on one, two
// and four boxes per column agrees, through an inflow face beside the walls.
void run_split_boxes ()
{
    Options options;
    options.inflow_x = true;
    options.plane = Real(0.05);
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    options.carrier_x = Real(0.5);
    options.profile = Profile::Smooth;
    auto whole = make_run(options);
    advance(*whole, HostIntegrator::CompressibleRK3, Real(0.02), 2);
    const Real scale = whole->spectrum().norm0(0, whole->spectrum().nComp(), IntVect(0));
    for (const int height : {4, 2}) {
        options.max_grid_z = height;
        auto split = make_run(options);
        advance(*split, HostIntegrator::CompressibleRK3, Real(0.02), 2);
        EXPECT_LE(difference(split->spectrum(), whole->spectrum()), Real(16.0) * eps() * scale)
            << "boxes of height " << height;
    }
}

// Into an empty domain at rest on the planar terrain z = zeta + s xi, one
// stage gives the edge cells of the inflow face its face-trace transfer,
// (8/3) K z_in dx_inv per unit area, and not the two-point exchange with the
// trace as a ghost centre. Its interpolation to the k faces of that column,
// -s / 2 of it on every interior k face, moves material down the column; the
// walls suppress their whole transfer, so only the bottom cell gains and the
// top cell loses that cross term. All other cells stay empty.
void run_inflow_face_trace ()
{
    Options options;
    options.inflow_x = true;
    options.hill = Real(0.0);
    options.plane = Real(0.1);
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.01);
    options.kv = Real(0.004);
    options.profile = Profile::Empty;
    auto run = make_run(options);
    const Real tau = Real(0.01);
    begin_step(*run, 0.0);
    stage(*run, HostIntegrator::AnelasticHeun, 0, 0.0, 0.0, tau, tau);

    const int nz = run->nz();
    const int ncomp = run->spectrum().nComp();
    std::vector<Real> host(static_cast<std::size_t>(ncomp));
    for (int n = 0; n < ncomp; ++n) {
        // Raw face flux; unit maps, density, Jacobian and side areas.
        host[static_cast<std::size_t>(n)] = (Real(8.0) / Real(3.0)) * options.kh *
                                            run->inflow_z[static_cast<std::size_t>(n)] *
                                            Real(options.nx);
    }
    amrex::Gpu::DeviceVector<Real> face_flux(host.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, host.begin(), host.end(), face_flux.begin());
    const Real* flux = face_flux.data();
    const Real s = options.plane;
    const int nx = options.nx;
    MultiFab expected(run->ba, run->dm, ncomp, 0);
    for (amrex::MFIter mfi(expected); mfi.isValid(); ++mfi) {
        const auto a = expected.array(mfi);
        amrex::ParallelFor(mfi.validbox(), ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                if (i != 0) {
                    a(i, j, k, n) = Real(0.0);
                    return;
                }
                const Real cross = Real(0.5) * s * flux[n] * Real(nz);
                const Real side = flux[n] * Real(nx);
                a(i, j, k, n) = tau * (side + (k == 0 ? cross : Real(0.0)) -
                                       (k == nz - 1 ? cross : Real(0.0)));
            });
    }
    const Real scale = *std::max_element(host.begin(), host.end()) * Real(nx) * tau;
    EXPECT_LE(difference(run->spectrum(), expected), Real(64.0) * eps() * scale);
}

// Diffusion alone carries material through the inflow face beside the walls,
// and the ledger records it: the inventory gains exactly what the xlo face
// recorded inward, and the outflow face, closed to mixing, records nothing.
void run_inflow_ledger ()
{
    Options options;
    options.inflow_x = true;
    options.plane = Real(0.05);
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.02);
    options.kv = Real(0.005);
    options.profile = Profile::Uniform;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    double t = 0.0;
    for (int step = 0; step < 3; ++step) {
        t = step_once(*run, HostIntegrator::CompressibleRK3, t, 0.02);
        std::string diagnostic;
        EXPECT_TRUE(erf_sbm::authoritative_state_admissible(run->spectrum(),
                                                            run->state->layout(), 0, &diagnostic))
            << "step " << step << ": " << diagnostic;
    }
    const auto& ledger = run->transport->boundary_ledger();
    const int ncomp = initial->nComp();
    // The ledger holds physical amounts: the inventory with the cell volume.
    const Real volume = Real(1.0) / Real(options.nx * options.ny * run->nz());
    for (int n = 0; n < ncomp; ++n) {
        const Real before =
            volume * mapped_inventory(*initial, run->transport->static_measure(0), n);
        const Real after =
            volume * mapped_inventory(run->spectrum(), run->transport->static_measure(0), n);
        const double xlo = ledger.outward[static_cast<std::size_t>(0 * ncomp + n)];
        const double xhi = ledger.outward[static_cast<std::size_t>(3 * ncomp + n)];
        EXPECT_NE(xlo, 0.0) << "component " << n;
        EXPECT_NEAR(after - before, -xlo, Real(256.0) * eps() * before) << "component " << n;
        EXPECT_EQ(xhi, 0.0) << "component " << n;
        EXPECT_EQ(ledger.gross[static_cast<std::size_t>(3 * ncomp + n)], 0.0) << "component " << n;
    }
}

// Where the carrier leaves through an inflow_outflow face, the face is closed
// to mixing before its side flux reaches the k faces, as a closed face is.
void run_inflow_outflow_masks_outward_faces ()
{
    Options options;
    options.inflow_x = true;
    options.plane = Real(0.05);
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.02);
    options.kv = Real(0.005);
    options.carrier_x = Real(0.5);
    options.carrier_ramp = true;
    options.profile = Profile::Smooth;
    options.xlo_diffusion = DiffusionFaceKind::DirichletWhereInward;
    auto inward_only = make_run(options);
    options.xlo_diffusion = DiffusionFaceKind::Closed;
    auto closed = make_run(options);
    options.xlo_diffusion = DiffusionFaceKind::Dirichlet;
    auto open = make_run(options);
    for (auto* run : {inward_only.get(), closed.get(), open.get()}) {
        advance(*run, HostIntegrator::CompressibleRK3, Real(0.02), 2);
    }
    EXPECT_EQ(difference(inward_only->spectrum(), closed->spectrum()), Real(0.0));
    EXPECT_GT(difference(open->spectrum(), closed->spectrum()),
              Real(1.0e-6) * closed->spectrum().norm0(0, closed->spectrum().nComp(), IntVect(0)));
}

// A step beyond the combined bound, whose terrain conductances keep their
// slope terms, fails before the stage changes any state.
void run_unsafe_step ()
{
    Options options;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.5);
    options.kv = Real(0.0);
    options.hill = Real(0.15);
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    begin_step(*run, 0.0);
    const std::string message = abort_message([&] {
        stage(*run, HostIntegrator::CompressibleRK3, 0, 0.0, 0.0, 0.1, 0.1);
    });
    EXPECT_NE(message.find("SBM combined advection-diffusion outgoing-demand condition failed"),
              std::string::npos) << message;
    EXPECT_EQ(difference(run->state->old_state(0), *initial), Real(0.0));
}

// The stage rejects terrain mixing it is not qualified for, before any state
// changes, and the transport rejects a terrain mesh outside qualification.
void run_terrain_envelope ()
{
    const auto layout = make_layout(erf_sbm::MomentMode::OneMoment, false);
    erf_sbm::SBMTransport transport(layout, 1, 16);
    erf_sbm::ParticleMixingConfig config;
    config.policy = ParticleMixingPolicy::HostPassiveScalar;
    EXPECT_THROW(transport.set_particle_mixing(config, erf_auxiliary::DiffusionMesh::VariableDz),
                 std::invalid_argument);
    EXPECT_NO_THROW(transport.set_particle_mixing(
        config, erf_auxiliary::DiffusionMesh::VariableDz, {true, true}));

    // Two cells in a column leave no third sample for the wall extension.
    Options shallow;
    shallow.nz = 2;
    auto run = make_run(shallow);
    begin_step(*run, 0.0);
    const std::string message = abort_message([&] {
        stage(*run, HostIntegrator::CompressibleRK3, 0, 0.0, 0.0, 0.001, 0.001);
    });
    EXPECT_NE(message.find("requires three cells in each column"), std::string::npos) << message;

    // The cell-centre heights are required.
    auto missing = make_run(Options{});
    missing->inputs.z_cc = nullptr;
    begin_step(*missing, 0.0);
    const std::string heights = abort_message([&] {
        stage(*missing, HostIntegrator::CompressibleRK3, 0, 0.0, 0.0, 0.001, 0.001);
    });
    EXPECT_NE(heights.find("cell-centre heights"), std::string::npos) << heights;
}

// Endpoint numbers of the last 2M bin, whose edges are 0.5 and 1:
// low = 2 N - 2 M and high = 2 M - N.
void bin1_endpoints (const Run& run, const MultiFab& u, MultiFab& endpoints)
{
    const auto& population = run.state->layout().populations().front();
    const int mass = population.mass_offset + 1;
    const int number = population.number_offset + 1;
    MultiFab::LinComb(endpoints, Real(2.0), u, number, Real(-2.0), u, mass, 0, 1, 0);
    MultiFab::LinComb(endpoints, Real(2.0), u, mass, Real(-1.0), u, number, 1, 1, 0);
}

// With K_v = 0 a spectrum of physical height alone is stationary on a slope,
// whose bottom, top and closed sides carry no physical flux either. Here A = 4,
// S = 1 and 32 cells span a wavelength; the step is 0.8 of the combined bound
// and the run reaches (2 pi / L)^2 K_h t = 2. Sharing one limiter with the dry
// low endpoint, the plain diagonal exchange mixes the wet high endpoint across
// height by over a tenth of its amplitude; matched to the native rate, the
// exchange leaves it within one percent (0.48% measured), about the native
// operator's own error.
void run_dry_endpoint ()
{
    Options options;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(0.01);
    options.kv = Real(0.0);
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.profile = Profile::DryEndpoint;
    options.inflow_x = true;
    options.xlo_diffusion = DiffusionFaceKind::Closed;
    options.hill = Real(0.0);
    options.plane = Real(0.25);
    options.nz = 32;
    options.max_grid_z = 16;
    const Real time = Real(0.05) / options.kh;
    Real wet_error[2] = {};
    for (const bool match : {true, false}) {
        options.match = match;
        auto run = make_run(options);
        const auto initial = copy_of(run->spectrum());
        erf_auxiliary::MappedFaceTransmissibility t;
        t.define(run->ba, run->dm, 1, 0);
        std::string diagnostic;
        ASSERT_TRUE(erf_sbm::build_particle_mixing_transmissibility(
            run->transport->particle_mixing(), run->geom, run->inputs.metrics, run->rho,
            Rho_comp, nullptr, t, diagnostic)) << diagnostic;
        erf_auxiliary::MappedFaceFluxRate at_rest;
        at_rest.define(run->ba, run->dm, 1, 0);
        at_rest.setVal(Real(0.0));
        Real rate = Real(0.0);
        ASSERT_TRUE(erf_auxiliary::ComputeMaxMappedCombinedOutgoingRate(
            at_rest, t, run->transport->static_measure(0), run->rho, Rho_comp,
            run->geom.InvCellSizeArray(), rate, diagnostic)) << diagnostic;
        const int steps = static_cast<int>(std::ceil(time * rate / Real(0.8)));
        advance(*run, HostIntegrator::CompressibleRK3, time / Real(steps), steps);
        EXPECT_TRUE(erf_sbm::authoritative_state_admissible(run->spectrum(),
                                                            run->state->layout(), 0, &diagnostic))
            << diagnostic;
        MultiFab before(run->ba, run->dm, 2, 0);
        MultiFab after(run->ba, run->dm, 2, 0);
        bin1_endpoints(*run, *initial, before);
        bin1_endpoints(*run, run->spectrum(), after);
        MultiFab::Subtract(after, before, 1, 1, 1, 0);
        // The wet endpoint's amplitude is 0.2 of its mean 0.01, at unit density.
        wet_error[match ? 0 : 1] = after.norm0(1) / (Real(0.2) * Real(0.01));
    }
    amrex::Print() << "SBM_M6B_DRY_ENDPOINT wet_peak_error matched=" << wet_error[0]
                   << " unmatched=" << wet_error[1] << '\n';
    EXPECT_LT(wet_error[0], Real(0.01));
    EXPECT_GT(wet_error[1], Real(0.1));
}

// Errors of a stripe run against physical vertical diffusion with K_v at
// flat, no-flux walls: the wet endpoint's peak error over its initial
// amplitude, the dry endpoint's L1 error over its initial inventory, and the
// vertical diffusivity ratio r_hat inferred from the wet fundamental.
struct StripeErrors
{
    Real wet_peak{0};
    Real dry_l1{0};
    Real r_hat{0};
    // With an attached property, its largest drift from the exact ratio.
    Real ratio_drift{0};
};

StripeErrors
stripe_errors (const Run& run, const MultiFab& initial, const Real kh, const Real kv,
               const Real time)
{
    // The dry endpoint's cosine series on [0, 1] with zero end gradients,
    // each mode decaying as exp(-K_v (n pi)^2 t).
    constexpr int modes = 128;
    constexpr int samples = 32768;
    std::vector<Real> coefficient(modes, Real(0.0));
    for (int q = 0; q < samples; ++q) {
        const Real z = (Real(q) + Real(0.5)) / Real(samples);
        const Real c = std::cos(two_pi * z);
        const Real c2 = c * c;
        const Real dry = c > Real(0.0) ? c2 * c2 * c2 * c2 : Real(0.0);
        for (int n = 0; n < modes; ++n) {
            coefficient[static_cast<std::size_t>(n)] +=
                (n == 0 ? Real(1.0) : Real(2.0)) * dry *
                std::cos(Real(0.5) * two_pi * Real(n) * z) / Real(samples);
        }
    }
    for (int n = 0; n < modes; ++n) {
        const Real wave = Real(0.5) * two_pi * Real(n);
        coefficient[static_cast<std::size_t>(n)] *= std::exp(-kv * wave * wave * time);
    }
    amrex::Gpu::DeviceVector<Real> series(coefficient.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, coefficient.begin(), coefficient.end(),
                     series.begin());
    const Real* c_dry = series.data();
    const Real wet_decay = std::exp(-kv * two_pi * two_pi * time);

    MultiFab before(run.ba, run.dm, 2, 0);
    MultiFab after(run.ba, run.dm, 2, 0);
    bin1_endpoints(run, initial, before);
    bin1_endpoints(run, run.spectrum(), after);
    // Per cell, with the measure omega: |dry - exact| omega, dry_0 omega,
    // |wet - exact|, and the wet fundamental's projections at t and 0.
    MultiFab terms(run.ba, run.dm, 6, 0);
    const MultiFab& measure = run.transport->static_measure(0);
    for (amrex::MFIter mfi(terms); mfi.isValid(); ++mfi) {
        const auto out = terms.array(mfi);
        const auto e0 = before.const_array(mfi);
        const auto e1 = after.const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto height = run.z_cc.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real z = height(i, j, k);
            Real dry = Real(0.0);
            for (int n = 0; n < modes; ++n) {
                dry += c_dry[n] * std::cos(Real(0.5) * two_pi * Real(n) * z);
            }
            const Real basis = std::cos(two_pi * z);
            const Real wet = Real(0.01) * (Real(1.0) + Real(0.2) * wet_decay * basis);
            const Real w = omega(i, j, k, 0);
            out(i, j, k, 0) = std::abs(e1(i, j, k, 0) - Real(0.01) * dry) * w;
            out(i, j, k, 1) = e0(i, j, k, 0) * w;
            out(i, j, k, 2) = std::abs(e1(i, j, k, 1) - wet);
            out(i, j, k, 3) = (e1(i, j, k, 1) - Real(0.01)) * basis * w;
            out(i, j, k, 4) = (e0(i, j, k, 1) - Real(0.01)) * basis * w;
            out(i, j, k, 5) = Real(0.0);
        });
    }
    StripeErrors errors;
    errors.dry_l1 = terms.sum(0) / terms.sum(1);
    errors.wet_peak = terms.norm0(2) / (Real(0.2) * Real(0.01));
    const Real star = two_pi * two_pi * kh * time;
    errors.r_hat = -std::log(std::abs(terms.sum(3)) / std::abs(terms.sum(4))) / star;
    return errors;
}

// A stripe run on the flat-walled wave terrain with A = Delta_x / Delta_Z = 10,
// cells_per_wave cells per vertical wavelength, a worst cell drop just under
// one thickness, K_h = 1 and K_v = r, to (2 pi / L)^2 K_h t = 2 at 0.8 of
// the combined bound.
// Runs a stripe configuration to (2 pi / L)^2 K_h t = 2, K_h the mean
// horizontal diffusivity, at 0.8 of the combined bound, and returns its
// errors against vertical diffusion with K_v.
StripeErrors run_height_only (const Options& options, const Real kh, const Real kv,
                              const std::string& label)
{
    auto run = make_run(options);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::admit_terrain_mixing_geometry(run->z_nd, run->ba, Real(1.0), diagnostic))
        << diagnostic;
    const auto initial = copy_of(run->spectrum());
    erf_auxiliary::MappedFaceTransmissibility t;
    t.define(run->ba, run->dm, 1, 0);
    EXPECT_TRUE(erf_sbm::build_particle_mixing_transmissibility(
        run->transport->particle_mixing(), run->geom, run->inputs.metrics, run->rho, Rho_comp,
        run->inputs.eddy_diffusivity, t, diagnostic)) << diagnostic;
    erf_auxiliary::MappedFaceFluxRate at_rest;
    at_rest.define(run->ba, run->dm, 1, 0);
    at_rest.setVal(Real(0.0));
    Real rate = Real(0.0);
    EXPECT_TRUE(erf_auxiliary::ComputeMaxMappedCombinedOutgoingRate(
        at_rest, t, run->transport->static_measure(0), run->rho, Rho_comp,
        run->geom.InvCellSizeArray(), rate, diagnostic)) << diagnostic;
    const Real time = Real(2.0) / (two_pi * two_pi * kh);
    const int steps = static_cast<int>(std::ceil(time * rate / Real(0.8)));
    advance(*run, HostIntegrator::CompressibleRK3, time / Real(steps), steps);
    EXPECT_TRUE(erf_sbm::authoritative_state_admissible(run->spectrum(), run->state->layout(), 0,
                                                        &diagnostic)) << diagnostic;
    auto errors = stripe_errors(*run, *initial, kh, kv, time);
    if (options.attached_property) {
        // The largest drift of the last bin's property from its exact ratio.
        const auto& layout = run->state->layout();
        const int number = layout.populations().front().number_offset + 1;
        const int property = layout.property_offset(0) + 1;
        MultiFab drift(run->ba, run->dm, 1, 0);
        MultiFab::Copy(drift, run->spectrum(), property, 0, 1, 0);
        MultiFab::Saxpy(drift, -options.property_ratio, run->spectrum(), number, 0, 1, 0);
        errors.ratio_drift = drift.norm0(0) /
                             (options.property_ratio * run->spectrum().norm0(number));
    }
    amrex::Print() << "SBM_M6B_STRIPE " << label << " matched=" << options.match
                   << " steps=" << steps << " wet_peak=" << errors.wet_peak
                   << " dry_l1=" << errors.dry_l1 << " r_hat=" << errors.r_hat
                   << " property_ratio_drift=" << errors.ratio_drift << '\n';
    return errors;
}

// The stripe on the flat-walled wave terrain with A = Delta_x / Delta_Z = 10,
// cells_per_wave cells per vertical wavelength, a worst cell drop just under
// one thickness, and ConstantReference K_h = 1, K_v = r.
Options stripe_options (const int cells_per_wave, const Real r, const bool match)
{
    Options options;
    options.terrain = Terrain::Wave;
    options.nz = cells_per_wave;
    options.max_grid_z = cells_per_wave;
    options.nx = 16;
    options.ny = 2;
    const Real dx = Real(10.0) / Real(cells_per_wave);
    options.lx = Real(options.nx) * dx;
    options.ly = Real(options.ny) * dx;
    // The worst drop is near hill 2 pi (nz / nx) / (1 - pi hill) thickness.
    options.hill = Real(0.85) * Real(options.nx) / (two_pi * Real(cells_per_wave));
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(1.0);
    options.kv = r;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.profile = Profile::Stripe;
    options.match = match;
    return options;
}

StripeErrors run_stripe (const int cells_per_wave, const Real r, const bool match)
{
    const auto options = stripe_options(cells_per_wave, r, match);
    return run_height_only(options, options.kh, options.kv,
                           "cells_per_wave=" + std::to_string(cells_per_wave) +
                           " r=" + std::to_string(r));
}

// The design's accuracy gate on moderate terrain, at 64 cells per vertical
// wavelength: wet endpoint peak error below 2% of its amplitude, dry-cohort
// L1 error below 2% of its inventory, and |r_hat - r| <= max(0.002, r / 10),
// for the closures' ratios r = 0 (K_v = 0), 0.01 (anisotropic Smagorinsky at
// A = 10) and 1 (isotropic).
void run_stripe_gate ()
{
    for (const Real r : {Real(0.0), Real(0.01), Real(1.0)}) {
        const auto errors = run_stripe(64, r, true);
        EXPECT_LT(errors.wet_peak, Real(0.02)) << "r = " << r;
        EXPECT_LT(errors.dry_l1, Real(0.02)) << "r = " << r;
        EXPECT_LE(std::abs(errors.r_hat - r), std::max(Real(0.002), Real(0.1) * r)) << "r = " << r;
    }
}

// Max and mean errors after a stage-integrated run on the tilted slab z = s x + w
// between walls w = 0 and 1: q = 1 + cos(a x + b w) / 2 with b = s a / (r + s^2)
// has zero mapped transfer through every w surface, the walls included, and
// decays as exp(-a^2 K_h r t / (r + s^2)). Next to a wall the native side
// stencil reads the quadratic extension in height, so the wall rule sets the
// error there.
struct SlabErrors
{
    Real max{0};
    Real mean{0};
};

SlabErrors tilted_slab_errors (const int n)
{
    const Real s = Real(0.5);
    const Real r = Real(0.1);
    const Real a = two_pi;
    const Real b = s * a / (r + s * s);
    Options options;
    options.nx = n;
    options.ny = 2;
    options.ly = Real(2.0) / Real(n);
    options.nz = n;
    options.max_grid_z = n / 2;
    options.hill = Real(0.0);
    options.plane = s;
    options.policy = ParticleMixingPolicy::ConstantReference;
    options.kh = Real(1.0);
    options.kv = r;
    options.profile = Profile::Uniform;
    auto run = make_run(options);
    const int ncomp = run->spectrum().nComp();
    const Real nx = Real(n);
    const auto fill = [&] (MultiFab& u, const Real amplitude) {
        for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
            const auto out = u.array(mfi);
            const auto height = run->z_cc.const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                    const Real x = (Real(i) + Real(0.5)) / nx;
                    const Real w = height(i, j, k) - s * x;
                    out(i, j, k, c) = Real(0.01) * Real(c + 1) *
                        (Real(1.0) + amplitude * std::cos(a * x + b * w));
                });
        }
    };
    fill(run->state->new_state_for_initialization(0), Real(0.5));
    const Real decay = a * a * options.kh * r / (r + s * s);
    const Real time = Real(0.5) / decay;
    erf_auxiliary::MappedFaceTransmissibility t;
    t.define(run->ba, run->dm, 1, 0);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::build_particle_mixing_transmissibility(
        run->transport->particle_mixing(), run->geom, run->inputs.metrics, run->rho, Rho_comp,
        nullptr, t, diagnostic)) << diagnostic;
    erf_auxiliary::MappedFaceFluxRate at_rest;
    at_rest.define(run->ba, run->dm, 1, 0);
    at_rest.setVal(Real(0.0));
    Real rate = Real(0.0);
    EXPECT_TRUE(erf_auxiliary::ComputeMaxMappedCombinedOutgoingRate(
        at_rest, t, run->transport->static_measure(0), run->rho, Rho_comp,
        run->geom.InvCellSizeArray(), rate, diagnostic)) << diagnostic;
    const int steps = static_cast<int>(std::ceil(time * rate / Real(0.8)));
    advance(*run, HostIntegrator::CompressibleRK3, time / Real(steps), steps);
    MultiFab exact(run->ba, run->dm, ncomp, 0);
    fill(exact, Real(0.5) * std::exp(-decay * time));
    // Relative to each component's amplitude 0.005 (c + 1).
    MultiFab::Subtract(exact, run->spectrum(), 0, 0, ncomp, 0);
    SlabErrors errors;
    const Real cells = Real(exact.boxArray().numPts());
    for (int c = 0; c < ncomp; ++c) {
        const Real amplitude = Real(0.005) * Real(c + 1);
        errors.max = std::max(errors.max, exact.norm0(c) / amplitude);
        errors.mean = std::max(errors.mean, exact.norm1(c) / (amplitude * cells));
    }
    return errors;
}

// The wall rule converges on the tilted slab. The interior error halves twice
// per refinement; the wall cells, whose first interior k face carries an
// O(h^2) error that the suppressed wall face does not balance, approach that
// rate more slowly (max-norm orders 2.0, 1.73, 1.86 from 16 to 128 cells per
// direction; mean-error orders 1.91, 1.92, 1.96).
void run_tilted_slab ()
{
    const auto e16 = tilted_slab_errors(16);
    const auto e32 = tilted_slab_errors(32);
    const auto e64 = tilted_slab_errors(64);
    amrex::Print() << "SBM_M6B_TILTED_SLAB max=" << e16.max << "," << e32.max << "," << e64.max
                   << " mean=" << e16.mean << "," << e32.mean << "," << e64.mean
                   << " mean_order=" << std::log2(e16.mean / e32.mean) << ","
                   << std::log2(e32.mean / e64.mean) << '\n';
    EXPECT_GE(std::log2(e32.mean / e64.mean), Real(1.8));
    EXPECT_GT(e32.max / e64.max, Real(3.0));
}

// Weighted, mean-free amplitude sqrt(sum C (q - qbar_C)^2), C = omega rho, of
// component 0 of a state.
Real weighted_amplitude (const Run& run, const MultiFab& u)
{
    const MultiFab& omega = run.transport->static_measure(0);
    MultiFab c(run.ba, run.dm, 1, 0);
    MultiFab::Copy(c, omega, 0, 0, 1, 0);
    MultiFab::Multiply(c, run.rho, Rho_comp, 0, 1, 0);
    MultiFab q(run.ba, run.dm, 1, 0);
    MultiFab::Copy(q, u, 0, 0, 1, 0);
    MultiFab::Divide(q, run.rho, Rho_comp, 0, 1, 0);
    const Real mean = MultiFab::Dot(c, 0, q, 0, 1, 0) / c.sum(0);
    q.plus(-mean, 0, 1, 0);
    MultiFab cq(run.ba, run.dm, 1, 0);
    MultiFab::Copy(cq, q, 0, 0, 1, 0);
    MultiFab::Multiply(cq, c, 0, 0, 1, 0);
    return std::sqrt(MultiFab::Dot(cq, 0, q, 0, 1, 0));
}

// The least eigenvector of the symmetric part of the C-weighted native
// operator, restricted to cells three or more layers from the walls, as an
// intensive field (zero outside the restriction): the numerical analysis's
// witness for instantaneous variance growth. Host-side; the grid is
// y-invariant and the vector is taken y-invariant too.
std::vector<Real> damping_witness (Run& run, Real& least_eigenvalue)
{
    const int nx = run.options.nx;
    const int nz = run.nz();
    std::vector<IntVect> cells;
    for (int k = 3; k <= nz - 4; ++k) {
        for (int i = 0; i < nx; ++i) { cells.emplace_back(i, 0, k); }
    }
    const int m = static_cast<int>(cells.size());
    std::vector<Real> weight(static_cast<std::size_t>(m));
    const MultiFab& omega = run.transport->static_measure(0);
    auto host_value = [&] (const MultiFab& mf, const int comp, const IntVect& cell) {
        Real value = Real(0.0);
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            if (mfi.validbox().contains(cell)) { value = mf.const_array(mfi)(cell, comp); }
        }
        amrex::ParallelAllReduce::Sum(value, amrex::ParallelContext::CommunicatorSub());
        return value;
    };
    for (int a = 0; a < m; ++a) {
        weight[static_cast<std::size_t>(a)] =
            host_value(omega, 0, cells[static_cast<std::size_t>(a)]) *
            host_value(run.rho, Rho_comp, cells[static_cast<std::size_t>(a)]);
    }
    // A = C^{-1/2} K C^{-1/2}, where K_ab = -omega_a T(e_b)_a and T(e_b) is the
    // tendency of U = rho e_b, so that the weighted variance changes by -q^T K q.
    std::vector<Real> matrix(static_cast<std::size_t>(m * m), Real(0.0));
    MultiFab u(run.ba, run.dm, run.spectrum().nComp(), 0);
    MultiFab tendency(run.ba, run.dm, run.spectrum().nComp(), 0);
    for (int b = 0; b < m; ++b) {
        const IntVect column = cells[static_cast<std::size_t>(b)];
        u.setVal(Real(0.0));
        for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
            const auto out = u.array(mfi);
            const auto r = run.rho.const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (i == column[0] && k == column[2]) { out(i, j, k, 0) = r(i, j, k, Rho_comp); }
            });
        }
        tendency.setVal(Real(0.0));
        native_tendency(run, u, tendency);
        for (int a = 0; a < m; ++a) {
            const Real value = host_value(tendency, 0, cells[static_cast<std::size_t>(a)]);
            matrix[static_cast<std::size_t>(a * m + b)] =
                -host_value(omega, 0, cells[static_cast<std::size_t>(a)]) * value /
                std::sqrt(weight[static_cast<std::size_t>(a)] * weight[static_cast<std::size_t>(b)]);
        }
    }
    // Its symmetric part S = (A + A^T) / 2, and a bound on its spectrum.
    std::vector<Real> sym(static_cast<std::size_t>(m * m));
    Real bound = Real(0.0);
    for (int a = 0; a < m; ++a) {
        Real row = Real(0.0);
        for (int b = 0; b < m; ++b) {
            const Real value = Real(0.5) * (matrix[static_cast<std::size_t>(a * m + b)] +
                                            matrix[static_cast<std::size_t>(b * m + a)]);
            sym[static_cast<std::size_t>(a * m + b)] = value;
            row += std::abs(value);
        }
        bound = std::max(bound, row);
    }
    // Power iteration on bound I - S for the least eigenvector of S.
    std::vector<Real> v(static_cast<std::size_t>(m)), w(static_cast<std::size_t>(m));
    for (int a = 0; a < m; ++a) {
        v[static_cast<std::size_t>(a)] = Real(1.0) + Real(0.1) * std::sin(Real(1.7) * Real(a));
    }
    Real eigen = Real(0.0);
    for (int iteration = 0; iteration < 200000; ++iteration) {
        Real norm = Real(0.0);
        for (int a = 0; a < m; ++a) {
            Real value = bound * v[static_cast<std::size_t>(a)];
            for (int b = 0; b < m; ++b) {
                value -= sym[static_cast<std::size_t>(a * m + b)] * v[static_cast<std::size_t>(b)];
            }
            w[static_cast<std::size_t>(a)] = value;
            norm += value * value;
        }
        norm = std::sqrt(norm);
        Real change = Real(0.0);
        for (int a = 0; a < m; ++a) {
            const Real next = w[static_cast<std::size_t>(a)] / norm;
            change = std::max(change, std::abs(next - v[static_cast<std::size_t>(a)]));
            v[static_cast<std::size_t>(a)] = next;
        }
        eigen = bound - norm;
        if (change < Real(1.0e-13)) { break; }
    }
    least_eigenvalue = eigen;
    std::vector<Real> field(static_cast<std::size_t>(nx * nz), Real(0.0));
    for (int a = 0; a < m; ++a) {
        const IntVect cell = cells[static_cast<std::size_t>(a)];
        field[static_cast<std::size_t>(cell[2] * nx + cell[0])] =
            v[static_cast<std::size_t>(a)] / std::sqrt(weight[static_cast<std::size_t>(a)]);
    }
    return field;
}

// Mean-free weighted amplitudes of a matched-stage run from an initial field
// 2 + f / max|f| (positive), at t K_h / L^2 = 0.001, 0.01, 0.1, 1, over the
// initial amplitude.
std::array<Real, 4> damping_trajectory (Run& run, const std::vector<Real>& f)
{
    const int nx = run.options.nx;
    Real peak = Real(0.0);
    for (const Real value : f) { peak = std::max(peak, std::abs(value)); }
    amrex::Gpu::DeviceVector<Real> field(f.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, f.begin(), f.end(), field.begin());
    const Real* fp = field.data();
    auto& u0 = run.state->new_state_for_initialization(0);
    const int ncomp = u0.nComp();
    for (amrex::MFIter mfi(u0); mfi.isValid(); ++mfi) {
        const auto out = u0.array(mfi);
        const auto r = run.rho.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                out(i, j, k, n) = r(i, j, k, Rho_comp) * Real(0.01) *
                                  (Real(2.0) + fp[k * nx + i] / peak);
            });
    }
    const Real initial = weighted_amplitude(run, run.spectrum());
    erf_auxiliary::MappedFaceTransmissibility t;
    t.define(run.ba, run.dm, 1, 0);
    std::string diagnostic;
    EXPECT_TRUE(erf_sbm::build_particle_mixing_transmissibility(
        run.transport->particle_mixing(), run.geom, run.inputs.metrics, run.rho, Rho_comp,
        nullptr, t, diagnostic)) << diagnostic;
    erf_auxiliary::MappedFaceFluxRate at_rest;
    at_rest.define(run.ba, run.dm, 1, 0);
    at_rest.setVal(Real(0.0));
    Real rate = Real(0.0);
    EXPECT_TRUE(erf_auxiliary::ComputeMaxMappedCombinedOutgoingRate(
        at_rest, t, run.transport->static_measure(0), run.rho, Rho_comp,
        run.geom.InvCellSizeArray(), rate, diagnostic)) << diagnostic;
    std::array<Real, 4> amplitude{};
    double now = 0.0;
    const std::array<double, 4> marks{0.001, 0.01, 0.1, 1.0};
    for (std::size_t s = 0; s < marks.size(); ++s) {
        const double span = marks[s] / run.options.kh - now;
        const int steps = static_cast<int>(std::ceil(span * rate / 0.8));
        const double dt = span / steps;
        for (int step = 0; step < steps; ++step) {
            now = step_once(run, HostIntegrator::CompressibleRK3, now, dt);
        }
        amplitude[s] = weighted_amplitude(run, run.spectrum()) / initial;
    }
    return amplitude;
}

// The weak-damping gate: from the curved witness, at K_v / K_h = 0, 1e-5,
// 0.01 and 1, the matched RK3 stage keeps the mean-free weighted amplitude
// within 2% of its initial value, and an ordinary resolved mode decays.
void run_weak_damping (const Real hill, const Real stretch, const int n,
                       const std::vector<Real>& ratios)
{
    for (const Real r : ratios) {
        Options options;
        options.terrain = Terrain::Witness;
        options.hill = hill;
        options.plane = stretch;
        options.nx = n;
        options.ny = 2;
        options.ly = Real(2.0) / Real(n);
        options.nz = n;
        options.max_grid_xy = n;
        options.max_grid_z = n;
        options.policy = ParticleMixingPolicy::ConstantReference;
        options.kh = Real(1.0);
        options.kv = r;
        options.profile = Profile::Uniform;
        auto run = make_run(options);
        std::string reason;
        const bool admitted =
            erf_sbm::admit_terrain_mixing_geometry(run->z_nd, run->ba, Real(1.0), reason);
        Real least = Real(0.0);
        const auto witness = damping_witness(*run, least);
        const auto amplitude = damping_trajectory(*run, witness);
        // An ordinary mode, cos(2 pi xi) cos(pi zeta) over the same grid.
        std::vector<Real> mode(static_cast<std::size_t>(n * n));
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                mode[static_cast<std::size_t>(k * n + i)] =
                    std::cos(two_pi * (Real(i) + Real(0.5)) / Real(n)) *
                    std::cos(Real(0.5) * two_pi * (Real(k) + Real(0.5)) / Real(n));
            }
        }
        auto ordinary = make_run(options);
        const auto decay = damping_trajectory(*ordinary, mode);
        amrex::Print() << "SBM_M6B_WEAK_DAMPING n=" << n << " hill=" << hill
                       << " stretch=" << stretch
                       << " moderate=" << admitted << " r=" << r
                       << " least_eigenvalue=" << least << " witness="
                       << amplitude[0] << "," << amplitude[1] << "," << amplitude[2] << ","
                       << amplitude[3] << " ordinary=" << decay[0] << "," << decay[1] << ","
                       << decay[2] << "," << decay[3] << '\n';
        for (const Real a : amplitude) { EXPECT_LE(a, Real(1.02)) << "r = " << r; }
        EXPECT_LT(decay[3], decay[0]) << "r = " << r;
    }
}

// The stripe gate with an attached property exactly proportional to its
// carrier, 0.37 of the number, in every cell. Its donor envelope has no width,
// so roundoff in that support can keep the limiter fully active and the
// accepted diffusion the matched low-order exchange; the gates still hold, and
// the ratio stays exact.
void run_proportional_property_gate ()
{
    for (const Real r : {Real(0.0), Real(0.01), Real(1.0)}) {
        auto options = stripe_options(64, r, true);
        options.attached_property = true;
        options.property_ratio = Real(0.37);
        const auto errors = run_height_only(options, options.kh, options.kv,
                                            "proportional_property r=" + std::to_string(r));
        EXPECT_LT(errors.wet_peak, Real(0.02)) << "r = " << r;
        EXPECT_LT(errors.dry_l1, Real(0.02)) << "r = " << r;
        EXPECT_LE(std::abs(errors.r_hat - r), std::max(Real(0.002), Real(0.1) * r)) << "r = " << r;
        EXPECT_LE(errors.ratio_drift, Real(64.0) * eps()) << "r = " << r;
    }
}

// The stripe with K_v = 0 on curved terrain sloping in x and y, with both
// signs, under HostPassiveScalar with a horizontal diffusivity varying in x
// and y. The spectrum is stationary, and the stripe gates hold.
void run_two_slope_gate ()
{
    auto options = stripe_options(64, Real(0.0), true);
    options.ny = 16;
    options.ly = options.lx;
    options.max_grid_xy = 16;
    // Equal slopes in x and y keep the worst combined drop under one thickness.
    options.hill /= std::sqrt(Real(2.0));
    options.hill_y = options.hill;
    options.policy = ParticleMixingPolicy::HostPassiveScalar;
    options.kh = Real(0.0);
    options.kv = Real(0.0);
    options.horizontal_eddy_only = true;
    options.eddy_scale = Real(1.0);
    const auto errors = run_height_only(options, Real(1.0), Real(0.0), "two_slope_varying_kh");
    EXPECT_LT(errors.wet_peak, Real(0.02));
    EXPECT_LT(errors.dry_l1, Real(0.02));
    EXPECT_LE(std::abs(errors.r_hat), Real(0.002));
}

// A nearly empty interior 2M bin whose mean mass lies within 1e-10 to 2e-9
// (relative) of its strict upper edge, with an attached property, advected and
// mixed over the hill with varying density: every group stays admissible,
// the strict edge's reserve included, every inventory is conserved, and the
// bin's mean stays below its upper edge in every cell.
void run_edge_bin ()
{
    Options options;
    options.mode = erf_sbm::MomentMode::TwoMoment;
    options.attached_property = true;
    options.profile = Profile::EdgeBin;
    options.carrier_x = Real(1.0);
    options.varying_density = true;
    auto run = make_run(options);
    const auto initial = copy_of(run->spectrum());
    const auto& population = run->state->layout().populations().front();
    const int mass = population.mass_offset;
    const int number = population.number_offset;
    double t = 0.0;
    for (int step = 0; step < 4; ++step) {
        t = step_once(*run, HostIntegrator::CompressibleRK3, t, 0.03);
        std::string diagnostic;
        EXPECT_TRUE(erf_sbm::authoritative_state_admissible(run->spectrum(),
                                                            run->state->layout(), 0, &diagnostic))
            << "step " << step << ": " << diagnostic;
        // M - 0.5 N < 0: the mean below the upper edge 0.5.
        MultiFab excess(run->ba, run->dm, 1, 0);
        MultiFab::LinComb(excess, Real(1.0), run->spectrum(), mass, Real(-0.5),
                          run->spectrum(), number, 0, 1, 0);
        EXPECT_LT(excess.max(0), Real(0.0)) << "step " << step;
    }
    MultiFab change(run->ba, run->dm, 1, 0);
    MultiFab::Copy(change, run->spectrum(), number, 0, 1, 0);
    MultiFab::Subtract(change, *initial, number, 0, 1, 0);
    EXPECT_GT(change.norm0(0), Real(1.0e-3) * initial->norm0(number));
    expect_conserved(*initial, *run);
}

// Qualification tool, disabled by default: the cost of matching against the
// unmatched comparator (theta = 1) on the same terrain stage, boundaries,
// chunks and timestep, with an advecting carrier and HostPassiveScalar
// diffusivities, for 2M populations of 32 and 128 bins on 16 boxes and on one.
// Prints the median wall time of an RK3 step after a warm-up, the ratio, and
// the terrain scratch bytes. Run it on one and several ranks and on a GPU.
void run_matching_cost ()
{
    constexpr int timed_steps = 5;
    for (const int nbins : {32, 128}) {
        for (const int max_grid_xy : {4, 16}) {
            Real median[2] = {};
            std::size_t scratch_bytes = 0;
            int ncomp = 0;
            for (const bool match : {true, false}) {
                Options options;
                options.nbins = nbins;
                options.mode = erf_sbm::MomentMode::TwoMoment;
                options.nx = 16;
                options.ny = 16;
                options.nz = 32;
                options.max_grid_xy = max_grid_xy;
                options.max_grid_z = 32;
                options.carrier_x = Real(0.5);
                options.hill = Real(0.05);
                options.match = match;
                auto run = make_run(options);
                ncomp = run->spectrum().nComp();
                const double dt = 2.0e-4;
                double t = 0.0;
                t = step_once(*run, HostIntegrator::CompressibleRK3, t, dt);
                std::vector<double> seconds;
                for (int step = 0; step < timed_steps; ++step) {
                    amrex::Gpu::streamSynchronize();
                    amrex::ParallelDescriptor::Barrier();
                    const auto start = std::chrono::steady_clock::now();
                    t = step_once(*run, HostIntegrator::CompressibleRK3, t, dt);
                    amrex::Gpu::streamSynchronize();
                    amrex::ParallelDescriptor::Barrier();
                    seconds.push_back(std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - start).count());
                }
                std::sort(seconds.begin(), seconds.end());
                median[match ? 0 : 1] = Real(seconds[seconds.size() / 2]);
                // Chunk-width terrain scratch: the stencil state with two ghost
                // cells, the raw face fluxes (side faces with one z layer), and
                // the two lateral open masks with one z layer.
                const int chunk = run->transport->max_chunk_components();
                scratch_bytes = 0;
                for (int b = 0; b < run->ba.size(); ++b) {
                    const IntVect len = run->ba[b].length();
                    const std::size_t cells = std::size_t(len[0] + 4) * (len[1] + 4) * (len[2] + 4);
                    const std::size_t faces =
                        std::size_t(len[0] + 1) * len[1] * (len[2] + 2) +
                        std::size_t(len[0]) * (len[1] + 1) * (len[2] + 2) +
                        std::size_t(len[0]) * len[1] * (len[2] + 1);
                    const std::size_t masks = std::size_t(len[0] + 1) * len[1] * (len[2] + 2) +
                                              std::size_t(len[0]) * (len[1] + 1) * (len[2] + 2);
                    scratch_bytes += sizeof(Real) * ((cells + faces) * chunk + masks);
                }
            }
            amrex::Print() << "SBM_M6B_MATCHING_COST bins=" << nbins << " components=" << ncomp
                           << " boxes=" << (max_grid_xy == 4 ? 16 : 1)
                           << " ranks=" << amrex::ParallelDescriptor::NProcs()
                           << " matched_step_s=" << median[0] << " unmatched_step_s=" << median[1]
                           << " ratio=" << median[0] / median[1]
                           << " terrain_scratch_bytes=" << scratch_bytes << '\n';
        }
    }
}

} // namespace

TEST(SBMTerrainMixing, InactiveLimiterRecoversNativeTerrainDiffusion) { run_recovery(); }
TEST(SBMTerrainMixing, WallsCloseTheMappedTransferAndGroupsStayAdmissible)
{
    run_conservation();
}
TEST(SBMTerrainMixing, ZeroCoefficientIsBitwiseAdvection) { run_zero_coefficient(); }
TEST(SBMTerrainMixing, ConstantRatioSurvivesCarrierTransport) { run_constant_ratio(); }
TEST(SBMTerrainMixing, VerticallySplitBoxesAgree) { run_split_boxes(); }
TEST(SBMTerrainMixing, InflowFaceTakesTheFaceTraceBesideTheWalls) { run_inflow_face_trace(); }
TEST(SBMTerrainMixing, DiffusionThroughAnInflowFaceClosesTheLedger) { run_inflow_ledger(); }
TEST(SBMTerrainMixing, InflowOutflowFaceClosesWhereTheCarrierLeaves)
{
    run_inflow_outflow_masks_outward_faces();
}
TEST(SBMTerrainMixing, UnsafeStepFailsBeforeAnyStateChange) { run_unsafe_step(); }
TEST(SBMTerrainMixing, StageRejectsUnqualifiedTerrainMixing) { run_terrain_envelope(); }
TEST(SBMTerrainMixing, MatchingKeepsTheWetEndpointOfADryGroupStationary)
{
    run_dry_endpoint();
}
TEST(SBMTerrainMixing, ModerateTerrainStripeAccuracyGate) { run_stripe_gate(); }
TEST(SBMTerrainMixing, WallRuleConvergesOnATiltedSlab) { run_tilted_slab(); }
TEST(SBMTerrainMixing, ProportionalPropertyStripeAccuracyGate) { run_proportional_property_gate(); }
TEST(SBMTerrainMixing, TwoSlopeCurvedDryHalfAccuracyGate) { run_two_slope_gate(); }
TEST(SBMTerrainMixing, NearlyEmptyBinAtItsStrictEdgeStaysCanonical) { run_edge_bin(); }
TEST(SBMTerrainMixing, DISABLED_MatchingCost) { run_matching_cost(); }
TEST(SBMTerrainMixing, DISABLED_WeakDampingWitness)
{
    // The numerical analysis's grid, then half its amplitudes, at 16^2; then
    // the weakest damping on the analysis's grid refined to 32^2.
    const std::vector<Real> all{Real(0.0), Real(1.0e-5), Real(0.01), Real(1.0)};
    run_weak_damping(Real(0.75) / two_pi, Real(0.32) / two_pi, 16, all);
    run_weak_damping(Real(0.375) / two_pi, Real(0.16) / two_pi, 16, all);
    run_weak_damping(Real(0.75) / two_pi, Real(0.32) / two_pi, 32, {Real(0.0), Real(1.0e-5)});
}
TEST(SBMTerrainMixing, DISABLED_StripeQualification)
{
    for (const int cells : {64, 128}) {
        for (const Real r : {Real(0.0), Real(0.01), Real(1.0)}) {
            for (const bool match : {true, false}) { run_stripe(cells, r, match); }
        }
    }
}
