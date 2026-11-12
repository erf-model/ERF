#include <gtest/gtest.h>

#include <AMReX_BCRec.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Reduce.H>

#include "ERF_AuxiliaryDiffusion.H"
#include "ERF_AuxiliaryMappedTransport.H"
#include "ERF_IndexDefines.H"
#include "ERF_SBMCanonicalIdentity.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMLayout.H"
#include "ERF_SBMParticleMixing.H"
#include "ERF_ScalarDiffusion.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

using amrex::Array4;
using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::GpuArray;
using amrex::IntVect;
using amrex::MultiFab;
using amrex::Real;
using namespace erf_auxiliary;

constexpr Real two_pi = Real(6.283185307179586476925286766559);
constexpr Real stretch_amplitude = Real(0.5);

Real eps () { return std::numeric_limits<Real>::epsilon(); }

// Physical height of computational level zeta in [0, 1]; dzp/dzeta is
// stationary at both walls.
Real stretched_height (const Real zeta)
{
    return zeta - stretch_amplitude / two_pi * std::sin(two_pi * zeta);
}

struct GridOptions
{
    IntVect cells{8, 6, 5};
    int max_grid{4};
    std::array<int, 3> periodic{1, 1, 1};
    bool stretched{false};
    bool variable_maps{true};
};

struct DiffusionGrid
{
    Box domain;
    Geometry geom;
    BoxArray ba;
    BoxArray ba2d;
    DistributionMapping dm;
    bool stretched{false};
    std::vector<Real> dz_host;
    amrex::Gpu::DeviceVector<Real> dz;
    MultiFab mf_ux, mf_uy, mf_vx, mf_vy, mf_mx, mf_my;
    MultiFab ax, ay, detj, omega;

    [[nodiscard]] MappedDiffusionMetrics metrics () const
    {
        MappedDiffusionMetrics m;
        m.mesh = stretched ? DiffusionMesh::StretchedDz : DiffusionMesh::ConstantDz;
        m.mf_ux = &mf_ux;
        m.mf_uy = &mf_uy;
        m.mf_vx = &mf_vx;
        m.mf_vy = &mf_vy;
        m.mf_mx = &mf_mx;
        m.mf_my = &mf_my;
        m.ax = &ax;
        m.ay = &ay;
        m.stretched_dz = stretched ? dz.data() : nullptr;
        return m;
    }

    [[nodiscard]] GpuArray<Real, AMREX_SPACEDIM> dx_inv () const { return geom.InvCellSizeArray(); }
};

void fill_map_factor (MultiFab& mf, const Real base, const Real amplitude,
                      const int phase, const IntVect cells)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box bx = mfi.fabbox();
        const auto a = mf.array(mfi);
        const Real nx = Real(cells[0]);
        const Real ny = Real(cells[1]);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
            a(i, j, 0, 0) = base + amplitude * std::sin(two_pi * Real(i + phase) / nx) +
                            Real(0.5) * amplitude * std::cos(two_pi * Real(j) / ny);
        });
    }
}

std::unique_ptr<DiffusionGrid> make_grid (const GridOptions& options)
{
    auto g = std::make_unique<DiffusionGrid>();
    g->domain = Box(IntVect(0), options.cells - IntVect(1));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    g->geom = Geometry(g->domain, &real_box, amrex::CoordSys::cartesian,
                       options.periodic.data());
    g->ba = BoxArray(g->domain);
    g->ba.maxSize(options.max_grid);
    g->dm = DistributionMapping(g->ba);
    amrex::BoxList boxes = g->ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    g->ba2d = BoxArray(std::move(boxes));
    g->stretched = options.stretched;

    g->mf_ux.define(amrex::convert(g->ba2d, IntVect(1, 0, 0)), g->dm, 1, 0);
    g->mf_uy.define(amrex::convert(g->ba2d, IntVect(1, 0, 0)), g->dm, 1, 0);
    g->mf_vx.define(amrex::convert(g->ba2d, IntVect(0, 1, 0)), g->dm, 1, 0);
    g->mf_vy.define(amrex::convert(g->ba2d, IntVect(0, 1, 0)), g->dm, 1, 0);
    g->mf_mx.define(g->ba2d, g->dm, 1, 0);
    g->mf_my.define(g->ba2d, g->dm, 1, 0);
    const Real amplitude = options.variable_maps ? Real(0.05) : Real(0.0);
    fill_map_factor(g->mf_ux, Real(1.10), amplitude, 0, options.cells);
    fill_map_factor(g->mf_uy, Real(0.92), amplitude, 1, options.cells);
    fill_map_factor(g->mf_vx, Real(1.04), amplitude, 2, options.cells);
    fill_map_factor(g->mf_vy, Real(0.97), amplitude, 3, options.cells);
    fill_map_factor(g->mf_mx, Real(1.07), amplitude, 4, options.cells);
    fill_map_factor(g->mf_my, Real(0.95), amplitude, 5, options.cells);
    if (!options.variable_maps) {
        for (auto* mf : {&g->mf_ux, &g->mf_uy, &g->mf_vx, &g->mf_vy, &g->mf_mx, &g->mf_my}) {
            mf->setVal(Real(1.0));
        }
    }

    const int nz = options.cells[2];
    const Real dzeta = Real(1.0) / Real(nz);
    g->dz_host.resize(static_cast<std::size_t>(nz));
    for (int k = 0; k < nz; ++k) {
        g->dz_host[static_cast<std::size_t>(k)] = options.stretched
            ? stretched_height(Real(k + 1) * dzeta) - stretched_height(Real(k) * dzeta)
            : dzeta;
    }
    g->dz.resize(g->dz_host.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, g->dz_host.begin(), g->dz_host.end(),
                     g->dz.begin());

    g->detj.define(g->ba, g->dm, 1, 0);
    g->ax.define(amrex::convert(g->ba, IntVect(1, 0, 0)), g->dm, 1, 0);
    g->ay.define(amrex::convert(g->ba, IntVect(0, 1, 0)), g->dm, 1, 0);
    const Real* dzp = g->dz.data();
    for (auto* mf : {&g->detj, &g->ax, &g->ay}) {
        for (amrex::MFIter mfi(*mf); mfi.isValid(); ++mfi) {
            const auto a = mf->array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                a(i, j, k, 0) = dzp[k] / dzeta;
            });
        }
    }
    g->omega.define(g->ba, g->dm, 1, 0);
    std::string diagnostic;
    EXPECT_TRUE(BuildMappedCellMeasure(g->omega, g->detj, g->mf_mx, g->mf_my, diagnostic))
        << diagnostic;
    return g;
}

struct DiffusionFields
{
    MultiFab rho;
    MultiFab eddy;
    MultiFab z;
};

template <class F>
void fill_cells (MultiFab& mf, const int comp, const Geometry& geom, F value)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, comp) = value(i, j, k);
        });
    }
    mf.FillBoundary(comp, 1, geom.periodicity());
}

std::unique_ptr<DiffusionFields> make_fields (const DiffusionGrid& g, const int ncomp)
{
    auto f = std::make_unique<DiffusionFields>();
    f->rho.define(g.ba, g.dm, 1, 1);
    f->eddy.define(g.ba, g.dm, EddyDiff::NumDiffs, 1);
    f->z.define(g.ba, g.dm, ncomp, 1);
    const Real nx = Real(g.domain.length(0));
    const Real ny = Real(g.domain.length(1));
    fill_cells(f->rho, 0, g.geom, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        return Real(1.1) + Real(0.15) * std::sin(two_pi * Real(i) / nx) +
               Real(0.08) * std::cos(two_pi * Real(j) / ny) + Real(0.03) * Real(k);
    });
    for (int c = 0; c < EddyDiff::NumDiffs; ++c) {
        fill_cells(f->eddy, c, g.geom, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            return Real(0.4) + Real(0.05) * Real(c) +
                   Real(0.1) * std::sin(two_pi * Real(i + c) / nx) +
                   Real(0.07) * std::cos(two_pi * Real(j) / ny) + Real(0.02) * Real(k);
        });
    }
    for (int c = 0; c < ncomp; ++c) {
        fill_cells(f->z, c, g.geom, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            return Real(0.5) + Real(0.3) * std::sin(two_pi * Real(i) / nx + Real(c)) *
                       std::cos(two_pi * Real(j) / ny) +
                   Real(0.1) * Real(k) + Real(0.05) * Real(c);
        });
    }
    return f;
}

MappedDiffusionCoefficients host_scalar_coefficients ()
{
    MappedDiffusionCoefficients c;
    c.horizontal = {Real(0.0), EddyDiff::Scalar_h, EddyDiff::Scalar_v};
    c.vertical = c.horizontal;
    c.mode = {false, true};
    return c;
}

MappedDiffusionCoefficients kinematic_coefficients (const Real kh, const Real kv)
{
    MappedDiffusionCoefficients c;
    c.horizontal.molecular_coeff = kh;
    c.vertical.molecular_coeff = kv;
    c.mode = {true, false};
    return c;
}

std::unique_ptr<MappedFaceTransmissibility>
build_transmissibility (const DiffusionGrid& g, const DiffusionFields& f,
                        const MappedDiffusionCoefficients& c, const bool expect_ok = true)
{
    auto t = std::make_unique<MappedFaceTransmissibility>();
    t->define(g.ba, g.dm, 1, 0);
    std::string diagnostic;
    const bool ok = BuildMappedDiffusionTransmissibility(
        *t, g.geom, g.metrics(), c, f.rho, 0, &f.eddy, diagnostic);
    EXPECT_EQ(ok, expect_ok) << diagnostic;
    return t;
}

std::unique_ptr<MappedFaceFluxRate>
build_rate (const DiffusionGrid& g, const MappedFaceTransmissibility& t,
            const MultiFab& z, const std::vector<int>& components)
{
    auto rate = std::make_unique<MappedFaceFluxRate>();
    rate->define(g.ba, g.dm, static_cast<int>(components.size()), 0);
    std::string diagnostic;
    EXPECT_TRUE(BuildLowOrderDiffusionRate(*rate, t, g.geom, z, components, diagnostic))
        << diagnostic;
    return rate;
}

// U' = U + tau dU/dt at fixed carrier for every rate slot n -> component n.
void explicit_update (const DiffusionGrid& g, const MappedFaceFluxRate& rate,
                      MultiFab& u, const Real tau)
{
    MultiFab tendency(g.ba, g.dm, rate.nComp(), 0);
    for (int n = 0; n < rate.nComp(); ++n) {
        ApplyMappedFluxTendency(rate, n, g.omega, 0, tendency, n, g.dx_inv());
    }
    MultiFab::Saxpy(u, tau, tendency, 0, 0, rate.nComp(), 0);
}

// Extensive state U = rho z of the selected components, with no ghost cells.
std::unique_ptr<MultiFab> extensive_state (const DiffusionGrid& g, const DiffusionFields& f)
{
    auto u = std::make_unique<MultiFab>(g.ba, g.dm, f.z.nComp(), 0);
    MultiFab::Copy(*u, f.z, 0, 0, f.z.nComp(), 0);
    for (int n = 0; n < f.z.nComp(); ++n) {
        MultiFab::Multiply(*u, f.rho, 0, n, 1, 0);
    }
    return u;
}

Real value_at (const MultiFab& mf, const int comp, const IntVect& iv)
{
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const Box point(iv, iv, mf.ixType());
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox() & point;
        if (!bx.ok()) { continue; }
        const auto a = mf.const_array(mfi);
        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> ReduceTuple {
                return {a(i, j, k, comp), Real(1.0)};
            });
    }
    auto result = reduce_data.value(reduce_op);
    Real sum = amrex::get<0>(result);
    Real count = amrex::get<1>(result);
    amrex::ParallelAllReduce::Sum(sum, amrex::ParallelContext::CommunicatorSub());
    amrex::ParallelAllReduce::Sum(count, amrex::ParallelContext::CommunicatorSub());
    return count > Real(0.0) ? sum / count : std::numeric_limits<Real>::quiet_NaN();
}

// Largest |a - b| over the faces that are not on a nonperiodic domain face.
template <typename TagA, typename TagB>
Real max_interior_face_difference (const MappedFaceField<TagA>& a, const int acomp,
                                   const MappedFaceField<TagB>& b, const int bcomp,
                                   const Geometry& geom)
{
    Real maximum = Real(0.0);
    const Box& domain = geom.Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto& fa = a.dir(dir);
        MultiFab difference(fa.boxArray(), fa.DistributionMap(), 1, 0);
        const int lo = domain.smallEnd(dir);
        const int hi = domain.bigEnd(dir);
        const bool periodic = geom.isPeriodic(dir);
        for (amrex::MFIter mfi(difference); mfi.isValid(); ++mfi) {
            const auto av = fa.const_array(mfi);
            const auto bv = b.dir(dir).const_array(mfi);
            const auto out = difference.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const IntVect cell(AMREX_D_DECL(i, j, k));
                const bool physical = !periodic && (cell[dir] == lo || cell[dir] == hi + 1);
                out(i, j, k, 0) = physical ? Real(0.0)
                    : amrex::Math::abs(av(i, j, k, acomp) - bv(i, j, k, bcomp));
            });
        }
        maximum = amrex::max(maximum, difference.norm0(0));
    }
    return maximum;
}

template <typename Tag>
Real max_abs_face (const MappedFaceField<Tag>& a, const int comp)
{
    Real maximum = Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        maximum = amrex::max(maximum, a.dir(dir).norm0(comp));
    }
    return maximum;
}

// The host N/S candidate mapped transfer of one intensive component.
std::unique_ptr<MappedFaceFluxRate>
host_mapped_diffusion (const DiffusionGrid& g, const DiffusionFields& f, const int comp,
                       const MappedDiffusionCoefficients& c)
{
    AMREX_ALWAYS_ASSERT(c.horizontal.molecular_coeff == c.vertical.molecular_coeff);
    ScalarDiffusionFluxPolicy policy;
    policy.coefficients = {c.horizontal.molecular_coeff, c.horizontal.eddy_h_comp,
                           c.vertical.eddy_v_comp};
    policy.coefficient_mode = c.mode;
    amrex::BCRec bc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int type = g.geom.isPeriodic(d) ? ERFBCType::int_dir : ERFBCType::foextrap;
        bc.setLo(d, type);
        bc.setHi(d, type);
    }
    amrex::Gpu::DeviceVector<amrex::BCRec> bcs(1);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, &bc, &bc + 1, bcs.begin());

    MappedFaceFluxRate raw;
    raw.define(g.ba, g.dm, 1, 0);
    auto mapped = std::make_unique<MappedFaceFluxRate>();
    mapped->define(g.ba, g.dm, 1, 0);
    const auto inv = g.dx_inv();
    const int klo = g.domain.smallEnd(2);
    const int khi = g.domain.bigEnd(2);
    for (amrex::MFIter mfi(f.rho); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        ScalarDiffusionFieldViews field;
        field.scalar = f.z.const_array(mfi);
        field.scalar_comp = comp;
        field.density = f.rho.const_array(mfi);
        field.rho_comp = 0;
        field.mu_turb = f.eddy.const_array(mfi);
        field.xflux = raw.dir(0).array(mfi);
        field.yflux = raw.dir(1).array(mfi);
        field.zflux = raw.dir(2).array(mfi);
        field.flux_comp = 0;
        if (g.stretched) {
            BuildScalarDiffusionFluxes_S(
                bx, g.domain, field, policy, inv[0], inv[1], g.dz.data(), klo, khi,
                g.mf_ux.const_array(mfi), g.mf_uy.const_array(mfi),
                g.mf_vx.const_array(mfi), g.mf_vy.const_array(mfi), bcs.data(), 0);
            BuildScalarDiffusionMappedTransfers_S(
                bx, raw.dir(0).const_array(mfi), 0, raw.dir(1).const_array(mfi), 0,
                raw.dir(2).const_array(mfi), 0, mapped->dir(0).array(mfi), 0,
                mapped->dir(1).array(mfi), 0, mapped->dir(2).array(mfi), 0,
                g.ax.const_array(mfi), g.ay.const_array(mfi),
                g.mf_mx.const_array(mfi), g.mf_my.const_array(mfi));
        } else {
            BuildScalarDiffusionFluxes_N(
                bx, g.domain, field, policy, inv, Array4<const Real>{}, Array4<const Real>{},
                g.mf_ux.const_array(mfi), g.mf_uy.const_array(mfi),
                g.mf_vx.const_array(mfi), g.mf_vy.const_array(mfi), bcs.data(), 0);
            BuildScalarDiffusionMappedTransfers_N(
                bx, raw.dir(0).const_array(mfi), 0, raw.dir(1).const_array(mfi), 0,
                raw.dir(2).const_array(mfi), 0, mapped->dir(0).array(mfi), 0,
                mapped->dir(1).array(mfi), 0, mapped->dir(2).array(mfi), 0,
                g.mf_mx.const_array(mfi), g.mf_my.const_array(mfi));
        }
    }
    amrex::Gpu::streamSynchronize();
    return mapped;
}

void expect_host_parity (const GridOptions& options)
{
    const auto g = make_grid(options);
    const auto f = make_fields(*g, 2);
    for (const auto& c : {host_scalar_coefficients(), kinematic_coefficients(Real(0.37), Real(0.37))}) {
        const auto t = build_transmissibility(*g, *f, c);
        const auto rate = build_rate(*g, *t, f->z, {1});
        const auto host = host_mapped_diffusion(*g, *f, 1, c);
        const Real scale = max_abs_face(*host, 0);
        ASSERT_GT(scale, Real(0.0));
        EXPECT_LE(max_interior_face_difference(*rate, 0, *host, 0, g->geom),
                  Real(64.0) * eps() * scale);
    }
}

void run_ngrid_parity ()
{
    GridOptions options;
    options.periodic = {1, 1, 1};
    expect_host_parity(options);
}

void run_stretched_parity ()
{
    GridOptions options;
    options.periodic = {1, 1, 0};
    options.stretched = true;
    expect_host_parity(options);
}

void run_physical_faces ()
{
    GridOptions options;
    options.cells = IntVect(6, 5, 4);
    options.max_grid = 3;
    options.periodic = {0, 1, 0};
    const auto g = make_grid(options);
    const auto f = make_fields(*g, 1);
    const auto t = build_transmissibility(*g, *f, host_scalar_coefficients());

    for (const int i : {0, 6}) {
        for (int j = 0; j < 5; ++j) {
            for (int k = 0; k < 4; ++k) {
                EXPECT_EQ(value_at(t->dir(0), 0, IntVect(i, j, k)), Real(0.0));
            }
        }
    }
    for (const int k : {0, 4}) {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 5; ++j) {
                EXPECT_EQ(value_at(t->dir(2), 0, IntVect(i, j, k)), Real(0.0));
            }
        }
    }
    EXPECT_GT(value_at(t->dir(0), 0, IntVect(3, 2, 1)), Real(0.0));
    EXPECT_GT(value_at(t->dir(1), 0, IntVect(3, 0, 1)), Real(0.0));
    EXPECT_GT(value_at(t->dir(2), 0, IntVect(3, 2, 2)), Real(0.0));

    // Ghosts beyond a physical face are never read.
    const Box& domain = g->domain;
    for (amrex::MFIter mfi(f->z); mfi.isValid(); ++mfi) {
        const auto z = f->z.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (i < domain.smallEnd(0) || i > domain.bigEnd(0) ||
                k < domain.smallEnd(2) || k > domain.bigEnd(2)) {
                z(i, j, k, 0) = std::numeric_limits<Real>::quiet_NaN();
            }
        });
    }
    const auto rate = build_rate(*g, *t, f->z, {0});
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        EXPECT_TRUE(rate->dir(dir).is_finite(0, 1, 0));
    }

    // A closed (walls plus periodic) domain conserves the extensive inventory.
    MultiFab tendency(g->ba, g->dm, 1, 0);
    ApplyMappedFluxTendency(*rate, 0, g->omega, 0, tendency, 0, g->dx_inv());
    MultiFab weighted(g->ba, g->dm, 1, 0);
    MultiFab::Copy(weighted, tendency, 0, 0, 1, 0);
    MultiFab::Multiply(weighted, g->omega, 0, 0, 1, 0);
    EXPECT_LE(std::abs(weighted.sum(0)), Real(64.0) * eps() * weighted.norm1(0));

    // A nonfinite ghost read across the periodic y face fails, on every rank.
    for (amrex::MFIter mfi(f->z); mfi.isValid(); ++mfi) {
        const auto z = f->z.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (i == 2 && j == -1 && k == 1) {
                z(i, j, k, 0) = std::numeric_limits<Real>::quiet_NaN();
            }
        });
    }
    MappedFaceFluxRate poisoned;
    poisoned.define(g->ba, g->dm, 1, 0);
    std::string diagnostic;
    EXPECT_FALSE(BuildLowOrderDiffusionRate(poisoned, *t, g->geom, f->z, {0}, diagnostic));
    EXPECT_NE(diagnostic.find("nonfinite"), std::string::npos) << diagnostic;
}

void run_zero_coefficient_identity ()
{
    const auto g = make_grid(GridOptions{});
    const auto f = make_fields(*g, 2);
    f->eddy.setVal(Real(0.0), EddyDiff::Scalar_h, 1, 1);
    f->eddy.setVal(Real(0.0), EddyDiff::Scalar_v, 1, 1);
    const auto t = build_transmissibility(*g, *f, host_scalar_coefficients());
    EXPECT_EQ(max_abs_face(*t, 0), Real(0.0));

    const auto rate = build_rate(*g, *t, f->z, {0, 1});
    EXPECT_EQ(max_abs_face(*rate, 0), Real(0.0));
    EXPECT_EQ(max_abs_face(*rate, 1), Real(0.0));
    const auto u = extensive_state(*g, *f);
    MultiFab before(g->ba, g->dm, u->nComp(), 0);
    MultiFab::Copy(before, *u, 0, 0, u->nComp(), 0);
    explicit_update(*g, *rate, *u, Real(0.7));
    MultiFab::Subtract(before, *u, 0, 0, u->nComp(), 0);
    EXPECT_EQ(before.norm0(0), Real(0.0));
    EXPECT_EQ(before.norm0(1), Real(0.0));

    MappedFaceFluxRate carrier;
    carrier.define(g->ba, g->dm, 1, 0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (amrex::MFIter mfi(carrier.dir(dir)); mfi.isValid(); ++mfi) {
            const auto a = carrier.dir(dir).array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                a(i, j, k, 0) = Real(0.3) * std::sin(Real(0.7) * Real(i + 2 * j) + Real(dir)) +
                                Real(0.05) * Real(k);
            });
        }
    }
    Real advective = Real(-1.0);
    Real combined = Real(-2.0);
    std::string diagnostic;
    ASSERT_TRUE(ComputeMaxMappedOutgoingRate(carrier, g->omega, f->rho, 0, g->dx_inv(),
                                             advective, diagnostic)) << diagnostic;
    ASSERT_TRUE(ComputeMaxMappedCombinedOutgoingRate(carrier, *t, g->omega, f->rho, 0,
                                                     g->dx_inv(), combined, diagnostic))
        << diagnostic;
    EXPECT_GT(advective, Real(0.0));
    EXPECT_EQ(combined, advective);
}

void run_conservation (const bool stretched)
{
    GridOptions options;
    options.periodic = stretched ? std::array<int, 3>{1, 1, 0} : std::array<int, 3>{1, 1, 1};
    options.stretched = stretched;
    const auto g = make_grid(options);
    const auto f = make_fields(*g, 3);
    const auto t = build_transmissibility(*g, *f, host_scalar_coefficients());
    const auto rate = build_rate(*g, *t, f->z, {0, 1, 2});
    MultiFab tendency(g->ba, g->dm, 3, 0);
    for (int n = 0; n < 3; ++n) {
        ApplyMappedFluxTendency(*rate, n, g->omega, 0, tendency, n, g->dx_inv());
        MultiFab::Multiply(tendency, g->omega, 0, n, 1, 0);
        EXPECT_GT(tendency.norm1(n), Real(0.0));
        EXPECT_LE(std::abs(tendency.sum(n)), Real(64.0) * eps() * tendency.norm1(n));
    }
}

// Two cells exchanging through one interior face, closed by walls in x.
struct PairCase
{
    std::unique_ptr<DiffusionGrid> grid;
    std::unique_ptr<DiffusionFields> fields;
    std::unique_ptr<MappedFaceTransmissibility> t;
};

PairCase make_pair_case (const Real rho0, const Real rho1, const Real z0, const Real z1)
{
    GridOptions options;
    options.cells = IntVect(2, 1, 1);
    options.max_grid = 2;
    options.periodic = {0, 1, 1};
    options.variable_maps = false;
    PairCase pair;
    pair.grid = make_grid(options);
    pair.fields = make_fields(*pair.grid, 1);
    fill_cells(pair.fields->rho, 0, pair.grid->geom,
               [=] AMREX_GPU_DEVICE(int i, int, int) noexcept { return i <= 0 ? rho0 : rho1; });
    fill_cells(pair.fields->z, 0, pair.grid->geom,
               [=] AMREX_GPU_DEVICE(int i, int, int) noexcept { return i <= 0 ? z0 : z1; });
    pair.fields->eddy.setVal(Real(0.8));
    pair.t = build_transmissibility(*pair.grid, *pair.fields, host_scalar_coefficients());
    return pair;
}

void run_positive_linear_combination ()
{
    constexpr Real rho0 = Real(1.0);
    constexpr Real rho1 = Real(2.5);
    constexpr Real z0 = Real(0.3);
    constexpr Real z1 = Real(1.7);
    auto pair = make_pair_case(rho0, rho1, z0, z1);
    const auto& g = *pair.grid;
    const Real dx_inv = g.dx_inv()[0];
    const Real face_t = value_at(pair.t->dir(0), 0, IntVect(1, 0, 0));
    ASSERT_GT(face_t, Real(0.0));
    // omega = 1, so C_i = rho_i per unit computational volume.
    const Real tau = Real(0.2);
    const Real b = tau * dx_inv * face_t;
    const auto rate = build_rate(g, *pair.t, pair.fields->z, {0});
    auto u = extensive_state(g, *pair.fields);
    explicit_update(g, *rate, *u, tau);
    const Real h0 = rho0 * z0;
    const Real h1 = rho1 * z1;
    const Real expected0 = (Real(1.0) - b / rho0) * h0 + (b / rho1) * h1;
    const Real expected1 = (b / rho0) * h0 + (Real(1.0) - b / rho1) * h1;
    EXPECT_NEAR(value_at(*u, 0, IntVect(0, 0, 0)), expected0, Real(16.0) * eps() * h1);
    EXPECT_NEAR(value_at(*u, 0, IntVect(1, 0, 0)), expected1, Real(16.0) * eps() * h1);

    // At fixed carrier the intensive update is convex: equal z stays equal.
    auto uniform = make_pair_case(rho0, rho1, Real(0.9), Real(0.9));
    const auto uniform_rate = build_rate(*uniform.grid, *uniform.t, uniform.fields->z, {0});
    auto uu = extensive_state(*uniform.grid, *uniform.fields);
    explicit_update(*uniform.grid, *uniform_rate, *uu, tau);
    EXPECT_NEAR(value_at(*uu, 0, IntVect(0, 0, 0)) / rho0, Real(0.9), Real(8.0) * eps());
    EXPECT_NEAR(value_at(*uu, 0, IntVect(1, 0, 0)) / rho1, Real(0.9), Real(8.0) * eps());

    // The amount update is not convex: equal H with unequal carrier changes.
    const Real h = Real(1.25);
    auto equal_h = make_pair_case(rho0, rho1, h / rho0, h / rho1);
    const auto equal_rate = build_rate(*equal_h.grid, *equal_h.t, equal_h.fields->z, {0});
    auto uh = extensive_state(*equal_h.grid, *equal_h.fields);
    explicit_update(*equal_h.grid, *equal_rate, *uh, tau);
    const Real expected_shift = b * h * (Real(1.0) / rho1 - Real(1.0) / rho0);
    EXPECT_NEAR(value_at(*uh, 0, IntVect(0, 0, 0)) - h, expected_shift, Real(16.0) * eps());
    EXPECT_GT(std::abs(expected_shift), Real(1.0e-3));

    // The row coefficient 1 - b/C_0 is nonnegative exactly up to b = C_0.
    auto pulse = make_pair_case(rho0, rho1, Real(1.0), Real(0.0));
    const auto pulse_rate = build_rate(*pulse.grid, *pulse.t, pulse.fields->z, {0});
    const Real tau_edge = rho0 / (dx_inv * face_t);
    auto at_edge = extensive_state(*pulse.grid, *pulse.fields);
    explicit_update(*pulse.grid, *pulse_rate, *at_edge, tau_edge);
    EXPECT_GE(value_at(*at_edge, 0, IntVect(0, 0, 0)), -Real(16.0) * eps());
    auto beyond = extensive_state(*pulse.grid, *pulse.fields);
    explicit_update(*pulse.grid, *pulse_rate, *beyond, Real(1.01) * tau_edge);
    EXPECT_LT(value_at(*beyond, 0, IntVect(0, 0, 0)), -Real(1.0e-3));
}

// A uniform, cubic-cell periodic grid with one nonzero cell.
struct UniformCase
{
    std::unique_ptr<DiffusionGrid> grid;
    std::unique_ptr<DiffusionFields> fields;
};

UniformCase make_uniform_pulse (const Real rho, const Real kh, const Real kv, const IntVect pulse)
{
    GridOptions options;
    options.cells = IntVect(4, 4, 4);
    options.max_grid = 2;
    options.periodic = {1, 1, 1};
    options.variable_maps = false;
    UniformCase uniform;
    uniform.grid = make_grid(options);
    uniform.fields = make_fields(*uniform.grid, 1);
    uniform.fields->rho.setVal(rho);
    uniform.fields->eddy.setVal(kh, EddyDiff::Scalar_h, 1, 1);
    uniform.fields->eddy.setVal(kv, EddyDiff::Scalar_v, 1, 1);
    fill_cells(uniform.fields->z, 0, uniform.grid->geom,
               [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   return (i == pulse[0] && j == pulse[1] && k == pulse[2]) ? Real(1.0)
                                                                          : Real(0.0);
               });
    return uniform;
}

Real combined_rate (const DiffusionGrid& g, const DiffusionFields& f,
                    const MappedFaceTransmissibility& t, const MappedFaceFluxRate& carrier)
{
    Real rate = Real(0.0);
    std::string diagnostic;
    EXPECT_TRUE(ComputeMaxMappedCombinedOutgoingRate(carrier, t, g.omega, f.rho, 0,
                                                     g.dx_inv(), rate, diagnostic))
        << diagnostic;
    return rate;
}

void run_sparse_pulse ()
{
    const IntVect pulse(2, 1, 3);
    auto uniform = make_uniform_pulse(Real(1.3), Real(0.6), Real(0.25), pulse);
    const auto& g = *uniform.grid;
    const auto& f = *uniform.fields;
    const auto t = build_transmissibility(g, f, host_scalar_coefficients());
    MappedFaceFluxRate no_carrier;
    no_carrier.define(g.ba, g.dm, 1, 0);
    no_carrier.setVal(Real(0.0));
    const Real bound = combined_rate(g, f, *t, no_carrier);
    ASSERT_GT(bound, Real(0.0));
    const auto rate = build_rate(g, *t, f.z, {0});

    auto at_bound = extensive_state(g, f);
    explicit_update(g, *rate, *at_bound, Real(1.0) / bound);
    EXPECT_GE(at_bound->min(0), -Real(64.0) * eps());

    auto beyond = extensive_state(g, f);
    explicit_update(g, *rate, *beyond, Real(1.01) / bound);
    EXPECT_LT(value_at(*beyond, 0, pulse), Real(0.0));
}

// Uniform +x carrier and uniform diffusion, each using `fraction` of one cell per unit time.
void run_combined_bound_case (const Real fraction, const bool expect_violation)
{
    const IntVect pulse(1, 2, 2);
    constexpr Real rho = Real(1.0);
    constexpr Real dx = Real(0.25);
    const Real k = fraction * rho * dx * dx / Real(6.0);
    auto uniform = make_uniform_pulse(rho, k, k, pulse);
    const auto& g = *uniform.grid;
    const auto& f = *uniform.fields;
    const auto t = build_transmissibility(g, f, host_scalar_coefficients());

    const Real carrier_value = fraction * rho * dx;
    MappedFaceFluxRate carrier;
    carrier.define(g.ba, g.dm, 1, 0);
    carrier.setVal(Real(0.0));
    carrier.dir(0).setVal(carrier_value);
    MappedFaceFluxRate no_carrier;
    no_carrier.define(g.ba, g.dm, 1, 0);
    no_carrier.setVal(Real(0.0));

    constexpr double tau = 1.0;
    Real advective = Real(0.0);
    std::string diagnostic;
    ASSERT_TRUE(ComputeMaxMappedOutgoingRate(carrier, g.omega, f.rho, 0, g.dx_inv(),
                                             advective, diagnostic)) << diagnostic;
    const Real diffusive = combined_rate(g, f, *t, no_carrier);
    const Real combined = combined_rate(g, f, *t, carrier);
    EXPECT_NEAR(advective, fraction, Real(64.0) * eps());
    EXPECT_NEAR(diffusive, fraction, Real(64.0) * eps());
    EXPECT_NEAR(combined, Real(2.0) * fraction, Real(64.0) * eps());
    double hard_limit = 0.0;
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(tau, advective, hard_limit));
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(tau, diffusive, hard_limit));
    EXPECT_EQ(FixedDtExceedsMappedDonorLimit(tau, combined, hard_limit), expect_violation);

    // Low-order advection (upwind donor) plus low-order diffusion from one anchor.
    auto rate = build_rate(g, *t, f.z, {0});
    for (amrex::MFIter mfi(rate->dir(0)); mfi.isValid(); ++mfi) {
        const auto out = rate->dir(0).array(mfi);
        const auto z = f.z.const_array(mfi);
        const Real carrier_rate = carrier_value;
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int kk) noexcept {
            out(i, j, kk, 0) += carrier_rate * z(i - 1, j, kk, 0);
        });
    }
    auto u = extensive_state(g, f);
    explicit_update(g, *rate, *u, Real(tau));
    const Real remaining = value_at(*u, 0, pulse) / rho;
    EXPECT_NEAR(remaining, Real(1.0) - Real(2.0) * fraction, Real(256.0) * eps());
    if (expect_violation) {
        EXPECT_LT(remaining, Real(0.0));
    } else {
        EXPECT_GE(u->min(0), Real(0.0));
    }
}

void run_dissipation ()
{
    const auto g = make_grid(GridOptions{});
    const auto f = make_fields(*g, 1);
    const auto t = build_transmissibility(*g, *f, host_scalar_coefficients());
    const auto rate = build_rate(*g, *t, f->z, {0});
    MultiFab tendency(g->ba, g->dm, 1, 0);
    ApplyMappedFluxTendency(*rate, 0, g->omega, 0, tendency, 0, g->dx_inv());

    // sum_i omega_i z_i dU_i/dt = -sum_f dx_inv_f T_f (z_R - z_L)^2.
    MultiFab lhs(g->ba, g->dm, 1, 0);
    MultiFab::Copy(lhs, tendency, 0, 0, 1, 0);
    MultiFab::Multiply(lhs, g->omega, 0, 0, 1, 0);
    MultiFab::Multiply(lhs, f->z, 0, 0, 1, 0);
    Real rhs = Real(0.0);
    const auto inv = g->dx_inv();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        MultiFab face_term(t->dir(dir).boxArray(), g->dm, 1, 0);
        const IntVect shift = IntVect::TheDimensionVector(dir);
        const Real scale = inv[dir];
        for (amrex::MFIter mfi(face_term); mfi.isValid(); ++mfi) {
            const auto out = face_term.array(mfi);
            const auto tf = t->dir(dir).const_array(mfi);
            const auto z = f->z.const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const IntVect cell(AMREX_D_DECL(i, j, k));
                const Real jump = z(cell, 0) - z(cell - shift, 0);
                out(i, j, k, 0) = scale * tf(i, j, k, 0) * jump * jump;
            });
        }
        rhs -= face_term.sum_unique(0, false, g->geom.periodicity());
    }
    EXPECT_LT(rhs, Real(0.0));
    EXPECT_NEAR(lhs.sum(0), rhs, Real(256.0) * eps() * std::abs(rhs));

    // Closed-domain energy decays monotonically within the positive bound.
    MappedFaceFluxRate no_carrier;
    no_carrier.define(g->ba, g->dm, 1, 0);
    no_carrier.setVal(Real(0.0));
    const Real tau = Real(0.5) / combined_rate(*g, *f, *t, no_carrier);
    auto energy = [&] () {
        MultiFab e(g->ba, g->dm, 1, 0);
        MultiFab::Copy(e, f->z, 0, 0, 1, 0);
        MultiFab::Multiply(e, f->z, 0, 0, 1, 0);
        MultiFab::Multiply(e, f->rho, 0, 0, 1, 0);
        MultiFab::Multiply(e, g->omega, 0, 0, 1, 0);
        return Real(0.5) * e.sum(0);
    };
    const Real initial = energy();
    Real previous = initial;
    for (int step = 0; step < 20; ++step) {
        const auto step_rate = build_rate(*g, *t, f->z, {0});
        auto u = extensive_state(*g, *f);
        explicit_update(*g, *step_rate, *u, tau);
        MultiFab::Copy(f->z, *u, 0, 0, 1, 0);
        MultiFab::Divide(f->z, f->rho, 0, 0, 1, 0);
        f->z.FillBoundary(g->geom.periodicity());
        const Real current = energy();
        EXPECT_LE(current, previous * (Real(1.0) + Real(16.0) * eps()));
        previous = current;
    }
    EXPECT_LT(previous, initial);
}

// Max error of the low-order divergence against the analytic div(K grad z).
Real manufactured_error (const int n, const bool stretched)
{
    GridOptions options;
    options.cells = IntVect(n, n, n);
    options.max_grid = 16;
    options.periodic = stretched ? std::array<int, 3>{1, 1, 0} : std::array<int, 3>{1, 1, 1};
    options.stretched = stretched;
    options.variable_maps = false;
    const auto g = make_grid(options);
    const auto f = make_fields(*g, 1);
    const Real h = Real(1.0) / Real(n);
    const Real* dzp = g->dz.data();
    const bool walls = stretched;
    auto height_at_center = [=] AMREX_GPU_DEVICE(int k) noexcept {
        if (!walls) { return (Real(k) + Real(0.5)) * h; }
        Real face = Real(0.0);
        for (int kk = 0; kk < k; ++kk) { face += dzp[kk]; }
        return face + Real(0.5) * dzp[k];
    };
    fill_cells(f->eddy, EddyDiff::Scalar_h, g->geom, [=] AMREX_GPU_DEVICE(int, int j, int) noexcept {
        return Real(0.5) + Real(0.2) * std::cos(two_pi * (Real(j) + Real(0.5)) * h);
    });
    fill_cells(f->eddy, EddyDiff::Scalar_v, g->geom, [=] AMREX_GPU_DEVICE(int, int j, int) noexcept {
        return Real(0.5) + Real(0.2) * std::cos(two_pi * (Real(j) + Real(0.5)) * h);
    });
    const Real vertical_wave = walls ? Real(0.5) * two_pi : two_pi;
    fill_cells(f->z, 0, g->geom, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real x = (Real(i) + Real(0.5)) * h;
        const Real y = (Real(j) + Real(0.5)) * h;
        const int kk = walls ? amrex::max(0, amrex::min(k, n - 1)) : k;
        return std::sin(two_pi * x) * std::cos(two_pi * y) +
               std::cos(vertical_wave * height_at_center(kk));
    });
    const auto t = build_transmissibility(*g, *f, host_scalar_coefficients());
    const auto rate = build_rate(*g, *t, f->z, {0});
    MultiFab tendency(g->ba, g->dm, 1, 0);
    ApplyMappedFluxTendency(*rate, 0, g->omega, 0, tendency, 0, g->dx_inv());
    for (amrex::MFIter mfi(tendency); mfi.isValid(); ++mfi) {
        const auto out = tendency.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real x = (Real(i) + Real(0.5)) * h;
            const Real y = (Real(j) + Real(0.5)) * h;
            const Real zc = height_at_center(k);
            const Real kappa = Real(0.5) + Real(0.2) * std::cos(two_pi * y);
            const Real dkappa = -Real(0.2) * two_pi * std::sin(two_pi * y);
            const Real sx = std::sin(two_pi * x);
            const Real horizontal =
                -Real(2.0) * kappa * two_pi * two_pi * sx * std::cos(two_pi * y) -
                dkappa * two_pi * sx * std::sin(two_pi * y);
            const Real vertical =
                -kappa * vertical_wave * vertical_wave * std::cos(vertical_wave * zc);
            out(i, j, k, 0) = amrex::Math::abs(out(i, j, k, 0) - horizontal - vertical);
        });
    }
    return tendency.norm0(0);
}

void run_refinement (const bool stretched)
{
    const Real e16 = manufactured_error(16, stretched);
    const Real e32 = manufactured_error(32, stretched);
    const Real e64 = manufactured_error(64, stretched);
    const Real order_coarse = std::log2(e16 / e32);
    const Real order_fine = std::log2(e32 / e64);
    amrex::Print() << "SBM_M5A_REFINEMENT mesh=" << (stretched ? "S" : "N")
                   << " e16=" << e16 << " e32=" << e32 << " e64=" << e64
                   << " order=" << order_coarse << "," << order_fine << '\n';
    EXPECT_GE(order_coarse, Real(1.8));
    EXPECT_GE(order_fine, Real(1.8));
}

// ERF's centred terrain-fitted stencil on affine terrain z = zeta + s x gives
// an empty cell a negative tendency from one diagonal neighbour.
void run_terrain_counterexample ()
{
    constexpr Real slope = Real(0.5);
    constexpr int n = 7;
    const Box domain(IntVect(0), IntVect(n - 1, 2, n - 1));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 3.0 / Real(n), 1.0});
    const std::array<int, 3> periodic{0, 1, 0};
    const Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, periodic.data());
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    amrex::BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    const BoxArray ba2d(std::move(boxes));
    const auto inv = geom.InvCellSizeArray();
    const Real h = Real(1.0) / Real(n);

    MultiFab z_nd(amrex::convert(ba, IntVect(1)), dm, 1, 1);
    for (amrex::MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
        const auto a = z_nd.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, 0) = Real(k) * h + slope * Real(i) * h;
        });
    }
    MultiFab z_cc(ba, dm, 1, 1);
    for (amrex::MFIter mfi(z_cc); mfi.isValid(); ++mfi) {
        const auto a = z_cc.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, 0) = (Real(k) + Real(0.5)) * h + slope * (Real(i) + Real(0.5)) * h;
        });
    }
    MultiFab scalar(ba, dm, 1, 2);
    MultiFab rho(ba, dm, 1, 2);
    MultiFab mu(ba, dm, EddyDiff::NumDiffs, 2);
    rho.setVal(Real(1.0));
    mu.setVal(Real(0.0));
    const IntVect center(3, 1, 3);
    for (amrex::MFIter mfi(scalar); mfi.isValid(); ++mfi) {
        const auto a = scalar.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, 0) = (i == center[0] + 1 && j == center[1] && k == center[2] + 1)
                ? Real(1.0) : Real(0.0);
        });
    }
    MultiFab one_x(amrex::convert(ba2d, IntVect(1, 0, 0)), dm, 1, 1);
    MultiFab one_y(amrex::convert(ba2d, IntVect(0, 1, 0)), dm, 1, 1);
    MultiFab one_m(ba2d, dm, 1, 1);
    MultiFab ax(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    MultiFab ay(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    MultiFab detj(ba, dm, 1, 0);
    for (auto* mf : {&one_x, &one_y, &one_m, &ax, &ay, &detj}) { mf->setVal(Real(1.0)); }
    MultiFab omega(ba, dm, 1, 0);
    std::string diagnostic;
    ASSERT_TRUE(BuildMappedCellMeasure(omega, detj, one_m, one_m, diagnostic)) << diagnostic;

    amrex::BCRec bc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int type = geom.isPeriodic(d) ? ERFBCType::int_dir : ERFBCType::foextrap;
        bc.setLo(d, type);
        bc.setHi(d, type);
    }
    amrex::Gpu::DeviceVector<amrex::BCRec> bcs(1);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, &bc, &bc + 1, bcs.begin());
    ScalarDiffusionFluxPolicy policy;
    policy.coefficients = {Real(1.0), EddyDiff::Scalar_h, EddyDiff::Scalar_v};
    policy.coefficient_mode = {false, false};

    MultiFab raw_x(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, IntVect(0, 0, 1));
    MultiFab raw_y(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, IntVect(0, 0, 1));
    MultiFab raw_z(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
    raw_x.setVal(Real(0.0));
    raw_y.setVal(Real(0.0));
    raw_z.setVal(Real(0.0));
    MappedFaceFluxRate mapped;
    mapped.define(ba, dm, 1, 0);
    for (amrex::MFIter mfi(scalar); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        ScalarDiffusionFieldViews field;
        field.scalar = scalar.const_array(mfi);
        field.density = rho.const_array(mfi);
        field.mu_turb = mu.const_array(mfi);
        field.xflux = raw_x.array(mfi);
        field.yflux = raw_y.array(mfi);
        field.zflux = raw_z.array(mfi);
        BuildScalarDiffusionFluxes_T(
            bx, domain, field, policy, false, z_nd.const_array(mfi), z_cc.const_array(mfi),
            inv, one_x.const_array(mfi), one_y.const_array(mfi), bcs.data(), 0);
        BuildScalarDiffusionMappedTransfers_T(
            bx, domain, raw_x.const_array(mfi), raw_y.const_array(mfi),
            raw_z.const_array(mfi), 0, mapped.dir(0).array(mfi), 0,
            mapped.dir(1).array(mfi), 0, mapped.dir(2).array(mfi), 0,
            z_nd.const_array(mfi), ax.const_array(mfi), ay.const_array(mfi), inv,
            one_m.const_array(mfi), one_x.const_array(mfi), one_m.const_array(mfi),
            one_y.const_array(mfi), false);
    }
    MultiFab tendency(ba, dm, 1, 0);
    ApplyMappedFluxTendency(mapped, 0, omega, 0, tendency, 0, inv);
    const Real center_tendency = value_at(tendency, 0, center);
    const Real pulse_tendency = value_at(tendency, 0, center + IntVect(1, 0, 1));
    amrex::Print() << "SBM_M5A_TERRAIN_COUNTEREXAMPLE slope=" << slope
                   << " empty_center_tendency=" << center_tendency
                   << " pulse_tendency=" << pulse_tendency << '\n';
    EXPECT_LT(pulse_tendency, Real(0.0));
    // The audited constant-coefficient value: -s/2 per unit cell area.
    EXPECT_NEAR(center_tendency, -Real(0.5) * slope / (h * h),
                Real(64.0) * eps() * std::abs(pulse_tendency));
}

erf_sbm::SBMLayout make_two_moment_layout ()
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
    population.moment_mode = erf_sbm::MomentMode::TwoMoment;
    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(std::move(population));
    spec.liquid_projection = {0, 1};
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
    return erf_sbm::SBMLayout(std::move(spec));
}

struct GroupState
{
    Real number;
    Real mean;
    Real property_fraction;
};

// Two walled cells; returns whether every group stays admissible after one
// update at the combined bound. A number-only transmissibility for the
// negative control breaks the shared-coefficient requirement.
bool groups_admissible_after_update (const bool shared_coefficient, std::string& failure)
{
    const auto layout = make_two_moment_layout();
    const auto groups = erf_sbm::make_constraint_groups(layout);
    const auto& pop = layout.populations().at(0);
    const int property_offset = layout.property_offset(0);
    const std::array<std::array<GroupState, 2>, 2> states{{
        {{{Real(1.0), Real(0.49), Real(0.8)}, {Real(2.0), Real(0.7), Real(0.6)}}},
        {{{Real(20.0), Real(0.11), Real(0.2)}, {Real(1.0), Real(0.95), Real(0.1)}}}}};

    GridOptions options;
    options.cells = IntVect(2, 1, 1);
    options.max_grid = 2;
    options.periodic = {0, 1, 1};
    options.variable_maps = false;
    const auto g = make_grid(options);
    const int ncomp = layout.ncomp();
    auto f = make_fields(*g, ncomp);
    f->rho.setVal(Real(1.0));
    f->eddy.setVal(Real(0.3));
    f->z.setVal(Real(0.0));
    for (int cell = 0; cell < 2; ++cell) {
        for (int bin = 0; bin < 2; ++bin) {
            const auto& s = states[static_cast<std::size_t>(cell)][static_cast<std::size_t>(bin)];
            const std::array<std::pair<int, Real>, 3> values{{
                {pop.number_offset + bin, s.number},
                {pop.mass_offset + bin, s.number * s.mean},
                {property_offset + bin, s.number * s.property_fraction * Real(0.5)}}};
            for (const auto& entry : values) {
                const int c = cell;
                const int component = entry.first;
                const Real v = entry.second;
                for (amrex::MFIter mfi(f->z); mfi.isValid(); ++mfi) {
                    const auto a = f->z.array(mfi);
                    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (i == c) { a(i, j, k, component) = v; }
                    });
                }
            }
        }
    }
    f->z.FillBoundary(g->geom.periodicity());
    for (int cell = 0; cell < 2; ++cell) {
        std::vector<Real> state(static_cast<std::size_t>(ncomp));
        for (int n = 0; n < ncomp; ++n) {
            state[static_cast<std::size_t>(n)] = value_at(f->z, n, IntVect(cell, 0, 0));
        }
        for (const auto& group : groups) {
            EXPECT_TRUE(group.admissible(state)) << "initial state is not admissible";
        }
    }

    const auto t = build_transmissibility(*g, *f, host_scalar_coefficients());
    MappedFaceFluxRate no_carrier;
    no_carrier.define(g->ba, g->dm, 1, 0);
    no_carrier.setVal(Real(0.0));
    const Real tau = Real(1.0) / combined_rate(*g, *f, *t, no_carrier);

    std::vector<int> all(static_cast<std::size_t>(ncomp));
    for (int n = 0; n < ncomp; ++n) { all[static_cast<std::size_t>(n)] = n; }
    auto rate = build_rate(*g, *t, f->z, all);
    if (!shared_coefficient) {
        for (int n = 0; n < ncomp; ++n) {
            const bool number = n >= pop.number_offset && n < pop.number_offset + 2;
            if (number) { continue; }
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                rate->dir(dir).setVal(Real(0.0), n, 1, 0);
            }
        }
    }
    auto u = extensive_state(*g, *f);
    explicit_update(*g, *rate, *u, tau);
    for (int cell = 0; cell < 2; ++cell) {
        std::vector<Real> state(static_cast<std::size_t>(ncomp));
        for (int n = 0; n < ncomp; ++n) {
            state[static_cast<std::size_t>(n)] = value_at(*u, n, IntVect(cell, 0, 0));
        }
        for (const auto& group : groups) {
            std::string failed;
            if (!group.admissible(state, nullptr, &failed)) {
                failure = "cell " + std::to_string(cell) + " group " + group.semantic_id + ": " + failed;
                return false;
            }
        }
    }
    return true;
}

void run_group_cone ()
{
    std::string failure;
    EXPECT_TRUE(groups_admissible_after_update(true, failure)) << failure;
    failure.clear();
    EXPECT_FALSE(groups_admissible_after_update(false, failure));
    amrex::Print() << "SBM_M5A_GROUP_CONE number_only_exchange_breaks=" << failure << '\n';
}

erf_sbm::ParticleMixingConfig policy_config (const erf_sbm::ParticleMixingPolicy policy,
                                             const Real kh = Real(0.0),
                                             const Real kv = Real(0.0))
{
    erf_sbm::ParticleMixingConfig config;
    config.policy = policy;
    config.reference_kh = kh;
    config.reference_kv = kv;
    return config;
}

bool build_policy (const erf_sbm::ParticleMixingConfig& config, const DiffusionGrid& g,
                   const MultiFab& rho, const MultiFab* eddy,
                   MappedFaceTransmissibility& t, std::string& diagnostic)
{
    t.define(g.ba, g.dm, 1, 0);
    return erf_sbm::build_particle_mixing_transmissibility(config, g.geom, g.metrics(), rho, 0,
                                                          eddy, t, diagnostic);
}

// Independent face oracle: T = 0.5 (a_L + a_R) w * geometry, with channel
// `hcomp` on x/y faces, `vcomp` on z faces and w the face density if requested.
std::unique_ptr<MappedFaceTransmissibility>
oracle (const DiffusionGrid& g, const MultiFab& coefficient, const int hcomp, const int vcomp,
        const MultiFab* density_weight, const Real factor = Real(1.0))
{
    auto t = std::make_unique<MappedFaceTransmissibility>();
    t->define(g.ba, g.dm, 1, 0);
    const auto inv = g.dx_inv();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const IntVect shift = IntVect::TheDimensionVector(dir);
        const int comp = dir < 2 ? hcomp : vcomp;
        const bool weighted = density_weight != nullptr;
        for (amrex::MFIter mfi(t->dir(dir)); mfi.isValid(); ++mfi) {
            const auto out = t->dir(dir).array(mfi);
            const auto a = coefficient.const_array(mfi);
            const auto w = weighted ? density_weight->const_array(mfi) : Array4<const Real>{};
            const auto num = dir == 0 ? g.mf_ux.const_array(mfi) : g.mf_vy.const_array(mfi);
            const auto den = dir == 0 ? g.mf_uy.const_array(mfi) : g.mf_vx.const_array(mfi);
            const auto mx = g.mf_mx.const_array(mfi);
            const auto my = g.mf_my.const_array(mfi);
            const Real spacing = inv[dir];
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const IntVect cell(AMREX_D_DECL(i, j, k));
                const IntVect left = cell - shift;
                Real coeff = Real(0.5) * (a(cell, comp) + a(left, comp));
                if (weighted) { coeff *= Real(0.5) * (w(cell, 0) + w(left, 0)); }
                const Real geometry = dir < 2 ? spacing * num(i, j, 0) / den(i, j, 0)
                                              : spacing / (mx(i, j, 0) * my(i, j, 0));
                out(i, j, k, 0) = factor * coeff * geometry;
            });
        }
    }
    return t;
}

void run_policy_names ()
{
    using erf_sbm::ParticleMixingPolicy;
    for (const auto policy : {ParticleMixingPolicy::None, ParticleMixingPolicy::ConstantReference,
                              ParticleMixingPolicy::HostPassiveScalar}) {
        ParticleMixingPolicy parsed = ParticleMixingPolicy::None;
        EXPECT_TRUE(erf_sbm::parse_particle_mixing_policy(erf_sbm::particle_mixing_name(policy),
                                                         parsed));
        EXPECT_EQ(parsed, policy);
    }
    ParticleMixingPolicy parsed = ParticleMixingPolicy::None;
    EXPECT_FALSE(erf_sbm::parse_particle_mixing_policy("hostpassivescalar", parsed));
    EXPECT_FALSE(erf_sbm::parse_particle_mixing_policy("HostMoisture", parsed));

    EXPECT_EQ(policy_config(ParticleMixingPolicy::None).identity(), "particle-mixing=none-v1");
    EXPECT_EQ(policy_config(ParticleMixingPolicy::HostPassiveScalar).identity(),
              "particle-mixing=host-passive-scalar-v1[Scalar_h,Scalar_v]");
    const auto reference = policy_config(ParticleMixingPolicy::ConstantReference,
                                         Real(0.25), Real(0.125));
    EXPECT_EQ(reference.identity(),
              "particle-mixing=constant-reference-v1[kh=" + erf_sbm::canonical_real(Real(0.25)) +
              ",kv=" + erf_sbm::canonical_real(Real(0.125)) + "]");
    EXPECT_NE(reference.identity(),
              policy_config(ParticleMixingPolicy::ConstantReference, Real(0.25),
                            std::nextafter(Real(0.125), Real(1.0))).identity());

    std::string reason;
    EXPECT_TRUE(erf_sbm::validate_particle_mixing_config(reference, reason)) << reason;
    EXPECT_FALSE(erf_sbm::validate_particle_mixing_config(
        policy_config(ParticleMixingPolicy::ConstantReference, Real(-1.0e-30), Real(0.1)), reason));
    EXPECT_FALSE(erf_sbm::validate_particle_mixing_config(
        policy_config(ParticleMixingPolicy::ConstantReference, Real(0.1),
                      std::numeric_limits<Real>::quiet_NaN()), reason));
    EXPECT_FALSE(erf_sbm::validate_particle_mixing_config(
        policy_config(ParticleMixingPolicy::HostPassiveScalar, Real(0.1)), reason));
    EXPECT_NE(reason.find("only to ConstantReference"), std::string::npos) << reason;
}

void run_host_channels ()
{
    const auto g = make_grid(GridOptions{});
    const auto f = make_fields(*g, 1);
    // The fixture gives every eddy channel a different field and rho != 1.
    MappedFaceTransmissibility t;
    std::string diagnostic;
    ASSERT_TRUE(build_policy(policy_config(erf_sbm::ParticleMixingPolicy::HostPassiveScalar),
                             *g, f->rho, &f->eddy, t, diagnostic)) << diagnostic;
    const auto expected = oracle(*g, f->eddy, EddyDiff::Scalar_h, EddyDiff::Scalar_v, nullptr);
    const Real scale = max_abs_face(*expected, 0);
    EXPECT_LE(max_interior_face_difference(t, 0, *expected, 0, g->geom), Real(8.0) * eps() * scale);

    // Each wrong reading differs materially, so the oracle above discriminates.
    const auto times_rho = oracle(*g, f->eddy, EddyDiff::Scalar_h, EddyDiff::Scalar_v, &f->rho);
    const auto times_schmidt = oracle(*g, f->eddy, EddyDiff::Scalar_h, EddyDiff::Scalar_v,
                                      nullptr, Real(3.0));
    const auto moisture = oracle(*g, f->eddy, EddyDiff::Q_h, EddyDiff::Q_v, nullptr);
    const auto momentum = oracle(*g, f->eddy, EddyDiff::Mom_h, EddyDiff::Mom_v, nullptr);
    const auto isotropic = oracle(*g, f->eddy, EddyDiff::Scalar_h, EddyDiff::Scalar_h, nullptr);
    for (const auto* wrong : {times_rho.get(), times_schmidt.get(), moisture.get(),
                              momentum.get(), isotropic.get()}) {
        EXPECT_GT(max_interior_face_difference(t, 0, *wrong, 0, g->geom), Real(1.0e-3) * scale);
    }
}

void run_channels_stay_separate ()
{
    const auto g = make_grid(GridOptions{});
    auto f = make_fields(*g, 1);
    MappedFaceTransmissibility t;
    std::string diagnostic;
    const auto config = policy_config(erf_sbm::ParticleMixingPolicy::HostPassiveScalar);
    f->eddy.setVal(Real(0.0), EddyDiff::Scalar_v, 1, 1);
    ASSERT_TRUE(build_policy(config, *g, f->rho, &f->eddy, t, diagnostic)) << diagnostic;
    EXPECT_GT(t.dir(0).min(0), Real(0.0));
    EXPECT_GT(t.dir(1).min(0), Real(0.0));
    EXPECT_EQ(t.dir(2).norm0(0), Real(0.0));

    f = make_fields(*g, 1);
    f->eddy.setVal(Real(0.0), EddyDiff::Scalar_h, 1, 1);
    ASSERT_TRUE(build_policy(config, *g, f->rho, &f->eddy, t, diagnostic)) << diagnostic;
    EXPECT_EQ(t.dir(0).norm0(0), Real(0.0));
    EXPECT_EQ(t.dir(1).norm0(0), Real(0.0));
    EXPECT_GT(t.dir(2).min(0), Real(0.0));
}

void run_constant_reference ()
{
    using erf_sbm::ParticleMixingPolicy;
    erf_auxiliary::MappedDiffusionCoefficients coefficients;
    std::string reason;
    ASSERT_TRUE(erf_sbm::particle_mixing_coefficients(
        policy_config(ParticleMixingPolicy::HostPassiveScalar), coefficients, reason)) << reason;
    EXPECT_EQ(coefficients.horizontal.molecular_coeff, Real(0.0));
    EXPECT_EQ(coefficients.vertical.molecular_coeff, Real(0.0));
    EXPECT_FALSE(coefficients.mode.molecular_is_kinematic);
    EXPECT_TRUE(coefficients.mode.add_turbulence);
    EXPECT_FALSE(erf_sbm::particle_mixing_coefficients(
        policy_config(ParticleMixingPolicy::None), coefficients, reason));

    const auto g = make_grid(GridOptions{});
    const auto f = make_fields(*g, 1);
    constexpr Real kh = Real(0.37);
    constexpr Real kv = Real(0.11);
    MappedFaceTransmissibility t;
    std::string diagnostic;
    ASSERT_TRUE(build_policy(policy_config(ParticleMixingPolicy::ConstantReference, kh, kv),
                             *g, f->rho, nullptr, t, diagnostic)) << diagnostic;
    MultiFab reference(g->ba, g->dm, 2, 1);
    reference.setVal(kh, 0, 1, 1);
    reference.setVal(kv, 1, 1, 1);
    const auto expected = oracle(*g, reference, 0, 1, &f->rho);
    EXPECT_LE(max_interior_face_difference(t, 0, *expected, 0, g->geom),
              Real(8.0) * eps() * max_abs_face(*expected, 0));
}

void run_invalid_coefficients ()
{
    using erf_sbm::ParticleMixingPolicy;
    GridOptions options;
    options.periodic = {1, 1, 0};
    const auto g = make_grid(options);
    const auto host = policy_config(ParticleMixingPolicy::HostPassiveScalar);
    MappedFaceTransmissibility t;
    std::string diagnostic;

    auto f = make_fields(*g, 1);
    EXPECT_FALSE(build_policy(host, *g, f->rho, nullptr, t, diagnostic));
    EXPECT_NE(diagnostic.find("eddy diffusivity"), std::string::npos) << diagnostic;

    MultiFab no_ghost(g->ba, g->dm, EddyDiff::NumDiffs, 0);
    no_ghost.setVal(Real(0.5));
    EXPECT_FALSE(build_policy(host, *g, f->rho, &no_ghost, t, diagnostic));

    auto poke = [&] (const IntVect& iv, const int comp, const Real value) {
        f = make_fields(*g, 1);
        for (amrex::MFIter mfi(f->eddy); mfi.isValid(); ++mfi) {
            const auto a = f->eddy.array(mfi);
            const Box point(iv, iv);
            const Box bx = mfi.fabbox() & point;
            if (!bx.ok()) { continue; }
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                a(i, j, k, comp) = value;
            });
        }
        return build_policy(host, *g, f->rho, &f->eddy, t, diagnostic);
    };
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real inf = std::numeric_limits<Real>::infinity();
    EXPECT_FALSE(poke(IntVect(3, 2, 1), EddyDiff::Scalar_h, nan));
    EXPECT_FALSE(poke(IntVect(3, 2, 1), EddyDiff::Scalar_v, inf));
    EXPECT_FALSE(poke(IntVect(3, 2, 1), EddyDiff::Scalar_v, -Real(1.0e-30)));
    // A periodic ghost cell is read by the face on the periodic boundary.
    EXPECT_FALSE(poke(IntVect(-1, 2, 1), EddyDiff::Scalar_h, -Real(1.0e-30)));
    EXPECT_NE(diagnostic.find("invalid reads"), std::string::npos) << diagnostic;
    // A ghost beyond a physical face, or an unused channel, is never read.
    EXPECT_TRUE(poke(IntVect(3, 2, -1), EddyDiff::Scalar_v, nan)) << diagnostic;
    EXPECT_TRUE(poke(IntVect(3, 2, 1), EddyDiff::Q_h, nan)) << diagnostic;

    f = make_fields(*g, 1);
    const auto reference = policy_config(ParticleMixingPolicy::ConstantReference, Real(0.2), Real(0.1));
    EXPECT_TRUE(build_policy(reference, *g, f->rho, nullptr, t, diagnostic)) << diagnostic;
    fill_cells(f->rho, 0, g->geom, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        return (i == 2 && j == 3 && k == 1) ? Real(0.0) : Real(1.0);
    });
    EXPECT_FALSE(build_policy(reference, *g, f->rho, nullptr, t, diagnostic));
    EXPECT_FALSE(build_policy(policy_config(ParticleMixingPolicy::ConstantReference, Real(-0.1),
                                            Real(0.1)), *g, f->rho, nullptr, t, diagnostic));
}

void run_geometry_rejections ()
{
    erf_auxiliary::DiffusionMesh mesh = erf_auxiliary::DiffusionMesh::ConstantDz;
    std::string reason;
    EXPECT_TRUE(erf_sbm::particle_mixing_mesh(MeshType::ConstantDz, mesh, reason));
    EXPECT_EQ(mesh, erf_auxiliary::DiffusionMesh::ConstantDz);
    EXPECT_TRUE(erf_sbm::particle_mixing_mesh(MeshType::StretchedDz, mesh, reason));
    EXPECT_EQ(mesh, erf_auxiliary::DiffusionMesh::StretchedDz);
    EXPECT_FALSE(erf_sbm::particle_mixing_mesh(MeshType::VariableDz, mesh, reason));
    EXPECT_NE(reason.find("terrain-fitted"), std::string::npos) << reason;

    GridOptions options;
    options.periodic = {1, 1, 1};
    options.stretched = true;
    const auto g = make_grid(options);
    const auto f = make_fields(*g, 1);
    MappedFaceTransmissibility t;
    std::string diagnostic;
    EXPECT_FALSE(build_policy(policy_config(erf_sbm::ParticleMixingPolicy::HostPassiveScalar),
                              *g, f->rho, &f->eddy, t, diagnostic));
    EXPECT_NE(diagnostic.find("periodic z"), std::string::npos) << diagnostic;
}

void run_none_policy ()
{
    const auto g = make_grid(GridOptions{});
    const auto f = make_fields(*g, 1);
    MappedFaceTransmissibility t;
    std::string diagnostic;
    ASSERT_TRUE(build_policy(policy_config(erf_sbm::ParticleMixingPolicy::None), *g, f->rho,
                             nullptr, t, diagnostic)) << diagnostic;
    EXPECT_EQ(max_abs_face(t, 0), Real(0.0));
}

} // namespace

TEST(AuxiliaryDiffusion, NGridTransmissibilityMatchesHostDiffusion) { run_ngrid_parity(); }
TEST(AuxiliaryDiffusion, StretchedGridTransmissibilityMatchesHostDiffusion) { run_stretched_parity(); }
TEST(AuxiliaryDiffusion, PhysicalFacesCarryNoDiffusiveTransfer) { run_physical_faces(); }
TEST(AuxiliaryDiffusion, ZeroCoefficientIsTheIdentity) { run_zero_coefficient_identity(); }
TEST(AuxiliaryDiffusion, ExchangeConservesInventory)
{
    run_conservation(false);
    run_conservation(true);
}
TEST(AuxiliaryDiffusion, UpdateIsAPositiveLinearCombination) { run_positive_linear_combination(); }
TEST(AuxiliaryDiffusion, SparsePulseStaysNonnegativeUpToTheBound) { run_sparse_pulse(); }
TEST(AuxiliaryDiffusion, CombinedBoundCatchesSeparatelyPassingTerms)
{
    run_combined_bound_case(Real(0.6), true);
    run_combined_bound_case(Real(0.4), false);
}
TEST(AuxiliaryDiffusion, DissipationIdentityAndEnergyDecay) { run_dissipation(); }
TEST(AuxiliaryDiffusion, NGridManufacturedRefinement) { run_refinement(false); }
TEST(AuxiliaryDiffusion, StretchedGridManufacturedRefinement) { run_refinement(true); }
TEST(AuxiliaryDiffusion, TerrainCentredStencilIsNotPositive) { run_terrain_counterexample(); }

TEST(SBMParticleMixing, PolicyNamesAndIdentity) { run_policy_names(); }
TEST(SBMParticleMixing, HostPassiveScalarUsesScalarChannelsAsGiven) { run_host_channels(); }
TEST(SBMParticleMixing, HorizontalAndVerticalChannelsStaySeparate) { run_channels_stay_separate(); }
TEST(SBMParticleMixing, ConstantReferenceAndNoInheritedMolecularDiffusion) { run_constant_reference(); }
TEST(SBMParticleMixing, InvalidCoefficientsFailClosed) { run_invalid_coefficients(); }
TEST(SBMParticleMixing, TerrainAndPeriodicStretchedGeometryAreRejected) { run_geometry_rejections(); }
TEST(SBMParticleMixing, NonePolicyHasNoTransmissibility) { run_none_policy(); }
TEST(SBMParticleMixing, SharedCoefficientPreservesEveryBinGroup) { run_group_cone(); }
