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
#include "ERF_SBMTransport.H"
#include "ERF_ScalarDiffusion.H"
#include "ERF_Utils.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
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
    bool terrain{false};
    std::vector<Real> dz_host;
    amrex::Gpu::DeviceVector<Real> dz;
    MultiFab mf_ux, mf_uy, mf_vx, mf_vy, mf_mx, mf_my;
    MultiFab ax, ay, detj, omega;
    MultiFab z_nd, z_cc;

    [[nodiscard]] MappedDiffusionMetrics metrics () const
    {
        MappedDiffusionMetrics m;
        m.mesh = terrain ? DiffusionMesh::VariableDz
               : stretched ? DiffusionMesh::StretchedDz : DiffusionMesh::ConstantDz;
        m.mf_ux = &mf_ux;
        m.mf_uy = &mf_uy;
        m.mf_vx = &mf_vx;
        m.mf_vy = &mf_vy;
        m.mf_mx = &mf_mx;
        m.mf_my = &mf_my;
        m.ax = &ax;
        m.ay = &ay;
        m.stretched_dz = stretched ? dz.data() : nullptr;
        m.z_nd = terrain ? &z_nd : nullptr;
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

std::unique_ptr<DiffusionFields> make_fields (const DiffusionGrid& g, const int ncomp,
                                              const int ngrow = 1)
{
    auto f = std::make_unique<DiffusionFields>();
    f->rho.define(g.ba, g.dm, 1, ngrow);
    f->eddy.define(g.ba, g.dm, EddyDiff::NumDiffs, ngrow);
    f->z.define(g.ba, g.dm, ncomp, ngrow);
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
                amrex::Array4<const amrex::Real>{}, amrex::Array4<const amrex::Real>{},
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

// The fields make_fields writes, at any index: an open face reads its ghost.
struct FieldFormulas
{
    Real nx, ny;
    Real rho (const int i, const int j, const int k) const
    {
        return Real(1.1) + Real(0.15) * std::sin(two_pi * Real(i) / nx) +
               Real(0.08) * std::cos(two_pi * Real(j) / ny) + Real(0.03) * Real(k);
    }
    Real eddy (const int c, const int i, const int j, const int k) const
    {
        return Real(0.4) + Real(0.05) * Real(c) + Real(0.1) * std::sin(two_pi * Real(i + c) / nx) +
               Real(0.07) * std::cos(two_pi * Real(j) / ny) + Real(0.02) * Real(k);
    }
    Real z (const int c, const int i, const int j, const int k) const
    {
        return Real(0.5) + Real(0.3) * std::sin(two_pi * Real(i) / nx + Real(c)) *
                   std::cos(two_pi * Real(j) / ny) +
               Real(0.1) * Real(k) + Real(0.05) * Real(c);
    }
};

// Positive (outward at the high face) on even cells of the face, negative on odd.
Real carrier_value (const int i, const int j, const int k)
{
    return (i + j + k) % 2 == 0 ? Real(0.3) : Real(-0.2);
}

// An open domain face exchanges over half the edge cell with the host face
// coefficient of the ghost and edge cells; an inward-only face closes where
// the carrier is outward. x faces on an N grid, z faces on a stretched one.
void expect_open_faces (const bool stretched)
{
    GridOptions options;
    options.cells = IntVect(6, 5, 4);
    options.max_grid = 3;
    options.periodic = stretched ? std::array<int, 3>{1, 1, 0} : std::array<int, 3>{0, 1, 1};
    options.stretched = stretched;
    const auto g = make_grid(options);
    const auto f = make_fields(*g, 1);
    const int dir = stretched ? 2 : 0;
    MappedFaceFluxRate carrier;
    carrier.define(g->ba, g->dm, 1, 0);
    for (amrex::MFIter mfi(carrier.dir(dir)); mfi.isValid(); ++mfi) {
        const auto a = carrier.dir(dir).array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, 0) = (i + j + k) % 2 == 0 ? Real(0.3) : Real(-0.2);
        });
    }
    MappedDiffusionBoundary boundary;
    boundary.kind[static_cast<std::size_t>(dir)] = DiffusionFaceKind::Dirichlet;
    boundary.kind[static_cast<std::size_t>(dir + AMREX_SPACEDIM)] =
        DiffusionFaceKind::DirichletWhereInward;
    boundary.carrier[static_cast<std::size_t>(dir)] = &carrier.dir(dir);

    const FieldFormulas formula{Real(options.cells[0]), Real(options.cells[1])};
    const int n = options.cells[dir];
    const auto inv = g->dx_inv();
    for (const bool kinematic : {false, true}) {
        const auto c = kinematic ? kinematic_coefficients(Real(0.37), Real(0.37))
                                 : host_scalar_coefficients();
        auto t = std::make_unique<MappedFaceTransmissibility>();
        t->define(g->ba, g->dm, 1, 0);
        std::string diagnostic;
        ASSERT_TRUE(BuildMappedDiffusionTransmissibility(*t, g->geom, g->metrics(), c, f->rho, 0,
                                                         &f->eddy, boundary, diagnostic))
            << diagnostic;
        const auto rate = build_rate(*g, *t, f->z, {0});
        const int eddy_comp = stretched ? EddyDiff::Scalar_v : EddyDiff::Scalar_h;
        int open_cells = 0;
        int closed_cells = 0;
        for (const int face : {0, n}) {
            for (int a = 0; a < options.cells[dir == 0 ? 1 : 0]; ++a) {
                for (int b = 0; b < options.cells[dir == 2 ? 1 : 2]; ++b) {
                    const IntVect iv = dir == 0 ? IntVect(face, a, b) : IntVect(a, b, face);
                    const IntVect edge = iv - (face == 0 ? IntVect(0) : IntVect::TheDimensionVector(dir));
                    const IntVect ghost = face == 0 ? iv - IntVect::TheDimensionVector(dir) : iv;
                    const auto avg = [&] (auto value) {
                        return Real(0.5) * (value(ghost) + value(edge));
                    };
                    const Real coeff = kinematic
                        ? Real(0.37) * avg([&] (const IntVect& v) { return formula.rho(v[0], v[1], v[2]); })
                        : avg([&] (const IntVect& v) { return formula.eddy(eddy_comp, v[0], v[1], v[2]); });
                    Real expected;
                    if (stretched) {
                        const Real dz = g->dz_host[static_cast<std::size_t>(edge[2])];
                        expected = coeff * Real(2.0) / dz /
                                   (value_at(g->mf_mx, 0, IntVect(a, b, 0)) *
                                    value_at(g->mf_my, 0, IntVect(a, b, 0)));
                    } else {
                        expected = coeff * Real(2.0) * inv[0] *
                                   value_at(g->mf_ux, 0, IntVect(face, a, 0)) /
                                   value_at(g->mf_uy, 0, IntVect(face, a, 0));
                    }
                    const bool open = face == 0 || carrier_value(iv[0], iv[1], iv[2]) <= Real(0.0);
                    const Real got = value_at(t->dir(dir), 0, iv);
                    if (!open) {
                        EXPECT_EQ(got, Real(0.0)) << iv;
                        EXPECT_EQ(value_at(rate->dir(dir), 0, iv), Real(0.0)) << iv;
                        ++closed_cells;
                        continue;
                    }
                    ++open_cells;
                    EXPECT_NEAR(got, expected, Real(16.0) * eps() * expected) << iv;
                    // The low-order rate exchanges with the ghost value, which is
                    // the value held on the face.
                    const Real dz_value = formula.z(0, iv[0], iv[1], iv[2]) -
                        formula.z(0, iv[0] - (dir == 0), iv[1], iv[2] - (dir == 2));
                    EXPECT_NEAR(value_at(rate->dir(dir), 0, iv), -got * dz_value,
                                Real(64.0) * eps() * std::abs(got));
                }
            }
        }
        EXPECT_GT(open_cells, 0);
        EXPECT_GT(closed_cells, 0);
    }
}

void run_open_faces () { expect_open_faces(false); expect_open_faces(true); }

// The validation reads ghost cells beyond open faces only, and an inward-only
// face cannot judge its carrier without one.
void run_open_face_validation ()
{
    GridOptions options;
    options.cells = IntVect(6, 5, 4);
    options.max_grid = 3;
    options.periodic = {0, 1, 1};
    const auto g = make_grid(options);
    const auto f = make_fields(*g, 1);
    for (amrex::MFIter mfi(f->eddy); mfi.isValid(); ++mfi) {
        const auto a = f->eddy.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (i == -1 && j == 2 && k == 1) {
                a(i, j, k, EddyDiff::Scalar_h) = std::numeric_limits<Real>::quiet_NaN();
            }
        });
    }
    MappedFaceTransmissibility t;
    t.define(g->ba, g->dm, 1, 0);
    std::string diagnostic;
    MappedDiffusionBoundary boundary;
    EXPECT_TRUE(BuildMappedDiffusionTransmissibility(t, g->geom, g->metrics(),
                                                     host_scalar_coefficients(), f->rho, 0,
                                                     &f->eddy, boundary, diagnostic))
        << diagnostic;
    boundary.kind[0] = DiffusionFaceKind::Dirichlet;
    EXPECT_FALSE(BuildMappedDiffusionTransmissibility(t, g->geom, g->metrics(),
                                                      host_scalar_coefficients(), f->rho, 0,
                                                      &f->eddy, boundary, diagnostic));
    EXPECT_NE(diagnostic.find("1 invalid reads"), std::string::npos) << diagnostic;

    boundary = {};
    boundary.kind[AMREX_SPACEDIM] = DiffusionFaceKind::DirichletWhereInward;
    EXPECT_FALSE(BuildMappedDiffusionTransmissibility(t, g->geom, g->metrics(),
                                                      host_scalar_coefficients(), f->rho, 0,
                                                      &f->eddy, boundary, diagnostic));
    EXPECT_NE(diagnostic.find("requires the dry-air carrier"), std::string::npos) << diagnostic;
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

// ERF's centered terrain-fitted stencil on affine terrain z = zeta + s x gives
// an empty cell a negative tendency from one diagonal neighbor.
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
            inv, amrex::Array4<const amrex::Real>{}, amrex::Array4<const amrex::Real>{},
            one_x.const_array(mfi), one_y.const_array(mfi), bcs.data(), 0);
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

// Static affine terrain z = height zeta + slope_x xi + slope_y eta on the unit
// computational cube: analytic nodal heights, with the cell-centre heights,
// Jacobian and side areas ERF derives from them. Every face then has
// h_xi = slope_x, h_eta = slope_y and J = height.
struct TerrainOptions
{
    IntVect cells{8, 6, 8};
    int max_grid{4};
    std::array<int, 3> periodic{0, 1, 0};
    Real height{Real(0.8)};
    Real slope_x{Real(0.6)};
    Real slope_y{Real(-0.35)};
    bool variable_maps{true};
    // Every x (y) map factor when the maps are not variable.
    Real map_x{Real(1.0)};
    Real map_y{Real(1.0)};
};

std::unique_ptr<DiffusionGrid> make_terrain_grid (const TerrainOptions& options)
{
    GridOptions grid;
    grid.cells = options.cells;
    grid.max_grid = options.max_grid;
    grid.periodic = options.periodic;
    grid.variable_maps = options.variable_maps;
    auto g = make_grid(grid);
    if (!options.variable_maps) {
        for (auto* mf : {&g->mf_ux, &g->mf_vx, &g->mf_mx}) { mf->setVal(options.map_x); }
        for (auto* mf : {&g->mf_uy, &g->mf_vy, &g->mf_my}) { mf->setVal(options.map_y); }
    }
    g->terrain = true;
    g->z_nd.define(amrex::convert(g->ba, IntVect(1)), g->dm, 1, 2);
    const Real nx = Real(options.cells[0]);
    const Real ny = Real(options.cells[1]);
    const Real nz = Real(options.cells[2]);
    const Real height = options.height;
    const Real slope_x = options.slope_x;
    const Real slope_y = options.slope_y;
    for (amrex::MFIter mfi(g->z_nd); mfi.isValid(); ++mfi) {
        const auto a = g->z_nd.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k, 0) = height * Real(k) / nz + slope_x * Real(i) / nx +
                            slope_y * Real(j) / ny;
        });
    }
    g->z_cc.define(g->ba, g->dm, 1, 2);
    make_zcc(g->geom, g->z_nd, g->z_cc);
    make_J(g->geom, g->z_nd, g->detj);
    MultiFab az(amrex::convert(g->ba, IntVect(0, 0, 1)), g->dm, 1, 0);
    make_areas(g->geom, g->z_nd, g->ax, g->ay, az);
    std::string diagnostic;
    EXPECT_TRUE(BuildMappedCellMeasure(g->omega, g->detj, g->mf_mx, g->mf_my, diagnostic))
        << diagnostic;
    return g;
}

// Constant density and no eddy diffusivity, with the two ghost cells the
// native terrain flux reads.
std::unique_ptr<DiffusionFields> make_uniform_terrain_fields (const DiffusionGrid& g,
                                                              const Real rho)
{
    auto f = std::make_unique<DiffusionFields>();
    f->rho.define(g.ba, g.dm, 1, 2);
    f->eddy.define(g.ba, g.dm, EddyDiff::NumDiffs, 2);
    f->z.define(g.ba, g.dm, 1, 2);
    f->rho.setVal(rho);
    f->eddy.setVal(Real(0.0));
    f->z.setVal(Real(0.0));
    return f;
}

// Largest |a - value| over the valid region.
Real max_deviation (const MultiFab& a, const int comp, const Real value)
{
    MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    MultiFab::Copy(d, a, comp, 0, 1, 0);
    d.plus(-value, 0, 1, 0);
    return d.norm0(0);
}

Real max_difference (const MultiFab& a, const MultiFab& b)
{
    MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    MultiFab::Copy(d, a, 0, 0, 1, 0);
    MultiFab::Subtract(d, b, 0, 0, 1, 0);
    return d.norm0(0);
}

// ERF's native terrain diffusion of component 0 of `scalar`: raw face fluxes,
// then one terrain mapping. A null `vertical` takes the one-coefficient
// builder, whose z faces read policy.coefficients.
struct TerrainFluxes
{
    MultiFab raw_x, raw_y, raw_z;
    MappedFaceFluxRate mapped;
};

std::unique_ptr<TerrainFluxes>
terrain_fluxes (const DiffusionGrid& g, const DiffusionFields& f,
                const ScalarDiffusionFluxPolicy& policy,
                const ScalarDiffusionCoefficients* vertical)
{
    amrex::BCRec bc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int type = g.geom.isPeriodic(d) ? ERFBCType::int_dir : ERFBCType::foextrap;
        bc.setLo(d, type);
        bc.setHi(d, type);
    }
    amrex::Gpu::DeviceVector<amrex::BCRec> bcs(1);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, &bc, &bc + 1, bcs.begin());

    auto out = std::make_unique<TerrainFluxes>();
    out->raw_x.define(amrex::convert(g.ba, IntVect(1, 0, 0)), g.dm, 1, IntVect(0, 0, 1));
    out->raw_y.define(amrex::convert(g.ba, IntVect(0, 1, 0)), g.dm, 1, IntVect(0, 0, 1));
    out->raw_z.define(amrex::convert(g.ba, IntVect(0, 0, 1)), g.dm, 1, 0);
    out->mapped.define(g.ba, g.dm, 1, 0);
    const auto inv = g.dx_inv();
    for (amrex::MFIter mfi(f.z); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        ScalarDiffusionFieldViews field;
        field.scalar = f.z.const_array(mfi);
        field.density = f.rho.const_array(mfi);
        field.mu_turb = f.eddy.const_array(mfi);
        field.xflux = out->raw_x.array(mfi);
        field.yflux = out->raw_y.array(mfi);
        field.zflux = out->raw_z.array(mfi);
        const auto z_nd = g.z_nd.const_array(mfi);
        const auto z_cc = g.z_cc.const_array(mfi);
        // No ext_dir_upwind face, so no velocity is read.
        const Array4<const Real> no_velocity{};
        if (vertical != nullptr) {
            BuildScalarDiffusionFluxes_T(bx, g.domain, field, policy, *vertical, false, z_nd,
                                         z_cc, inv, no_velocity, no_velocity,
                                         g.mf_ux.const_array(mfi), g.mf_vy.const_array(mfi),
                                         bcs.data(), 0);
        } else {
            BuildScalarDiffusionFluxes_T(bx, g.domain, field, policy, false, z_nd, z_cc, inv,
                                         no_velocity, no_velocity, g.mf_ux.const_array(mfi),
                                         g.mf_vy.const_array(mfi), bcs.data(), 0);
        }
        BuildScalarDiffusionMappedTransfers_T(
            bx, g.domain, out->raw_x.const_array(mfi), out->raw_y.const_array(mfi),
            out->raw_z.const_array(mfi), 0, out->mapped.dir(0).array(mfi), 0,
            out->mapped.dir(1).array(mfi), 0, out->mapped.dir(2).array(mfi), 0, z_nd,
            g.ax.const_array(mfi), g.ay.const_array(mfi), inv, g.mf_mx.const_array(mfi),
            g.mf_uy.const_array(mfi), g.mf_my.const_array(mfi), g.mf_vx.const_array(mfi), false);
    }
    amrex::Gpu::streamSynchronize();
    return out;
}

// Each terrain face matches an oracle written from the affine slopes, the
// Jacobian and the host face averages: a (m_num / m_den) mu_h / dxi on the
// sides, half the edge cell at the open xlo face, zero at the closed xhi and
// z faces, and (mu_v + mu_h |m s|^2) / (m_x m_y J dzeta) on interior k faces.
void expect_terrain_conductances (const Real sign)
{
    TerrainOptions options;
    options.slope_x = sign * Real(0.6);
    options.slope_y = -sign * Real(0.35);
    const auto g = make_terrain_grid(options);
    const auto f = make_fields(*g, 1);
    MappedDiffusionBoundary boundary;
    boundary.kind[0] = DiffusionFaceKind::Dirichlet;
    const Real nx = Real(options.cells[0]);
    const Real ny = Real(options.cells[1]);
    const Real nz = Real(options.cells[2]);
    const int x_last = options.cells[0];
    const int z_last = options.cells[2];
    const Real height = options.height;
    const Real slope_x = options.slope_x;
    const Real slope_y = options.slope_y;

    for (const auto& c : {host_scalar_coefficients(),
                          kinematic_coefficients(Real(0.37), Real(0.11)),
                          kinematic_coefficients(Real(0.37), Real(0.0))}) {
        MappedFaceTransmissibility t;
        t.define(g->ba, g->dm, 1, 0);
        std::string diagnostic;
        ASSERT_TRUE(BuildMappedDiffusionTransmissibility(t, g->geom, g->metrics(), c, f->rho, 0,
                                                         &f->eddy, boundary, diagnostic))
            << diagnostic;

        const bool eddy = c.mode.add_turbulence;
        const Real kh = c.horizontal.molecular_coeff;
        const Real kv = c.vertical.molecular_coeff;
        const int hcomp = c.horizontal.eddy_h_comp;
        const int vcomp = c.vertical.eddy_v_comp;
        MappedFaceTransmissibility expected;
        expected.define(g->ba, g->dm, 1, 0);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const MultiFab& num = dir == 0 ? g->mf_ux : dir == 1 ? g->mf_vy : g->mf_mx;
            const MultiFab& den = dir == 0 ? g->mf_uy : dir == 1 ? g->mf_vx : g->mf_my;
            for (amrex::MFIter mfi(expected.dir(dir)); mfi.isValid(); ++mfi) {
                const auto e = expected.dir(dir).array(mfi);
                const auto rho = f->rho.const_array(mfi);
                const auto mu = f->eddy.const_array(mfi);
                const auto m_num = num.const_array(mfi);
                const auto m_den = den.const_array(mfi);
                amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const int di = dir == 0 ? 1 : 0;
                    const int dj = dir == 1 ? 1 : 0;
                    const int dk = dir == 2 ? 1 : 0;
                    const Real rho_face =
                        Real(0.5) * (rho(i, j, k, 0) + rho(i - di, j - dj, k - dk, 0));
                    const Real mu_h = eddy
                        ? Real(0.5) * (mu(i, j, k, hcomp) + mu(i - di, j - dj, k - dk, hcomp))
                        : kh * rho_face;
                    const Real mu_v = eddy
                        ? Real(0.5) * (mu(i, j, k, vcomp) + mu(i - di, j - dj, k - dk, vcomp))
                        : kv * rho_face;
                    Real value = Real(0.0);
                    if (dir == 0) {
                        if (i < x_last) {
                            const Real spacing = i == 0 ? Real(0.5) / nx : Real(1.0) / nx;
                            value = height * m_num(i, j, 0) / m_den(i, j, 0) * mu_h / spacing;
                        }
                    } else if (dir == 1) {
                        value = height * m_num(i, j, 0) / m_den(i, j, 0) * mu_h * ny;
                    } else if (k > 0 && k < z_last) {
                        const Real sx = m_num(i, j, 0) * slope_x;
                        const Real sy = m_den(i, j, 0) * slope_y;
                        value = (mu_v + mu_h * (sx * sx + sy * sy)) /
                                (m_num(i, j, 0) * m_den(i, j, 0) * height / nz);
                    }
                    e(i, j, k, 0) = value;
                });
            }
            const Real scale = expected.dir(dir).norm0(0);
            ASSERT_GT(scale, Real(0.0));
            EXPECT_LE(max_difference(t.dir(dir), expected.dir(dir)), Real(64.0) * eps() * scale)
                << "dir " << dir;
        }
        // K_v = 0: the slope alone makes every interior k face exchange.
        if (!eddy && kv == Real(0.0)) {
            EXPECT_GT(value_at(t.dir(2), 0, IntVect(3, 2, 4)), Real(0.0));
        }
    }

    // On flat terrain K_v = 0 leaves no vertical exchange, and zero
    // coefficients none anywhere, on the slope too.
    TerrainOptions flat_options = options;
    flat_options.slope_x = Real(0.0);
    flat_options.slope_y = Real(0.0);
    const auto flat = make_terrain_grid(flat_options);
    const auto flat_t = build_transmissibility(*flat, *make_fields(*flat, 1),
                                               kinematic_coefficients(Real(0.37), Real(0.0)));
    EXPECT_EQ(flat_t->dir(2).norm0(0), Real(0.0));
    EXPECT_GT(flat_t->dir(0).norm0(0), Real(0.0));
    const auto zero = build_transmissibility(*g, *f, kinematic_coefficients(Real(0.0), Real(0.0)));
    EXPECT_EQ(max_abs_face(*zero, 0), Real(0.0));
}

void run_terrain_conductances ()
{
    expect_terrain_conductances(Real(1.0));
    expect_terrain_conductances(Real(-1.0));
}

// On affine terrain ERF's mapped transfer of a field linear in physical space
// is the physical flux -rho (K_h q_x, K_h q_y, K_v q_z) through each face;
// a k face's normal is (-m_x h_xi, -m_y h_eta, 1) / (m_x m_y). With K_v = 0
// the raw vertical flux vanishes, but horizontal diffusion along the slope
// still crosses the k faces, reversing with the slope. One pass giving every
// raw flux K_v loses that transfer.
void expect_horizontal_diffusion_crosses_k_faces (const Real sign)
{
    TerrainOptions options;
    options.periodic = {0, 0, 0};
    options.slope_x = sign * Real(0.6);
    options.slope_y = -sign * Real(0.35);
    options.variable_maps = false;
    options.map_x = Real(1.1);
    options.map_y = Real(0.9);
    const auto g = make_terrain_grid(options);
    const Real rho0 = Real(1.3);
    const auto f = make_uniform_terrain_fields(*g, rho0);
    const Real kh = Real(0.7);
    const Real alpha = Real(0.9);
    const Real beta = Real(-0.4);
    const Real gamma = Real(0.5);
    const Real mx = options.map_x;
    const Real my = options.map_y;
    const Real nx = Real(options.cells[0]);
    const Real ny = Real(options.cells[1]);
    const Real nz = Real(options.cells[2]);
    const Real height = options.height;
    const Real slope_x = options.slope_x;
    const Real slope_y = options.slope_y;
    // q = alpha x + beta y + gamma z with physical x = xi / m_x, y = eta / m_y.
    for (amrex::MFIter mfi(f->z); mfi.isValid(); ++mfi) {
        const auto q = f->z.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real xi = (Real(i) + Real(0.5)) / nx;
            const Real eta = (Real(j) + Real(0.5)) / ny;
            const Real z = height * (Real(k) + Real(0.5)) / nz + slope_x * xi + slope_y * eta;
            q(i, j, k, 0) = alpha * xi / mx + beta * eta / my + gamma * z;
        });
    }
    const Real along_slope = alpha * mx * slope_x + beta * my * slope_y;
    ASSERT_GT(std::abs(along_slope), Real(0.5));

    for (const Real kv : {Real(0.0), Real(0.2)}) {
        ScalarDiffusionFluxPolicy policy;
        policy.coefficients = {kh, EddyDiff::Scalar_h, EddyDiff::Scalar_v};
        policy.coefficient_mode = {true, false};
        const ScalarDiffusionCoefficients vertical{kv, EddyDiff::Scalar_h, EddyDiff::Scalar_v};
        const auto high = terrain_fluxes(*g, *f, policy, &vertical);

        const Real raw_z = -rho0 * kv * gamma;
        const Real mapped_z = (raw_z + rho0 * kh * along_slope) / (mx * my);
        const Real mapped_x = -rho0 * kh * alpha * height / my;
        const Real tol = Real(256.0) * eps() * rho0 * kh;
        if (kv == Real(0.0)) {
            EXPECT_EQ(high->raw_z.norm0(0), Real(0.0));
        }
        EXPECT_LE(max_deviation(high->raw_z, 0, raw_z), tol);
        EXPECT_LE(max_deviation(high->mapped.dir(2), 0, mapped_z), tol);
        EXPECT_LE(max_deviation(high->mapped.dir(0), 0, mapped_x), tol);

        ScalarDiffusionFluxPolicy vertical_only = policy;
        vertical_only.coefficients = vertical;
        const auto single = terrain_fluxes(*g, *f, vertical_only, nullptr);
        const Real single_z = (raw_z + rho0 * kv * along_slope) / (mx * my);
        EXPECT_LE(max_deviation(single->mapped.dir(2), 0, single_z), tol);
        EXPECT_GT(std::abs(mapped_z - single_z), Real(0.25));

        // The low-order vertical exchange keeps the slope term at K_v = 0.
        const auto t = build_transmissibility(*g, *f, kinematic_coefficients(kh, kv));
        const Real slope_sq = mx * mx * slope_x * slope_x + my * my * slope_y * slope_y;
        const Real tz = rho0 * (kv + kh * slope_sq) * nz / (mx * my * height);
        EXPECT_NEAR(value_at(t->dir(2), 0, IntVect(3, 2, 4)), tz, tol * nz);
        EXPECT_EQ(value_at(t->dir(2), 0, IntVect(3, 2, 0)), Real(0.0));
    }
}

void run_horizontal_diffusion_crosses_k_faces ()
{
    expect_horizontal_diffusion_crosses_k_faces(Real(1.0));
    expect_horizontal_diffusion_crosses_k_faces(Real(-1.0));
}

// The directional builder reproduces the one-coefficient builder bit for bit
// when both coefficient sets agree, and a distinct vertical set changes only
// the raw z flux.
void run_terrain_directional_builder ()
{
    const auto g = make_terrain_grid(TerrainOptions{});
    const auto f = make_fields(*g, 1, 2);

    ScalarDiffusionFluxPolicy policy;
    policy.coefficients = {Real(0.0), EddyDiff::Scalar_h, EddyDiff::Scalar_v};
    policy.coefficient_mode = {false, true};
    const auto shared = terrain_fluxes(*g, *f, policy, nullptr);
    const auto same = terrain_fluxes(*g, *f, policy, &policy.coefficients);
    EXPECT_EQ(max_difference(shared->raw_x, same->raw_x), Real(0.0));
    EXPECT_EQ(max_difference(shared->raw_y, same->raw_y), Real(0.0));
    EXPECT_EQ(max_difference(shared->raw_z, same->raw_z), Real(0.0));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        EXPECT_EQ(max_difference(shared->mapped.dir(dir), same->mapped.dir(dir)), Real(0.0));
    }

    const ScalarDiffusionCoefficients other{Real(0.0), EddyDiff::Scalar_h, EddyDiff::Theta_v};
    const auto split = terrain_fluxes(*g, *f, policy, &other);
    EXPECT_EQ(max_difference(shared->raw_x, split->raw_x), Real(0.0));
    EXPECT_EQ(max_difference(shared->raw_y, split->raw_y), Real(0.0));
    EXPECT_GT(max_difference(shared->raw_z, split->raw_z), Real(0.0));
}

// The design's affine witness: z = zeta / 4 + xi with equal computational
// spacing in x and zeta, a rise of four cells per cell, K_h = 1 and K_v = 0.
std::unique_ptr<DiffusionGrid> make_witness_grid (const Real slope)
{
    TerrainOptions options;
    options.cells = IntVect(8, 4, 8);
    options.periodic = {1, 1, 0};
    options.height = Real(0.25);
    options.slope_x = slope;
    options.slope_y = Real(0.0);
    options.variable_maps = false;
    return make_terrain_grid(options);
}

// A step the slope-blind combined bound admits multiplies the native terrain
// mode q = cos(pi k / 2) by less than -1; the terrain bound rejects that step
// and admits one that damps the mode.
void run_terrain_bound_includes_slope ()
{
    const Real slope = Real(1.0);
    const Real height = Real(0.25);
    const Real n = Real(8.0);
    const Real ny = Real(4.0);
    const auto sloped = make_witness_grid(slope);
    const auto flat = make_witness_grid(Real(0.0));
    const auto c = kinematic_coefficients(Real(1.0), Real(0.0));
    const auto f = make_uniform_terrain_fields(*sloped, Real(1.0));
    const auto flat_f = make_uniform_terrain_fields(*flat, Real(1.0));
    MappedFaceFluxRate no_carrier;
    no_carrier.define(sloped->ba, sloped->dm, 1, 0);
    no_carrier.setVal(Real(0.0));
    const Real rate = combined_rate(*sloped, *f, *build_transmissibility(*sloped, *f, c),
                                    no_carrier);
    const Real flat_rate = combined_rate(*flat, *flat_f,
                                         *build_transmissibility(*flat, *flat_f, c), no_carrier);
    const Real side_rate = Real(2.0) * n * n + Real(2.0) * ny * ny;
    const Real slope_rate = Real(2.0) * (slope / height) * (slope / height) * n * n;
    EXPECT_NEAR(flat_rate, side_rate, Real(64.0) * eps() * side_rate);
    EXPECT_NEAR(rate, side_rate + slope_rate, Real(64.0) * eps() * (side_rate + slope_rate));

    const double flat_step = 0.8 / static_cast<double>(flat_rate);
    double hard_limit = 0.0;
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(flat_step, flat_rate, hard_limit));
    EXPECT_TRUE(FixedDtExceedsMappedDonorLimit(flat_step, rate, hard_limit));

    const Real half_pi = Real(0.25) * two_pi;
    for (amrex::MFIter mfi(f->z); mfi.isValid(); ++mfi) {
        const auto q = f->z.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            q(i, j, k, 0) = std::cos(half_pi * Real(k));
        });
    }
    ScalarDiffusionFluxPolicy policy;
    policy.coefficients = c.horizontal;
    policy.coefficient_mode = c.mode;
    const auto high = terrain_fluxes(*sloped, *f, policy, &c.vertical);
    MultiFab tendency(sloped->ba, sloped->dm, 1, 0);
    ApplyMappedFluxTendency(high->mapped, 0, sloped->omega, 0, tendency, 0, sloped->dx_inv());
    // q = 1 two cells clear of the walls, where the stencil sees no boundary.
    const Real decay = -value_at(tendency, 0, IntVect(3, 1, 4));
    const Real expected_decay = (slope / height) * (slope / height) * n * n;
    EXPECT_NEAR(decay, expected_decay, Real(1.0e-10) * expected_decay);
    const Real flat_multiplier = Real(1.0) - Real(flat_step) * decay;
    const Real multiplier = Real(1.0) - decay / rate;
    amrex::Print() << "SBM_M6B_TERRAIN_BOUND slope_cells=" << slope / height
                   << " flat_rate=" << flat_rate << " terrain_rate=" << rate
                   << " flat_step_multiplier=" << flat_multiplier
                   << " terrain_step_multiplier=" << multiplier << '\n';
    EXPECT_LT(flat_multiplier, Real(-1.0));
    EXPECT_GE(multiplier, Real(-1.0));
    EXPECT_LE(multiplier, Real(1.0));
}

// The low-order terrain update of an isolated pulse stays nonnegative at the
// combined bound and goes negative just beyond it, with K_v = 0.
void run_terrain_low_order_bound ()
{
    const auto g = make_witness_grid(Real(1.0));
    const auto f = make_uniform_terrain_fields(*g, Real(1.0));
    const auto t = build_transmissibility(*g, *f, kinematic_coefficients(Real(1.0), Real(0.0)));
    MappedFaceFluxRate no_carrier;
    no_carrier.define(g->ba, g->dm, 1, 0);
    no_carrier.setVal(Real(0.0));
    const Real bound = combined_rate(*g, *f, *t, no_carrier);
    ASSERT_GT(bound, Real(0.0));
    const IntVect pulse(3, 1, 4);
    for (amrex::MFIter mfi(f->z); mfi.isValid(); ++mfi) {
        const auto q = f->z.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            q(i, j, k, 0) = (i == pulse[0] && j == pulse[1] && k == pulse[2]) ? Real(1.0)
                                                                             : Real(0.0);
        });
    }
    f->z.FillBoundary(g->geom.periodicity());
    const auto rate = build_rate(*g, *t, f->z, {0});

    auto at_bound = extensive_state(*g, *f);
    explicit_update(*g, *rate, *at_bound, Real(1.0) / bound);
    EXPECT_GE(at_bound->min(0), -Real(64.0) * eps());
    EXPECT_NEAR(value_at(*at_bound, 0, pulse), Real(0.0), Real(64.0) * eps());

    auto beyond = extensive_state(*g, *f);
    explicit_update(*g, *rate, *beyond, Real(1.01) / bound);
    EXPECT_LT(value_at(*beyond, 0, pulse), Real(0.0));
}

// Max interior error of ERF's native terrain diffusion, horizontal and
// vertical eddy coefficients mapped once, against the physical-space
// divergence d_x(mu_h q_x) + d_y(mu_h q_y) + d_z(mu_v q_z), on the curved
// terrain z = zeta + (1 - zeta) h (sin 2 pi xi + cos(2 pi eta) / 2) with map
// factors 1.1 in x and 0.9 in y. Cells two or more from a wall are measured.
Real curved_manufactured_error (const int n)
{
    TerrainOptions options;
    options.cells = IntVect(n);
    options.max_grid = n / 2;
    options.periodic = {1, 1, 0};
    options.variable_maps = false;
    options.map_x = Real(1.1);
    options.map_y = Real(0.9);
    auto g = make_terrain_grid(options);
    const Real hill = Real(0.06);
    const Real nn = Real(n);
    for (amrex::MFIter mfi(g->z_nd); mfi.isValid(); ++mfi) {
        const auto a = g->z_nd.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real xi = Real(i) / nn;
            const Real eta = Real(j) / nn;
            const Real zeta = Real(k) / nn;
            a(i, j, k, 0) = zeta + (Real(1.0) - zeta) * hill *
                                       (std::sin(two_pi * xi) + Real(0.5) * std::cos(two_pi * eta));
        });
    }
    make_zcc(g->geom, g->z_nd, g->z_cc);
    make_J(g->geom, g->z_nd, g->detj);
    MultiFab az(amrex::convert(g->ba, IntVect(0, 0, 1)), g->dm, 1, 0);
    make_areas(g->geom, g->z_nd, g->ax, g->ay, az);
    std::string diagnostic;
    EXPECT_TRUE(BuildMappedCellMeasure(g->omega, g->detj, g->mf_mx, g->mf_my, diagnostic))
        << diagnostic;

    // Physical x = xi / m_x and y = eta / m_y, so one period in each.
    const Real kx = two_pi * options.map_x;
    const Real ky = two_pi * options.map_y;
    const Real kz = Real(0.5) * two_pi;
    auto f = make_uniform_terrain_fields(*g, Real(1.0));
    for (amrex::MFIter mfi(f->z); mfi.isValid(); ++mfi) {
        const auto q = f->z.array(mfi);
        const auto mu = f->eddy.array(mfi);
        const auto zc = g->z_cc.const_array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real x = (Real(i) + Real(0.5)) / nn / Real(1.1);
            const Real y = (Real(j) + Real(0.5)) / nn / Real(0.9);
            const Real z = zc(i, j, k);
            q(i, j, k, 0) = std::sin(kx * x) * std::cos(ky * y) * std::cos(kz * z) +
                            Real(0.5) * z * z;
            mu(i, j, k, EddyDiff::Scalar_h) = Real(1.0) + Real(0.3) * std::cos(kx * x) +
                                              Real(0.2) * z;
            mu(i, j, k, EddyDiff::Scalar_v) = Real(0.5) + Real(0.2) * std::sin(ky * y) +
                                              Real(0.1) * z;
        });
    }
    ScalarDiffusionFluxPolicy policy;
    policy.coefficients = {Real(0.0), EddyDiff::Scalar_h, EddyDiff::Scalar_v};
    policy.coefficient_mode = {false, true};
    const auto high = terrain_fluxes(*g, *f, policy, &policy.coefficients);
    MultiFab tendency(g->ba, g->dm, 1, 0);
    ApplyMappedFluxTendency(high->mapped, 0, g->omega, 0, tendency, 0, g->dx_inv());

    MultiFab error(g->ba, g->dm, 1, 0);
    error.setVal(Real(0.0));
    for (amrex::MFIter mfi(error); mfi.isValid(); ++mfi) {
        const auto out = error.array(mfi);
        const auto got = tendency.const_array(mfi);
        const auto zc = g->z_cc.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (k < 2 || k > n - 3) return;
            const Real x = (Real(i) + Real(0.5)) / nn / Real(1.1);
            const Real y = (Real(j) + Real(0.5)) / nn / Real(0.9);
            const Real z = zc(i, j, k);
            const Real sx = std::sin(kx * x), cx = std::cos(kx * x);
            const Real sy = std::sin(ky * y), cy = std::cos(ky * y);
            const Real sz = std::sin(kz * z), cz = std::cos(kz * z);
            const Real wave = sx * cy * cz;
            const Real mu_h = Real(1.0) + Real(0.3) * cx + Real(0.2) * z;
            const Real mu_v = Real(0.5) + Real(0.2) * sy + Real(0.1) * z;
            const Real q_x = kx * cx * cy * cz;
            const Real q_z = -kz * sx * cy * sz + z;
            const Real d_x = -Real(0.3) * kx * sx * q_x - mu_h * kx * kx * wave;
            const Real d_y = -mu_h * ky * ky * wave;
            const Real d_z = Real(0.1) * q_z + mu_v * (-kz * kz * wave + Real(1.0));
            out(i, j, k, 0) = std::abs(got(i, j, k, 0) - (d_x + d_y + d_z));
        });
    }
    return error.norm0(0);
}

// The native terrain operator, the high-order candidate, converges at second
// order on the curved interior with distinct horizontal and vertical
// coefficients and nonunit maps: the design's gate is order 1.8 or better.
void run_curved_manufactured_order ()
{
    const Real e16 = curved_manufactured_error(16);
    const Real e32 = curved_manufactured_error(32);
    const Real e64 = curved_manufactured_error(64);
    const Real coarse = std::log2(e16 / e32);
    const Real fine = std::log2(e32 / e64);
    amrex::Print() << "SBM_M6B_CURVED_MANUFACTURED e16=" << e16 << " e32=" << e32 << " e64=" << e64
                   << " order=" << coarse << "," << fine << '\n';
    EXPECT_GE(fine, Real(1.8));
}

// Terrain conductances fail closed on missing or unusable metrics and on
// vertical faces they do not support.
void run_terrain_validation ()
{
    TerrainOptions options;
    const auto g = make_terrain_grid(options);
    const auto f = make_fields(*g, 1);
    const auto c = host_scalar_coefficients();
    std::string diagnostic;
    const auto build = [&] (const DiffusionGrid& grid, const DiffusionFields& fields,
                            const MappedDiffusionMetrics& metrics,
                            const MappedDiffusionBoundary& boundary) {
        MappedFaceTransmissibility out;
        out.define(grid.ba, grid.dm, 1, 0);
        return BuildMappedDiffusionTransmissibility(out, grid.geom, metrics, c, fields.rho, 0,
                                                    &fields.eddy, boundary, diagnostic);
    };
    EXPECT_TRUE(build(*g, *f, g->metrics(), {})) << diagnostic;

    auto no_heights = g->metrics();
    no_heights.z_nd = nullptr;
    EXPECT_FALSE(build(*g, *f, no_heights, {}));
    EXPECT_NE(diagnostic.find("nodal heights"), std::string::npos) << diagnostic;

    for (const int face : {2, 2 + AMREX_SPACEDIM}) {
        MappedDiffusionBoundary open;
        open.kind[static_cast<std::size_t>(face)] = DiffusionFaceKind::Dirichlet;
        EXPECT_FALSE(build(*g, *f, g->metrics(), open));
        EXPECT_NE(diagnostic.find("closed bottom and top"), std::string::npos) << diagnostic;
    }

    TerrainOptions periodic_options = options;
    periodic_options.periodic = {0, 1, 1};
    const auto periodic = make_terrain_grid(periodic_options);
    EXPECT_FALSE(build(*periodic, *make_fields(*periodic, 1), periodic->metrics(), {}));
    EXPECT_NE(diagnostic.find("closed bottom and top"), std::string::npos) << diagnostic;

    // One nonfinite side area fails on every rank.
    const auto set_ax = [&] (const Real value) {
        for (amrex::MFIter mfi(g->ax); mfi.isValid(); ++mfi) {
            const auto a = g->ax.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (i == 2 && j == 1 && k == 3) { a(i, j, k, 0) = value; }
            });
        }
    };
    set_ax(std::numeric_limits<Real>::quiet_NaN());
    EXPECT_FALSE(build(*g, *f, g->metrics(), {}));
    EXPECT_NE(diagnostic.find("(1 invalid samples)"), std::string::npos) << diagnostic;
    set_ax(options.height);
    EXPECT_TRUE(build(*g, *f, g->metrics(), {})) << diagnostic;

    // An inverted column has J < 0 on every interior k face.
    g->z_nd.mult(Real(-1.0), 0, 1, g->z_nd.nGrow());
    EXPECT_FALSE(build(*g, *f, g->metrics(), {}));
    EXPECT_NE(diagnostic.find("invalid samples"), std::string::npos) << diagnostic;
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

    // The transport refuses a terrain mesh too until its stage operator is qualified.
    const auto layout = make_two_moment_layout();
    erf_sbm::SBMTransport transport(layout, 1);
    EXPECT_THROW(transport.set_particle_mixing(
                     policy_config(erf_sbm::ParticleMixingPolicy::HostPassiveScalar),
                     erf_auxiliary::DiffusionMesh::VariableDz),
                 std::invalid_argument);
    EXPECT_NO_THROW(transport.set_particle_mixing(
        policy_config(erf_sbm::ParticleMixingPolicy::None),
        erf_auxiliary::DiffusionMesh::VariableDz));

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

// The terrain matching rule: theta keeps theta * low between zero and high,
// without forming high / low where it could overflow, and fails on NaN.
void run_matched_fraction ()
{
    using erf_sbm::matched_diffusion_fraction;
    EXPECT_EQ(matched_diffusion_fraction(Real(0.0), Real(0.0)), Real(1.0));
    EXPECT_EQ(matched_diffusion_fraction(Real(3.0), Real(0.0)), Real(1.0));
    EXPECT_EQ(matched_diffusion_fraction(Real(0.0), Real(2.0)), Real(0.0));
    EXPECT_EQ(matched_diffusion_fraction(Real(-1.0), Real(2.0)), Real(0.0));
    EXPECT_EQ(matched_diffusion_fraction(Real(1.0), Real(-2.0)), Real(0.0));
    EXPECT_EQ(matched_diffusion_fraction(Real(1.0), Real(2.0)), Real(0.5));
    EXPECT_EQ(matched_diffusion_fraction(Real(-1.0), Real(-4.0)), Real(0.25));
    EXPECT_EQ(matched_diffusion_fraction(Real(5.0), Real(2.0)), Real(1.0));
    EXPECT_EQ(matched_diffusion_fraction(Real(2.0), Real(2.0)), Real(1.0));
    const Real huge = std::numeric_limits<Real>::max();
    const Real tiny = std::numeric_limits<Real>::denorm_min();
    EXPECT_EQ(matched_diffusion_fraction(huge, tiny), Real(1.0));
    const Real small = matched_diffusion_fraction(tiny, huge);
    EXPECT_GE(small, Real(0.0));
    EXPECT_LT(small, eps());
    for (const Real bad : {std::numeric_limits<Real>::quiet_NaN(),
                           std::numeric_limits<Real>::infinity()}) {
        EXPECT_TRUE(std::isnan(matched_diffusion_fraction(bad, Real(1.0))));
        EXPECT_TRUE(std::isnan(matched_diffusion_fraction(Real(1.0), bad)));
    }
}

// The restart schema pins the terrain operator and boundary policy for a
// terrain mesh only, leaving every other mesh's schema unchanged.
void run_terrain_schema ()
{
    erf_sbm::SBMPhysicalBoundaryConfig boundaries;
    boundaries.kind[2] = erf_sbm::AdvectionFaceKind::ImpermeableWall;
    boundaries.kind[5] = erf_sbm::AdvectionFaceKind::ImpermeableWall;
    const auto mixing = policy_config(erf_sbm::ParticleMixingPolicy::HostPassiveScalar);
    const std::string base = erf_sbm::particle_mixing_schema(mixing, boundaries);
    for (const auto mesh : {erf_auxiliary::DiffusionMesh::ConstantDz,
                            erf_auxiliary::DiffusionMesh::StretchedDz}) {
        EXPECT_EQ(erf_sbm::particle_mixing_schema(mixing, boundaries, mesh), base);
    }
    const std::string terrain = erf_sbm::particle_mixing_schema(
        mixing, boundaries, erf_auxiliary::DiffusionMesh::VariableDz);
    EXPECT_EQ(terrain.rfind(base, 0), 0U);
    EXPECT_NE(terrain.find("\nterrain-diffusion="), std::string::npos) << terrain;
    EXPECT_TRUE(erf_sbm::particle_mixing_schema(policy_config(erf_sbm::ParticleMixingPolicy::None),
                                                boundaries,
                                                erf_auxiliary::DiffusionMesh::VariableDz)
                    .empty());
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

erf_sbm::ParticleMixingConfig parsed_config (const std::string& name,
                                             const Real kh = Real(0.0),
                                             const Real kv = Real(0.0))
{
    erf_sbm::ParticleMixingConfig config;
    std::string reason;
    EXPECT_TRUE(erf_sbm::make_particle_mixing_config(name, kh, kv, config, reason)) << reason;
    return config;
}

TurbChoice les (const LESType type)
{
    TurbChoice turbulence;
    turbulence.les_type = type;
    return turbulence;
}

void run_config_parsing ()
{
    const auto reference = parsed_config("ConstantReference", Real(0.25), Real(0.125));
    EXPECT_EQ(reference.policy, erf_sbm::ParticleMixingPolicy::ConstantReference);
    EXPECT_EQ(reference.reference_kh, Real(0.25));
    EXPECT_EQ(reference.reference_kv, Real(0.125));
    EXPECT_EQ(parsed_config("HostPassiveScalar").policy,
              erf_sbm::ParticleMixingPolicy::HostPassiveScalar);

    erf_sbm::ParticleMixingConfig config;
    std::string reason;
    EXPECT_FALSE(erf_sbm::make_particle_mixing_config("Eddy", Real(0.0), Real(0.0), config, reason));
    EXPECT_NE(reason.find("not 'Eddy'"), std::string::npos) << reason;
    // Coefficients given to any other policy are a configuration error, not ignored.
    for (const char* name : {"None", "HostPassiveScalar"}) {
        EXPECT_FALSE(erf_sbm::make_particle_mixing_config(name, Real(0.0), Real(1.0e-3),
                                                          config, reason)) << name;
        EXPECT_NE(reason.find("only to ConstantReference"), std::string::npos) << reason;
    }
    EXPECT_FALSE(erf_sbm::make_particle_mixing_config("ConstantReference", Real(-1.0), Real(0.0),
                                                      config, reason));
    EXPECT_NE(reason.find("finite and nonnegative"), std::string::npos) << reason;
}

void run_host_validation ()
{
    const auto none = parsed_config("None");
    const auto host = parsed_config("HostPassiveScalar");
    const auto reference = parsed_config("ConstantReference", Real(0.1), Real(0.05));
    std::string reason;
    const auto accepts = [&] (const erf_sbm::ParticleMixingConfig& config,
                              const std::vector<TurbChoice>& turbulence,
                              const MeshType mesh = MeshType::ConstantDz,
                              const MolecDiffType molecular = MolecDiffType::None) {
        return erf_sbm::validate_particle_mixing_host(config, turbulence, molecular, mesh,
                                                      reason);
    };
    const auto rejects = [&] (const bool accepted, const std::string& fragment) {
        EXPECT_FALSE(accepted) << fragment;
        EXPECT_NE(reason.find(fragment), std::string::npos) << reason;
    };
    const std::vector<TurbChoice> quiet{les(LESType::None)};
    const std::vector<TurbChoice> smagorinsky{les(LESType::Smagorinsky)};
    const std::vector<TurbChoice> deardorff{les(LESType::Deardorff)};

    EXPECT_TRUE(accepts(none, quiet)) << reason;
    // None adds nothing, so it constrains neither the mesh nor the levels.
    EXPECT_TRUE(accepts(none, {les(LESType::None), les(LESType::None)},
                        MeshType::VariableDz)) << reason;
    EXPECT_TRUE(accepts(host, smagorinsky)) << reason;
    EXPECT_TRUE(accepts(host, deardorff, MeshType::StretchedDz)) << reason;
    EXPECT_TRUE(accepts(reference, quiet, MeshType::StretchedDz)) << reason;

    rejects(accepts(none, smagorinsky),
            "SBM particle mixing None does not mix the spectrum with ERF's LES closure");
    rejects(accepts(reference, deardorff),
            "SBM particle mixing ConstantReference does not mix the spectrum");
    rejects(accepts(host, quiet), "requires erf.les_type = Smagorinsky or Deardorff");
    auto smag2d = les(LESType::Smagorinsky);
    smag2d.smag2d = true;
    rejects(accepts(host, {smag2d}), "does not support erf.les_type = Smagorinsky2D (level 0)");
    rejects(accepts(host, {les(LESType::Smagorinsky), les(LESType::None)}),
            "Smagorinsky or Deardorff on every level (level 1)");
    rejects(accepts(host, {les(LESType::Smagorinsky), les(LESType::Deardorff)}),
            "the same erf.les_type on every level");

    auto pbl = les(LESType::None);
    pbl.pbl_type = PBLType::MYNN25;
    for (const auto& config : {none, host, reference}) {
        rejects(accepts(config, {pbl}), "SBM does not support a PBL closure (erf.pbl_type = MYNN25 (level 0))");
        rejects(accepts(config, quiet, MeshType::ConstantDz, MolecDiffType::Constant),
                "SBM does not support molecular diffusion");
    }
    auto rans = les(LESType::None);
    rans.rans_type = RANSType::kEqn;
    rejects(accepts(none, {rans}), "SBM does not support a RANS closure");

    rejects(accepts(host, smagorinsky, MeshType::VariableDz), "terrain-fitted (VariableDz)");
    rejects(accepts(reference, quiet, MeshType::VariableDz), "terrain-fitted (VariableDz)");
    // Two levels: SBM's own AMR envelope bounds the levels, not the mixing policy.
    EXPECT_TRUE(accepts(host, {les(LESType::Smagorinsky), les(LESType::Smagorinsky)})) << reason;
    EXPECT_TRUE(accepts(reference, {les(LESType::None), les(LESType::None)})) << reason;

    // Terrain qualification admits a terrain mesh on one level only.
    const auto qualifies = [&] (const erf_sbm::ParticleMixingConfig& config,
                                const std::vector<TurbChoice>& turbulence) {
        return erf_sbm::validate_particle_mixing_host(config, turbulence, MolecDiffType::None,
                                                      MeshType::VariableDz, true, reason);
    };
    EXPECT_TRUE(qualifies(host, smagorinsky)) << reason;
    EXPECT_TRUE(qualifies(reference, quiet)) << reason;
    rejects(qualifies(host, {les(LESType::Smagorinsky), les(LESType::Smagorinsky)}),
            "qualified on one level only");
    erf_auxiliary::DiffusionMesh mesh = erf_auxiliary::DiffusionMesh::ConstantDz;
    EXPECT_TRUE(erf_sbm::particle_mixing_mesh(MeshType::VariableDz, mesh, reason, true));
    EXPECT_EQ(mesh, erf_auxiliary::DiffusionMesh::VariableDz);
}

// The terrain envelope: every cell drops at most one thickness across, by
// ERF's unclamped cell drops; the worst cell is named otherwise.
void run_terrain_admission ()
{
    const Box domain(IntVect(0), IntVect(7, 3, 5));
    BoxArray ba(domain);
    ba.maxSize(4);
    const DistributionMapping dm(ba);
    MultiFab z_nd(amrex::convert(ba, IntVect(1)), dm, 1, 0);
    // Each cell is 1/6 thick; a plane rising `slope` per cell in x, plus a
    // bump of `bump` on the nodes i = 5, 6 of j = 3.
    const auto fill = [&] (const Real slope, const Real bump) {
        for (amrex::MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
            const auto a = z_nd.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                a(i, j, k, 0) = Real(k) / Real(6) + slope * Real(i) +
                                ((i == 5 || i == 6) && j == 3 ? bump : Real(0.0));
            });
        }
    };
    std::string reason;
    fill(Real(1.0) / Real(6.0), Real(0.0));
    EXPECT_TRUE(erf_sbm::admit_terrain_mixing_geometry(z_nd, ba, Real(1.0) + Real(8.0) * eps(),
                                                       reason)) << reason;
    fill(Real(0.1), Real(0.4));
    EXPECT_FALSE(erf_sbm::admit_terrain_mixing_geometry(z_nd, ba, Real(1.0), reason));
    // Cells (5,2,k) are worst: all four y edges rise 0.4 and the x edges 0.1,
    // so sqrt(0.1^2 + 0.4^2) * 6 = 2.474. Roundoff in the thickness picks the k.
    EXPECT_NE(reason.find("cell (5,2,"), std::string::npos) << reason;
    EXPECT_NE(reason.find("drops 2.47386"), std::string::npos) << reason;
}

void run_boundary_validation ()
{
    using erf_auxiliary::DiffusionFaceKind;
    using erf_sbm::AdvectionFaceKind;
    const auto none = parsed_config("None");
    const auto host = parsed_config("HostPassiveScalar");
    std::string reason;
    erf_sbm::SBMPhysicalBoundaryConfig faces;
    EXPECT_TRUE(erf_sbm::validate_particle_mixing_boundaries(host, faces, reason)) << reason;
    faces.kind[2] = AdvectionFaceKind::ImpermeableWall;
    faces.kind[5] = AdvectionFaceKind::ImpermeableWall;
    EXPECT_TRUE(erf_sbm::validate_particle_mixing_boundaries(host, faces, reason)) << reason;

    // Outflow and inflow faces, open or closed to mixing, on x.
    auto open = faces;
    open.kind[0] = AdvectionFaceKind::PrescribedInflow;
    open.kind[3] = AdvectionFaceKind::AdvectiveOutflow;
    for (const auto diffusion : {DiffusionFaceKind::Closed, DiffusionFaceKind::Dirichlet,
                                 DiffusionFaceKind::DirichletWhereInward}) {
        open.diffusion[0] = diffusion;
        EXPECT_TRUE(erf_sbm::validate_particle_mixing_boundaries(none, open, reason)) << reason;
        EXPECT_TRUE(erf_sbm::validate_particle_mixing_boundaries(host, open, reason)) << reason;
    }

    // Only an inflow face holds a spectrum to exchange with.
    auto outflow_open = open;
    outflow_open.diffusion[3] = DiffusionFaceKind::Dirichlet;
    EXPECT_FALSE(erf_sbm::validate_particle_mixing_boundaries(none, outflow_open, reason));
    EXPECT_NE(reason.find("only through an inflow face, not the outflow face xhi"),
              std::string::npos) << reason;

    // ERF's vertical scalar diffusion has no inward-only rule.
    auto vertical = faces;
    vertical.kind[2] = AdvectionFaceKind::PrescribedInflow;
    vertical.diffusion[2] = DiffusionFaceKind::DirichletWhereInward;
    EXPECT_TRUE(erf_sbm::validate_particle_mixing_boundaries(none, vertical, reason)) << reason;
    EXPECT_FALSE(erf_sbm::validate_particle_mixing_boundaries(host, vertical, reason));
    EXPECT_NE(reason.find("does not support an inflow_outflow zlo face"), std::string::npos)
        << reason;
    vertical.diffusion[2] = DiffusionFaceKind::Dirichlet;
    EXPECT_TRUE(erf_sbm::validate_particle_mixing_boundaries(host, vertical, reason)) << reason;

    // ERF's boundary types map onto the faces' mixing semantics.
    EXPECT_EQ(erf_sbm::map_erf_diffusion(ERF_BC::inflow), DiffusionFaceKind::Dirichlet);
    EXPECT_EQ(erf_sbm::map_erf_diffusion(ERF_BC::inflow_outflow),
              DiffusionFaceKind::DirichletWhereInward);
    for (const auto bc : {ERF_BC::outflow, ERF_BC::ho_outflow, ERF_BC::slip_wall,
                          ERF_BC::no_slip_wall, ERF_BC::surface_layer, ERF_BC::periodic}) {
        EXPECT_EQ(erf_sbm::map_erf_diffusion(bc), DiffusionFaceKind::Closed);
    }
}

} // namespace

TEST(AuxiliaryDiffusion, NGridTransmissibilityMatchesHostDiffusion) { run_ngrid_parity(); }
TEST(AuxiliaryDiffusion, StretchedGridTransmissibilityMatchesHostDiffusion) { run_stretched_parity(); }
TEST(AuxiliaryDiffusion, PhysicalFacesCarryNoDiffusiveTransfer) { run_physical_faces(); }
TEST(AuxiliaryDiffusion, OpenFacesExchangeOverHalfTheEdgeCell) { run_open_faces(); }
TEST(AuxiliaryDiffusion, OpenFacesValidateTheirGhostsAndCarrier) { run_open_face_validation(); }
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
TEST(AuxiliaryDiffusion, TerrainConductancesMatchAffineSlopeOracle) { run_terrain_conductances(); }
TEST(AuxiliaryDiffusion, TerrainHorizontalDiffusionCrossesKFaces)
{
    run_horizontal_diffusion_crosses_k_faces();
}
TEST(AuxiliaryDiffusion, TerrainDirectionalBuilderForwardsSharedCoefficients)
{
    run_terrain_directional_builder();
}
TEST(AuxiliaryDiffusion, TerrainBoundRejectsStepSlopeBlindBoundAdmits)
{
    run_terrain_bound_includes_slope();
}
TEST(AuxiliaryDiffusion, TerrainLowOrderPulseStaysNonnegativeUpToTheBound)
{
    run_terrain_low_order_bound();
}
TEST(AuxiliaryDiffusion, TerrainMetricsFailClosed) { run_terrain_validation(); }
TEST(AuxiliaryDiffusion, TerrainNativeOperatorConvergesOnCurvedTerrain)
{
    run_curved_manufactured_order();
}

TEST(SBMParticleMixing, PolicyNamesAndIdentity) { run_policy_names(); }
TEST(SBMParticleMixing, HostPassiveScalarUsesScalarChannelsAsGiven) { run_host_channels(); }
TEST(SBMParticleMixing, HorizontalAndVerticalChannelsStaySeparate) { run_channels_stay_separate(); }
TEST(SBMParticleMixing, ConstantReferenceAndNoInheritedMolecularDiffusion) { run_constant_reference(); }
TEST(SBMParticleMixing, InvalidCoefficientsFailClosed) { run_invalid_coefficients(); }
TEST(SBMParticleMixing, TerrainAndPeriodicStretchedGeometryAreRejected) { run_geometry_rejections(); }
TEST(SBMParticleMixing, NonePolicyHasNoTransmissibility) { run_none_policy(); }
TEST(SBMParticleMixing, TerrainMatchingFractionStaysBetweenZeroAndTheNativeRate)
{
    run_matched_fraction();
}
TEST(SBMParticleMixing, RestartSchemaPinsTheTerrainOperator) { run_terrain_schema(); }
TEST(SBMParticleMixing, TerrainGeometryAdmissionNamesTheWorstCell) { run_terrain_admission(); }
TEST(SBMParticleMixing, SharedCoefficientPreservesEveryBinGroup) { run_group_cone(); }
TEST(SBMParticleMixing, ConfigurationParsingFailsClosed) { run_config_parsing(); }
TEST(SBMParticleMixing, HostClosureAndMeshValidation) { run_host_validation(); }
TEST(SBMParticleMixing, MixingFacesFollowERFScalarDiffusion) { run_boundary_validation(); }
