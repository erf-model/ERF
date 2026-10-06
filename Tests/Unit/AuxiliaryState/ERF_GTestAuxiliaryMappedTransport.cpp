#include <AMReX_BCRec.H>

#include "ERF_AuxiliaryInertTracer.H"
#include "ERF_AuxiliaryMappedTransport.H"
#include "ERF_AuxiliaryStage.H"
#include "ERF_ScalarDiffusion.H"
#include "ERF_AdvectionSrcForScalars.H"
#include "ERF_IndexDefines.H"
#include "ERF_TerrainMetrics.H"

#include <AMReX_Gpu.H>
#include <AMReX_Math.H>

#include <gtest/gtest.h>

#include <cmath>
#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>

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

static_assert(!std::is_same_v<MappedFaceFluxRate, IntegratedMappedFaceFlux>);

Geometry make_geometry (const Box& domain, const bool periodic = true)
{
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int is_periodic[AMREX_SPACEDIM] = {
        periodic ? 1 : 0, periodic ? 1 : 0, periodic ? 1 : 0};
    return Geometry(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
}

BoxArray project_to_xy (const BoxArray& cell_ba)
{
    amrex::BoxList boxes = cell_ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    return BoxArray(std::move(boxes));
}

struct TestGrid {
    Box domain;
    Geometry geom;
    BoxArray ba;
    BoxArray map_ba;
    DistributionMapping dm;
    MultiFab detj;
    MultiFab mx;
    MultiFab my;
    MultiFab omega;

    explicit TestGrid (const int nx = 4, const int ny = 3, const int nz = 2)
        : domain(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1)),
          geom(make_geometry(domain)),
          ba(domain), map_ba(project_to_xy(ba)), dm(ba),
          detj(ba, dm, 1, 0), mx(map_ba, dm, 1, 0), my(map_ba, dm, 1, 0),
          omega(ba, dm, 1, 0)
    {
        fill_metrics();
    }

    void fill_metrics ()
    {
        for (amrex::MFIter mfi(detj); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto d = detj.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                d(i, j, k, 0) = Real(1.31) + Real(0.037) * i + Real(0.019) * j + Real(0.011) * k;
            });
        }
        for (amrex::MFIter mfi(mx); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto x = mx.array(mfi);
            const auto y = my.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
                x(i, j, 0, 0) = Real(1.17) + Real(0.023) * i;
                y(i, j, 0, 0) = Real(0.83) + Real(0.017) * j + Real(0.006) * i;
            });
        }
    }

    bool build_measure ()
    {
        std::string diagnostic;
        const bool ok = BuildMappedCellMeasure(omega, detj, mx, my, diagnostic);
        EXPECT_TRUE(ok) << diagnostic;
        return ok;
    }
};

void fill_rate (MappedFaceFluxRate& rate, const Box& domain,
                const int comp, const bool periodic = false,
                const Real factor = Real(1.0))
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto& face = rate.dir(dir);
        for (amrex::MFIter mfi(face); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto out = face.array(mfi);
            const int ncomp = face.nComp();
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                Real ii = Real(i);
                Real jj = Real(j);
                Real kk = Real(k);
                if (periodic) {
                    if (dir == 0 && i == domain.bigEnd(0) + 1) ii = Real(domain.smallEnd(0));
                    if (dir == 1 && j == domain.bigEnd(1) + 1) jj = Real(domain.smallEnd(1));
                    if (dir == 2 && k == domain.bigEnd(2) + 1) kk = Real(domain.smallEnd(2));
                }
                const Real value = Real(0.9) + Real(0.31) * ii - Real(0.27) * jj +
                    Real(0.19) * kk + Real(0.043) * ii * jj +
                    Real(0.017) * jj * kk + Real(0.011) * dir;
                out(i, j, k, n) = n == comp ? factor * value : Real(91.0) + Real(n);
            });
        }
    }
}

void fill_constant_rate (MappedFaceFluxRate& rate, const Real value)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto& face = rate.dir(dir);
        for (amrex::MFIter mfi(face); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto out = face.array(mfi);
            const int ncomp = face.nComp();
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                out(i, j, k, n) = n == 0 ? value : Real(41.0) + Real(n);
            });
        }
    }
}

void fill_componentwise_constant_rate (MappedFaceFluxRate& rate,
                                       const std::array<Real, 3>& values)
{
    AMREX_ALWAYS_ASSERT(rate.nComp() == static_cast<int>(values.size()));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto& face = rate.dir(dir);
        for (int comp = 0; comp < rate.nComp(); ++comp) {
            face.setVal(values[static_cast<std::size_t>(comp)], comp, 1, 0);
        }
    }
}

Real max_face_component_error (const IntegratedMappedFaceFlux& flux,
                               const int component,
                               const Real expected)
{
    Real maximum = Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto& face = flux.dir(dir);
        MultiFab difference(face.boxArray(), face.DistributionMap(), 1, 0);
        for (amrex::MFIter mfi(face); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto in = face.const_array(mfi);
            const auto out = difference.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                out(i, j, k, 0) = amrex::Math::abs(in(i, j, k, component) - expected);
            });
        }
        maximum = amrex::max(maximum, difference.norm0(0));
    }
    return maximum;
}

void fill_diffusion_raw (MultiFab& face, const int dir, const int comp)
{
    for (amrex::MFIter mfi(face); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto out = face.array(mfi);
        const int ncomp = face.nComp();
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            const Real value = Real(0.42) + Real(0.14) * i - Real(0.09) * j +
                Real(0.12) * k + Real(0.025) * i * j + Real(0.013) * dir;
            out(i, j, k, n) = n == comp ? value : Real(300.0) + Real(10 * dir + n);
        });
    }
}

Real max_component_difference (const MultiFab& a, const int acomp,
                               const MultiFab& b, const int bcomp)
{
    MultiFab difference(a.boxArray(), a.DistributionMap(), 1, 0);
    for (amrex::MFIter mfi(a); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto av = a.const_array(mfi);
        const auto bv = b.const_array(mfi);
        const auto out = difference.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = amrex::Math::abs(av(i, j, k, acomp) - bv(i, j, k, bcomp));
        });
    }
    return difference.norm0(0);
}

Real mapped_diffusion_parity (const char kind,
                              Real* raw_vertical_difference = nullptr,
                              Real* omitted_cross_difference = nullptr)
{
    TestGrid g;
    g.build_measure();
    constexpr int raw_comp = 1;
    constexpr int mapped_comp = 2;
    constexpr int rhs_comp = 1;
    MultiFab raw_x(amrex::convert(g.ba, IntVect::TheDimensionVector(0)), g.dm, 3, 1);
    MultiFab raw_y(amrex::convert(g.ba, IntVect::TheDimensionVector(1)), g.dm, 3, 1);
    MultiFab raw_z(amrex::convert(g.ba, IntVect::TheDimensionVector(2)), g.dm, 3, 1);
    fill_diffusion_raw(raw_x, 0, raw_comp);
    fill_diffusion_raw(raw_y, 1, raw_comp);
    fill_diffusion_raw(raw_z, 2, raw_comp);

    MappedFaceFluxRate mapped;
    mapped.define(g.ba, g.dm, 3, 0);
    mapped.setVal(Real(-71.0));
    MultiFab ax(amrex::convert(g.ba, IntVect::TheDimensionVector(0)), g.dm, 1, 0);
    MultiFab ay(amrex::convert(g.ba, IntVect::TheDimensionVector(1)), g.dm, 1, 0);
    MultiFab mf_uy(amrex::convert(g.ba, IntVect::TheDimensionVector(0)), g.dm, 1, 0);
    MultiFab mf_vx(amrex::convert(g.ba, IntVect::TheDimensionVector(1)), g.dm, 1, 0);
    MultiFab z_nd(amrex::convert(g.ba, IntVect::TheNodeVector()), g.dm, 1, 0);
    ax.setVal(Real(0.79)); ay.setVal(Real(0.86));
    mf_uy.setVal(Real(1.43)); mf_vx.setVal(Real(1.28));
    for (amrex::MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto z = z_nd.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            z(i, j, k, 0) = Real(0.37) * i - Real(0.23) * j + Real(0.19) * k +
                             Real(0.014) * i * j;
        });
    }

    const auto inv = g.geom.InvCellSizeArray();
    for (amrex::MFIter mfi(g.omega); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        if (kind == 'N') {
            BuildScalarDiffusionMappedTransfers_N(
                bx, raw_x.const_array(mfi), raw_comp, raw_y.const_array(mfi), raw_comp,
                raw_z.const_array(mfi), raw_comp,
                mapped.dir(0).array(mfi), mapped_comp,
                mapped.dir(1).array(mfi), mapped_comp,
                mapped.dir(2).array(mfi), mapped_comp,
                g.mx.const_array(mfi), g.my.const_array(mfi));
        } else if (kind == 'S') {
            BuildScalarDiffusionMappedTransfers_S(
                bx, raw_x.const_array(mfi), raw_comp, raw_y.const_array(mfi), raw_comp,
                raw_z.const_array(mfi), raw_comp,
                mapped.dir(0).array(mfi), mapped_comp,
                mapped.dir(1).array(mfi), mapped_comp,
                mapped.dir(2).array(mfi), mapped_comp,
                ax.const_array(mfi), ay.const_array(mfi),
                g.mx.const_array(mfi), g.my.const_array(mfi));
        } else {
            BuildScalarDiffusionMappedTransfers_T(
                bx, g.domain, raw_x.const_array(mfi), raw_y.const_array(mfi),
                raw_z.const_array(mfi), raw_comp,
                mapped.dir(0).array(mfi), mapped_comp,
                mapped.dir(1).array(mfi), mapped_comp,
                mapped.dir(2).array(mfi), mapped_comp,
                z_nd.const_array(mfi), ax.const_array(mfi), ay.const_array(mfi), inv,
                g.mx.const_array(mfi), mf_uy.const_array(mfi), g.my.const_array(mfi),
                mf_vx.const_array(mfi), false);
        }
    }

    MultiFab native(g.ba, g.dm, 3, 0);
    MultiFab generic(g.ba, g.dm, 3, 0);
    native.setVal(Real(0.0)); generic.setVal(Real(0.0));
    for (amrex::MFIter mfi(native); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        ApplyScalarMappedFluxDivergence(
            bx, mapped.dir(0).const_array(mfi), mapped_comp,
            mapped.dir(1).const_array(mfi), mapped_comp,
            mapped.dir(2).const_array(mfi), mapped_comp,
            native.array(mfi), rhs_comp, g.detj.const_array(mfi), inv,
            g.mx.const_array(mfi), g.my.const_array(mfi));
    }
    ApplyMappedFluxTendency(mapped, mapped_comp, g.omega, 0, generic, rhs_comp, inv);

    if (kind == 'T') {
        MappedFaceFluxRate raw_fz;
        MappedFaceFluxRate no_cross;
        raw_fz.define(g.ba, g.dm, 3, 0);
        no_cross.define(g.ba, g.dm, 3, 0);
        raw_fz.setVal(Real(0.0));
        no_cross.setVal(Real(0.0));
        for (int dir = 0; dir < 2; ++dir) {
            MultiFab::Copy(raw_fz.dir(dir), mapped.dir(dir), mapped_comp, mapped_comp, 1, 0);
            MultiFab::Copy(no_cross.dir(dir), mapped.dir(dir), mapped_comp, mapped_comp, 1, 0);
        }
        for (amrex::MFIter mfi(g.omega); mfi.isValid(); ++mfi) {
            const Box zbx = amrex::surroundingNodes(mfi.tilebox(), 2);
            const auto rawz = raw_z.const_array(mfi);
            const auto mx = g.mx.const_array(mfi);
            const auto my = g.my.const_array(mfi);
            const auto wrong_raw = raw_fz.dir(2).array(mfi);
            const auto wrong_cross = no_cross.dir(2).array(mfi);
            amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                wrong_raw(i, j, k, mapped_comp) = rawz(i, j, k, raw_comp);
                wrong_cross(i, j, k, mapped_comp) = rawz(i, j, k, raw_comp) /
                    (mx(i, j, 0) * my(i, j, 0));
            });
        }
        MultiFab wrong_raw_rhs(g.ba, g.dm, 3, 0), wrong_cross_rhs(g.ba, g.dm, 3, 0);
        wrong_raw_rhs.setVal(Real(0.0)); wrong_cross_rhs.setVal(Real(0.0));
        ApplyMappedFluxTendency(raw_fz, mapped_comp, g.omega, 0,
                                wrong_raw_rhs, rhs_comp, inv);
        ApplyMappedFluxTendency(no_cross, mapped_comp, g.omega, 0,
                                wrong_cross_rhs, rhs_comp, inv);
        if (raw_vertical_difference) {
            *raw_vertical_difference = max_component_difference(native, rhs_comp,
                                                                 wrong_raw_rhs, rhs_comp);
        }
        if (omitted_cross_difference) {
            *omitted_cross_difference = max_component_difference(native, rhs_comp,
                                                                  wrong_cross_rhs, rhs_comp);
        }
    }
    return max_component_difference(native, rhs_comp, generic, rhs_comp);
}

void run_auxiliary_mapped_transport_MappedMeasureIsExplicitAndRejectsInvalidMetrics ()
{
    TestGrid g;
    std::string diagnostic;
    ASSERT_TRUE(BuildMappedCellMeasure(g.omega, g.detj, g.mx, g.my, diagnostic)) << diagnostic;
    MultiFab expected(g.ba, g.dm, 1, 0);
    for (amrex::MFIter mfi(expected); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto d = g.detj.const_array(mfi);
        const auto x = g.mx.const_array(mfi);
        const auto y = g.my.const_array(mfi);
        const auto e = expected.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            e(i, j, k, 0) = d(i, j, k, 0) / (x(i, j, 0) * y(i, j, 0));
        });
    }
    EXPECT_LT(max_component_difference(g.omega, 0, expected, 0), Real(32.0) *
              std::numeric_limits<Real>::epsilon());

    AuxiliaryInertTracer tracer(1);
    tracer.define(0, g.ba, g.dm);
    EXPECT_FALSE(tracer.measure_is_ready(0));
    ASSERT_TRUE(tracer.rebuild_static_measure(0, g.detj, g.mx, g.my, diagnostic))
        << diagnostic;
    EXPECT_TRUE(tracer.measure_is_ready(0));
    EXPECT_LT(max_component_difference(tracer.static_measure(0), 0, expected, 0),
              Real(32.0) * std::numeric_limits<Real>::epsilon());

    MultiFab previous_measure(g.ba, g.dm, 1, 0);
    MultiFab::Copy(previous_measure, tracer.static_measure(0), 0, 0, 1, 0);
    g.detj.mult(Real(1.23), 0, 1, 0);
    g.mx.mult(Real(0.91), 0, 1, 0);
    g.my.mult(Real(1.07), 0, 1, 0);
    ASSERT_TRUE(tracer.rebuild_static_measure(0, g.detj, g.mx, g.my, diagnostic))
        << diagnostic;
    MultiFab rebuilt_expected(g.ba, g.dm, 1, 0);
    for (amrex::MFIter mfi(rebuilt_expected); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto d = g.detj.const_array(mfi);
        const auto x = g.mx.const_array(mfi);
        const auto y = g.my.const_array(mfi);
        const auto e = rebuilt_expected.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            e(i, j, k, 0) = d(i, j, k, 0) / (x(i, j, 0) * y(i, j, 0));
        });
    }
    EXPECT_LT(max_component_difference(tracer.static_measure(0), 0, rebuilt_expected, 0),
              Real(32.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_GT(max_component_difference(tracer.static_measure(0), 0, previous_measure, 0),
              Real(1.0e-3));

    g.detj.setVal(Real(0.0));
    EXPECT_FALSE(tracer.rebuild_static_measure(0, g.detj, g.mx, g.my, diagnostic));
    EXPECT_FALSE(tracer.measure_is_ready(0));

    for (int bad_kind = 0; bad_kind < 8; ++bad_kind) {
        g.detj.setVal(Real(1.2)); g.mx.setVal(Real(1.1)); g.my.setVal(Real(0.9));
        if (bad_kind == 0) g.detj.setVal(std::numeric_limits<Real>::quiet_NaN());
        if (bad_kind == 1) g.detj.setVal(std::numeric_limits<Real>::infinity());
        if (bad_kind == 2) g.detj.setVal(Real(0.0));
        if (bad_kind == 3) g.detj.setVal(Real(-1.0));
        if (bad_kind == 4) g.mx.setVal(Real(0.0));
        if (bad_kind == 5) g.my.setVal(Real(-0.2));
        if (bad_kind == 6) g.mx.setVal(std::numeric_limits<Real>::infinity());
        if (bad_kind == 7) g.my.setVal(std::numeric_limits<Real>::quiet_NaN());
        EXPECT_FALSE(BuildMappedCellMeasure(g.omega, g.detj, g.mx, g.my, diagnostic))
            << "invalid metric case " << bad_kind;
    }
}

void run_auxiliary_mapped_transport_SharedLayoutPredicatesRejectIncompatibleLayouts ()
{
    TestGrid g;
    MultiFab compressed(g.map_ba, g.dm, 1, 0);
    EXPECT_FALSE(SameCellLayout(compressed, g.omega));

    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 1, 0);
    EXPECT_TRUE(MappedFaceLayoutMatchesCellLayout(rate, g.omega));
    EXPECT_FALSE(MappedFaceLayoutMatchesCellLayout(rate, compressed));

    Box shifted_domain = g.domain;
    shifted_domain.shift(0, 17);
    const BoxArray shifted_ba(shifted_domain);
    const DistributionMapping shifted_dm(shifted_ba);
    MappedFaceFluxRate shifted_rate;
    shifted_rate.define(shifted_ba, shifted_dm, 1, 0);
    EXPECT_FALSE(MappedFaceLayoutMatchesCellLayout(shifted_rate, g.omega));
    EXPECT_FALSE(SameMappedFaceLayout(rate, shifted_rate));
}

void run_auxiliary_mapped_transport_LedgerRejectsMismatchedFaceLayout ()
{
    TestGrid g;
    CompletedStepFluxLedger ledger;
    ledger.define(g.ba, g.dm);

    Box shifted_domain = g.domain;
    shifted_domain.shift(1, 13);
    const BoxArray shifted_ba(shifted_domain);
    const DistributionMapping shifted_dm(shifted_ba);
    MappedFaceFluxRate mismatched_rate;
    mismatched_rate.define(shifted_ba, shifted_dm, 1, 0);
    AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 0,
                                         0.25, recipe, diagnostic)) << diagnostic;
    EXPECT_FALSE(ledger.accept_stage(HostIntegrator::CompressibleRK3, 0, 0.0,
                                     recipe, mismatched_rate, diagnostic));
    EXPECT_FALSE(diagnostic.empty());
    EXPECT_FALSE(ledger.step_active());
    EXPECT_EQ(ledger.next_stage(), 0);
}

void run_auxiliary_mapped_transport_FluxRateAndIntegratedFluxHaveDistinctLayoutsAndTypes ()
{
    TestGrid g;
    MappedFaceFluxRate rate;
    IntegratedMappedFaceFlux integral;
    rate.define(g.ba, g.dm, 3, 0);
    integral.define(g.ba, g.dm, 3, 0);
    EXPECT_TRUE(rate.is_defined());
    EXPECT_TRUE(integral.is_defined());
    EXPECT_EQ(rate.nComp(), 3);
    EXPECT_EQ(integral.nComp(), 3);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        EXPECT_EQ(rate.dir(dir).boxArray(),
                  amrex::convert(g.ba, IntVect::TheDimensionVector(dir)));
        EXPECT_EQ(integral.dir(dir).boxArray(),
                  amrex::convert(g.ba, IntVect::TheDimensionVector(dir)));
        EXPECT_EQ(rate.dir(dir).DistributionMap(), g.dm);
    }
}

void run_auxiliary_mapped_transport_ComputationalMappedDivergenceMatchesIndependentArithmetic ()
{
    TestGrid g;
    ASSERT_TRUE(g.build_measure());
    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 3, 0);
    fill_rate(rate, g.domain, 2);
    MultiFab actual(g.ba, g.dm, 4, 0), oracle(g.ba, g.dm, 4, 0);
    actual.setVal(Real(-9.0)); oracle.setVal(Real(-9.0));
    const auto inv = g.geom.InvCellSizeArray();
    ApplyMappedFluxTendency(rate, 2, g.omega, 0, actual, 3, inv);
    for (amrex::MFIter mfi(actual); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto fx = rate.dir(0).const_array(mfi);
        const auto fy = rate.dir(1).const_array(mfi);
        const auto fz = rate.dir(2).const_array(mfi);
        const auto omega = g.omega.const_array(mfi);
        const auto out = oracle.array(mfi);
        const Real dx = inv[0], dy = inv[1], dz = inv[2];
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real arithmetic =
                (fx(i + 1, j, k, 2) - fx(i, j, k, 2)) * dx +
                (fy(i, j + 1, k, 2) - fy(i, j, k, 2)) * dy +
                (fz(i, j, k + 1, 2) - fz(i, j, k, 2)) * dz;
            out(i, j, k, 3) = -arithmetic / omega(i, j, k, 0);
        });
    }
    EXPECT_LT(max_component_difference(actual, 3, oracle, 3),
              Real(32.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_DOUBLE_EQ(ComputationalMappedDivergence(7, 2, 5, 1, 11, 3,
                                                    Real(0.5), Real(0.25), Real(2.0)),
                     Real(5.0) * Real(0.5) + Real(4.0) * Real(0.25) + Real(8.0) * Real(2.0));
}

void run_auxiliary_mapped_transport_PeriodicArbitraryMappedFluxTelescopes ()
{
    TestGrid g;
    ASSERT_TRUE(g.build_measure());
    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 1, 0);
    fill_rate(rate, g.domain, 0, true);
    MultiFab tendency(g.ba, g.dm, 1, 0), state(g.ba, g.dm, 1, 0);
    MultiFab h0(g.ba, g.dm, 1, 0), h1(g.ba, g.dm, 1, 0);
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto u = state.array(mfi);
        const auto omega = g.omega.const_array(mfi);
        const auto a = h0.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            u(i, j, k, 0) = Real(1.3) + Real(0.1) * i - Real(0.04) * j + Real(0.03) * k;
            a(i, j, k, 0) = omega(i, j, k, 0) * u(i, j, k, 0);
        });
    }
    ApplyMappedFluxTendency(rate, 0, g.omega, 0, tendency, 0,
                            g.geom.InvCellSizeArray());
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto u = state.array(mfi);
        const auto tendency_arr = tendency.const_array(mfi);
        const auto omega = g.omega.const_array(mfi);
        const auto b = h1.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            u(i, j, k, 0) += Real(0.071) * tendency_arr(i, j, k, 0);
            b(i, j, k, 0) = omega(i, j, k, 0) * u(i, j, k, 0);
        });
    }
    const Real total_change = h1.sum(0) - h0.sum(0);
    EXPECT_NEAR(total_change, Real(0.0),
                Real(256.0) * std::numeric_limits<Real>::epsilon() *
                std::max(Real(1.0), amrex::Math::abs(h0.sum(0))));
}

void run_auxiliary_mapped_transport_NativeScalarAdvectionParityAndMetricNegativeControls ()
{
    TestGrid g;
    ASSERT_TRUE(g.build_measure());
    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 3, 0);
    fill_rate(rate, g.domain, 2);
    MultiFab native(g.ba, g.dm, 4, 0), generic(g.ba, g.dm, 4, 0);
    MultiFab wrong_measure(g.ba, g.dm, 4, 0), wrong_map(g.ba, g.dm, 4, 0);
    native.setVal(Real(-99.0)); generic.setVal(Real(-99.0));
    wrong_measure.setVal(Real(-99.0)); wrong_map.setVal(Real(-99.0));
    MultiFab omega_without_detj(g.ba, g.dm, 1, 0);
    for (amrex::MFIter mfi(g.omega); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto x = g.mx.const_array(mfi);
        const auto y = g.my.const_array(mfi);
        const auto out = omega_without_detj.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = Real(1.0) / (x(i, j, 0) * y(i, j, 0));
        });
    }
    const auto inv = g.geom.InvCellSizeArray();
    for (amrex::MFIter mfi(native); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const GpuArray<const Array4<Real>, AMREX_SPACEDIM> flux{{
            rate.dir(0).array(mfi), rate.dir(1).array(mfi),
            rate.dir(2).array(mfi)}};
        ApplyScalarAdvectionFluxDivergence(
            bx, flux, 2, native.array(mfi), 1, g.detj.const_array(mfi), inv,
            g.mx.const_array(mfi), g.my.const_array(mfi));
    }
    ApplyMappedFluxTendency(rate, 2, omega_without_detj, 0, wrong_measure, 1, inv);
    ApplyMappedFluxTendency(rate, 2, g.omega, 0, generic, 1, inv);
    EXPECT_LT(max_component_difference(native, 1, generic, 1),
              Real(64.0) * std::numeric_limits<Real>::epsilon());
    const Real omitted_detj_discrepancy =
        max_component_difference(native, 1, wrong_measure, 1);
    EXPECT_GT(omitted_detj_discrepancy, Real(1.0e-3));

    for (amrex::MFIter mfi(native); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto correct = generic.const_array(mfi);
        const auto mx = g.mx.const_array(mfi);
        const auto out = wrong_map.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 1) = correct(i, j, k, 1) * mx(i, j, 0);
        });
    }
    const Real doubled_map_discrepancy = max_component_difference(native, 1, wrong_map, 1);
    EXPECT_GT(doubled_map_discrepancy, Real(1.0e-3));
    amrex::Print() << "AUX_M2_METRIC_CONTROL omitted_detJ=" << omitted_detj_discrepancy
                   << " doubled_map_factor=" << doubled_map_discrepancy << std::endl;
}

void run_auxiliary_mapped_transport_CanonicalNGridDiffusionTransferParity ()
{
    EXPECT_LT(mapped_diffusion_parity('N'), Real(64.0) * std::numeric_limits<Real>::epsilon());
}

void run_auxiliary_mapped_transport_CanonicalStretchedGridDiffusionTransferParity ()
{
    EXPECT_LT(mapped_diffusion_parity('S'), Real(64.0) * std::numeric_limits<Real>::epsilon());
}

void run_auxiliary_mapped_transport_CanonicalStaticTerrainParityAndMetricNegativeControls ()
{
    Real raw_fz_discrepancy = 0.0;
    Real missing_cross_discrepancy = 0.0;
    const Real terrain_parity = mapped_diffusion_parity('T', &raw_fz_discrepancy,
                                                         &missing_cross_discrepancy);
    EXPECT_LT(terrain_parity,
              Real(128.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_GT(raw_fz_discrepancy, Real(1.0e-3));
    EXPECT_GT(missing_cross_discrepancy, Real(1.0e-3));
    amrex::Print() << "AUX_M2_METRIC_CONTROL terrain_parity=" << terrain_parity
                   << " raw_terrain_Fz=" << raw_fz_discrepancy
                   << " omitted_terrain_cross_terms=" << missing_cross_discrepancy
                   << std::endl;
}

void run_auxiliary_mapped_transport_CompressibleRK3RecipeUsesAuditedStageCoefficients ()
{
    constexpr double dt = 0.37;
    const double trial[] = {dt / 3.0, dt / 2.0, dt};
    const double face[] = {dt / 3.0, dt / 2.0, dt};
    const double ledger[] = {0.0, 0.0, dt};
    AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    for (int stage = 0; stage < 3; ++stage) {
        ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3,
                                             stage, trial[stage], recipe, diagnostic))
            << diagnostic;
        EXPECT_DOUBLE_EQ(recipe.anchor_weight, 1.0);
        EXPECT_DOUBLE_EQ(recipe.input_weight, 0.0);
        EXPECT_EQ(recipe.limiter_trial_base, LimiterTrialBase::Anchor);
        EXPECT_DOUBLE_EQ(recipe.limiter_trial_interval, trial[stage]);
        EXPECT_DOUBLE_EQ(recipe.face_rate_time_coefficient, face[stage]);
        EXPECT_DOUBLE_EQ(recipe.completed_ledger_time, ledger[stage]);
    }
}

void run_auxiliary_mapped_transport_HeunRecipeSeparatesTrialAndWeightedFaceTime ()
{
    constexpr double dt = 0.61;
    AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 0, dt,
                                         recipe, diagnostic)) << diagnostic;
    EXPECT_DOUBLE_EQ(recipe.anchor_weight, 1.0);
    EXPECT_DOUBLE_EQ(recipe.input_weight, 0.0);
    EXPECT_EQ(recipe.limiter_trial_base, LimiterTrialBase::Anchor);
    EXPECT_DOUBLE_EQ(recipe.limiter_trial_interval, dt);
    EXPECT_DOUBLE_EQ(recipe.face_rate_time_coefficient, dt);
    EXPECT_DOUBLE_EQ(recipe.completed_ledger_time, 0.5 * dt);

    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 1, dt,
                                         recipe, diagnostic)) << diagnostic;
    EXPECT_DOUBLE_EQ(recipe.anchor_weight, 0.5);
    EXPECT_DOUBLE_EQ(recipe.input_weight, 0.5);
    EXPECT_EQ(recipe.limiter_trial_base, LimiterTrialBase::Input);
    EXPECT_DOUBLE_EQ(recipe.limiter_trial_interval, dt);
    EXPECT_DOUBLE_EQ(recipe.face_rate_time_coefficient, 0.5 * dt);
    EXPECT_DOUBLE_EQ(recipe.completed_ledger_time, 0.5 * dt);
    EXPECT_NE(recipe.limiter_trial_interval, recipe.face_rate_time_coefficient);
}

void run_auxiliary_mapped_transport_TimedInputViewsBuildIntensiveStateFromRhoInput ()
{
    TestGrid g;
    MultiFab u(g.ba, g.dm, 3, 0), rho_anchor(g.ba, g.dm, 3, 0);
    MultiFab rho_input(g.ba, g.dm, 3, 0), rho_target(g.ba, g.dm, 3, 0);
    MultiFab intensive(g.ba, g.dm, 1, 1), fast(g.ba, g.dm, 1, 1);
    MultiFab wrong(g.ba, g.dm, 1, 0), expected(g.ba, g.dm, 1, 0);
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto ua = u.array(mfi);
        const auto ra = rho_anchor.array(mfi);
        const auto ri = rho_input.array(mfi);
        const auto rt = rho_target.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int kidx) noexcept {
            const Real a = Real(0.8) + Real(0.01) * i + Real(0.02) * j;
            const Real in = Real(1.2) + Real(0.03) * i + Real(0.01) * kidx;
            const Real target = Real(1.7) + Real(0.02) * j + Real(0.04) * kidx;
            ra(i, j, kidx, 0) = a;
            ri(i, j, kidx, 0) = in;
            rt(i, j, kidx, 0) = target;
            ua(i, j, kidx, 0) = Real(0.2) + Real(0.07) * i +
                                Real(0.03) * j + Real(0.02) * kidx;
            ua(i, j, kidx, 1) = Real(120.0) + target;
            ua(i, j, kidx, 2) = Real(-31.0) + a;
        });
    }
    std::string diagnostic;
    ASSERT_TRUE(BuildAuxiliaryIntensiveState(
        {&u, 0, 1.25}, {&rho_input, 0, 1.25}, 1.25, intensive,
        AuxiliaryFieldValidationPolicy::Global, diagnostic)) << diagnostic;
    ASSERT_TRUE(BuildAuxiliaryIntensiveState(
        {&u, 0, 1.25}, {&rho_input, 0, 1.25}, 1.25, fast,
        AuxiliaryFieldValidationPolicy::AssumeValid, diagnostic)) << diagnostic;
    for (amrex::MFIter mfi(expected); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = u.const_array(mfi);
        const auto rho_in = rho_input.const_array(mfi);
        const auto rho_wrong = rho_target.const_array(mfi);
        const auto out = expected.array(mfi);
        const auto w = wrong.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int kidx) noexcept {
            out(i, j, kidx, 0) = state(i, j, kidx, 0) / rho_in(i, j, kidx, 0);
            w(i, j, kidx, 0) = state(i, j, kidx, 0) / rho_wrong(i, j, kidx, 0);
        });
    }
    EXPECT_LT(max_component_difference(intensive, 0, expected, 0),
              Real(32.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_LT(max_component_difference(fast, 0, expected, 0),
              Real(32.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_LT(max_component_difference(intensive, 0, fast, 0),
              Real(32.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_GT(max_component_difference(intensive, 0, wrong, 0), Real(0.05));

    EXPECT_FALSE(BuildAuxiliaryIntensiveState(
        {&u, 0, 0.0}, {&rho_input, 0, 1.0}, 1.0, intensive,
        AuxiliaryFieldValidationPolicy::Global, diagnostic));
    const std::array<Real, 3> invalid_density{
        Real(0.0), std::numeric_limits<Real>::quiet_NaN(),
        std::numeric_limits<Real>::infinity()};
    for (const Real invalid : invalid_density) {
        rho_input.setVal(Real(1.0));
        rho_input.setVal(invalid, 0, 1, 0);
        EXPECT_FALSE(BuildAuxiliaryIntensiveState(
            {&u, 0, 1.25}, {&rho_input, 0, 1.25}, 1.25, intensive,
            AuxiliaryFieldValidationPolicy::Global, diagnostic)) << diagnostic;
    }
}

void run_auxiliary_mapped_transport_ConstantConstituentRatioSurvivesVariableDensityStages ()
{
    auto run = [](const HostIntegrator method) {
        TestGrid g;
        ASSERT_TRUE(g.build_measure());
        constexpr Real k = Real(0.43);
        constexpr double full_dt = 0.0017;
        MultiFab rho_anchor(g.ba, g.dm, 1, 0), rho_input(g.ba, g.dm, 1, 0);
        MultiFab rho_target(g.ba, g.dm, 1, 0);
        MultiFab u_anchor(g.ba, g.dm, 1, 0), u_input(g.ba, g.dm, 1, 0);
        MultiFab u_target(g.ba, g.dm, 1, 0), intensive(g.ba, g.dm, 1, 1);
        MultiFab ratio(g.ba, g.dm, 1, 0);
        MappedFaceFluxRate carrier, aux_rate;
        carrier.define(g.ba, g.dm, 1, 0);
        aux_rate.define(g.ba, g.dm, 1, 0);
        MultiFab avg_x(carrier.dir(0).boxArray(), g.dm, 1, 0);
        MultiFab avg_y(carrier.dir(1).boxArray(), g.dm, 1, 0);
        MultiFab avg_z(carrier.dir(2).boxArray(), g.dm, 1, 0);

        for (amrex::MFIter mfi(rho_anchor); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto rho = rho_anchor.array(mfi);
            const auto u = u_anchor.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int kidx) noexcept {
                const Real value = Real(1.4) + Real(0.021) * i + Real(0.017) * j + Real(0.013) * kidx;
                rho(i, j, kidx, 0) = value;
                u(i, j, kidx, 0) = k * value;
            });
        }
        amrex::MultiFab::Copy(rho_input, rho_anchor, 0, 0, 1, 0);
        amrex::MultiFab::Copy(u_input, u_anchor, 0, 0, 1, 0);
        std::string diagnostic;
        const auto inv = g.geom.InvCellSizeArray();
        const int stages = method == HostIntegrator::CompressibleRK3 ? 3 : 2;
        double input_time = 0.0;
        for (int stage = 0; stage < stages; ++stage) {
            const double interval = method == HostIntegrator::CompressibleRK3 ?
                (stage == 0 ? full_dt / 3.0 : (stage == 1 ? full_dt / 2.0 : full_dt)) : full_dt;
            const double target_time = method == HostIntegrator::CompressibleRK3 ?
                (stage == 0 ? full_dt / 3.0 : (stage == 1 ? full_dt / 2.0 : full_dt)) : full_dt;
            const Real carrier_factor = Real(1.0) - Real(0.13) * stage;
            fill_rate(carrier, g.domain, 0, false, carrier_factor);
            amrex::MultiFab::Copy(avg_x, carrier.dir(0), 0, 0, 1, 0);
            amrex::MultiFab::Copy(avg_y, carrier.dir(1), 0, 0, 1, 0);
            amrex::MultiFab::Copy(avg_z, carrier.dir(2), 0, 0, 1, 0);
            if (!BuildAuxiliaryIntensiveState({&u_input, 0, input_time},
                                              {&rho_input, 0, input_time},
                                              input_time, intensive,
                                              AuxiliaryFieldValidationPolicy::Global,
                                              diagnostic)) {
                ADD_FAILURE() << diagnostic;
                return;
            }
            intensive.FillBoundary(g.geom.periodicity());
            aux_rate.setVal(Real(0.0));
            for (amrex::MFIter mfi(intensive); mfi.isValid(); ++mfi) {
                const Box bx = mfi.validbox();
                GpuArray<const Array4<Real>, AMREX_SPACEDIM> flux{{
                    aux_rate.dir(0).array(mfi), aux_rate.dir(1).array(mfi),
                    aux_rate.dir(2).array(mfi)}};
                BuildScalarAdvectionFluxes(
                    bx, intensive.const_array(mfi), 0, flux, 0,
                    avg_x.const_array(mfi), avg_y.const_array(mfi), avg_z.const_array(mfi),
                    AdvType::Centered_2nd, AdvType::Centered_2nd, Real(1.0), Real(1.0));
            }
            AuxiliaryStageRecipe recipe;
            ASSERT_TRUE(MakeAuxiliaryStageRecipe(method, stage, interval,
                                                 recipe, diagnostic)) << diagnostic;
            AuxiliaryStageContext density_context;
            density_context.method = method;
            density_context.level = 0;
            density_context.stage = stage;
            density_context.step_old_time = 0.0;
            density_context.input_time = input_time;
            density_context.target_time = target_time;
            density_context.recurrence = recipe;
            density_context.state_anchor = {&rho_anchor, 0, 0.0};
            density_context.state_input = {&rho_input, 0, input_time};
            density_context.state_target = {&rho_target, 0, target_time};
            density_context.rho_anchor = {&rho_anchor, 0, 0.0};
            density_context.rho_input = {&rho_input, 0, input_time};
            density_context.rho_target = {&rho_target, 0, target_time};
            density_context.measure_anchor = {&g.omega, 0, 0.0};
            density_context.measure_input = {&g.omega, 0, input_time};
            density_context.measure_target = {&g.omega, 0, target_time};
            density_context.carrier = {&avg_x, &avg_y, &avg_z};
            ApplyAuxiliaryMappedStage(density_context, carrier, inv);

            AuxiliaryStageContext auxiliary_context = density_context;
            auxiliary_context.state_anchor = {&u_anchor, 0, 0.0};
            auxiliary_context.state_input = {&u_input, 0, input_time};
            auxiliary_context.state_target = {&u_target, 0, target_time};
            ApplyAuxiliaryMappedStage(auxiliary_context, aux_rate, inv);

            for (amrex::MFIter mfi(ratio); mfi.isValid(); ++mfi) {
                const Box bx = mfi.validbox();
                const auto rho = rho_target.const_array(mfi);
                const auto u = u_target.const_array(mfi);
                const auto r = ratio.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int kidx) noexcept {
                    r(i, j, kidx, 0) = u(i, j, kidx, 0) / rho(i, j, kidx, 0);
                });
            }
            MultiFab expected(g.ba, g.dm, 1, 0);
            expected.setVal(k);
            EXPECT_LT(max_component_difference(ratio, 0, expected, 0),
                      Real(128.0) * std::numeric_limits<Real>::epsilon());
            EXPECT_GT(max_component_difference(rho_anchor, 0, rho_target, 0), Real(1.0e-9));
            EXPECT_GT(max_component_difference(rho_input, 0, rho_target, 0), Real(1.0e-9));

            amrex::MultiFab::Copy(rho_input, rho_target, 0, 0, 1, 0);
            amrex::MultiFab::Copy(u_input, u_target, 0, 0, 1, 0);
            input_time = target_time;
        }
    };
    run(HostIntegrator::CompressibleRK3);
    run(HostIntegrator::AnelasticHeun);
}

void run_auxiliary_mapped_transport_CompletedLedgerUsesExactHostTemporalWeights ()
{
    TestGrid g;
    constexpr double dt = 0.41;
    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 1, 0);

    CompletedStepFluxLedger compressible;
    compressible.define(g.ba, g.dm);
    AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    fill_constant_rate(rate, Real(1.0));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 0,
                                         dt / 3.0, recipe, diagnostic));
    ASSERT_TRUE(compressible.accept_stage(HostIntegrator::CompressibleRK3, 0,
                                          0.0, recipe, rate, diagnostic)) << diagnostic;
    fill_constant_rate(rate, Real(3.0));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 1,
                                         dt / 2.0, recipe, diagnostic));
    ASSERT_TRUE(compressible.accept_stage(HostIntegrator::CompressibleRK3, 1,
                                          0.0, recipe, rate, diagnostic)) << diagnostic;
    fill_constant_rate(rate, Real(7.0));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 2,
                                         dt, recipe, diagnostic));
    ASSERT_TRUE(compressible.accept_stage(HostIntegrator::CompressibleRK3, 2,
                                          0.0, recipe, rate, diagnostic)) << diagnostic;
    ASSERT_TRUE(compressible.step_complete());
    const Real compressible_expected = Real(dt) * Real(7.0);
    const Real compressible_wrong = Real(dt / 3.0) * Real(1.0) +
                                    Real(dt / 2.0) * Real(3.0) + Real(dt) * Real(7.0);
    EXPECT_NEAR(compressible.integrated_flux().dir(0).norm0(0), compressible_expected,
                Real(32.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_GT(amrex::Math::abs(compressible_expected - compressible_wrong), Real(0.7));

    CompletedStepFluxLedger heun;
    heun.define(g.ba, g.dm);
    fill_constant_rate(rate, Real(2.0));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 0,
                                         dt, recipe, diagnostic));
    ASSERT_TRUE(heun.accept_stage(HostIntegrator::AnelasticHeun, 0,
                                  0.0, recipe, rate, diagnostic)) << diagnostic;
    fill_constant_rate(rate, Real(9.0));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 1,
                                         dt, recipe, diagnostic));
    ASSERT_TRUE(heun.accept_stage(HostIntegrator::AnelasticHeun, 1,
                                  0.0, recipe, rate, diagnostic)) << diagnostic;
    ASSERT_TRUE(heun.step_complete());
    const Real heun_expected = Real(0.5 * dt) * Real(2.0) + Real(0.5 * dt) * Real(9.0);
    EXPECT_NEAR(heun.integrated_flux().dir(0).norm0(0), heun_expected,
                Real(32.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_GT(amrex::Math::abs(heun_expected - Real(dt) * Real(9.0)), Real(0.1));
    EXPECT_GT(amrex::Math::abs(heun_expected - Real(dt) * Real(11.0)), Real(0.1));
}

void run_auxiliary_mapped_transport_CompletedLedgerAccumulatesEveryComponent ()
{
    TestGrid g;
    constexpr int ncomp = 3;
    constexpr double dt = 0.41;
    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, ncomp, 0);
    AuxiliaryStageRecipe recipe;
    std::string diagnostic;

    CompletedStepFluxLedger compressible;
    compressible.define(g.ba, g.dm, ncomp);
    fill_componentwise_constant_rate(rate, {Real(1.0), Real(2.0), Real(3.0)});
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 0,
                                         dt / 3.0, recipe, diagnostic)) << diagnostic;
    ASSERT_TRUE(compressible.accept_stage(HostIntegrator::CompressibleRK3, 0,
                                          0.0, recipe, rate, diagnostic)) << diagnostic;
    fill_componentwise_constant_rate(rate, {Real(4.0), Real(5.0), Real(6.0)});
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 1,
                                         dt / 2.0, recipe, diagnostic)) << diagnostic;
    ASSERT_TRUE(compressible.accept_stage(HostIntegrator::CompressibleRK3, 1,
                                          0.0, recipe, rate, diagnostic)) << diagnostic;
    fill_componentwise_constant_rate(rate, {Real(7.0), Real(11.0), Real(13.0)});
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 2,
                                         dt, recipe, diagnostic)) << diagnostic;
    ASSERT_TRUE(compressible.accept_stage(HostIntegrator::CompressibleRK3, 2,
                                          0.0, recipe, rate, diagnostic)) << diagnostic;
    ASSERT_TRUE(compressible.step_complete());

    const std::array<Real, ncomp> compressible_expected{
        Real(dt) * Real(7.0), Real(dt) * Real(11.0), Real(dt) * Real(13.0)};
    const Real tolerance = Real(32.0) * std::numeric_limits<Real>::epsilon();
    for (int comp = 0; comp < ncomp; ++comp) {
        EXPECT_NEAR(max_face_component_error(compressible.integrated_flux(), comp,
                                             compressible_expected[static_cast<std::size_t>(comp)]),
                    Real(0.0), tolerance) << "compressible component " << comp;
    }

    CompletedStepFluxLedger heun;
    heun.define(g.ba, g.dm, ncomp);
    fill_componentwise_constant_rate(rate, {Real(2.0), Real(4.0), Real(6.0)});
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 0,
                                         dt, recipe, diagnostic)) << diagnostic;
    ASSERT_TRUE(heun.accept_stage(HostIntegrator::AnelasticHeun, 0,
                                  0.0, recipe, rate, diagnostic)) << diagnostic;
    fill_componentwise_constant_rate(rate, {Real(10.0), Real(20.0), Real(30.0)});
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 1,
                                         dt, recipe, diagnostic)) << diagnostic;
    ASSERT_TRUE(heun.accept_stage(HostIntegrator::AnelasticHeun, 1,
                                  0.0, recipe, rate, diagnostic)) << diagnostic;
    ASSERT_TRUE(heun.step_complete());

    const std::array<Real, ncomp> heun_expected{
        Real(0.5 * dt) * Real(12.0),
        Real(0.5 * dt) * Real(24.0),
        Real(0.5 * dt) * Real(36.0)};
    for (int comp = 0; comp < ncomp; ++comp) {
        EXPECT_NEAR(max_face_component_error(heun.integrated_flux(), comp,
                                             heun_expected[static_cast<std::size_t>(comp)]),
                    Real(0.0), tolerance) << "Heun component " << comp;
    }
}

void run_auxiliary_mapped_transport_StageTargetMustBeDisjoint ()
{
    TestGrid g;
    MultiFab anchor(g.ba, g.dm, 1, 0);
    MultiFab input(g.ba, g.dm, 1, 0);
    MultiFab target(g.ba, g.dm, 1, 0);
    AuxiliaryStageContext context;
    context.state_anchor = {&anchor, 0, 0.0};
    context.state_input = {&input, 0, 1.0};
    context.state_target = {&target, 0, 1.0};
    EXPECT_TRUE(AuxiliaryStageTargetIsDisjoint(context));

    context.state_target.field = &input;
    EXPECT_FALSE(AuxiliaryStageTargetIsDisjoint(context));
    context.state_target.field = &anchor;
    EXPECT_FALSE(AuxiliaryStageTargetIsDisjoint(context));

    MultiFab input_alias(input, amrex::make_alias, 0, 1);
    context.state_target.field = &input_alias;
    EXPECT_FALSE(AuxiliaryStageTargetIsDisjoint(context));

    // Read-only anchor and input views may intentionally alias.
    context.state_input.field = &anchor;
    context.state_target.field = &target;
    EXPECT_TRUE(AuxiliaryStageTargetIsDisjoint(context));
}

void run_auxiliary_mapped_transport_ZeroRateStillAppliesHostAnchorInputRecurrence ()
{
    TestGrid g;
    MappedFaceFluxRate zero_rate;
    zero_rate.define(g.ba, g.dm, 1, 0);
    zero_rate.setVal(Real(0.0));
    MultiFab anchor(g.ba, g.dm, 1, 0), input(g.ba, g.dm, 1, 0), target(g.ba, g.dm, 1, 0);
    MultiFab rho_anchor(g.ba, g.dm, 1, 0), rho_input(g.ba, g.dm, 1, 0), rho_target(g.ba, g.dm, 1, 0);
    MultiFab omega_anchor(g.ba, g.dm, 1, 0), omega_input(g.ba, g.dm, 1, 0), omega_target(g.ba, g.dm, 1, 0);
    MultiFab expected(g.ba, g.dm, 1, 0);
    for (amrex::MFIter mfi(anchor); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto a = anchor.array(mfi); const auto i = input.array(mfi);
        const auto ra = rho_anchor.array(mfi); const auto ri = rho_input.array(mfi);
        const auto rt = rho_target.array(mfi);
        const auto oa = omega_anchor.array(mfi); const auto oi = omega_input.array(mfi);
        const auto omega_target_array = omega_target.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int x, int y, int z) noexcept {
            a(x, y, z, 0) = Real(1.1) + Real(0.1) * x;
            i(x, y, z, 0) = Real(2.4) - Real(0.07) * y + Real(0.02) * z;
            ra(x, y, z, 0) = Real(1.0) + Real(0.03) * x;
            ri(x, y, z, 0) = Real(1.2) + Real(0.04) * y;
            rt(x, y, z, 0) = Real(1.4) + Real(0.05) * z;
            oa(x, y, z, 0) = Real(0.9) + Real(0.02) * x;
            oi(x, y, z, 0) = Real(1.3) + Real(0.03) * y;
            omega_target_array(x, y, z, 0) = Real(1.7) + Real(0.01) * z;
        });
    }
    auto make_context = [&](const double aweight, const double iweight) {
        AuxiliaryStageContext context;
        context.level = 0; context.stage = 1;
        context.step_old_time = 0.0; context.input_time = 1.0; context.target_time = 1.0;
        context.recurrence.anchor_weight = aweight;
        context.recurrence.input_weight = iweight;
        context.recurrence.face_rate_time_coefficient = 0.5;
        context.state_anchor = {&anchor, 0, 0.0};
        context.state_input = {&input, 0, 1.0};
        context.state_target = {&target, 0, 1.0};
        context.rho_anchor = {&rho_anchor, 0, 0.0};
        context.rho_input = {&rho_input, 0, 1.0};
        context.rho_target = {&rho_target, 0, 1.0};
        context.measure_anchor = {&omega_anchor, 0, 0.0};
        context.measure_input = {&omega_input, 0, 1.0};
        context.measure_target = {&omega_target, 0, 1.0};
        context.carrier = {&rho_anchor, &rho_input, &rho_target};
        return context;
    };
    auto comp = make_context(1.0, 0.0);
    ApplyAuxiliaryMappedStage(comp, zero_rate, g.geom.InvCellSizeArray());
    for (amrex::MFIter mfi(expected); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto a = anchor.const_array(mfi); const auto oa = omega_anchor.const_array(mfi);
        const auto omega_target_array = omega_target.const_array(mfi);
        const auto e = expected.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            e(i, j, k, 0) = oa(i, j, k, 0) * a(i, j, k, 0) /
                            omega_target_array(i, j, k, 0);
        });
    }
    EXPECT_LT(max_component_difference(target, 0, expected, 0),
              Real(32.0) * std::numeric_limits<Real>::epsilon());

    auto heun = make_context(0.5, 0.5);
    ApplyAuxiliaryMappedStage(heun, zero_rate, g.geom.InvCellSizeArray());
    for (amrex::MFIter mfi(expected); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto a = anchor.const_array(mfi); const auto i = input.const_array(mfi);
        const auto oa = omega_anchor.const_array(mfi); const auto oi = omega_input.const_array(mfi);
        const auto omega_target_array = omega_target.const_array(mfi);
        const auto e = expected.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int x, int y, int z) noexcept {
            e(x, y, z, 0) = (Real(0.5) * oa(x, y, z, 0) * a(x, y, z, 0) +
                             Real(0.5) * oi(x, y, z, 0) * i(x, y, z, 0)) /
                            omega_target_array(x, y, z, 0);
        });
    }
    EXPECT_LT(max_component_difference(target, 0, expected, 0),
              Real(32.0) * std::numeric_limits<Real>::epsilon());
}

void run_auxiliary_mapped_transport_StageSequenceFailsClosed ()
{
    TestGrid g;
    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 1, 0);
    fill_constant_rate(rate, Real(0.0));
    AuxiliaryStageRecipe recipe;
    std::string diagnostic;
    EXPECT_FALSE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 3, 1.0,
                                          recipe, diagnostic));
    EXPECT_FALSE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 2, 1.0,
                                          recipe, diagnostic));
    EXPECT_FALSE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 0, -1.0,
                                          recipe, diagnostic));
    EXPECT_FALSE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 0,
                                          std::numeric_limits<double>::quiet_NaN(),
                                          recipe, diagnostic));
    EXPECT_FALSE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 0,
                                          std::numeric_limits<double>::infinity(),
                                          recipe, diagnostic));
    EXPECT_STREQ(HostIntegratorName(HostIntegrator::AnelasticMidPoint), "AnelasticMidPoint");
    EXPECT_FALSE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticMidPoint, 0, 1.0,
                                          recipe, diagnostic));

    CompletedStepFluxLedger before_zero;
    before_zero.define(g.ba, g.dm);
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 1, 1.0,
                                         recipe, diagnostic));
    EXPECT_FALSE(before_zero.accept_stage(HostIntegrator::AnelasticHeun, 1,
                                          0.0, recipe, rate, diagnostic));

    CompletedStepFluxLedger sequence;
    sequence.define(g.ba, g.dm);
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 0,
                                         1.0 / 3.0, recipe, diagnostic));
    ASSERT_TRUE(sequence.accept_stage(HostIntegrator::CompressibleRK3, 0,
                                     0.0, recipe, rate, diagnostic));
    EXPECT_FALSE(sequence.accept_stage(HostIntegrator::CompressibleRK3, 0,
                                       0.0, recipe, rate, diagnostic));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 2,
                                         1.0, recipe, diagnostic));
    EXPECT_FALSE(sequence.accept_stage(HostIntegrator::CompressibleRK3, 2,
                                       0.0, recipe, rate, diagnostic));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 1,
                                         1.0, recipe, diagnostic));
    EXPECT_FALSE(sequence.accept_stage(HostIntegrator::AnelasticHeun, 1,
                                       0.0, recipe, rate, diagnostic));

    CompletedStepFluxLedger method_change;
    method_change.define(g.ba, g.dm);
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::CompressibleRK3, 0,
                                         1.0 / 3.0, recipe, diagnostic));
    ASSERT_TRUE(method_change.accept_stage(HostIntegrator::CompressibleRK3, 0,
                                           0.0, recipe, rate, diagnostic));
    ASSERT_TRUE(MakeAuxiliaryStageRecipe(HostIntegrator::AnelasticHeun, 1,
                                         1.0, recipe, diagnostic));
    EXPECT_FALSE(method_change.accept_stage(HostIntegrator::AnelasticHeun, 1,
                                            0.0, recipe, rate, diagnostic));
}

// Keep GPU kernels in free functions: NVCC rejects extended device lambdas nested
// in GoogleTest's private TestBody().
TEST(AuxiliaryMappedTransport, MappedMeasureIsExplicitAndRejectsInvalidMetrics)
{
    run_auxiliary_mapped_transport_MappedMeasureIsExplicitAndRejectsInvalidMetrics();
}
TEST(AuxiliaryMappedTransport, SharedLayoutPredicatesRejectIncompatibleLayouts)
{
    run_auxiliary_mapped_transport_SharedLayoutPredicatesRejectIncompatibleLayouts();
}
TEST(AuxiliaryMappedTransport, LedgerRejectsMismatchedFaceLayout)
{
    run_auxiliary_mapped_transport_LedgerRejectsMismatchedFaceLayout();
}
TEST(AuxiliaryMappedTransport, FluxRateAndIntegratedFluxHaveDistinctLayoutsAndTypes)
{
    run_auxiliary_mapped_transport_FluxRateAndIntegratedFluxHaveDistinctLayoutsAndTypes();
}
TEST(AuxiliaryMappedTransport, ComputationalMappedDivergenceMatchesIndependentArithmetic)
{
    run_auxiliary_mapped_transport_ComputationalMappedDivergenceMatchesIndependentArithmetic();
}
TEST(AuxiliaryMappedTransport, PeriodicArbitraryMappedFluxTelescopes)
{
    run_auxiliary_mapped_transport_PeriodicArbitraryMappedFluxTelescopes();
}
TEST(AuxiliaryMappedTransport, NativeScalarAdvectionParityAndMetricNegativeControls)
{
    run_auxiliary_mapped_transport_NativeScalarAdvectionParityAndMetricNegativeControls();
}
TEST(AuxiliaryMappedTransport, CanonicalNGridDiffusionTransferParity)
{
    run_auxiliary_mapped_transport_CanonicalNGridDiffusionTransferParity();
}
TEST(AuxiliaryMappedTransport, CanonicalStretchedGridDiffusionTransferParity)
{
    run_auxiliary_mapped_transport_CanonicalStretchedGridDiffusionTransferParity();
}
TEST(AuxiliaryMappedTransport, CanonicalStaticTerrainParityAndMetricNegativeControls)
{
    run_auxiliary_mapped_transport_CanonicalStaticTerrainParityAndMetricNegativeControls();
}
TEST(AuxiliaryMappedTransport, CompressibleRK3RecipeUsesAuditedStageCoefficients)
{
    run_auxiliary_mapped_transport_CompressibleRK3RecipeUsesAuditedStageCoefficients();
}
TEST(AuxiliaryMappedTransport, HeunRecipeSeparatesTrialAndWeightedFaceTime)
{
    run_auxiliary_mapped_transport_HeunRecipeSeparatesTrialAndWeightedFaceTime();
}
TEST(AuxiliaryMappedTransport, TimedInputViewsBuildIntensiveStateFromRhoInput)
{
    run_auxiliary_mapped_transport_TimedInputViewsBuildIntensiveStateFromRhoInput();
}
TEST(AuxiliaryMappedTransport, ConstantConstituentRatioSurvivesVariableDensityStages)
{
    run_auxiliary_mapped_transport_ConstantConstituentRatioSurvivesVariableDensityStages();
}
TEST(AuxiliaryMappedTransport, CompletedLedgerUsesExactHostTemporalWeights)
{
    run_auxiliary_mapped_transport_CompletedLedgerUsesExactHostTemporalWeights();
}
TEST(AuxiliaryMappedTransport, CompletedLedgerAccumulatesEveryComponent)
{
    run_auxiliary_mapped_transport_CompletedLedgerAccumulatesEveryComponent();
}
TEST(AuxiliaryMappedTransport, StageTargetMustBeDisjoint)
{
    run_auxiliary_mapped_transport_StageTargetMustBeDisjoint();
}
TEST(AuxiliaryMappedTransport, ZeroRateStillAppliesHostAnchorInputRecurrence)
{
    run_auxiliary_mapped_transport_ZeroRateStillAppliesHostAnchorInputRecurrence();
}
TEST(AuxiliaryMappedTransport, StageSequenceFailsClosed)
{
    run_auxiliary_mapped_transport_StageSequenceFailsClosed();
}

void run_mapped_donor_rate_cartesian_sum_test ()
{
    TestGrid g;
    g.detj.setVal(Real(1.0));
    g.mx.setVal(Real(1.0));
    g.my.setVal(Real(1.0));
    ASSERT_TRUE(g.build_measure());

    const BoxArray xb = amrex::convert(g.ba, IntVect::TheDimensionVector(0));
    const BoxArray yb = amrex::convert(g.ba, IntVect::TheDimensionVector(1));
    const BoxArray zb = amrex::convert(g.ba, IntVect::TheDimensionVector(2));
    MultiFab rho_u(xb, g.dm, 1, 0), rho_v(yb, g.dm, 1, 0);
    MultiFab omega(zb, g.dm, 1, 0);
    MultiFab ax(xb, g.dm, 1, 0), ay(yb, g.dm, 1, 0), az(zb, g.dm, 1, 0);
    MultiFab mf_uy(project_to_xy(xb), g.dm, 1, 0);
    MultiFab mf_vx(project_to_xy(yb), g.dm, 1, 0);
    ax.setVal(Real(1.0)); ay.setVal(Real(1.0)); az.setVal(Real(1.0));
    mf_uy.setVal(Real(1.0)); mf_vx.setVal(Real(1.0));
    const auto dx_inv = g.geom.InvCellSizeArray();
    rho_u.setVal(Real(1.0) / dx_inv[0]);
    rho_v.setVal(Real(1.0) / dx_inv[1]);
    omega.setVal(Real(1.0) / dx_inv[2]);
    MultiFab density(g.ba, g.dm, 1, 0);
    density.setVal(Real(1.0));

    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 1, 0);
    std::string diagnostic;
    ASSERT_TRUE(BuildMappedDryAirCarrierFluxRate(
        rate, rho_u, rho_v, omega, ax, ay, az, mf_uy, mf_vx, g.mx, g.my,
        diagnostic)) << diagnostic;
    Real max_rate = Real(0.0);
    ASSERT_TRUE(ComputeMaxMappedOutgoingRate(rate, g.omega, density, 0,
                                              dx_inv, max_rate, diagnostic))
        << diagnostic;

    // Each positive direction contributes one unit of outgoing rate.  The
    // host directional-max estimator would return 1, while donor demand sums
    // the three mapped directions and must return 3.
    EXPECT_NEAR(max_rate, Real(3.0), Real(32.0) *
                                      std::numeric_limits<Real>::epsilon());
}

void run_mapped_donor_rate_metric_oracle_test ()
{
    TestGrid g;
    ASSERT_TRUE(g.build_measure());
    const BoxArray xb = amrex::convert(g.ba, IntVect::TheDimensionVector(0));
    const BoxArray yb = amrex::convert(g.ba, IntVect::TheDimensionVector(1));
    const BoxArray zb = amrex::convert(g.ba, IntVect::TheDimensionVector(2));
    MultiFab rho_u(xb, g.dm, 1, 0), rho_v(yb, g.dm, 1, 0);
    MultiFab omega(zb, g.dm, 1, 0);
    MultiFab ax(xb, g.dm, 1, 0), ay(yb, g.dm, 1, 0), az(zb, g.dm, 1, 0);
    MultiFab mf_uy(project_to_xy(xb), g.dm, 1, 0);
    MultiFab mf_vx(project_to_xy(yb), g.dm, 1, 0);
    MultiFab density(g.ba, g.dm, 1, 0);
    for (amrex::MFIter mfi(rho_u); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto ru = rho_u.array(mfi);
        const auto area = ax.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            ru(i,j,k,0) = Real(0.3) + Real(0.02)*i + Real(0.01)*j + Real(0.005)*k;
            area(i,j,k,0) = Real(0.7) + Real(0.01)*i + Real(0.02)*k;
        });
    }
    for (amrex::MFIter mfi(rho_v); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto rv = rho_v.array(mfi);
        const auto area = ay.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            rv(i,j,k,0) = Real(-0.25) + Real(0.01)*j + Real(0.02)*k;
            area(i,j,k,0) = Real(0.8) + Real(0.02)*i + Real(0.01)*j;
        });
    }
    for (amrex::MFIter mfi(omega); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto o = omega.array(mfi);
        const auto area = az.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            o(i,j,k,0) = Real(0.15) + Real(0.03)*i - Real(0.02)*j + Real(0.01)*k;
            area(i,j,k,0) = Real(0.9) + Real(0.01)*i + Real(0.02)*j;
        });
    }
    for (amrex::MFIter mfi(mf_uy); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto m = mf_uy.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
            m(i,j,0,0) = Real(1.2) + Real(0.03)*j;
        });
    }
    for (amrex::MFIter mfi(mf_vx); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto m = mf_vx.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
            m(i,j,0,0) = Real(1.4) + Real(0.02)*i;
        });
    }
    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto rho = density.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            rho(i,j,k,0) = Real(1.1) + Real(0.05)*i + Real(0.03)*j + Real(0.02)*k;
        });
    }

    MappedFaceFluxRate rate;
    rate.define(g.ba, g.dm, 1, 0);
    std::string diagnostic;
    ASSERT_TRUE(BuildMappedDryAirCarrierFluxRate(
        rate, rho_u, rho_v, omega, ax, ay, az, mf_uy, mf_vx, g.mx, g.my,
        diagnostic)) << diagnostic;
    const auto dx_inv = g.geom.InvCellSizeArray();
    Real max_rate = Real(0.0);
    ASSERT_TRUE(ComputeMaxMappedOutgoingRate(rate, g.omega, density, 0,
                                              dx_inv, max_rate, diagnostic))
        << diagnostic;

    // Independent host oracle evaluates the face mapping, mapped volume,
    // density denominator, and directional outgoing sum from their formulas.
    Real expected = Real(0.0);
    for (int k = g.domain.smallEnd(2); k <= g.domain.bigEnd(2); ++k) {
        for (int j = g.domain.smallEnd(1); j <= g.domain.bigEnd(1); ++j) {
            for (int i = g.domain.smallEnd(0); i <= g.domain.bigEnd(0); ++i) {
                auto fx = [=](int fi) {
                    return (Real(0.7) + Real(0.01)*fi + Real(0.02)*k) *
                           (Real(0.3) + Real(0.02)*fi + Real(0.01)*j + Real(0.005)*k) /
                           (Real(1.2) + Real(0.03)*j);
                };
                auto fy = [=](int fj) {
                    return (Real(0.8) + Real(0.02)*i + Real(0.01)*fj) *
                           (Real(-0.25) + Real(0.01)*fj + Real(0.02)*k) /
                           (Real(1.4) + Real(0.02)*i);
                };
                auto fz = [=](int fk) {
                    return (Real(0.9) + Real(0.01)*i + Real(0.02)*j) *
                           (Real(0.15) + Real(0.03)*i - Real(0.02)*j + Real(0.01)*fk) /
                           ((Real(1.17) + Real(0.023)*i) *
                            (Real(0.83) + Real(0.017)*j + Real(0.006)*i));
                };
                const Real measure =
                    (Real(1.31) + Real(0.037)*i + Real(0.019)*j + Real(0.011)*k) /
                    ((Real(1.17) + Real(0.023)*i) *
                     (Real(0.83) + Real(0.017)*j + Real(0.006)*i));
                const Real rho = Real(1.1) + Real(0.05)*i + Real(0.03)*j + Real(0.02)*k;
                const Real outward =
                    (std::max(fx(i+1), Real(0.0)) + std::max(-fx(i), Real(0.0))) * dx_inv[0] +
                    (std::max(fy(j+1), Real(0.0)) + std::max(-fy(j), Real(0.0))) * dx_inv[1] +
                    (std::max(fz(k+1), Real(0.0)) + std::max(-fz(k), Real(0.0))) * dx_inv[2];
                expected = std::max(expected, outward / (measure * rho));
            }
        }
    }
    EXPECT_NEAR(max_rate, expected, Real(128.0) *
                                   std::numeric_limits<Real>::epsilon() * expected);
}

void run_mapped_donor_fixed_dt_bound_test ()
{
    double hard_limit = 0.0;
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(-1.0, Real(10.0), hard_limit));
    EXPECT_DOUBLE_EQ(hard_limit, 0.1);
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(0.1, Real(10.0), hard_limit));
    EXPECT_DOUBLE_EQ(hard_limit, 0.1);
    EXPECT_TRUE(FixedDtExceedsMappedDonorLimit(0.2, Real(10.0), hard_limit));
    EXPECT_DOUBLE_EQ(hard_limit, 0.1);
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(0.2, Real(0.0), hard_limit));
    EXPECT_TRUE(std::isinf(hard_limit));
}

void run_native_mapped_carrier_sloping_terrain_counterexample_test ()
{
    TestGrid g(8, 8, 8);
    g.detj.setVal(Real(1.0));
    g.mx.setVal(Real(1.0));
    g.my.setVal(Real(1.0));
    ASSERT_TRUE(g.build_measure());

    const BoxArray xb = amrex::convert(g.ba, IntVect::TheDimensionVector(0));
    const BoxArray yb = amrex::convert(g.ba, IntVect::TheDimensionVector(1));
    const BoxArray zb = amrex::convert(g.ba, IntVect::TheDimensionVector(2));
    MultiFab rho_u(xb, g.dm, 1, 2), rho_v(yb, g.dm, 1, 2);
    MultiFab rho_w(zb, g.dm, 1, 0);
    rho_u.setVal(Real(0.5), 0, 1, 2);
    rho_v.setVal(Real(0.3), 0, 1, 2);
    rho_w.setVal(Real(0.45));
    const int top = g.geom.Domain().bigEnd(2) + 1;
    for (amrex::MFIter mfi(rho_w); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto w = rho_w.array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (k == top) { w(i, j, k, 0) = Real(0.01); }
            });
    }

    MultiFab mf_ux(project_to_xy(xb), g.dm, 1, 0);
    MultiFab mf_vy(project_to_xy(yb), g.dm, 1, 0);
    mf_ux.setVal(Real(1.0));
    mf_vy.setVal(Real(1.0));
    MultiFab z_nd(amrex::convert(g.ba, IntVect::TheNodeVector()), g.dm, 1, 2);
    constexpr Real slope_x = Real(0.5);
    constexpr Real slope_y = Real(0.4);
    const auto dx = g.geom.CellSizeArray();
    for (amrex::MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
        const Box bx = mfi.fabbox();
        const auto z = z_nd.array(mfi);
        const auto cell_size = dx;
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                z(i, j, k, 0) = Real(k) * cell_size[2] +
                    slope_x * Real(i) * cell_size[0] +
                    slope_y * Real(j) * cell_size[1];
            });
    }

    MappedFaceFluxRate native_rate;
    native_rate.define(g.ba, g.dm, 1, 0);
    std::string diagnostic;
    ASSERT_TRUE(CopyNativeMappedDryAirCarrierFluxRate(
        native_rate, rho_u, rho_v, rho_w, diagnostic)) << diagnostic;
    EXPECT_NEAR(rho_w.min(0), Real(0.01), Real(1.0e-6));
    EXPECT_NEAR(rho_w.max(0), Real(0.45), Real(1.0e-6));
    EXPECT_EQ(max_component_difference(native_rate.dir(0), 0, rho_u, 0),
              Real(0.0));
    EXPECT_EQ(max_component_difference(native_rate.dir(1), 0, rho_v, 0),
              Real(0.0));
    EXPECT_EQ(max_component_difference(native_rate.dir(2), 0, rho_w, 0),
              Real(0.0));

    const int i = 3;
    const int j = 3;
    const int k = 3;
    const auto inv_dx = g.geom.InvCellSizeArray();
    MultiFab omega_result(g.ba, g.dm, 1, 0);
    omega_result.setVal(Real(0.0));
    const IntVect omega_cell(i, j, k);
    // OmegaFromW is device-only; evaluate it on the execution backend.
    for (amrex::MFIter mfi(omega_result); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(omega_cell)) {
            const Box point_box(omega_cell, omega_cell);
            const auto result = omega_result.array(mfi);
            const auto w = rho_w.const_array(mfi);
            const auto u = rho_u.const_array(mfi);
            const auto v = rho_v.const_array(mfi);
            const auto map_u = mf_ux.const_array(mfi);
            const auto map_v = mf_vy.const_array(mfi);
            const auto z = z_nd.const_array(mfi);
            amrex::ParallelFor(
                point_box, [=] AMREX_GPU_DEVICE(int ii, int jj, int kk) noexcept {
                    int oi = ii;
                    int oj = jj;
                    int ok = kk;
                    result(ii, jj, kk, 0) = OmegaFromW(
                        oi, oj, ok, w(oi, oj, ok, 0), u, v, map_u, map_v, z,
                        inv_dx);
                });
        }
    }
    const Real omega = omega_result.sum(0);
    EXPECT_NEAR(omega, Real(0.08), Real(1.0e-5));
    EXPECT_GT(amrex::Math::abs(Real(0.45) - omega), Real(0.2));

    // Negative control: reproduce the fitted-terrain carrier construction
    // from 0045b83acef6e7f5d8c711c7ebba8b2d0cca80de. That implementation kept
    // the native carrier at the top face, zeroed the bottom face, and used
    // OmegaFromW on every interior face.
    MappedFaceFluxRate historical_rate;
    historical_rate.define(g.ba, g.dm, 1, 0);
    MultiFab::Copy(historical_rate.dir(0), rho_u, 0, 0, 1, 0);
    MultiFab::Copy(historical_rate.dir(1), rho_v, 0, 0, 1, 0);
    const int bottom = g.geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(historical_rate.dir(2)); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto rw = rho_w.const_array(mfi);
        const auto ru = rho_u.const_array(mfi);
        const auto rv = rho_v.const_array(mfi);
        const auto ux = mf_ux.const_array(mfi);
        const auto vy = mf_vy.const_array(mfi);
        const auto z = z_nd.const_array(mfi);
        const auto out = historical_rate.dir(2).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int ii, int jj, int kk) noexcept {
                if (kk == bottom) {
                    out(ii, jj, kk, 0) = Real(0.0);
                } else if (kk == top) {
                    out(ii, jj, kk, 0) = rw(ii, jj, kk, 0);
                } else {
                    int oi = ii;
                    int oj = jj;
                    int ok = kk;
                    out(ii, jj, kk, 0) = OmegaFromW(
                        oi, oj, ok, rw(oi, oj, ok, 0), ru, rv, ux, vy, z,
                        inv_dx);
                }
            });
    }
    EXPECT_GT(max_component_difference(native_rate.dir(2), 0,
                                       historical_rate.dir(2), 0),
              Real(0.2));

    MultiFab density(g.ba, g.dm, 1, 0);
    density.setVal(Real(1.0));
    Real native_max_rate = Real(0.0);
    Real historical_max_rate = Real(0.0);
    ASSERT_TRUE(ComputeMaxMappedOutgoingRate(
        native_rate, g.omega, density, 0, inv_dx, native_max_rate,
        diagnostic)) << diagnostic;
    ASSERT_TRUE(ComputeMaxMappedOutgoingRate(
        historical_rate, g.omega, density, 0, inv_dx, historical_max_rate,
        diagnostic)) << diagnostic;
    ASSERT_TRUE(std::isfinite(native_max_rate));
    ASSERT_TRUE(std::isfinite(historical_max_rate));
    ASSERT_GT(native_max_rate, Real(0.0));
    ASSERT_GT(historical_max_rate, Real(0.0));
    EXPECT_NEAR(native_max_rate, Real(10.0), Real(1.0e-5));
    EXPECT_NEAR(historical_max_rate, Real(7.04), Real(1.0e-5));
    EXPECT_GT(amrex::Math::abs(native_max_rate - historical_max_rate) /
                  amrex::max(native_max_rate, historical_max_rate),
              Real(0.1));

    double native_hard_limit = 0.0;
    double historical_hard_limit = 0.0;
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(
        1.0e-12, native_max_rate, native_hard_limit));
    EXPECT_FALSE(FixedDtExceedsMappedDonorLimit(
        1.0e-12, historical_max_rate, historical_hard_limit));
    ASSERT_NE(native_hard_limit, historical_hard_limit);
    EXPECT_NEAR(native_hard_limit, 0.1, 1.0e-5);
    EXPECT_NEAR(historical_hard_limit, 1.0 / 7.04, 1.0e-5);
    const double between_limits =
        0.5 * (native_hard_limit + historical_hard_limit);
    const bool native_unsafe = FixedDtExceedsMappedDonorLimit(
        between_limits, native_max_rate, native_hard_limit);
    const bool historical_unsafe = FixedDtExceedsMappedDonorLimit(
        between_limits, historical_max_rate, historical_hard_limit);
    EXPECT_TRUE(native_unsafe);
    EXPECT_FALSE(historical_unsafe);

    MultiFab wrong_layout(yb, g.dm, 1, 0);
    EXPECT_FALSE(CopyNativeMappedDryAirCarrierFluxRate(
        native_rate, rho_u, rho_v, wrong_layout, diagnostic));
    EXPECT_FALSE(diagnostic.empty());

    rho_u.setVal(std::numeric_limits<Real>::quiet_NaN());
    EXPECT_FALSE(CopyNativeMappedDryAirCarrierFluxRate(
        native_rate, rho_u, rho_v, rho_w, diagnostic));
    EXPECT_NE(diagnostic.find("nonfinite"), std::string::npos);
}

TEST(AuxiliaryMappedTransport, MappedDonorRateSumsCartesianDirections)
{
    run_mapped_donor_rate_cartesian_sum_test();
}
TEST(AuxiliaryMappedTransport, MappedDonorRateMatchesMappedMetricOracle)
{
    run_mapped_donor_rate_metric_oracle_test();
}
TEST(AuxiliaryMappedTransport, FixedDtUsesHardMappedDonorLimit)
{
    run_mapped_donor_fixed_dt_bound_test();
}
TEST(AuxiliaryMappedTransport, NativeCarrierPreservesRhoWOnSlopingTerrain)
{
    run_native_mapped_carrier_sloping_terrain_counterexample_test();
}

} // namespace
