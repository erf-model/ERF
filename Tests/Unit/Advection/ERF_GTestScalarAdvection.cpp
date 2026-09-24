#include <algorithm>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_Advection.H>
#include <ERF_AdvectionSrcForScalars.H>

#include <gtest/gtest.h>

#include "../Utils/Interpolation/ERF_GTestInterpolationCommon.H"

using namespace interpolation_test;

namespace {

using amrex::Array4;
using amrex::Box;
using amrex::FArrayBox;
using amrex::GpuArray;
using amrex::IntVect;
using amrex::Real;

constexpr int kScalarComp = 2;
constexpr int kFluxComp = 1;
constexpr int kRhsComp = 4;
constexpr int kNumComponents = 6;
constexpr int kNumFluxComponents = 3;

Box test_cell_box ()
{
    return Box(IntVect(0,0,0), IntVect(2,2,2));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real scalar_value (const int i, const int j, const int k) noexcept
{
    return Real(0.7) + Real(0.13)*i - Real(0.09)*j + Real(0.11)*k
         + Real(0.017)*i*j - Real(0.012)*i*k + Real(0.023)*j*k
         + Real(0.002)*i*j*k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real momentum_value (const int i, const int j, const int k, const int dir) noexcept
{
    if (dir == 0) {
        return Real(-0.71) + Real(0.11)*i + Real(0.08)*j - Real(0.04)*k;
    } else if (dir == 1) {
        return Real(-0.31) + Real(0.06)*i + Real(0.07)*j - Real(0.025)*k;
    }
    return Real(0.45) - Real(0.04)*i + Real(0.035)*j + Real(0.09)*k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real component_sentinel (const int comp) noexcept
{
    return Real(120.0) + Real(31.0)*comp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real flux_sentinel (const int comp) noexcept
{
    return Real(5000.0) + Real(137.0)*comp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rhs_sentinel (const int comp) noexcept
{
    return Real(-3000.0) - Real(53.0)*comp;
}

void fill_scalar_components (FArrayBox& fab, const int scalar_comp)
{
    auto arr = fab.array();
    amrex::ParallelFor(fab.box(), fab.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        arr(i,j,k,n) = (n == scalar_comp) ? scalar_value(i,j,k) : component_sentinel(n);
    });
}

void fill_flux_sentinels (FArrayBox& fab)
{
    auto arr = fab.array();
    amrex::ParallelFor(fab.box(), fab.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        arr(i,j,k,n) = flux_sentinel(n);
    });
}

void fill_rhs_sentinels (FArrayBox& fab)
{
    auto arr = fab.array();
    amrex::ParallelFor(fab.box(), fab.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        arr(i,j,k,n) = rhs_sentinel(n);
    });
}

void fill_momentum (FArrayBox& fab, const int dir)
{
    auto arr = fab.array();
    amrex::ParallelFor(fab.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        arr(i,j,k,0) = momentum_value(i,j,k,dir);
    });
}

void fill_metrics (FArrayBox& detJ, FArrayBox& mf_mx, FArrayBox& mf_my,
                   const bool include_nonpositive_detJ)
{
    auto det_arr = detJ.array();
    auto mx_arr = mf_mx.array();
    auto my_arr = mf_my.array();
    amrex::ParallelFor(detJ.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        det_arr(i,j,k,0) = include_nonpositive_detJ && i == 2 && j == 2 && k == 2
                         ? Real(0.0)
                         : Real(1.4) + Real(0.06)*i + Real(0.04)*j + Real(0.03)*k;
        mx_arr(i,j,k,0) = Real(0.91) + Real(0.03)*i + Real(0.02)*j;
        my_arr(i,j,k,0) = Real(1.17) + Real(0.025)*j + Real(0.015)*k;
    });
}

void copy_to_host (const FArrayBox& device_fab, FArrayBox& host_fab)
{
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_fab.dataPtr(0), device_fab.dataPtr(0) + device_fab.size(),
                     host_fab.dataPtr(0));
}

GpuArray<const Array4<Real>, AMREX_SPACEDIM>
make_flux_views (FArrayBox& xflux, FArrayBox& yflux, FArrayBox& zflux)
{
    return {xflux.array(), yflux.array(), zflux.array()};
}

void expect_close (const Real actual, const Real expected)
{
    EXPECT_LE(normalized_error(actual, expected, kDeviceRelTol), Real(1.0))
        << "actual=" << static_cast<double>(actual)
        << " expected=" << static_cast<double>(expected);
}

} // namespace

// Motivation: The reusable scalar operator must not infer the source or
// destination component from ERF's native cons layout. Deliberately unrelated
// source, flux, and RHS component numbers, together with spatially varying
// non-unit map factors and detJ, make an accidental rhs_comp-1 lookup,
// component-0 flux write, or omitted mapped-divergence factor observable.
TEST(ScalarAdvectionPrimitives, ExplicitComponentsAndMappedDivergenceAreIndependent)
{
    const Box bx = test_cell_box();
    Box scalar_box = bx;
    scalar_box.grow(1);
    FArrayBox scalar(scalar_box, kNumComponents);
    FArrayBox xmom(amrex::surroundingNodes(bx,0), 1);
    FArrayBox ymom(amrex::surroundingNodes(bx,1), 1);
    FArrayBox zmom(amrex::surroundingNodes(bx,2), 1);
    FArrayBox xflux(amrex::surroundingNodes(bx,0), kNumFluxComponents);
    FArrayBox yflux(amrex::surroundingNodes(bx,1), kNumFluxComponents);
    FArrayBox zflux(amrex::surroundingNodes(bx,2), kNumFluxComponents);
    FArrayBox rhs(bx, kNumComponents);
    FArrayBox detJ(bx, 1);
    FArrayBox mf_mx(bx, 1);
    FArrayBox mf_my(bx, 1);

    fill_scalar_components(scalar, kScalarComp);
    fill_momentum(xmom, 0);
    fill_momentum(ymom, 1);
    fill_momentum(zmom, 2);
    fill_flux_sentinels(xflux);
    fill_flux_sentinels(yflux);
    fill_flux_sentinels(zflux);
    fill_rhs_sentinels(rhs);
    fill_metrics(detJ, mf_mx, mf_my, true);

    const auto flux_views = make_flux_views(xflux, yflux, zflux);
    const GpuArray<Real, AMREX_SPACEDIM> cellSizeInv{{Real(0.5), Real(0.8), Real(1.2)}};
    BuildScalarAdvectionFluxesCentered2(bx, scalar.const_array(), kScalarComp, flux_views,
                                        kFluxComp, xmom.const_array(), ymom.const_array(),
                                        zmom.const_array());
    ApplyScalarAdvectionFluxDivergence(bx, flux_views, kFluxComp, rhs.array(), kRhsComp,
                                       detJ.const_array(), cellSizeInv,
                                       mf_mx.const_array(), mf_my.const_array());
    gpu_sync();

    FArrayBox host_xflux(xflux.box(), xflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_yflux(yflux.box(), yflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_zflux(zflux.box(), zflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_rhs(rhs.box(), rhs.nComp(), amrex::The_Pinned_Arena());
    copy_to_host(xflux, host_xflux);
    copy_to_host(yflux, host_yflux);
    copy_to_host(zflux, host_zflux);
    copy_to_host(rhs, host_rhs);

    const auto hx = host_xflux.const_array();
    const auto hy = host_yflux.const_array();
    const auto hz = host_zflux.const_array();
    const auto hrhs = host_rhs.const_array();
    const Real dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    for (int k = host_xflux.box().smallEnd(2); k <= host_xflux.box().bigEnd(2); ++k) {
        for (int j = host_xflux.box().smallEnd(1); j <= host_xflux.box().bigEnd(1); ++j) {
            for (int i = host_xflux.box().smallEnd(0); i <= host_xflux.box().bigEnd(0); ++i) {
                const Real expected = momentum_value(i,j,k,0) * Real(0.5) *
                    (scalar_value(i,j,k) + scalar_value(i-1,j,k));
                expect_close(hx(i,j,k,kFluxComp), expected);
                for (int n = 0; n < kNumFluxComponents; ++n) {
                    if (n != kFluxComp) { EXPECT_EQ(hx(i,j,k,n), flux_sentinel(n)); }
                }
            }
        }
    }
    for (int k = host_yflux.box().smallEnd(2); k <= host_yflux.box().bigEnd(2); ++k) {
        for (int j = host_yflux.box().smallEnd(1); j <= host_yflux.box().bigEnd(1); ++j) {
            for (int i = host_yflux.box().smallEnd(0); i <= host_yflux.box().bigEnd(0); ++i) {
                const Real expected = momentum_value(i,j,k,1) * Real(0.5) *
                    (scalar_value(i,j,k) + scalar_value(i,j-1,k));
                expect_close(hy(i,j,k,kFluxComp), expected);
                for (int n = 0; n < kNumFluxComponents; ++n) {
                    if (n != kFluxComp) { EXPECT_EQ(hy(i,j,k,n), flux_sentinel(n)); }
                }
            }
        }
    }
    for (int k = host_zflux.box().smallEnd(2); k <= host_zflux.box().bigEnd(2); ++k) {
        for (int j = host_zflux.box().smallEnd(1); j <= host_zflux.box().bigEnd(1); ++j) {
            for (int i = host_zflux.box().smallEnd(0); i <= host_zflux.box().bigEnd(0); ++i) {
                const Real expected = momentum_value(i,j,k,2) * Real(0.5) *
                    (scalar_value(i,j,k) + scalar_value(i,j,k-1));
                expect_close(hz(i,j,k,kFluxComp), expected);
                for (int n = 0; n < kNumFluxComponents; ++n) {
                    if (n != kFluxComp) { EXPECT_EQ(hz(i,j,k,n), flux_sentinel(n)); }
                }
            }
        }
    }

    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                const Real det = (i == 2 && j == 2 && k == 2)
                               ? Real(0.0)
                               : Real(1.4) + Real(0.06)*i + Real(0.04)*j + Real(0.03)*k;
                Real expected = Real(0.0);
                if (det > Real(0.0)) {
                    const Real mx = Real(0.91) + Real(0.03)*i + Real(0.02)*j;
                    const Real my = Real(1.17) + Real(0.025)*j;
                    const Real div =
                        (momentum_value(i+1,j,k,0) * Real(0.5) * (scalar_value(i+1,j,k)+scalar_value(i,j,k)) -
                         momentum_value(i,j,k,0) * Real(0.5) * (scalar_value(i,j,k)+scalar_value(i-1,j,k))) * dxInv +
                        (momentum_value(i,j+1,k,1) * Real(0.5) * (scalar_value(i,j+1,k)+scalar_value(i,j,k)) -
                         momentum_value(i,j,k,1) * Real(0.5) * (scalar_value(i,j,k)+scalar_value(i,j-1,k))) * dyInv +
                        (momentum_value(i,j,k+1,2) * Real(0.5) * (scalar_value(i,j,k+1)+scalar_value(i,j,k)) -
                         momentum_value(i,j,k,2) * Real(0.5) * (scalar_value(i,j,k)+scalar_value(i,j,k-1))) * dzInv;
                    expected = -(Real(1.0) / det) * (mx * my) * div;
                }
                expect_close(hrhs(i,j,k,kRhsComp), expected);
                for (int n = 0; n < kNumComponents; ++n) {
                    if (n != kRhsComp) { EXPECT_EQ(hrhs(i,j,k,n), rhs_sentinel(n)); }
                }
            }
        }
    }
}

// Motivation: The higher-order scalar path previously derived its primitive
// component from a conserved-state index. Distinct neighboring components
// verify that reconstruction reads only the explicitly requested scalar
// component and writes only the explicitly requested flux component.
TEST(ScalarAdvectionPrimitives, HigherOrderFluxUsesRequestedScalarAndFluxComponents)
{
    const Box bx = test_cell_box();
    Box scalar_box = bx;
    scalar_box.grow(2);
    FArrayBox scalar(scalar_box, kNumComponents);
    FArrayBox xmom(amrex::surroundingNodes(bx,0), 1);
    FArrayBox ymom(amrex::surroundingNodes(bx,1), 1);
    FArrayBox zmom(amrex::surroundingNodes(bx,2), 1);
    FArrayBox xflux(amrex::surroundingNodes(bx,0), kNumFluxComponents);
    FArrayBox yflux(amrex::surroundingNodes(bx,1), kNumFluxComponents);
    FArrayBox zflux(amrex::surroundingNodes(bx,2), kNumFluxComponents);
    fill_scalar_components(scalar, kScalarComp);
    fill_momentum(xmom, 0);
    fill_momentum(ymom, 1);
    fill_momentum(zmom, 2);
    fill_flux_sentinels(xflux);
    fill_flux_sentinels(yflux);
    fill_flux_sentinels(zflux);

    const auto flux_views = make_flux_views(xflux, yflux, zflux);
    BuildScalarAdvectionFluxes(bx, scalar.const_array(), kScalarComp, flux_views, kFluxComp,
                               xmom.const_array(), ymom.const_array(), zmom.const_array(),
                               AdvType::Centered_4th, AdvType::Centered_4th,
                               Real(0.0), Real(0.0));
    gpu_sync();

    FArrayBox host_xflux(xflux.box(), xflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_yflux(yflux.box(), yflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_zflux(zflux.box(), zflux.nComp(), amrex::The_Pinned_Arena());
    copy_to_host(xflux, host_xflux);
    copy_to_host(yflux, host_yflux);
    copy_to_host(zflux, host_zflux);

    const auto check_centered4 = [](const FArrayBox& flux, const int dir) {
        const auto arr = flux.const_array();
        for (int k = flux.box().smallEnd(2); k <= flux.box().bigEnd(2); ++k) {
            for (int j = flux.box().smallEnd(1); j <= flux.box().bigEnd(1); ++j) {
                for (int i = flux.box().smallEnd(0); i <= flux.box().bigEnd(0); ++i) {
                    int im = i, jm = j, km = k;
                    int ip = i, jp = j, kp = k;
                    if (dir == 0) { --im; ++ip; }
                    if (dir == 1) { --jm; ++jp; }
                    if (dir == 2) { --km; ++kp; }
                    int imm = im, jmm = jm, kmm = km;
                    if (dir == 0) { --imm; }
                    if (dir == 1) { --jmm; }
                    if (dir == 2) { --kmm; }
                    const Real q_face = (Real(7.0)/Real(12.0)) *
                        (scalar_value(i,j,k) + scalar_value(im,jm,km)) -
                        (Real(1.0)/Real(12.0)) *
                        (scalar_value(ip,jp,kp) + scalar_value(imm,jmm,kmm));
                    const Real expected = momentum_value(i,j,k,dir) * q_face;
                    expect_close(arr(i,j,k,kFluxComp), expected);
                    for (int n = 0; n < kNumFluxComponents; ++n) {
                        if (n != kFluxComp) { EXPECT_EQ(arr(i,j,k,n), flux_sentinel(n)); }
                    }
                }
            }
        }
    };
    check_centered4(host_xflux, 0);
    check_centered4(host_yflux, 1);
    check_centered4(host_zflux, 2);
}

// Motivation: Native ERF stores primitive scalar n one component below its
// conserved counterpart, and stores the face flux in the component matching
// the conserved counterpart, because the flux registers reflux state
// component n from flux component n. The public native adapter must preserve
// both mappings while delegating numerical work to the component-explicit
// scalar primitives, and must not clobber the flux components belonging to
// other conserved variables (in particular component 0, which carries the
// density flux built just before the scalars in the pre-RHS).
TEST(ScalarAdvectionPrimitives, NativeAdapterMatchesExplicitPrimitiveMapping)
{
    const Box bx = test_cell_box();
    Box scalar_box = bx;
    scalar_box.grow(1);
    FArrayBox scalar(scalar_box, kNumComponents + 1);
    FArrayBox xmom(amrex::surroundingNodes(bx,0), 1);
    FArrayBox ymom(amrex::surroundingNodes(bx,1), 1);
    FArrayBox zmom(amrex::surroundingNodes(bx,2), 1);
    FArrayBox native_xflux(amrex::surroundingNodes(bx,0), kNumComponents);
    FArrayBox native_yflux(amrex::surroundingNodes(bx,1), kNumComponents);
    FArrayBox native_zflux(amrex::surroundingNodes(bx,2), kNumComponents);
    FArrayBox explicit_xflux(amrex::surroundingNodes(bx,0), kNumComponents);
    FArrayBox explicit_yflux(amrex::surroundingNodes(bx,1), kNumComponents);
    FArrayBox explicit_zflux(amrex::surroundingNodes(bx,2), kNumComponents);
    FArrayBox native_rhs(bx, kNumComponents + 1);
    FArrayBox explicit_rhs(bx, kNumComponents + 1);
    FArrayBox detJ(bx, 1);
    FArrayBox mf_mx(bx, 1);
    FArrayBox mf_my(bx, 1);

    fill_scalar_components(scalar, kScalarComp);
    fill_momentum(xmom, 0);
    fill_momentum(ymom, 1);
    fill_momentum(zmom, 2);
    fill_flux_sentinels(native_xflux);
    fill_flux_sentinels(native_yflux);
    fill_flux_sentinels(native_zflux);
    fill_flux_sentinels(explicit_xflux);
    fill_flux_sentinels(explicit_yflux);
    fill_flux_sentinels(explicit_zflux);
    fill_rhs_sentinels(native_rhs);
    fill_rhs_sentinels(explicit_rhs);
    fill_metrics(detJ, mf_mx, mf_my, false);

    const int cons_comp = kScalarComp + 1;
    GpuArray<const Array4<Real>, AMREX_SPACEDIM> native_flux =
        make_flux_views(native_xflux, native_yflux, native_zflux);
    GpuArray<const Array4<Real>, AMREX_SPACEDIM> explicit_flux =
        make_flux_views(explicit_xflux, explicit_yflux, explicit_zflux);
    GpuArray<amrex::BCRec, BCVars::NumTypes> closed_bcs{};
    for (auto& bc : closed_bcs) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bc.setLo(dir, ERFBCType::foextrap);
            bc.setHi(dir, ERFBCType::foextrap);
        }
    }
    const GpuArray<Real, AMREX_SPACEDIM> cellSizeInv{{Real(0.5), Real(0.8), Real(1.2)}};

    AdvectionSrcForScalars(bx, cons_comp, 1,
                           xmom.const_array(), ymom.const_array(), zmom.const_array(),
                           scalar.const_array(), native_rhs.array(), detJ.const_array(),
                           cellSizeInv, mf_mx.const_array(), mf_my.const_array(),
                           AdvType::Centered_2nd, AdvType::Centered_2nd,
                           Real(0.0), Real(0.0), native_flux, bx, closed_bcs.data());
    BuildScalarAdvectionFluxesCentered2(bx, scalar.const_array(), kScalarComp, explicit_flux,
                                        cons_comp, xmom.const_array(), ymom.const_array(), zmom.const_array());
    ApplyScalarAdvectionFluxDivergence(bx, explicit_flux, cons_comp, explicit_rhs.array(), cons_comp,
                                       detJ.const_array(), cellSizeInv,
                                       mf_mx.const_array(), mf_my.const_array());
    gpu_sync();

    FArrayBox host_native_x(native_xflux.box(), native_xflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_native_y(native_yflux.box(), native_yflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_native_z(native_zflux.box(), native_zflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_explicit_x(explicit_xflux.box(), explicit_xflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_explicit_y(explicit_yflux.box(), explicit_yflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_explicit_z(explicit_zflux.box(), explicit_zflux.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_native_rhs(native_rhs.box(), native_rhs.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_explicit_rhs(explicit_rhs.box(), explicit_rhs.nComp(), amrex::The_Pinned_Arena());
    copy_to_host(native_xflux, host_native_x);
    copy_to_host(native_yflux, host_native_y);
    copy_to_host(native_zflux, host_native_z);
    copy_to_host(explicit_xflux, host_explicit_x);
    copy_to_host(explicit_yflux, host_explicit_y);
    copy_to_host(explicit_zflux, host_explicit_z);
    copy_to_host(native_rhs, host_native_rhs);
    copy_to_host(explicit_rhs, host_explicit_rhs);

    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                expect_close(host_native_rhs.const_array()(i,j,k,cons_comp),
                             host_explicit_rhs.const_array()(i,j,k,cons_comp));
            }
        }
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const FArrayBox* native_flux_fab[] = {&host_native_x, &host_native_y, &host_native_z};
        const FArrayBox* explicit_flux_fab[] = {&host_explicit_x, &host_explicit_y, &host_explicit_z};
        const auto native_arr = native_flux_fab[dir]->const_array();
        const auto explicit_arr = explicit_flux_fab[dir]->const_array();
        const Box flux_box = native_flux_fab[dir]->box();
        for (int k = flux_box.smallEnd(2); k <= flux_box.bigEnd(2); ++k) {
            for (int j = flux_box.smallEnd(1); j <= flux_box.bigEnd(1); ++j) {
                for (int i = flux_box.smallEnd(0); i <= flux_box.bigEnd(0); ++i) {
                    expect_close(native_arr(i,j,k,cons_comp), explicit_arr(i,j,k,cons_comp));
                    // Every other flux component -- component 0 above all --
                    // must be left exactly as the caller supplied it.
                    for (int n = 0; n < native_flux_fab[dir]->nComp(); ++n) {
                        if (n == cons_comp) { continue; }
                        expect_close(native_arr(i,j,k,n), flux_sentinel(n));
                    }
                }
            }
        }
    }
}

// Motivation: Lateral open-boundary scalar advection previously recovered the
// primitive component as cons_index-1. A reusable scalar view must preserve
// the existing boundary formula while reading and writing independently
// specified scalar and RHS components. The native adapter must also preserve
// ERF's cons-to-primitive mapping when it delegates to that generic operation.
TEST(ScalarAdvectionPrimitives, OpenBoundaryUsesExplicitScalarAndRhsComponents)
{
    const Box domain = test_cell_box();
    const Box bx(IntVect(0,1,1), IntVect(0,1,1));
    Box scalar_box = domain;
    scalar_box.grow(1);
    FArrayBox scalar(scalar_box, kNumComponents);
    FArrayBox xmom(amrex::surroundingNodes(domain,0), 1);
    FArrayBox ymom(amrex::surroundingNodes(domain,1), 1);
    FArrayBox zmom(amrex::surroundingNodes(domain,2), 1);
    FArrayBox rhs(domain, kNumComponents);
    FArrayBox native_rhs(domain, kNumComponents);
    FArrayBox detJ(domain, 1);
    FArrayBox dummy_mx(domain, 1);
    FArrayBox dummy_my(domain, 1);

    fill_scalar_components(scalar, kScalarComp);
    fill_momentum(xmom, 0);
    fill_momentum(ymom, 1);
    fill_momentum(zmom, 2);
    fill_rhs_sentinels(rhs);
    fill_rhs_sentinels(native_rhs);
    fill_metrics(detJ, dummy_mx, dummy_my, false);

    const GpuArray<Real, AMREX_SPACEDIM> cellSizeInv{{Real(0.45), Real(0.8), Real(1.15)}};
    AdvectionSrcForOpenBC_Tangent_Scalars(bx, OpenSide::lo, OpenSide::none,
                                          kScalarComp, kRhsComp, 1,
                                          rhs.array(), scalar.const_array(),
                                          xmom.const_array(), ymom.const_array(), zmom.const_array(),
                                          detJ.const_array(), cellSizeInv);
    const int cons_comp = kScalarComp + 1;
    AdvectionSrcForOpenBC_Tangent_Cons(bx, OpenSide::lo, OpenSide::none,
                                       cons_comp, 1,
                                       native_rhs.array(), scalar.const_array(),
                                       xmom.const_array(), ymom.const_array(), zmom.const_array(),
                                       detJ.const_array(), cellSizeInv);
    gpu_sync();

    FArrayBox host_rhs(rhs.box(), rhs.nComp(), amrex::The_Pinned_Arena());
    FArrayBox host_native_rhs(native_rhs.box(), native_rhs.nComp(), amrex::The_Pinned_Arena());
    copy_to_host(rhs, host_rhs);
    copy_to_host(native_rhs, host_native_rhs);
    const auto actual_rhs = host_rhs.const_array();
    const auto actual_native_rhs = host_native_rhs.const_array();
    const int i = 0, j = 1, k = 1;
    const Real dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    const Real mom_hi = momentum_value(i+1,j,k,0);
    const Real mom_lo = momentum_value(i,j,k,0);
    const Real mom_cc = Real(0.5) * (mom_hi + mom_lo);
    const Real mom_star = -std::max(-mom_cc, Real(0.0));
    const Real mom_grad = (mom_hi - mom_lo) * dxInv;
    const Real scalar_grad = (scalar_value(i+1,j,k) - scalar_value(i,j,k)) * dxInv;
    const Real x_src = mom_star * scalar_grad + scalar_value(i,j,k) * mom_grad;

    const Real y_flux_lo = momentum_value(i,j,k,1) * Real(0.5) *
        (scalar_value(i,j,k) + scalar_value(i,j-1,k));
    const Real y_flux_hi = momentum_value(i,j+1,k,1) * Real(0.5) *
        (scalar_value(i,j+1,k) + scalar_value(i,j,k));
    const Real z_flux_lo = momentum_value(i,j,k,2) * Real(0.5) *
        (scalar_value(i,j,k) + scalar_value(i,j,k-1));
    const Real z_flux_hi = momentum_value(i,j,k+1,2) * Real(0.5) *
        (scalar_value(i,j,k+1) + scalar_value(i,j,k));
    const Real y_src = (y_flux_hi - y_flux_lo) * dyInv;
    const Real z_src = (z_flux_hi - z_flux_lo) * dzInv;
    const Real det = Real(1.4) + Real(0.06)*i + Real(0.04)*j + Real(0.03)*k;
    const Real expected = -(x_src + y_src + z_src) / det;

    expect_close(actual_rhs(i,j,k,kRhsComp), expected);
    expect_close(actual_native_rhs(i,j,k,cons_comp), actual_rhs(i,j,k,kRhsComp));
    for (int n = 0; n < kNumComponents; ++n) {
        if (n != kRhsComp) { EXPECT_EQ(actual_rhs(i,j,k,n), rhs_sentinel(n)); }
        if (n != cons_comp) { EXPECT_EQ(actual_native_rhs(i,j,k,n), rhs_sentinel(n)); }
    }
}
