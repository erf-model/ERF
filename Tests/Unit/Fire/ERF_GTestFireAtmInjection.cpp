#include <gtest/gtest.h>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Math.H>
#include <AMReX_Vector.H>
#include <cmath>

#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include "ERF_FireAtmCoupling.H"
#include "ERF_FireSmokeEmission.H"

/**
 * @file ERF_GTestFireAtmInjection.cpp
 * @brief What the fire puts into the atmospheric column: with
 *        erf.fire.heat_tendency_exner the rho theta tendency carries 1 / Pi, so
 *        the enthalpy c_p Pi d(rho theta)/dt integrated over the column is the
 *        surface flux (without it the column is Pi times short, 4.5 % at 850
 *        hPa); the heat on a stretched column still integrates to the flux's
 *        share of the profile, because each cell's tendency is divided by its
 *        own thickness; the smoke source divides the mass flux by the thickness
 *        of the k = 0 cell from its nodal heights, not by the centre-to-centre
 *        spacing.
 */

using namespace amrex;

namespace {

constexpr int  NX = 4, NZ = 50;
constexpr Real DXY = 100.0, DZ = 10.0;
/// Relative round-off tolerance of the build precision
constexpr double REL = (sizeof(Real) == 8) ? 1.0e-9 : 1.0e-4;

struct Column
{
    BoxArray            ba;
    DistributionMapping dm;
    Geometry            geom;
    Column ()
    {
        Box domain(IntVect(0, 0, 0), IntVect(NX - 1, NX - 1, NZ - 1));
        ba = BoxArray(domain);
        dm = DistributionMapping(ba);
        geom = Geometry(domain, RealBox(0.0, 0.0, 0.0, NX * DXY, NX * DXY, NZ * DZ), CoordSys::cartesian, {1, 1, 0});
    }
};

} // namespace

// nvcc refuses an extended __device__ lambda whose enclosing function is
// private, as gtest's TestBody is: the bodies that launch kernels are
// namespace-scope functions, and each TEST calls its own.
namespace {
void ExnerTendency_TheColumnEnthalpyIsTheSurfaceFlux ()
{
    Column c;
    const Real q = 5.0e4_rt, alfg = 10.0_rt, p = 8.5e4_rt;
    MultiFab src_off(c.ba, c.dm, NVAR_max, 0), src_on(c.ba, c.dm, NVAR_max, 0);
    MultiFab Q(c.ba, c.dm, 1, 0), zcc(c.ba, c.dm, 1, 1), S(c.ba, c.dm, NVAR_max, 0);
    src_off.setVal(0.0_rt); src_on.setVal(0.0_rt);
    Q.setVal(0.0_rt);
    Q.setVal(q, makeSlab(c.geom.Domain(), 2, 0), 0, 1);
    S.setVal(0.0_rt);
    S.setVal(1.0_rt, Rho_comp, 1);
    S.setVal(getRhoThetagivenP(p), RhoTheta_comp, 1);
    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        auto z = zcc.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { z(i, j, k) = (k + 0.5_rt) * DZ; });
    }
    // d(rho theta)/dt = -(1/c_p) dQ/dz without the density factor; alfg = 10 m
    // puts the whole profile inside the 500 m column (1 - e^{-50})
    apply_fire_tendency_to_cc_source(src_off, Q, nullptr, zcc, nullptr, S, c.geom, alfg, 1.0_rt, false,
                                     false, false, nullptr, nullptr, false, false, false);
    apply_fire_tendency_to_cc_source(src_on,  Q, nullptr, zcc, nullptr, S, c.geom, alfg, 1.0_rt, false,
                                     false, false, nullptr, nullptr, false, false, true);

    const Real exner = std::pow(p / p_0, R_d / Cp_d);
    constexpr Real Cp_fire = Cp_d;   // the same dry-air c_p used by ERF
    EXPECT_NEAR(exner, 0.9546, 5.0e-4);
    Real err_ratio = 0.0, sum_off = 0.0, sum_on = 0.0;
    for (MFIter mfi(src_on); mfi.isValid(); ++mfi) {
        auto a = src_off.const_array(mfi, RhoTheta_comp);
        auto b = src_on.const_array(mfi, RhoTheta_comp);
        for (int k = 0; k < NZ; ++k) {
            const Real toff = a(1, 2, k), ton = b(1, 2, k);
            if (toff > 0.0_rt) { err_ratio = amrex::max(err_ratio, std::abs(ton / toff - 1.0_rt / exner)); }
            sum_off += Cp_fire * exner * toff * DZ;   // the enthalpy the plain tendency carries
            sum_on  += Cp_fire * exner * ton * DZ;
        }
    }
    EXPECT_LT(err_ratio, REL) << "every cell's tendency is divided by Pi";
    EXPECT_NEAR(sum_on,  q, 1.0e3 * REL * q) << "c_p Pi d(rho theta)/dt integrates to the flux";
    EXPECT_NEAR(sum_off, exner * q, 1.0e3 * REL * q) << "without the factor the column enthalpy is Pi q";
}
}  // namespace

TEST(ExnerTendency, TheColumnEnthalpyIsTheSurfaceFlux)
{
    ExnerTendency_TheColumnEnthalpyIsTheSurfaceFlux();
}

namespace {
void HeatInjection_UsesERFDryAirHeatCapacity ()
{
    // Q is a physical flux in W/m2. With heat_tendency_exner on, ERF's
    // c_p * Pi times its rho_d*theta_d tendency must recover that flux.
    // With the option off, c_p times the tendency represents the
    // potential-temperature-equivalent input, not physical enthalpy.
    Column c;
    const Real heat_flux = 2.0e4_rt, alfg = 45.0_rt;
    MultiFab src(c.ba, c.dm, NVAR_max, 0), Q(c.ba, c.dm, 1, 0);
    MultiFab zcc(c.ba, c.dm, 1, 1), S(c.ba, c.dm, NVAR_max, 0);
    Q.setVal(0.0_rt);
    Q.setVal(heat_flux, makeSlab(c.geom.Domain(), 2, 0), 0, 1);
    S.setVal(0.0_rt);
    S.setVal(1.0_rt, Rho_comp, 1);
    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        auto z = zcc.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = (Real(k) + 0.5_rt) * DZ;
        });
    }
    const Real retained = 1.0_rt - std::exp(-Real(NZ) * DZ / alfg);
    const Real expected = Real(NX * NX) * heat_flux * retained;
    const Real tolerance = (sizeof(Real) == 8) ? Real(3.0e-8) : Real(4.0e-5);

    for (Real pressure : {8.5e4_rt, 1.0e5_rt}) {
        S.setVal(getRhoThetagivenP(pressure), RhoTheta_comp, 1);
        const Real exner = std::pow(pressure / p_0, R_d / Cp_d);
        for (bool with_exner : {false, true}) {
            SCOPED_TRACE(::testing::Message() << "p=" << pressure
                         << " exner_enabled=" << with_exner);
            src.setVal(0.0_rt);
            apply_fire_tendency_to_cc_source(src, Q, nullptr, zcc, nullptr, S,
                c.geom, alfg, 1.0_rt, false, false, false,
                nullptr, nullptr, false, false, with_exner);
            // MultiFab::sum is a CPU/GPU and MPI collective reduction.
            const Real deposited = Cp_d * (with_exner ? exner : Real(1.0))
                                   * src.sum(RhoTheta_comp) * DZ;
            EXPECT_NEAR(deposited, expected, tolerance * expected)
                << "the prescribed physical heat flux must use ERF's Cp_d";
        }
    }
}
} // namespace

TEST(HeatInjection, UsesERFDryAirHeatCapacity)
{
    HeatInjection_UsesERFDryAirHeatCapacity();
}

namespace {
void HeatInjection_TheColumnGetsTheFluxOnAStretchedGrid (Real r)
{
    // nodal heights stretched by r per cell up to 500 m: r = 1.05 gives a
    // 2.39 m first cell, r = 1.1 a 0.43 m one, below a tenth of the mean
    // spacing (10 m). The centre-to-centre spacing used to divide the
    // tendency put 2 / (1 + r) of each cell's share into the column; a floor
    // of 0.1 dz_avg on the thickness put 0.98 of the share in at r = 1.1.
    Column c;
    const Real q = 5.0e4_rt, alfg = 100.0_rt;
    const Real d0 = 500.0_rt * (r - 1.0_rt) / (std::pow(r, Real(NZ)) - 1.0_rt);
    MultiFab src(c.ba, c.dm, NVAR_max, 0), Q(c.ba, c.dm, 1, 0), S(c.ba, c.dm, NVAR_max, 0);
    MultiFab znd(convert(c.ba, IntVect(1, 1, 1)), c.dm, 1, 1), zcc(c.ba, c.dm, 1, 1);
    src.setVal(0.0_rt); S.setVal(0.0_rt); S.setVal(1.0_rt, Rho_comp, 1);
    S.setVal(getRhoThetagivenP(1.0e5_rt), RhoTheta_comp, 1);
    Q.setVal(0.0_rt);
    Q.setVal(q, makeSlab(c.geom.Domain(), 2, 0), 0, 1);
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        auto z = znd.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = d0 * (std::pow(r, Real(k)) - 1.0_rt) / (r - 1.0_rt);
        });
    }
    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        auto z = zcc.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real lo = d0 * (std::pow(r, Real(k)) - 1.0_rt) / (r - 1.0_rt);
            const Real hi = d0 * (std::pow(r, Real(k + 1)) - 1.0_rt) / (r - 1.0_rt);
            z(i, j, k) = 0.5_rt * (lo + hi);
        });
    }
    apply_fire_tendency_to_cc_source(src, Q, nullptr, zcc, &znd, S, c.geom, alfg, 1.0_rt, false,
                                     false, false, nullptr, nullptr, false, false, false);
    constexpr Real Cp_fire = Cp_d;   // the same dry-air c_p used by ERF
    // the profile's share inside the column: z from the first cell centre, the
    // top cell closed one geometry dz above its centre (10 m)
    auto zc = [=] (int k) {
        return 0.5_rt * d0 * ((std::pow(r, Real(k)) - 1.0_rt) + (std::pow(r, Real(k + 1)) - 1.0_rt)) / (r - 1.0_rt);
    };
    const Real share = 1.0_rt - std::exp(-((zc(NZ - 1) - zc(0)) + DZ) / alfg);
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        auto t = src.const_array(mfi, RhoTheta_comp);
        Real sum = 0.0_rt;
        for (int k = 0; k < NZ; ++k) {
            const Real dz_cell = d0 * std::pow(r, Real(k));
            sum += Cp_fire * t(1, 2, k) * dz_cell;
        }
        // 1e-6 q in double precision, 1e-3 q in single (the 2 % and 4 % errors
        // of the old divisors fail both)
        const Real tol = (sizeof(Real) == 8) ? 1.0e3 * REL * q : 1.0e-3 * q;
        EXPECT_NEAR(sum, q * share, tol) << "c_p sum_k d(rho theta)/dt dz_cell is the flux's share, r = " << r;
    }
}
}  // namespace

TEST(HeatInjection, TheColumnGetsTheFluxOnAStretchedGrid)
{
    HeatInjection_TheColumnGetsTheFluxOnAStretchedGrid(1.05_rt);
    HeatInjection_TheColumnGetsTheFluxOnAStretchedGrid(1.1_rt);
}

namespace {
// erf.fire.heat_profile = bottom_cell: the whole flux, sensible and latent,
// enters the lowest cell of the column (the LSFIRE coupling), on a stretched
// grid whose first cell is thinner than the mean spacing
void HeatInjection_BottomCellTakesTheWholeFlux ()
{
    Column c;
    const Real r = 1.05_rt, q = 5.0e4_rt, ql = 2.0e4_rt, alfg = 100.0_rt;
    const Real d0 = 500.0_rt * (r - 1.0_rt) / (std::pow(r, Real(NZ)) - 1.0_rt);
    MultiFab src(c.ba, c.dm, NVAR_max, 0), Q(c.ba, c.dm, 1, 0), QL(c.ba, c.dm, 1, 0), S(c.ba, c.dm, NVAR_max, 0);
    MultiFab znd(convert(c.ba, IntVect(1, 1, 1)), c.dm, 1, 1), zcc(c.ba, c.dm, 1, 1);
    src.setVal(0.0_rt); S.setVal(0.0_rt); S.setVal(1.0_rt, Rho_comp, 1);
    S.setVal(getRhoThetagivenP(1.0e5_rt), RhoTheta_comp, 1);
    Q.setVal(0.0_rt);  Q.setVal(q,  makeSlab(c.geom.Domain(), 2, 0), 0, 1);
    QL.setVal(0.0_rt); QL.setVal(ql, makeSlab(c.geom.Domain(), 2, 0), 0, 1);
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        auto z = znd.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = d0 * (std::pow(r, Real(k)) - 1.0_rt) / (r - 1.0_rt);
        });
    }
    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        auto z = zcc.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real lo = d0 * (std::pow(r, Real(k)) - 1.0_rt) / (r - 1.0_rt);
            const Real hi = d0 * (std::pow(r, Real(k + 1)) - 1.0_rt) / (r - 1.0_rt);
            z(i, j, k) = 0.5_rt * (lo + hi);
        });
    }
    apply_fire_tendency_to_cc_source(src, Q, &QL, zcc, &znd, S, c.geom, alfg, 1.0_rt, true,
                                     false, false, nullptr, nullptr, false, false, false,
                                     /*bottom_cell=*/true);
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        auto t  = src.const_array(mfi, RhoTheta_comp);
        auto qv = src.const_array(mfi, RhoQ1_comp);
        EXPECT_NEAR(Cp_d * t(1, 2, 0) * d0, q, 1.0e-9 * q) << "the sensible flux, all in the first cell";
        EXPECT_NEAR(L_v * qv(1, 2, 0) * d0, ql, 1.0e-9 * ql) << "the vapour, latent / L_v, all in the first cell";
        for (int k = 1; k < NZ; ++k) {
            EXPECT_EQ(t(1, 2, k), 0.0_rt) << "k = " << k;
            EXPECT_EQ(qv(1, 2, k), 0.0_rt) << "k = " << k;
        }
    }
}
}  // namespace

TEST(HeatInjection, BottomCellTakesTheWholeFlux)
{
    HeatInjection_BottomCellTakesTheWholeFlux();
}

namespace {
void SmokeInjection_SourceDividesByTheCellThickness ()
{
    Column c;
    const Real heat = 1.0e5_rt, ef = 0.02_rt, hoc = 1.8e7_rt;
    MultiFab src(c.ba, c.dm, 1, 0), Q(c.ba, c.dm, 1, 0);
    Q.setVal(0.0_rt);
    Q.setVal(heat, makeSlab(c.geom.Domain(), 2, 0), 0, 1);
    // nodal heights: z(i, j, k) = 0.1 i + k (8 + 0.2 i + 0.2 j), so the k = 0
    // cell is 8 + 0.2 (i + 1/2) + 0.2 (j + 1/2) thick (the mean of its four
    // node pairs) while the uniform dz is 10 and the centre-to-centre spacing
    // of a stretched column would differ again
    MultiFab znd(convert(c.ba, IntVect(1, 1, 1)), c.dm, 1, 1);
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        auto z = znd.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = 0.1_rt * i + k * (8.0_rt + 0.2_rt * i + 0.2_rt * j);
        });
    }
    const Real flux = ef * heat / hoc;
    src.setVal(0.0_rt);
    inject_smoke_from_fire(src, Q, &znd, c.geom, ef, hoc, 0, false, 0);
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        auto s = src.const_array(mfi);
        for (int j = 0; j < NX; ++j) {
            for (int i = 0; i < NX; ++i) {
                const Real dz_cell = 8.0_rt + 0.2_rt * (i + 0.5_rt) + 0.2_rt * (j + 0.5_rt);
                EXPECT_NEAR(s(i, j, 0), flux / dz_cell, REL * flux / dz_cell) << "column " << i << "," << j;
                EXPECT_EQ(s(i, j, 1), 0.0_rt) << "only k = 0 takes the source";
            }
        }
    }
    // a first cell of 0.5 + 0.01 (i + 1/2) m, thinner than a tenth of the
    // uniform dz: the source still divides by the cell's own thickness (a
    // floor of 0.1 dz put about half of the mass in)
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        auto z = znd.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = 0.01_rt * i * k + 0.5_rt * k;
        });
    }
    src.setVal(0.0_rt);
    inject_smoke_from_fire(src, Q, &znd, c.geom, ef, hoc, 0, false, 0);
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        auto s = src.const_array(mfi);
        for (int j = 0; j < NX; ++j) {
            for (int i = 0; i < NX; ++i) {
                const Real dz_cell = 0.5_rt + 0.01_rt * (i + 0.5_rt);
                EXPECT_NEAR(s(i, j, 0), flux / dz_cell, REL * flux / dz_cell) << "thin first cell, column " << i << "," << j;
            }
        }
    }
    // without terrain: the uniform dz
    src.setVal(0.0_rt);
    inject_smoke_from_fire(src, Q, nullptr, c.geom, ef, hoc, 0, false, 0);
    EXPECT_NEAR(src.max(0), flux / DZ, REL * flux / DZ);
    EXPECT_NEAR(src.min(0), 0.0, 1.0e-30);
}
}  // namespace

TEST(SmokeInjection, SourceDividesByTheCellThickness)
{
    SmokeInjection_SourceDividesByTheCellThickness();
}

namespace {

// Check the production source kernel against independent analytic layer
// weights and column budgets. AMReX reductions work with device data and
// across MPI ranks; do not read GPU-resident MultiFabs on the host.
void check_exner_heat_and_water(Real pressure, Real density, Real qv,
                                bool open_weights, bool additive,
                                bool legacy_density)
{
    SCOPED_TRACE(::testing::Message()
                 << "p=" << pressure << ", rho_d=" << density << ", qv=" << qv
                 << ", open=" << open_weights << ", add=" << additive
                 << ", legacy_density=" << legacy_density);

    Column c;
    // Split into full-height boxes; this also permits a multi-rank run.
    c.ba.maxSize(IntVect(2, 2, NZ));
    c.dm = DistributionMapping(c.ba);
    BoxArray slabs = c.ba;
    slabs.coarsen(IntVect(1, 1, NZ));

    constexpr Real sensible = 5.0e4_rt, latent = 3.0e3_rt;
    // The coupling's conversion: ERF's Cp_d (hgopalan/ERF#521)
    constexpr Real cp = Cp_d, lv = L_v, depth = 100.0_rt;
    // Roof halfway through k=2, not just aligned to a cell face.
    constexpr Real open_fraction = 0.25_rt, roof = 2.5_rt * DZ;
    constexpr Real background_heat = 0.125_rt, background_water = 1.0e-8_rt;
    constexpr Real background_other = 0.25_rt;

    MultiFab off(c.ba, c.dm, NVAR_max, 0), on(c.ba, c.dm, NVAR_max, 0);
    MultiFab state(c.ba, c.dm, NVAR_max, 0);
    MultiFab qs(slabs, c.dm, 1, 0), ql(slabs, c.dm, 1, 0);
    MultiFab open(slabs, c.dm, 1, 1), roofs(slabs, c.dm, 1, 1);
    MultiFab zcc(c.ba, c.dm, 1, 1);
    MultiFab znd(convert(c.ba, IntVect(1, 1, 1)), c.dm, 1, 1);
    for (MultiFab* source : {&off, &on}) {
        source->setVal(background_other);
        source->setVal(background_heat, RhoTheta_comp, 1);
        source->setVal(background_water, RhoQ1_comp, 1);
    }
    state.setVal(0.0_rt);
    state.setVal(density, Rho_comp, 1);
    state.setVal(getRhoThetagivenP(pressure, qv), RhoTheta_comp, 1);
    state.setVal(density * qv, RhoQ1_comp, 1);
    qs.setVal(sensible); ql.setVal(latent);
    open.setVal(open_fraction); roofs.setVal(roof);

    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        const auto z = zcc.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = (Real(k) + 0.5_rt) * DZ;
        });
    }
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        const auto z = znd.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = Real(k) * DZ;
        });
    }
    apply_fire_tendency_to_cc_source(off, qs, &ql, zcc, &znd, state, c.geom,
        depth, 1.0_rt, true, false, additive, &open, &roofs,
        open_weights, legacy_density, false);
    apply_fire_tendency_to_cc_source(on, qs, &ql, zcc, &znd, state, c.geom,
        depth, 1.0_rt, true, false, additive, &open, &roofs,
        open_weights, legacy_density, true);

    // Closed-form weights, including the cell half-covered by the roof.
    const Real pi = std::pow(pressure / p_0, R_d / Cp_d);
    const Real decay = -std::expm1(-DZ / depth);
    const Real retained = -std::expm1(-Real(NZ) * DZ / depth);
    const Real below_two = -std::expm1(-2.0_rt * DZ / depth);
    const Real third_weight = std::exp(-2.0_rt * DZ / depth) * decay;
    const Real third_share = open_fraction + 0.5_rt * (1.0_rt - open_fraction);
    const Real normalization = retained - (1.0_rt - open_fraction) * below_two
                                        - (1.0_rt - third_share) * third_weight;
    const Real scale = legacy_density ? density : 1.0_rt;
    const Real bh = additive ? background_heat : 0.0_rt;
    const Real bw = additive ? background_water : 0.0_rt;

    // Residuals: heat off/on, water off/on, exact Exner water invariance,
    // and every unrelated source component. No host-device staging required.
    MultiFab residual(c.ba, c.dm, 6, 0);
    for (MFIter mfi(residual); mfi.isValid(); ++mfi) {
        const auto a = off.const_array(mfi);
        const auto b = on.const_array(mfi);
        const auto e = residual.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real weight = std::exp(-Real(k) * DZ / depth) * decay;
            if (open_weights) {
                const Real share = k < 2 ? open_fraction : (k == 2 ? third_share : 1.0_rt);
                weight *= share * retained / normalization;
            }
            const Real heat = scale * sensible * weight / (cp * DZ);
            const Real water = scale * latent * weight / (lv * DZ);
            e(i, j, k, 0) = a(i, j, k, RhoTheta_comp) - bh - heat;
            e(i, j, k, 1) = b(i, j, k, RhoTheta_comp) - bh - heat / pi;
            e(i, j, k, 2) = a(i, j, k, RhoQ1_comp) - bw - water;
            e(i, j, k, 3) = b(i, j, k, RhoQ1_comp) - bw - water;
            e(i, j, k, 4) = b(i, j, k, RhoQ1_comp) - a(i, j, k, RhoQ1_comp);
            Real other = 0.0_rt;
            for (int n = 0; n < NVAR_max; ++n) {
                if (n == RhoTheta_comp || n == RhoQ1_comp) { continue; }
                other = amrex::max(other, amrex::max(
                    amrex::Math::abs(a(i, j, k, n) - background_other),
                    amrex::Math::abs(b(i, j, k, n) - background_other)));
            }
            e(i, j, k, 5) = other;
        });
    }

    EXPECT_FALSE(off.contains_nan(0, NVAR_max));
    EXPECT_FALSE(on.contains_nan(0, NVAR_max));
    EXPECT_FALSE(off.contains_inf(0, NVAR_max));
    EXPECT_FALSE(on.contains_inf(0, NVAR_max));

    const Real tolerance = (sizeof(Real) == 8) ? Real(1.0e-10) : Real(8.0e-5);
    const Real peak_heat = scale * sensible * decay / (cp * DZ);
    const Real peak_water = scale * latent * decay / (lv * DZ);
    const auto err = residual.norm0(Vector<int>{0, 1, 2, 3, 4, 5});
    EXPECT_LE(err[0], tolerance * peak_heat) << "heat with Exner off";
    EXPECT_LE(err[1], tolerance * peak_heat / pi) << "heat with Exner on";
    EXPECT_LE(err[2], tolerance * peak_water) << "water with Exner off";
    EXPECT_LE(err[3], tolerance * peak_water) << "water without 1/Pi";
    EXPECT_EQ(err[4], 0.0_rt) << "Exner must leave vapor deposition identical cell by cell";
    EXPECT_EQ(err[5], 0.0_rt) << "unrelated sources must be unchanged";

    const Real volume = DXY * DXY * DZ;
    const Real area = Real(NX * NX) * DXY * DXY;
    const Real cells = Real(c.geom.Domain().numPts());
    const Real heat_budget = scale * sensible * retained * area;
    const Real water_budget = scale * (latent / lv) * retained * area;
    const Real heat_off = cp * pi * (off.sum(RhoTheta_comp) - bh * cells) * volume;
    const Real heat_on = cp * pi * (on.sum(RhoTheta_comp) - bh * cells) * volume;
    const Real water_off = (off.sum(RhoQ1_comp) - bw * cells) * volume;
    const Real water_on = (on.sum(RhoQ1_comp) - bw * cells) * volume;
    EXPECT_NEAR(heat_off, pi * heat_budget, tolerance * heat_budget);
    EXPECT_NEAR(heat_on, heat_budget, tolerance * heat_budget);
    EXPECT_NEAR(water_off, water_budget, tolerance * water_budget);
    EXPECT_NEAR(water_on, water_budget, tolerance * water_budget);
}

} // namespace

TEST(ExnerTendency, WaterMassIndependentOfExner)
{
    check_exner_heat_and_water(7.0e4_rt, 0.75_rt, 0.0_rt,   false, false, false);
    check_exner_heat_and_water(8.5e4_rt, 1.25_rt, 0.015_rt, false, false, false);
    check_exner_heat_and_water(1.0e5_rt, 1.0_rt,  0.015_rt, false, false, false);
}

TEST(ExnerTendency, WaterMassIndependentOfExnerWithOpenFraction)
{
    check_exner_heat_and_water(8.5e4_rt, 1.25_rt, 0.015_rt, true, false, false);
}

TEST(ExnerTendency, HeatAndWaterPreserveAdditiveBackground)
{
    check_exner_heat_and_water(8.5e4_rt, 0.75_rt, 0.015_rt, false, true, false);
    check_exner_heat_and_water(8.5e4_rt, 1.25_rt, 0.015_rt, true,  true, false);
}

TEST(ExnerTendency, LegacyDensityFactorDoesNotApplyExnerToWater)
{
    // Compatibility only, not endorsement of legacy density scaling for W/m2.
    check_exner_heat_and_water(8.5e4_rt, 1.25_rt, 0.015_rt, false, true, true);
    check_exner_heat_and_water(7.0e4_rt, 0.75_rt, 0.015_rt, true,  true, true);
}
