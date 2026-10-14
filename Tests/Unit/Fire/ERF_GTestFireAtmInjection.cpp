#include <gtest/gtest.h>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Arena.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <algorithm>
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
    constexpr Real Cp_fire = 1004.64;   // the coupling's c_p
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
    constexpr Real Cp_fire = 1004.64;   // the coupling's c_p
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

// These tests exercise the production kernel, not a copy of its arithmetic.
// Keep the PR's profile, constants, density option and defaults fixed: only
// the Exner conversion of heat versus water is under test here.
void check_exner_heat_and_water(Real pressure, Real density, Real qv,
                               bool open_weights, bool additive,
                               bool legacy_density)
{
    SCOPED_TRACE(::testing::Message()
                 << "p=" << pressure << ", rho_d=" << density << ", qv=" << qv
                 << ", open=" << open_weights << ", add=" << additive
                 << ", legacy_density=" << legacy_density);

    Column c;
    // Several full-height boxes: a rank must never inspect another rank's FAB.
    c.ba.maxSize(IntVect(2, 2, NZ));
    c.dm = DistributionMapping(c.ba);
    BoxArray slabs = c.ba;
    slabs.coarsen(IntVect(1, 1, NZ));
    constexpr Real sensible = 5.0e4_rt, latent = 3.0e3_rt;
    constexpr Real cp = 1004.64_rt, lv = 2.501e6_rt, depth = 100.0_rt;
    constexpr Real open_fraction = 0.25_rt, roof = 2.0_rt * DZ;
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
    Gpu::synchronize();

    // Closed-form column weights, independently of the production helpers.
    const Real pi = std::pow(pressure / p_0, R_d / Cp_d);
    const Real retained = -std::expm1(-Real(NZ) * DZ / depth);
    const Real below = -std::expm1(-roof / depth);
    const Real normalization = retained - (1.0_rt - open_fraction) * below;
    const Real scale = legacy_density ? density : 1.0_rt;
    const Real bh = additive ? background_heat : 0.0_rt;
    const Real bw = additive ? background_water : 0.0_rt;
    const Real volume = DXY * DXY * DZ;
    Real error[5] = {};  // heat off/on, water off/on, unrelated components
    Real totals[4] = {}; // enthalpy off/on, water mass rate off/on
    int nonfinite = 0;
    Long inspected = 0;
    for (MFIter mfi(on); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        // Explicit host staging works with either device or managed arenas.
        // No host dereference of a device-resident source array is permitted.
        FArrayBox off_host(bx, NVAR_max, The_Pinned_Arena());
        FArrayBox on_host(bx, NVAR_max, The_Pinned_Arena());
        const Long count = bx.numPts() * NVAR_max;
        Gpu::copy(Gpu::deviceToHost, off[mfi].dataPtr(),
                  off[mfi].dataPtr() + count, off_host.dataPtr());
        Gpu::copy(Gpu::deviceToHost, on[mfi].dataPtr(),
                  on[mfi].dataPtr() + count, on_host.dataPtr());
        const auto a = off_host.const_array();
        const auto b = on_host.const_array();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            Real weight = std::exp(-Real(k) * DZ / depth)
                          * (-std::expm1(-DZ / depth));
            if (open_weights) {
                weight *= (k < 2 ? open_fraction : 1.0_rt) * retained / normalization;
            }
            const Real heat = scale * sensible * weight / (cp * DZ);
            const Real water = scale * latent * weight / (lv * DZ);
            const Real expected[4] = {heat, heat / pi, water, water};
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const Real actual[4] = {
                        a(i, j, k, RhoTheta_comp) - bh,
                        b(i, j, k, RhoTheta_comp) - bh,
                        a(i, j, k, RhoQ1_comp) - bw,
                        b(i, j, k, RhoQ1_comp) - bw};
                    for (int n = 0; n < 4; ++n) {
                        if (!std::isfinite(actual[n])) { nonfinite = 1; }
                        else {
                            error[n] = std::max(error[n], std::abs(actual[n] - expected[n]) / expected[n]);
                        }
                    }
                    for (int n = 0; n < NVAR_max; ++n) {
                        if (n == RhoTheta_comp || n == RhoQ1_comp) { continue; }
                        if (!std::isfinite(a(i, j, k, n)) || !std::isfinite(b(i, j, k, n))) {
                            nonfinite = 1;
                        } else {
                            error[4] = std::max(error[4], std::max(
                                std::abs(a(i, j, k, n) - background_other),
                                std::abs(b(i, j, k, n) - background_other)));
                        }
                    }
                    totals[0] += cp * pi * actual[0] * volume;
                    totals[1] += cp * pi * actual[1] * volume;
                    totals[2] += actual[2] * volume;
                    totals[3] += actual[3] * volume;
                    ++inspected;
                }
            }
        }
    }
    ParallelDescriptor::ReduceRealMax(error, 5);
    ParallelDescriptor::ReduceRealSum(totals, 4);
    ParallelDescriptor::ReduceIntMax(nonfinite);
    ParallelDescriptor::ReduceLongSum(inspected);
    const Real tolerance = (sizeof(Real) == 8) ? Real(1.0e-10) : Real(8.0e-5);
    EXPECT_EQ(nonfinite, 0);
    EXPECT_EQ(inspected, c.geom.Domain().numPts());
    EXPECT_LE(error[0], tolerance) << "heat with Exner off";
    EXPECT_LE(error[1], tolerance) << "heat with Exner on";
    EXPECT_LE(error[2], tolerance) << "water with Exner off";
    EXPECT_LE(error[3], tolerance) << "water must not acquire 1/Pi";
    EXPECT_EQ(error[4], 0.0_rt) << "unrelated sources must be unchanged";
    const Real area = Real(NX * NX) * DXY * DXY;
    const Real heat_budget = scale * sensible * retained * area;
    const Real water_budget = scale * (latent / lv) * retained * area;
    EXPECT_NEAR(totals[0], pi * heat_budget, tolerance * heat_budget);
    EXPECT_NEAR(totals[1], heat_budget, tolerance * heat_budget);
    EXPECT_NEAR(totals[2], water_budget, tolerance * water_budget);
    EXPECT_NEAR(totals[3], water_budget, tolerance * water_budget);
}

void check_exner_cases(bool open_weights, bool additive, bool legacy_density)
{
    for (Real pressure : {7.0e4_rt, 8.5e4_rt, 1.0e5_rt}) {
        for (Real density : {0.75_rt, 1.25_rt}) {
            for (Real qv : {0.0_rt, 0.015_rt}) {
                check_exner_heat_and_water(pressure, density, qv,
                                          open_weights, additive, legacy_density);
            }
        }
    }
}

} // namespace

TEST(ExnerTendency, WaterMassIndependentOfExner)
{
    check_exner_cases(false, false, false);
}

TEST(ExnerTendency, WaterMassIndependentOfExnerWithOpenFraction)
{
    check_exner_cases(true, false, false);
}

TEST(ExnerTendency, HeatAndWaterPreserveAdditiveBackground)
{
    check_exner_cases(false, true, false);
    check_exner_cases(true, true, false);
}

TEST(ExnerTendency, LegacyDensityFactorDoesNotApplyExnerToWater)
{
    // Compatibility coverage only, not an endorsement of density scaling as
    // the physical conversion of the prescribed W/m2 inputs.
    check_exner_cases(false, true, true);
    check_exner_cases(true, true, true);
}
