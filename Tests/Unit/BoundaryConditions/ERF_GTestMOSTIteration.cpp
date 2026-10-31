// The MOST flux kernels must return a bounded surface state with the right sign when the
// prescribed-flux problem has no solution (stable flux, wind below the fold), when the
// log-law bracket ln(zref/z0) - psi would change sign (strongly unstable, rough or calm
// points), and must leave ordinary points as they were.
#include <cmath>
#include <vector>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "BoundaryConditions/ERF_MOSTStress.H"

namespace
{

using amrex::Real;

struct Point
{
    Real zref;
    Real z0;
    Real umm;
    Real theta;
    Real t_surf; // only read by the prescribed-temperature kernels
};

struct Result
{
    Real ustar;
    Real tstar;
    Real t_surf;
    Real olen;
    Real qstar;
};

enum Slot { ZREF, Z0, UMM, TM, TVM, QVM, USTAR, TSTAR, QSTAR, TSURF, QSURF, OLEN, PBLH,
            HWAVE, LWAVE, ETA, NSLOT };

// Run one cold-start iterate_flux on a single surface point (dry air, no w*)
template <typename Kernel>
Result run (const Kernel& kernel, const Point& p)
{
    amrex::Gpu::DeviceVector<Real> buf(NSLOT, Real(0.0));
    Real* d = buf.data();
    amrex::single_task([=] AMREX_GPU_DEVICE () noexcept {
        auto slot = [d] (int s) {
            return amrex::Array4<Real>(d + s, amrex::Dim3{0,0,0}, amrex::Dim3{1,1,1}, 1);
        };
        d[ZREF]  = p.zref;
        d[Z0]    = p.z0;
        d[UMM]   = p.umm;
        d[TM]    = p.theta;
        d[TVM]   = p.theta;
        d[TSURF] = p.t_surf;
        d[USTAR] = bogus_large_value;
        d[OLEN]  = bogus_large_value;
        d[PBLH]  = Real(1000.0);
        kernel.iterate_flux(0, 0, 0, 100,
                            slot(ZREF), slot(Z0), slot(UMM), slot(TM), slot(TVM), slot(QVM),
                            slot(USTAR), amrex::Array4<Real>{}, slot(TSTAR), slot(QSTAR),
                            slot(TSURF), slot(QSURF), slot(OLEN), slot(PBLH),
                            slot(HWAVE), slot(LWAVE), slot(ETA));
    });
    std::vector<Real> h(NSLOT);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, buf.begin(), buf.end(), h.begin());
    amrex::Gpu::streamSynchronize();
    return Result{h[USTAR], h[TSTAR], h[TSURF], h[OLEN], h[QSTAR]};
}

Result run_flux (const Point& p, Real flux, Real qflux = Real(0.0))
{
    return run(surface_flux(flux, qflux, true), p);
}

Result run_temp (const Point& p)
{
    return run(surface_temp(Real(0.0), Real(0.0), false, 2, true), p);
}

} // namespace

// Motivation: the floor is 1 where the bracket is large, never exceeds the neutral value C
// for C >= 1/2, and stays at 1/2 when z0 approaches or exceeds zref (C <= 0), so u* is
// positive and finite there.
TEST(MOSTIteration, BracketFloorIsOneForLargeCAndBelowNeutralForSmallC)
{
    EXPECT_EQ(most_bracket(Real(3.0), Real(5.0)), Real(1.0));
    EXPECT_EQ(most_bracket(Real(3.0), Real(1.0)), Real(2.0));
    EXPECT_EQ(most_bracket(Real(1.4), Real(2.0)), Real(0.7));
    EXPECT_EQ(most_bracket(Real(0.5), Real(2.0)), Real(0.5));
    EXPECT_EQ(most_bracket(Real(0.5), Real(0.0)), Real(0.5));
    EXPECT_EQ(most_bracket(Real(0.0), Real(0.0)), Real(0.5));
    EXPECT_EQ(most_bracket(Real(-0.2), Real(-1.0)), Real(0.8));
    EXPECT_EQ(most_bracket(Real(-0.2), Real(0.0)), Real(0.5));
}

// Motivation: with the Businger-Dyer stable arm, the wind u* (C + beta_m z/L) / kappa that a
// fixed stable flux supports is smallest at z/L = C / (2 beta_m).
TEST(MOSTIteration, FoldIsTheMinimumOfTheSupportedWind)
{
    const similarity_funs sfuns{};
    const Real C = std::log(Real(5.0) / Real(0.1));
    EXPECT_NEAR(sfuns.calc_zeta_fold(C), C / Real(10.0), Real(1.0e-12));

    // U(u*) for F = -0.004 K m/s at zref 5 m: its minimum sits where z/L reaches the fold
    const Real a = KAPPA * CONST_GRAV * Real(0.004) * Real(5.0) / Real(300.0);
    Real best_u = Real(1.0e30);
    Real best_zeta = Real(0.0);
    for (int n = 1; n < 4000; ++n) {
        const Real ustar = Real(1.0e-4) * n;
        const Real zeta  = a / (ustar * ustar * ustar);
        const Real U     = ustar * (C - sfuns.calc_psi_m(zeta)) / KAPPA;
        if (U < best_u) { best_u = U; best_zeta = zeta; }
    }
    EXPECT_NEAR(best_zeta, sfuns.calc_zeta_fold(C), Real(5.0e-3));
}

// Motivation (review of PR #4866): with z0 >= zref, C <= 0 and C / (2 beta_m) <= 0 would pull
// every stable z/L onto the unstable side. With the bracket floor the supported wind
// u* most_bracket(C, psi_m) / kappa still has a minimum, at the kink where the floor binds.
TEST(MOSTIteration, FoldIsPositiveWhenZ0ReachesZref)
{
    const similarity_funs sfuns{};
    for (const Real z0 : {Real(1.0), Real(1.2), Real(0.85)}) {
        const Real C = std::log(Real(1.0) / z0);
        const Real zf = sfuns.calc_zeta_fold(C);
        EXPECT_GT(zf, Real(0.0)) << "z0 = " << z0;

        const Real a = KAPPA * CONST_GRAV * Real(0.004) * Real(1.0) / Real(300.0);
        Real best_u = Real(1.0e30);
        Real best_zeta = Real(0.0);
        for (int n = 1; n < 20000; ++n) {
            const Real ustar = Real(1.0e-5) * n;
            const Real zeta  = a / (ustar * ustar * ustar);
            const Real U     = ustar * most_bracket(C, sfuns.calc_psi_m(zeta)) / KAPPA;
            if (U < best_u) { best_u = U; best_zeta = zeta; }
        }
        EXPECT_NEAR(best_zeta / zf, Real(1.0), Real(1.0e-2)) << "z0 = " << z0;
    }
}

// Motivation: the fold flux scale must only act at the fold, and never raise the flux.
TEST(MOSTIteration, FoldFluxScaleActsOnlyAtTheFold)
{
    const Real zf = Real(0.39);
    auto scale = [=] (bool at_fold, Real ustar, Real ts, Real tq, bool spec_q) {
        return most_fold_flux_scale(at_fold, zf, ustar, Real(300.0), ts, tq, spec_q,
                                    KAPPA, CONST_GRAV, Real(5.0));
    };
    EXPECT_EQ(scale(false, Real(0.09), Real(-0.004), Real(0.0), true).temp, Real(1.0));
    EXPECT_EQ(scale(true,  Real(0.09), Real( 0.004), Real(0.0), true).temp, Real(1.0));
    const most_flux_scale s = scale(true, Real(0.08), Real(-0.004), Real(0.0), true);
    EXPECT_GT(s.temp, Real(0.0));
    EXPECT_LT(s.temp, Real(1.0));
}

// Motivation (review of PR #4866): the fold limits the net virtual flux, so only the
// prescribed components that are downward may be reduced, and the net must reach the flux
// the wind carries at the fold.
TEST(MOSTIteration, FoldFluxScaleReducesOnlyTheDownwardPrescribedFlux)
{
    const Real zf = Real(0.39);
    const Real ustar = Real(0.03);
    const Real tflux_fold = -ustar * ustar * ustar * Real(300.0) * zf / (KAPPA * CONST_GRAV * Real(5.0));
    auto scale = [=] (Real ts, Real tq, bool spec_q) {
        return most_fold_flux_scale(true, zf, ustar, Real(300.0), ts, tq, spec_q,
                                    KAPPA, CONST_GRAV, Real(5.0));
    };

    // Cooling with evaporation: only the sensible part is reduced
    for (const bool spec_q : {true, false}) {
        const most_flux_scale s = scale(Real(-0.004), Real(0.0035), spec_q);
        EXPECT_EQ(s.moist, Real(1.0));
        EXPECT_NEAR(s.temp * Real(-0.004) + Real(0.0035), tflux_fold, Real(1.0e-12));
        EXPECT_LT(s.temp * Real(-0.004) + Real(0.0035), Real(0.0));
    }

    // Heating with prescribed condensation: only the moisture part is reduced
    {
        const most_flux_scale s = scale(Real(0.001), Real(-0.004), true);
        EXPECT_EQ(s.temp, Real(1.0));
        EXPECT_NEAR(Real(0.001) + s.moist * Real(-0.004), tflux_fold, Real(1.0e-12));
    }

    // Heating with diagnosed condensation: no prescribed flux drives the fold
    {
        const most_flux_scale s = scale(Real(0.001), Real(-0.004), false);
        EXPECT_EQ(s.temp, Real(1.0));
        EXPECT_EQ(s.moist, Real(1.0));
    }

    // Both downward and prescribed: a common factor
    {
        const most_flux_scale s = scale(Real(-0.002), Real(-0.002), true);
        EXPECT_EQ(s.temp, s.moist);
        EXPECT_NEAR(s.temp * Real(-0.004), tflux_fold, Real(1.0e-12));
    }
}

// Motivation (issue #4857, canopy column): the stable flux -0.004 K m/s at zref 5 m needs a
// wind of at least 1.2621 m/s (kappa = 0.41). The iteration hit its cap or collapsed onto
// u* = 0.001 with t_surf thousands of K below the air; it must return the fold instead. The
// issue's wind of 1.26178 m/s is only 3e-4 below the fold, so 0.8 m/s checks a point well
// below it.
TEST(MOSTIteration, StableFluxBelowTheFoldReturnsTheFold)
{
    for (const Real U : {Real(1.26178), Real(0.8)}) {
        const Point p{Real(5.0), Real(0.1), U, Real(299.866), Real(0.0)};
        const Real F = Real(-0.004);
        const Result r = run_flux(p, F);
        const Real C = std::log(p.zref / p.z0);

        EXPECT_NEAR(p.zref / r.olen, C / Real(10.0), Real(1.0e-10)) << "U = " << U;
        EXPECT_NEAR(r.ustar, Real(2.0) * KAPPA * p.umm / (Real(3.0) * C), Real(1.0e-4)) << "U = " << U;
        const Real applied = -r.ustar * r.tstar;
        EXPECT_LT(applied, Real(0.0)) << "U = " << U;
        EXPECT_GT(applied, F) << "U = " << U;    // limited, never stronger than prescribed
        EXPECT_LT(r.t_surf, p.theta) << "U = " << U;
        EXPECT_GT(r.t_surf, p.theta - Real(1.0)) << "U = " << U;
    }
}

// Motivation (review of PR #4866): with evaporation offsetting part of a prescribed cooling
// flux, the fold must reduce the cooling only, leave the prescribed moisture flux as given,
// and keep the net virtual flux downward, as the stored Obukhov length says.
TEST(MOSTIteration, StableFluxBelowTheFoldKeepsThePrescribedEvaporation)
{
    const Point p{Real(5.0), Real(0.1), Real(0.4), Real(300.0), Real(0.0)};
    const Real F = Real(-0.004);
    const Real Q = Real(0.0035) / (epsv * p.theta);
    const Result r = run_flux(p, F, Q);
    const Real C = std::log(p.zref / p.z0);

    EXPECT_NEAR(p.zref / r.olen, C / Real(10.0), Real(1.0e-10));
    EXPECT_NEAR(-r.ustar * r.qstar, Q, Real(1.0e-12));
    const Real applied = -r.ustar * r.tstar;
    EXPECT_LT(applied, Real(0.0));
    EXPECT_GT(applied, F);
    EXPECT_LT(applied + Q * epsv * p.theta, Real(0.0));
}

// Motivation (review of PR #4866): a canopy cell with z0 = 1.2 m above zref = 1 m must stay
// on the stable side, with u* > 0 and the prescribed cooling applied no more strongly than
// given, both below its 0.15 m/s fold wind (at the wind floor) and well above it.
TEST(MOSTIteration, StableFluxWithZ0AboveZrefStaysStable)
{
    for (const Real U : {Real(0.1), Real(3.0)}) {
        const Point p{Real(1.0), Real(1.2), U, Real(300.0), Real(0.0)};
        const Real F = Real(-0.02);
        const Result r = run_flux(p, F);
        EXPECT_GT(r.ustar, Real(0.0)) << "U = " << U;
        EXPECT_GT(r.olen, Real(0.0)) << "U = " << U;
        const Real applied = -r.ustar * r.tstar;
        EXPECT_LT(applied, Real(0.0)) << "U = " << U;
        EXPECT_GE(applied, F - Real(1.0e-12)) << "U = " << U;
        EXPECT_LT(r.t_surf, p.theta) << "U = " << U;
    }
}

// Motivation (issue #4857, 3D terrain case, forest cell): 1.01425 m/s is 5e-4 below the
// 1.01475 m/s fold wind of this cell.
TEST(MOSTIteration, StableFluxBelowTheFoldAtAForestCell)
{
    const Point p{Real(4.42926), Real(0.1), Real(1.01425), Real(299.942), Real(0.0)};
    const Result r = run_flux(p, Real(-0.0025));
    EXPECT_NEAR(r.ustar, Real(0.07313), Real(1.0e-4));
    EXPECT_LT(r.t_surf, p.theta);
    EXPECT_GT(r.t_surf, p.theta - Real(1.0));
}

// Motivation: above the fold the stable point keeps its solution and the full prescribed flux.
TEST(MOSTIteration, StableFluxAboveTheFoldIsUnchanged)
{
    const Point p{Real(5.0), Real(0.1), Real(3.0), Real(299.866), Real(0.0)};
    const Real F = Real(-0.004);
    const Result r = run_flux(p, F);
    EXPECT_NEAR(r.ustar, Real(0.31088), Real(1.0e-4));
    EXPECT_NEAR(-r.ustar * r.tstar, F, Real(1.0e-12));
}

// Motivation (issue #4857, z0 = 0.5 m cell at the wind floor): ln(zref/z0) - psi_m changed
// sign and the iterate went negative; it must keep u* > 0 and t_surf above the air.
TEST(MOSTIteration, UnstableFluxAtARoughCalmCellKeepsItsSign)
{
    const Point p{Real(4.9531), Real(0.5), Real(0.11076), Real(300.848), Real(0.0)};
    const Result r = run_flux(p, Real(0.018));
    EXPECT_NEAR(r.ustar, Real(0.04541), Real(1.0e-4));
    EXPECT_GT(r.t_surf, p.theta);
}

// Motivation (PR #4764 reproducer, and zref/z0 = 20 which the z/L >= -5 bound left open):
// the diagnosed surface must stay warmer than the air under an upward flux.
TEST(MOSTIteration, UnstableWeakWindFluxKeepsTheSurfaceWarmer)
{
    for (const Real z0 : {Real(0.05), Real(0.1)}) {
        const Point p{Real(2.0), z0, Real(0.1), Real(300.0), Real(0.0)};
        const Result r = run_flux(p, Real(0.1408));
        EXPECT_NEAR(r.ustar, Real(0.04099), Real(1.0e-4)) << "z0 = " << z0;
        EXPECT_GT(r.t_surf, p.theta) << "z0 = " << z0;
    }
}

// Motivation: the Charnock kernel shares the iteration and must converge with the right sign.
TEST(MOSTIteration, CharnockUnstableWeakWindKeepsItsSign)
{
    const Point p{Real(10.0), Real(1.0e-4), Real(0.1), Real(300.0), Real(0.0)};
    const Result r = run(surface_flux_charnock(Real(0.1), Real(0.0), Real(0.0185), true, true), p);
    EXPECT_GT(r.ustar, Real(0.0));
    EXPECT_GT(r.t_surf, p.theta);
}

// Motivation: with zref < e z0 the old floor of 1 exceeded ln(zref/z0) and lowered u* in
// neutral air; neutral u* must follow the log law down to zref = e^(1/2) z0.
TEST(MOSTIteration, NeutralSurfaceTempFollowsTheLogLawAtSmallC)
{
    const Point p{Real(2.0), Real(1.0), Real(5.0), Real(300.0), Real(300.0)};
    const Result r = run_temp(p);
    const Real loglaw = KAPPA * p.umm / std::log(p.zref / p.z0);
    // The old floor gave 0.69 here; the ~2% left comes from the kernel's absolute 1e-3
    // tolerance on z/L, which stops a cold start (z/L = zref/1000) just short of 0
    EXPECT_NEAR(r.ustar / loglaw, Real(1.0), Real(0.05));
}

// Motivation: strongly unstable, rough, low reference height: u* > 0 and theta* < 0.
TEST(MOSTIteration, UnstableSurfaceTempAtSmallCKeepsItsSign)
{
    const Point p{Real(1.6), Real(0.5), Real(0.3), Real(300.0), Real(310.0)};
    const Result r = run_temp(p);
    EXPECT_GT(r.ustar, Real(0.0));
    EXPECT_LT(r.tstar, Real(0.0));
}

// Motivation: in stable air with z0 at or above zref the relaxed z/L iteration of
// surface_temp crept through a near-touch of its residual and hit the iteration cap
// (an abort). The bracketed fallback must return a root of zeta = Rib F(zeta) instead.
TEST(MOSTIteration, StableSurfaceTempWithZ0AtZrefConverges)
{
    const std::vector<Point> pts{
        {Real(1.0), Real(1.0), Real(0.8), Real(300.0), Real(295.0)},
        {Real(5.0), Real(5.0), Real(0.8), Real(300.0), Real(299.0)},
    };
    similarity_funs sfuns;
    for (const auto& p : pts) {
        const Result r = run_temp(p);
        const Real C    = std::log(p.zref / p.z0);
        const Real Rib  = CONST_GRAV * p.zref / p.theta * (p.theta - p.t_surf) / (p.umm * p.umm);
        const Real zeta = p.zref / r.olen;
        const Real CmPsim = most_bracket(C, sfuns.calc_psi_m2(zeta));
        const Real CmPsih = most_bracket(C, sfuns.calc_psi_h2(zeta));
        EXPECT_GT(zeta, Real(0.0));
        EXPECT_NEAR(zeta, Rib * CmPsim * CmPsim / CmPsih, Real(2.0e-3));
        EXPECT_GT(r.ustar, Real(0.0));
        EXPECT_GT(r.tstar, Real(0.0));
    }
}

// Motivation: the adiabatic kernel divided by ln(zref/z0) directly, so u* was infinite at
// z0 = zref and negative above it; it must use the bounded bracket like the others.
TEST(MOSTIteration, AdiabaticUstarIsBoundedWhenZ0ReachesZref)
{
    for (const Real z0 : {Real(1.0), Real(2.0)}) {
        const Point p{Real(1.0), z0, Real(2.0), Real(300.0), Real(300.0)};
        const Result r = run(adiabatic(Real(0.0), Real(0.0)), p);
        EXPECT_TRUE(std::isfinite(r.ustar));
        EXPECT_NEAR(r.ustar, KAPPA * p.umm / most_bracket_floor(std::log(p.zref / z0)), Real(1.0e-12));
    }
}
