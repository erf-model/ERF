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
    return Result{h[USTAR], h[TSTAR], h[TSURF], h[OLEN]};
}

Result run_flux (const Point& p, Real flux)
{
    return run(surface_flux(flux, Real(0.0), true), p);
}

Result run_temp (const Point& p)
{
    return run(surface_temp(Real(0.0), Real(0.0), false, 2, true), p);
}

} // namespace

// Motivation: the floor is 1 where the bracket is large and never exceeds the neutral value C.
TEST(MOSTIteration, BracketFloorIsOneForLargeCAndBelowNeutralForSmallC)
{
    EXPECT_EQ(most_bracket(Real(3.0), Real(5.0)), Real(1.0));
    EXPECT_EQ(most_bracket(Real(3.0), Real(1.0)), Real(2.0));
    EXPECT_EQ(most_bracket(Real(0.5), Real(2.0)), Real(0.25));
    EXPECT_EQ(most_bracket(Real(0.5), Real(0.0)), Real(0.5));
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

// Motivation: the fold flux scale must only act at the fold, and never raise the flux.
TEST(MOSTIteration, FoldFluxScaleActsOnlyAtTheFold)
{
    const Real zf = Real(0.39);
    EXPECT_EQ(most_fold_flux_scale(Real(0.2), zf, Real(0.09), Real(300.0), Real(-0.004),
                                   KAPPA, CONST_GRAV, Real(5.0)), Real(1.0));
    EXPECT_EQ(most_fold_flux_scale(zf, zf, Real(0.09), Real(300.0), Real(0.004),
                                   KAPPA, CONST_GRAV, Real(5.0)), Real(1.0));
    const Real s = most_fold_flux_scale(zf, zf, Real(0.08), Real(300.0), Real(-0.004),
                                        KAPPA, CONST_GRAV, Real(5.0));
    EXPECT_GT(s, Real(0.0));
    EXPECT_LT(s, Real(1.0));
}

// Motivation (issue #4857, canopy column): the wind 1.26 m/s is below the 1.28 m/s the
// stable flux -0.004 K m/s needs at zref 5 m. The iteration hit its cap or collapsed onto
// u* = 0.001 with t_surf thousands of K below the air; it must return the fold instead.
TEST(MOSTIteration, StableFluxBelowTheFoldReturnsTheFold)
{
    const Point p{Real(5.0), Real(0.1), Real(1.26178), Real(299.866), Real(0.0)};
    const Real F = Real(-0.004);
    const Result r = run_flux(p, F);
    const Real C = std::log(p.zref / p.z0);

    EXPECT_NEAR(p.zref / r.olen, C / Real(10.0), Real(1.0e-10));
    EXPECT_NEAR(r.ustar, Real(2.0) * KAPPA * p.umm / (Real(3.0) * C), Real(1.0e-4));
    const Real applied = -r.ustar * r.tstar;
    EXPECT_LT(applied, Real(0.0));
    EXPECT_GT(applied, F);                       // limited, never stronger than prescribed
    EXPECT_LT(r.t_surf, p.theta);
    EXPECT_GT(r.t_surf, p.theta - Real(1.0));
}

// Motivation (issue #4857, 3D terrain case, forest cell).
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
// neutral air by half; neutral u* must follow the log law.
TEST(MOSTIteration, NeutralSurfaceTempFollowsTheLogLawAtSmallC)
{
    const Point p{Real(1.6), Real(1.0), Real(5.0), Real(300.0), Real(300.0)};
    const Result r = run_temp(p);
    const Real loglaw = KAPPA * p.umm / std::log(p.zref / p.z0);
    // The old floor gave 0.47 here; the ~2% left comes from the kernel's absolute 1e-3
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
