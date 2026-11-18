// The stable side of the MOST similarity functions can be the QNSE form
// phi = (1 + a zeta)/(1 + b zeta), psi = -((a - b)/b) ln(1 + b zeta) (erf.most.stable_functions =
// qnse) instead of Businger-Dyer. The functions must match their closed forms, phi must stay the
// derivative companion of psi, the unstable side and the default must be untouched, the fold of
// a prescribed stable flux must be the first minimum of the wind the flux can be carried by, and
// the prescribed-flux kernel must solve the QNSE profile once given the functions.
#include <cmath>
#include <vector>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "BoundaryConditions/ERF_MOSTStress.H"

namespace
{

using amrex::Real;

constexpr Real tol  = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-5);
constexpr Real ftol = (sizeof(Real) == 8) ? Real(1.0e-9)  : Real(1.0e-4);
constexpr Real xtol = (sizeof(Real) == 8) ? Real(1.0e-10) : Real(1.0e-7); // identities up to rounding

similarity_funs qnse_default ()
{
    similarity_funs sf;
    sf.set_stable_functions(true, Real(2.5), Real(0.2), Real(2.5), Real(0.2));
    return sf;
}

// Wind that a friction velocity supports for a prescribed stable flux: with
// zeta = A / u*^3 (A = kappa g zref |F| / thetav), U = u* most_bracket(C, psi_m(zeta)) / kappa
Real supported_wind (const similarity_funs& sf, Real C, Real A, Real ustar)
{
    const Real zeta = A / (ustar * ustar * ustar);
    return ustar * most_bracket(C, sf.calc_psi_m(zeta)) / KAPPA;
}

// z/L of the smallest supported wind over the sampled u* (zeta in (zeta_min, zeta_max])
Real brute_force_fold (const similarity_funs& sf, Real C, Real A, Real zeta_min, Real zeta_max, int n)
{
    Real best_u = Real(1.0e30);
    Real best_zeta = Real(0.0);
    const Real ustar_lo = std::cbrt(A / zeta_max);
    const Real ustar_hi = std::cbrt(A / zeta_min);
    for (int i = 0; i <= n; ++i) {
        const Real ustar = ustar_lo + (ustar_hi - ustar_lo) * Real(i) / Real(n);
        const Real U = supported_wind(sf, C, A, ustar);
        if (U < best_u) { best_u = U; best_zeta = A / (ustar * ustar * ustar); }
    }
    return best_zeta;
}

struct Point
{
    Real zref;
    Real z0;
    Real umm;
    Real theta;
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

// One cold-start iterate_flux of the prescribed-flux kernel on a single dry point (no w*)
Result run_flux (const Point& p, Real flux, const similarity_funs& sf)
{
    surface_flux kernel(flux, Real(0.0), true);
    kernel.set_similarity_funs(sf);
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

// u* of the physical branch (u* above the fold's) that carries the wind U with the flux,
// by bisection of supported_wind, independent of the kernel's relaxed iteration
Real bisect_ustar (const similarity_funs& sf, Real C, Real A, Real U)
{
    Real lo = std::cbrt(A / sf.calc_zeta_fold(C));
    Real hi = Real(10.0);
    EXPECT_LE(supported_wind(sf, C, A, lo), U);
    for (int n = 0; n < 200; ++n) {
        const Real m = Real(0.5) * (lo + hi);
        if (supported_wind(sf, C, A, m) < U) { lo = m; } else { hi = m; }
    }
    return Real(0.5) * (lo + hi);
}

} // namespace

// Motivation: a default-constructed similarity_funs must be exactly Businger-Dyer, so every
// kernel, diagnostic and PBL gradient that is not given the surface layer's functions keeps
// its answer bit for bit.
TEST(MOSTStableFunctions, DefaultIsBusingerDyer)
{
    const similarity_funs sf{};
    EXPECT_FALSE(sf.uses_qnse());
    EXPECT_EQ(sf.calc_psi_m(Real(0.5)), Real(-2.5));
    EXPECT_EQ(sf.calc_psi_h(Real(0.5)), Real(-2.5));
    EXPECT_EQ(sf.calc_phi_m(Real(0.5)), Real(3.5));
    EXPECT_EQ(sf.calc_phi_h(Real(0.5)), Real(3.5));
    EXPECT_EQ(sf.calc_zeta_fold(Real(3.0)), Real(0.3));
}

// Motivation: psi = -((a - b)/b) ln(1 + b zeta) and phi = (1 + a zeta)/(1 + b zeta), with the
// momentum and heat coefficients kept apart. Values from a Python evaluation of the closed forms.
TEST(MOSTStableFunctions, QnseMatchesTheClosedForms)
{
    similarity_funs sf;
    sf.set_stable_functions(true, Real(2.5), Real(0.2), Real(3.0), Real(0.5));
    EXPECT_TRUE(sf.uses_qnse());
    const Real zeta[4]  = {Real(0.1), Real(1.0), Real(10.0), Real(100.0)};
    const Real psi_m[4] = {Real(-0.227730213906067), Real(-2.09669790313048),
                           Real(-12.6340413196833), Real(-35.0120080338194)};
    const Real phi_m[4] = {Real(1.22549019607843), Real(2.91666666666667),
                           Real(8.66666666666667), Real(11.952380952381)};
    const Real psi_h[4] = {Real(-0.24395082084716), Real(-2.02732554054082),
                           Real(-8.95879734614027), Real(-19.6591281636216)};
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(sf.calc_psi_m(zeta[n]), psi_m[n], tol * std::abs(psi_m[n])) << "zeta = " << zeta[n];
        EXPECT_NEAR(sf.calc_phi_m(zeta[n]), phi_m[n], tol * phi_m[n]) << "zeta = " << zeta[n];
        EXPECT_NEAR(sf.calc_psi_h(zeta[n]), psi_h[n], tol * std::abs(psi_h[n])) << "zeta = " << zeta[n];
        // phi_h with a_h = 3, b_h = 0.5: (1 + 3 zeta)/(1 + 0.5 zeta)
        EXPECT_NEAR(sf.calc_phi_h(zeta[n]), (Real(1.0) + Real(3.0) * zeta[n]) / (Real(1.0) + Real(0.5) * zeta[n]),
                    tol * Real(10.0)) << "zeta = " << zeta[n];
    }
    // phi tends to a/b in very stable air
    EXPECT_NEAR(sf.calc_phi_m(Real(1.0e6)), Real(12.5), Real(1.0e-4));
    EXPECT_NEAR(sf.calc_phi_h(Real(1.0e6)), Real(6.0), Real(1.0e-4));
}

// Motivation: the profiles integrate phi, so psi must satisfy zeta dpsi/dzeta = 1 - phi on the
// stable side for both pairs, as the Businger-Dyer pair does.
TEST(MOSTStableFunctions, QnsePhiIsTheDerivativeCompanionOfPsi)
{
    similarity_funs sf;
    sf.set_stable_functions(true, Real(2.5), Real(0.2), Real(2.0), Real(0.25));
    // central difference: truncation error h^2, rounding error eps / h
    const Real dtol = (sizeof(Real) == 8) ? Real(1.0e-7) : Real(2.0e-3);
    const Real hfac = (sizeof(Real) == 8) ? Real(1.0e-4) : Real(1.0e-2);
    for (const Real zeta : {Real(0.05), Real(0.5), Real(2.0), Real(20.0)}) {
        const Real h = hfac * zeta;
        const Real dpsi_m = (sf.calc_psi_m(zeta + h) - sf.calc_psi_m(zeta - h)) / (Real(2.0) * h);
        const Real dpsi_h = (sf.calc_psi_h(zeta + h) - sf.calc_psi_h(zeta - h)) / (Real(2.0) * h);
        EXPECT_NEAR(zeta * dpsi_m, Real(1.0) - sf.calc_phi_m(zeta), dtol * sf.calc_phi_m(zeta)) << "zeta = " << zeta;
        EXPECT_NEAR(zeta * dpsi_h, Real(1.0) - sf.calc_phi_h(zeta), dtol * sf.calc_phi_h(zeta)) << "zeta = " << zeta;
    }
}

// Motivation: the option is the stable side only; the unstable side must stay bit for bit
// Businger-Dyer, and the Jimenez pair must not change either.
TEST(MOSTStableFunctions, UnstableSideAndJimenezPairUnchanged)
{
    const similarity_funs bd{};
    const similarity_funs qn = qnse_default();
    for (const Real zeta : {Real(-0.01), Real(-1.0), Real(-30.0)}) {
        EXPECT_EQ(qn.calc_psi_m(zeta), bd.calc_psi_m(zeta)) << "zeta = " << zeta;
        EXPECT_EQ(qn.calc_psi_h(zeta), bd.calc_psi_h(zeta)) << "zeta = " << zeta;
        EXPECT_EQ(qn.calc_phi_m(zeta), bd.calc_phi_m(zeta)) << "zeta = " << zeta;
        EXPECT_EQ(qn.calc_phi_h(zeta), bd.calc_phi_h(zeta)) << "zeta = " << zeta;
    }
    for (const Real zeta : {Real(-1.0), Real(0.0), Real(2.0)}) {
        EXPECT_EQ(qn.calc_psi_m2(zeta), bd.calc_psi_m2(zeta)) << "zeta = " << zeta;
        EXPECT_EQ(qn.calc_psi_h2(zeta), bd.calc_psi_h2(zeta)) << "zeta = " << zeta;
    }
    EXPECT_EQ(qn.calc_psi_m(Real(0.0)), Real(0.0));
    EXPECT_EQ(qn.calc_phi_m(Real(0.0)), Real(1.0));
}

// Motivation: with the QNSE psi_m the wind a prescribed stable flux can be carried by has its
// first minimum where f = C + c ln(1 + b zeta) - 3 c b zeta / (1 + b zeta) crosses zero, or at
// the bracket floor's kink when z0 approaches or exceeds zref. The pinned values are the roots
// of f from a Python bisection and the closed-form kink.
TEST(MOSTStableFunctions, QnseFoldIsTheFirstMinimumOfTheSupportedWind)
{
    const similarity_funs sf = qnse_default();

    // C = ln(5/0.1): root of f
    {
        const Real C = std::log(Real(5.0) / Real(0.1));
        const Real zf = sf.calc_zeta_fold(C);
        EXPECT_NEAR(zf, Real(1.0918333752866505), ftol);
        const Real A = KAPPA * CONST_GRAV * Real(0.004) * Real(5.0) / Real(300.0);
        // the first minimum lies below 2/b = 10 (beyond it the supported wind rises again,
        // then falls to zero as u* does)
        EXPECT_NEAR(brute_force_fold(sf, C, A, Real(1.0e-3), Real(10.0), 400000) / zf, Real(1.0), Real(1.0e-2));
        EXPECT_GT(supported_wind(sf, C, A, std::cbrt(A / zf) * Real(1.05)), supported_wind(sf, C, A, std::cbrt(A / zf)));
        EXPECT_GT(supported_wind(sf, C, A, std::cbrt(A / zf) * Real(0.95)), supported_wind(sf, C, A, std::cbrt(A / zf)));
    }

    // z0 above zref (C < 0): the kink, positive
    {
        const Real C = std::log(Real(1.0) / Real(1.2));
        const Real zf = sf.calc_zeta_fold(C);
        EXPECT_NEAR(zf, Real(0.30563902390190356), ftol);
        const Real A = KAPPA * CONST_GRAV * Real(0.004) * Real(1.0) / Real(300.0);
        EXPECT_NEAR(brute_force_fold(sf, C, A, Real(1.0e-3), Real(10.0), 400000) / zf, Real(1.0), Real(1.0e-2));
    }

    // Large C (z0 far below zref): no fold; the supported wind falls with z/L everywhere
    {
        const Real C = Real(12.0);
        EXPECT_EQ(sf.calc_zeta_fold(C), most_zeta_lim);
        const Real A = KAPPA * CONST_GRAV * Real(0.004) * Real(10.0) / Real(300.0);
        Real last = Real(1.0e30);
        for (int n = 1; n <= 1000; ++n) {
            const Real zeta = Real(0.01) * Real(n);
            const Real U = supported_wind(sf, C, A, std::cbrt(A / zeta));
            EXPECT_LT(U, last) << "zeta = " << zeta;
            last = U;
        }
    }

    // a = b: psi_m = 0 on the stable side, no fold
    {
        similarity_funs flat;
        flat.set_stable_functions(true, Real(0.2), Real(0.2), Real(2.5), Real(0.2));
        EXPECT_EQ(flat.calc_psi_m(Real(3.0)), Real(0.0));
        EXPECT_EQ(flat.calc_zeta_fold(Real(3.0)), most_zeta_lim);
    }
}

// Motivation: near neutral the QNSE correction has slope a - b (2.3 with the defaults, against
// Businger-Dyer's 5), so the option changes the weakly stable regime too, not only the strongly
// stable tail; a user who wants the Businger-Dyer slope with a bounded phi can take a = 5.2,
// b = 0.2 (asymptote 26). The defaults also give phi_h = phi_m, a turbulent Prandtl number of 1.
TEST(MOSTStableFunctions, QnseNearNeutralSlopeIsAMinusB)
{
    const similarity_funs sf = qnse_default();
    const similarity_funs bd{};
    const Real z = Real(1.0e-4);
    const Real stol = (sizeof(Real) == 8) ? Real(1.0e-6) : Real(1.0e-3); // (phi - 1)/z at z = 1e-4 loses digits in single precision
    EXPECT_NEAR((sf.calc_phi_m(z) - Real(1.0)) / z, Real(2.3), Real(1.0e-3));
    EXPECT_NEAR((bd.calc_phi_m(z) - Real(1.0)) / z, Real(5.0), stol);
    EXPECT_NEAR(sf.calc_phi_m(Real(0.1)) - Real(1.0), Real(0.5) * (bd.calc_phi_m(Real(0.1)) - Real(1.0)), Real(0.03));
    EXPECT_EQ(sf.calc_phi_h(Real(0.7)), sf.calc_phi_m(Real(0.7)));
    similarity_funs steep;
    steep.set_stable_functions(true, Real(5.2), Real(0.2), Real(5.2), Real(0.2));
    EXPECT_NEAR((steep.calc_phi_m(z) - Real(1.0)) / z, Real(5.0), Real(2.0e-3));
    EXPECT_NEAR(steep.calc_phi_m(Real(1.0e6)), Real(26.0), Real(1.0e-3));
}

// Motivation: the Newton solve of the fold must give the root of f to its tolerance over the
// whole range of C that has one (a 200-pass bisection of the same f is the reference), with
// the kink and the no-fold cases unchanged.
TEST(MOSTStableFunctions, QnseFoldNewtonMatchesABisectionOfF)
{
    const similarity_funs sf = qnse_default();
    const Real a = Real(2.5), b = Real(0.2), c = (a - b) / b;
    auto f = [&] (Real C, Real z) { return C + c * std::log(Real(1.0) + b * z) - Real(3.0) * c * b * z / (Real(1.0) + b * z); };
    const Real ntol = (sizeof(Real) == 8) ? Real(1.0e-9) : Real(1.0e-5); // in single precision f is flat to rounding near its root for large C
    int n_root = 0;
    for (int i = 0; i <= 400; ++i) {
        const Real C = Real(-1.0) + Real(0.03) * Real(i);   // -1 .. 11
        const Real bf = most_bracket_floor(C);
        Real lo = (C < bf) ? (std::exp((bf - C) / c) - Real(1.0)) / b : Real(0.0);
        Real hi = Real(2.0) / b;
        const Real zf = sf.calc_zeta_fold(C);
        if (f(C, lo) <= Real(0.0)) { EXPECT_EQ(zf, lo) << "C = " << C; continue; }
        if (lo >= hi || f(C, hi) >= Real(0.0)) { EXPECT_EQ(zf, most_zeta_lim) << "C = " << C; continue; }
        for (int k = 0; k < 200; ++k) {
            const Real m = Real(0.5) * (lo + hi);
            if (f(C, m) > Real(0.0)) { lo = m; } else { hi = m; }
        }
        EXPECT_NEAR(zf / lo, Real(1.0), ntol) << "C = " << C;
        ++n_root;
    }
    EXPECT_GT(n_root, 300);
}

// Motivation: once given the functions, the prescribed-flux kernel must converge to the u* of
// the QNSE profile (checked against an independent bisection of U(u*)), which differs from the
// Businger-Dyer u* at the same point; a kernel that ignored set_similarity_funs would return
// the Businger-Dyer value.
TEST(MOSTStableFunctions, FluxKernelSolvesTheQnseProfile)
{
    const Point p{Real(5.0), Real(0.1), Real(3.0), Real(300.0)};
    const Real F = Real(-0.02);   // z/L about 0.09: the two stable arms differ by 5 % in u*
    const Real C = std::log(p.zref / p.z0);
    const Real A = KAPPA * CONST_GRAV * (-F) * p.zref / p.theta;
    const similarity_funs qn = qnse_default();
    const similarity_funs bd{};
    const Result rq = run_flux(p, F, qn);
    const Result rb = run_flux(p, F, bd);
    const Real uq = bisect_ustar(qn, C, A, p.umm);
    const Real ub = bisect_ustar(bd, C, A, p.umm);
    const Real ktol = (sizeof(Real) == 8) ? Real(2.0e-4) : Real(2.0e-3);
    EXPECT_NEAR(rq.ustar, uq, ktol);
    EXPECT_NEAR(rb.ustar, ub, ktol);
    EXPECT_GT(std::abs(uq - ub), Real(1.0e-2));        // the two profiles are told apart
    EXPECT_NEAR(p.zref / rq.olen, A / (rq.ustar * rq.ustar * rq.ustar), Real(1.0e-3) * A / (rq.ustar * rq.ustar * rq.ustar));
    EXPECT_NEAR(-rq.ustar * rq.tstar, F, xtol);   // above the fold the full flux is applied
    EXPECT_LT(rq.t_surf, p.theta);
}

// Motivation: QNSE's fold sits at a lower wind than Businger-Dyer's (0.53 against 0.70 m/s at
// zref 4.43 m, z0 0.5 m, flux -0.0025 K m/s: the forest cell that stopped the WINSENT runs), so
// at 0.6 m/s the QNSE kernel carries the full prescribed flux while the Businger-Dyer kernel is
// held at its fold with a weaker flux; at 0.5 m/s both are held at their own folds.
TEST(MOSTStableFunctions, QnseCarriesAStableFluxAtWindsBelowTheBusingerDyerFold)
{
    const Real F = Real(-0.0025);
    const similarity_funs qn = qnse_default();
    const similarity_funs bd{};
    const Real C = std::log(Real(4.43) / Real(0.5));
    EXPECT_NEAR(qn.calc_zeta_fold(C), Real(0.539472968792851), ftol);
    EXPECT_NEAR(bd.calc_zeta_fold(C), C / Real(10.0), ftol);

    const Point p6{Real(4.43), Real(0.5), Real(0.6), Real(300.0)};
    const Result q6 = run_flux(p6, F, qn);
    const Result b6 = run_flux(p6, F, bd);
    EXPECT_NEAR(-q6.ustar * q6.tstar, F, xtol);
    EXPECT_LT(p6.zref / q6.olen, qn.calc_zeta_fold(C));
    EXPECT_NEAR(p6.zref / b6.olen, bd.calc_zeta_fold(C), xtol);
    EXPECT_GT(-b6.ustar * b6.tstar, F);   // weaker than prescribed
    EXPECT_LT(-b6.ustar * b6.tstar, Real(0.0));

    const Point p5{Real(4.43), Real(0.5), Real(0.5), Real(300.0)};
    const Result q5 = run_flux(p5, F, qn);
    const Result b5 = run_flux(p5, F, bd);
    EXPECT_NEAR(p5.zref / q5.olen, qn.calc_zeta_fold(C), xtol);
    EXPECT_NEAR(p5.zref / b5.olen, bd.calc_zeta_fold(C), xtol);
    EXPECT_GT(-q5.ustar * q5.tstar, F);
    EXPECT_LT(-q5.ustar * q5.tstar, Real(0.0));
    // held at the fold, u* is kappa U / bracket(C, psi_m(zeta_fold)) for either pair
    const Real ktol = (sizeof(Real) == 8) ? Real(2.0e-4) : Real(2.0e-3);
    EXPECT_NEAR(q5.ustar, KAPPA * p5.umm / most_bracket(C, qn.calc_psi_m(qn.calc_zeta_fold(C))), ktol);
    EXPECT_NEAR(b5.ustar, KAPPA * p5.umm / most_bracket(C, bd.calc_psi_m(bd.calc_zeta_fold(C))), ktol);
}
