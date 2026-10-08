// Contract of the immersed-forcing wall law's stability bounds (erf-model/ERF#4016). Always on: a
// derived Obukhov length is held at |L| >= 1.5 dz / 100 and every zeta at |zeta| <= 100, as the
// flat-ground surface layer does (a prescribed erf.if_Olen is used as given), and u* stays within
// [0, 2] m/s on every branch. A calm cell above the wall used to give u* = 0, L = 0 and zeta = z / 0;
// under a cooling flux zeta = +inf, psi_m = -inf and the momentum and temperature targets became
// 0 * inf = NaN. Opt-in, as in WRF's revised surface layer: erf.if_stability_wind_floor (default 0)
// floors the wind of the stability estimate, and erf.if_psi_cap_factor (default 1) caps psi_m in the
// momentum law at factor ln(z / z0), also before it forms u*. The temperature forcing keeps caps of
// ln(z / z0) (legacy law kept as it is, erf-model/ERF#4134). Inside the bounds and with the defaults,
// nothing changes.

#include <cmath>
#include <limits>

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_TableData.H>

#include <gtest/gtest.h>

#include "ERF_DataStruct.H"
#include "ERF_ImmersedForcing.H"
#include "ERF_ImmersedWallStability.H"
#include "ERF_IndexDefines.H"
#include "ERF_MOSTUtils.H"

using amrex::Real;

TEST(ImmersedWallStability, WindFloorIsOffByDefaultAndFloorsWhenSet)
{
    const SolverChoice sc{};
    EXPECT_EQ(sc.if_stability_wind_floor, Real(0.0));
    const Real f = Real(0.1);  // WRF's and the flat-ground surface layer's floor
    EXPECT_EQ(ib_stability::floored_wind(Real(0.0), f), Real(0.1));
    EXPECT_EQ(ib_stability::floored_wind(Real(0.05), f), Real(0.1));
    // above the floor the speed is untouched, bit for bit, and a floor of 0 never changes it
    for (Real s : {Real(0.1), Real(0.1000001), Real(3.7), Real(25.0)}) {
        EXPECT_EQ(ib_stability::floored_wind(s, f), s);
    }
    for (Real s : {Real(0.0), Real(0.05), Real(3.7)}) {
        EXPECT_EQ(ib_stability::floored_wind(s, Real(0.0)), s);
    }
}

TEST(ImmersedWallStability, ObukhovLengthIsBoundedWithItsSign)
{
    const Real z = Real(15.0);  // 1.5 dz for dz = 10 m
    const Real Lmin = z / ib_stability::zeta_max();
    // L = 0 from a zero friction velocity: -0 is the unstable bound, +0 the stable one
    EXPECT_EQ(ib_stability::bounded_obukhov_length(Real(-0.0), z), -Lmin);
    EXPECT_EQ(ib_stability::bounded_obukhov_length(Real(+0.0), z),  Lmin);
    EXPECT_EQ(ib_stability::bounded_obukhov_length(Real(-1.e-6), z), -Lmin);
    EXPECT_EQ(ib_stability::bounded_obukhov_length(Real( 1.e-6), z),  Lmin);
    // zeta at every height up to z stays within the bound
    for (Real L : {Real(-0.0), Real(0.0), Real(-0.01), Real(0.01)}) {
        const Real Lb = ib_stability::bounded_obukhov_length(L, z);
        for (Real zz : {Real(5.0), Real(15.0)}) {
            EXPECT_TRUE(std::isfinite(zz / Lb));
            EXPECT_LE(std::abs(zz / Lb), ib_stability::zeta_max());
        }
    }
    // inside the bound the length is untouched, bit for bit
    for (Real L : {Lmin, -Lmin, Real(0.5), Real(-0.5), Real(150.0), Real(-186.0), Real(1.e30)}) {
        EXPECT_EQ(ib_stability::bounded_obukhov_length(L, z), L);
    }
}

TEST(ImmersedWallStability, ZetaIsBounded)
{
    const Real zmax = ib_stability::zeta_max();
    EXPECT_EQ(ib_stability::bounded_zeta(Real(5.0), Real(-0.01)), -zmax);
    EXPECT_EQ(ib_stability::bounded_zeta(Real(5.0), Real(0.01)), zmax);
    EXPECT_EQ(ib_stability::bounded_zeta(Real(5.0), Real(1.e-30)), zmax);
    // inside the bound, z / L bit for bit
    for (Real L : {Real(-0.2), Real(0.5), Real(-186.0), Real(1.e30)}) {
        EXPECT_EQ(ib_stability::bounded_zeta(Real(5.0), L), Real(5.0) / L);
    }
    // a NaN stays NaN rather than becoming a bound
    EXPECT_TRUE(std::isnan(ib_stability::bounded_zeta(Real(5.0), std::numeric_limits<Real>::quiet_NaN())));
}

TEST(ImmersedWallStability, PsiCapIsWrfsAndNeverLooserThanTheLog)
{
    const SolverChoice sc{};
    EXPECT_EQ(sc.if_psi_cap_factor, Real(1.0));  // off by default
    const Real f = Real(0.9);                    // WRF's factor
    // z > z0: f ln(z / z0), so ln(z / z0) - psi >= (1 - f) ln(z / z0) > 0
    for (Real z : {Real(5.0), Real(15.0)}) {
        const Real ln = std::log(z / Real(0.1));
        const Real cap = ib_stability::psi_cap(z, Real(0.1), f);
        EXPECT_NEAR(cap, f * ln, Real(1.e-5) * cap);
        EXPECT_GT(ln - cap, Real(0.0));
        // a factor of 1 is the cap of ln(z / z0), bit for bit
        EXPECT_EQ(ib_stability::psi_cap(z, Real(0.1), Real(1.0)), ln);
    }
    // z <= z0, where ln(z / z0) <= 0 and f ln would be looser: the cap stays ln(z / z0)
    for (Real z : {Real(0.5), Real(1.0)}) {
        const Real ln = std::log(z / Real(1.0));
        EXPECT_EQ(ib_stability::psi_cap(z, Real(1.0), f), ln);
        EXPECT_LE(ib_stability::psi_cap(z, Real(1.0), f), f * ln);
    }
}

TEST(ImmersedWallStability, PsiMIsCappedBeforeItFormsUstar)
{
    const Real z0 = Real(0.5);
    const Real z = Real(15.0);
    const Real ln = std::log(z / z0);
    // factor < 1, z > z0: the denominator ln(z / z0) - psi_m stays at least (1 - f) ln(z / z0)
    const Real capped = ib_stability::psi_m_for_ustar(Real(3.9), z, z0, Real(0.9));
    EXPECT_NEAR(capped, Real(0.9) * ln, Real(1.e-5));
    EXPECT_GE(ln - capped, Real(0.1) * ln * (Real(1.0) - Real(1.e-5)));
    // below the cap, untouched
    EXPECT_EQ(ib_stability::psi_m_for_ustar(Real(1.0), z, z0, Real(0.9)), Real(1.0));
    // factor 1 and z <= z0: no cap on the u* side (development's law)
    EXPECT_EQ(ib_stability::psi_m_for_ustar(Real(3.9), z, z0, Real(1.0)), Real(3.9));
    EXPECT_EQ(ib_stability::psi_m_for_ustar(Real(3.9), Real(0.4), z0, Real(0.9)), Real(3.9));
}

TEST(ImmersedWallStability, FrictionVelocityIsClamped)
{
    EXPECT_EQ(ib_stability::clamped_ustar(Real(0.3)), Real(0.3));
    EXPECT_EQ(ib_stability::clamped_ustar(Real(-0.3)), Real(0.0));
    EXPECT_EQ(ib_stability::clamped_ustar(Real(91.0)), Real(2.0));
    EXPECT_EQ(ib_stability::clamped_ustar(std::numeric_limits<Real>::infinity()), Real(2.0));
    EXPECT_EQ(ib_stability::clamped_ustar(-std::numeric_limits<Real>::infinity()), Real(0.0));
    EXPECT_EQ(ib_stability::clamped_ustar(std::numeric_limits<Real>::quiet_NaN()), Real(0.0));
}

namespace {

// A column of 4 x 4 x 6 cells, dz = 10 m: solid below k = 1, a half-blanked wall cell at k = 1,
// fluid above, and no wind anywhere (the calm case of erf-model/ERF#4016).
struct CalmTerrainColumn {
    amrex::Box bx{amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 5)};
    amrex::Box gbx = amrex::grow(bx, 2);
    amrex::Geometry geom;
    amrex::FArrayBox u, v, w, cell, blank, src;

    CalmTerrainColumn ()
        : u(gbx, 1, amrex::The_Managed_Arena()), v(gbx, 1, amrex::The_Managed_Arena()),
          w(gbx, 1, amrex::The_Managed_Arena()), cell(gbx, NDRY, amrex::The_Managed_Arena()),
          blank(gbx, 1, amrex::The_Managed_Arena()), src(gbx, NDRY, amrex::The_Managed_Arena())
    {
        const amrex::RealBox rb({0.0, 0.0, 0.0}, {40.0, 40.0, 60.0});
        const amrex::Array<int, AMREX_SPACEDIM> periodic{1, 1, 0};
        geom.define(bx, rb, amrex::CoordSys::cartesian, periodic);

        u.setVal<amrex::RunOn::Host>(0.0);
        v.setVal<amrex::RunOn::Host>(0.0);
        w.setVal<amrex::RunOn::Host>(0.0);
        src.setVal<amrex::RunOn::Host>(0.0);
        const auto c = cell.array();
        const auto b = blank.array();
        amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) {
            c(i, j, k, Rho_comp)      = Real(1.2);
            c(i, j, k, RhoTheta_comp) = Real(1.2) * Real(300.0);
            c(i, j, k, 2)             = Real(0.0);
            b(i, j, k) = (k < 1) ? Real(1.0) : ((k == 1) ? Real(0.5) : Real(0.0));
        });
    }

    SolverChoice choice (Real tflux) const
    {
        SolverChoice sc{};
        sc.if_use_most       = true;
        sc.if_z0             = Real(0.1);
        sc.if_Cd_momentum    = Real(50.0);
        sc.if_Cd_scalar      = Real(5.0);
        sc.if_surf_temp_flux = tflux;
        return sc;
    }

    // whether every component of the source on box b is finite
    bool all_finite (const amrex::Box& b, int ncomp) const
    {
        amrex::Gpu::streamSynchronize();
        const auto s = src.const_array();
        bool ok = true;
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            for (int n = 0; n < ncomp; ++n) { if (!std::isfinite(s(i, j, k, n))) { ok = false; } }
        });
        return ok;
    }
};

} // namespace

TEST(ImmersedWallStability, CalmTerrainMomentumTargetIsFinite)
{
    // cooling (stable: zeta = +inf on development) and heating (unstable: zeta = -inf), with the
    // default inputs: the bounds on L and zeta alone keep the target finite
    for (Real tflux : {Real(-0.05), Real(0.05)}) {
        CalmTerrainColumn col;
        const SolverChoice sc = col.choice(tflux);
        const amrex::Box xbx = amrex::surroundingNodes(col.bx, 0);
        ImmersedForcingTerrain_Xmom(xbx, col.u.const_array(), col.v.const_array(), col.w.const_array(),
                                    col.cell.const_array(), col.blank.const_array(),
                                    amrex::Array4<const Real>{}, amrex::Array4<const Real>{},
                                    col.src.array(), col.geom, sc, Real(1.0));
        // every face the kernel writes, the last plane i = 4 included
        EXPECT_TRUE(col.all_finite(xbx, 1)) << "x-momentum, tflux = " << tflux;

        col.src.setVal<amrex::RunOn::Host>(0.0);
        const amrex::Box ybx = amrex::surroundingNodes(col.bx, 1);
        ImmersedForcingTerrain_Ymom(ybx, col.u.const_array(), col.v.const_array(), col.w.const_array(),
                                    col.cell.const_array(), col.blank.const_array(),
                                    amrex::Array4<const Real>{}, amrex::Array4<const Real>{},
                                    col.src.array(), col.geom, sc, Real(1.0));
        EXPECT_TRUE(col.all_finite(ybx, 1)) << "y-momentum, tflux = " << tflux;
    }
}

TEST(ImmersedWallStability, CalmTerrainHeatFluxTargetIsFinite)
{
    // default inputs: finite under cooling and heating (development returned 0 / 0)
    for (Real tflux : {Real(-0.05), Real(0.05)}) {
        CalmTerrainColumn col;
        const SolverChoice sc = col.choice(tflux);
        ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                      col.cell.const_array(), col.blank.const_array(),
                                      amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                      amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
        EXPECT_TRUE(col.all_finite(col.bx, 2)) << "rho theta, tflux = " << tflux;
    }
    // With the 0.1 m/s floor the cooling flux is carried by the calm cell, at a sane size: the
    // target stays within 20 K of the air above (theta* = -q / u* with this L, the formula of
    // 41e39f9c4, put it about 1.7e5 K away)
    CalmTerrainColumn col;
    SolverChoice sc = col.choice(Real(-0.05));
    sc.if_stability_wind_floor = Real(0.1);
    ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                  col.cell.const_array(), col.blank.const_array(),
                                  amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                  amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
    ASSERT_TRUE(col.all_finite(col.bx, 2));
    amrex::Gpu::streamSynchronize();
    const Real src = col.src.const_array()(1, 1, 1, RhoTheta_comp);
    const Real drag_rho = Real(5.0) / Real(10.0) * Real(1.2);
    EXPECT_GT(std::abs(src), Real(1.e-6));
    EXPECT_LT(std::abs(src) / drag_rho, Real(20.0)) << "target " << src / drag_rho << " K from the air";
}

TEST(ImmersedWallStability, PrescribedObukhovLengthAcrossTheUstarCancellation)
{
    // erf.if_Olen with z0 = 0.5 m on 10 m cells, 5 m/s: psi_m at 0.5 dz crosses ln(1.5 dz / z0) = ln 30
    // near L = -0.167 m, where the u* denominator cancels and then turns negative (u* clamped to
    // [0, 2] m/s). Across the sweep both psi_h caps bind (psi_h >= psi_m in unstable air), so the
    // log brackets close and the target is the air above: every source is finite and round-off.
    // Without the caps the target was up to 160 K above the air (u* at the clamp, theta* ~ -1.5e3 K).
    for (int n = 0; n <= 400; ++n) {
        const Real L = Real(-0.25) + Real(0.13) * Real(n) / Real(400);
        CalmTerrainColumn col;
        col.u.setVal<amrex::RunOn::Host>(5.0);
        SolverChoice sc = col.choice(Real(1.e-8));
        sc.if_z0 = Real(0.5);
        sc.if_Olen_in = L;
        ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                      col.cell.const_array(), col.blank.const_array(),
                                      amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                      amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
        ASSERT_TRUE(col.all_finite(col.bx, 2)) << "L = " << L;
        amrex::Gpu::streamSynchronize();
        EXPECT_LT(std::abs(col.src.const_array()(1, 1, 1, RhoTheta_comp)), Real(1.e-9)) << "L = " << L;
    }
}

TEST(ImmersedWallStability, PsiMCapKeepsAVelocityTargetInUnstableAir)
{
    // z0 = 0.5 m on 10 m cells, 1 m/s wind and a heating flux of 0.12 K m/s: L is about -1 m and
    // psi_m at 1.5 dz is about 2.8, past ln(0.5 dz / z0) = 2.30 and below ln(1.5 dz / z0) = 3.40
    // (so u* > 0). A cap of ln 10 makes the velocity target zero, and the wall face is relaxed to
    // rest; the 0.9 cap keeps a log-law target, so the source toward it is weaker. Sources at the
    // wall face (i, j, k) = (1, 1, 1), u = 1 m/s everywhere.
    Real src[2] = {0.0, 0.0};
    const Real factors[2] = {Real(0.9), Real(1.0)};
    for (int n = 0; n < 2; ++n) {
        CalmTerrainColumn col;
        col.u.setVal<amrex::RunOn::Host>(1.0);
        SolverChoice sc = col.choice(Real(0.12));
        sc.if_z0 = Real(0.5);
        sc.if_psi_cap_factor = factors[n];
        const amrex::Box xbx = amrex::surroundingNodes(col.bx, 0);
        ImmersedForcingTerrain_Xmom(xbx, col.u.const_array(), col.v.const_array(), col.w.const_array(),
                                    col.cell.const_array(), col.blank.const_array(),
                                    amrex::Array4<const Real>{}, amrex::Array4<const Real>{},
                                    col.src.array(), col.geom, sc, Real(1.0));
        ASSERT_TRUE(col.all_finite(xbx, 1));
        src[n] = col.src.const_array()(1, 1, 1, 0);
    }
    // both pull the face toward a target below 1 m/s; the target is zero only with the ln cap
    EXPECT_LT(src[1], Real(0.0));
    EXPECT_GT(src[0], src[1]) << "0.9 cap " << src[0] << ", ln cap " << src[1];
}

TEST(ImmersedWallStability, PsiMCapBeforeUstarKeepsTheWallLawInStrongInstability)
{
    // z0 = 0.5 m on 10 m cells and 1 m/s wind: psi_m passes ln(1.5 dz / z0) = ln 30, the u*
    // denominator turns negative and the u* clamp sets u* = 0. With a factor of 1 (no cap before
    // u*) the velocity target is zero; the 0.9 cap before u*, as WRF caps PSIM before UST, keeps
    // u* > 0 and the target. The temperature forcing does not take the factor (erf-model/ERF#4134):
    // its source is the same at both.
    Real xsrc[2] = {0.0, 0.0};
    Real tsrc[2] = {0.0, 0.0};
    const Real factors[2] = {Real(0.9), Real(1.0)};
    for (int n = 0; n < 2; ++n) {
        {   // momentum: 0.5 K m/s heating, L about -0.25 m, zeta at 1.5 dz about -60
            CalmTerrainColumn col;
            col.u.setVal<amrex::RunOn::Host>(1.0);
            SolverChoice sc = col.choice(Real(0.5));
            sc.if_z0 = Real(0.5);
            sc.if_psi_cap_factor = factors[n];
            const amrex::Box xbx = amrex::surroundingNodes(col.bx, 0);
            ImmersedForcingTerrain_Xmom(xbx, col.u.const_array(), col.v.const_array(), col.w.const_array(),
                                        col.cell.const_array(), col.blank.const_array(),
                                        amrex::Array4<const Real>{}, amrex::Array4<const Real>{},
                                        col.src.array(), col.geom, sc, Real(1.0));
            ASSERT_TRUE(col.all_finite(xbx, 1));
            xsrc[n] = col.src.const_array()(1, 1, 1, 0);
        }
        {   // heat: 2 K m/s heating, L at its bound -0.15 m, zeta at 0.5 dz about -33
            CalmTerrainColumn col;
            col.u.setVal<amrex::RunOn::Host>(1.0);
            SolverChoice sc = col.choice(Real(2.0));
            sc.if_z0 = Real(0.5);
            sc.if_psi_cap_factor = factors[n];
            ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                          col.cell.const_array(), col.blank.const_array(),
                                          amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                          amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
            ASSERT_TRUE(col.all_finite(col.bx, 2));
            tsrc[n] = col.src.const_array()(1, 1, 1, RhoTheta_comp);
        }
    }
    // factor 1: zero velocity target (u = 1 m/s relaxed to rest)
    EXPECT_LT(xsrc[1], Real(0.0));
    // factor 0.9: a positive target, so a weaker pull
    EXPECT_GT(xsrc[0], xsrc[1]) << "0.9 " << xsrc[0] << ", 1 " << xsrc[1];
    // the temperature forcing is development's at both factors, bit for bit
    EXPECT_EQ(tsrc[0], tsrc[1]) << "0.9 " << tsrc[0] << ", 1 " << tsrc[1];
}

TEST(ImmersedWallStability, PrescribedObukhovLengthIsUsedAsGiven)
{
    // A prescribed erf.if_Olen is not clamped to 1.5 dz / 100 (a bound that moves with the local
    // cell size); only its zeta is held within +-100. On 10 m cells that bound is 0.15 m, so
    // L = 0.12 and 0.14 m would both run as 0.15 m if the length were clamped; used as given,
    // zeta at 0.5 dz (42 and 36) and theta* differ, and so does the wall-cell heat source. (Stable
    // lengths: the psi_h caps never bind there; unstable ones this short close both brackets.)
    Real src[2] = {0.0, 0.0};
    const Real lengths[2] = {Real(0.12), Real(0.14)};
    for (int n = 0; n < 2; ++n) {
        CalmTerrainColumn col;
        col.u.setVal<amrex::RunOn::Host>(5.0);
        SolverChoice sc = col.choice(Real(1.e-8));
        sc.if_Olen_in = lengths[n];
        ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                      col.cell.const_array(), col.blank.const_array(),
                                      amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                      amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
        ASSERT_TRUE(col.all_finite(col.bx, 2));
        src[n] = col.src.const_array()(1, 1, 1, RhoTheta_comp);
    }
    EXPECT_GT(std::abs(src[0]), Real(1.e-6)) << "L = 0.12: " << src[0];
    EXPECT_NE(src[0], src[1]) << "L = 0.12: " << src[0] << ", L = 0.14: " << src[1];
}

TEST(ImmersedWallStability, PrescribedObukhovLengthTargetStaysNearTheAir)
{
    // erf.if_Olen = -0.2 m, z0 = 0.5 m, 10 m cells, 5 m/s: u* reaches the 2 m/s clamp and
    // theta* = theta u*^2 / (kappa g L) is about -1.5e3 K. With psi_h left uncapped
    // (psi_h = 4.7 and 5.8 against ln 10 and ln 30) the target was 160 K above the air; with the
    // ln(z / z0) caps of the other temperature branches both log brackets close and the target is
    // the air above.
    CalmTerrainColumn col;
    col.u.setVal<amrex::RunOn::Host>(5.0);
    SolverChoice sc = col.choice(Real(1.e-8));
    sc.if_z0 = Real(0.5);
    sc.if_Olen_in = Real(-0.2);
    ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                  col.cell.const_array(), col.blank.const_array(),
                                  amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                  amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
    ASSERT_TRUE(col.all_finite(col.bx, 2));
    amrex::Gpu::streamSynchronize();
    const Real drag_rho = Real(5.0) / Real(10.0) * Real(1.2);
    const Real dT = col.src.const_array()(1, 1, 1, RhoTheta_comp) / drag_rho;
    EXPECT_LT(std::abs(dT), Real(1.0)) << "target " << dT << " K from the air";
}
