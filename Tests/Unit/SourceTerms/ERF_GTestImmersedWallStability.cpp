// Contract of the immersed-forcing wall law's stability bounds (erf-model/ERF#4016): the friction
// velocity behind the Obukhov length uses the wind floored at 0.1 m/s, and the Obukhov length is
// held at |L| >= 1.5 dz / 100, as the flat-ground surface layer does; psi_h is capped at
// 0.9 ln(z / z0), as WRF's revised surface layer caps it. A calm cell above the wall used to give
// u* = 0, L = 0 and zeta = z / 0; under a cooling flux zeta = +inf, psi_m = -inf and the momentum
// and temperature targets became 0 * inf = NaN. Above the floor, inside the bound and below the
// cap, nothing changes.

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

using amrex::Real;

TEST(ImmersedWallStability, WindFloorIsTheFlatGroundOne)
{
    EXPECT_EQ(ib_stability::wind_floor(), Real(0.1));
    EXPECT_EQ(ib_stability::floored_wind(Real(0.0)), Real(0.1));
    EXPECT_EQ(ib_stability::floored_wind(Real(0.05)), Real(0.1));
    // above the floor the speed is untouched, bit for bit
    for (Real s : {Real(0.1), Real(0.1000001), Real(3.7), Real(25.0)}) {
        EXPECT_EQ(ib_stability::floored_wind(s), s);
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

TEST(ImmersedWallStability, PsiHCapIsWrfs)
{
    // 0.9 ln(z / z0), so ln(z / z0) - psi_h >= 0.1 ln(z / z0) > 0 for z > z0
    for (Real z : {Real(5.0), Real(15.0)}) {
        const Real cap = ib_stability::psi_h_cap(z, Real(0.1));
        EXPECT_NEAR(cap, Real(0.9) * std::log(z / Real(0.1)), Real(1.e-5) * cap);
        EXPECT_GT(std::log(z / Real(0.1)) - cap, Real(0.0));
    }
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

    // whether every component of the source is finite, and whether the wall cells got any forcing
    bool all_finite (int ncomp) const
    {
        amrex::Gpu::streamSynchronize();
        const auto s = src.const_array();
        bool ok = true;
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            for (int n = 0; n < ncomp; ++n) { if (!std::isfinite(s(i, j, k, n))) { ok = false; } }
        });
        return ok;
    }
};

} // namespace

TEST(ImmersedWallStability, CalmTerrainMomentumTargetIsFinite)
{
    // cooling (stable: zeta = +inf on development) and heating (unstable: zeta = -inf)
    for (Real tflux : {Real(-0.05), Real(0.05)}) {
        CalmTerrainColumn col;
        const SolverChoice sc = col.choice(tflux);
        const amrex::Box xbx = amrex::surroundingNodes(col.bx, 0);
        ImmersedForcingTerrain_Xmom(xbx, col.u.const_array(), col.v.const_array(), col.w.const_array(),
                                    col.cell.const_array(), col.blank.const_array(),
                                    amrex::Array4<const Real>{}, amrex::Array4<const Real>{},
                                    col.src.array(), col.geom, sc, Real(1.0));
        EXPECT_TRUE(col.all_finite(1)) << "x-momentum, tflux = " << tflux;

        col.src.setVal<amrex::RunOn::Host>(0.0);
        const amrex::Box ybx = amrex::surroundingNodes(col.bx, 1);
        ImmersedForcingTerrain_Ymom(ybx, col.u.const_array(), col.v.const_array(), col.w.const_array(),
                                    col.cell.const_array(), col.blank.const_array(),
                                    amrex::Array4<const Real>{}, amrex::Array4<const Real>{},
                                    col.src.array(), col.geom, sc, Real(1.0));
        EXPECT_TRUE(col.all_finite(1)) << "y-momentum, tflux = " << tflux;
    }
}

TEST(ImmersedWallStability, CalmTerrainHeatFluxTargetIsFinite)
{
    for (Real tflux : {Real(-0.05), Real(0.05)}) {
        CalmTerrainColumn col;
        const SolverChoice sc = col.choice(tflux);
        ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                      col.cell.const_array(), col.blank.const_array(),
                                      amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                      amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
        EXPECT_TRUE(col.all_finite(2)) << "rho theta, tflux = " << tflux;
        // and the wall cell is forced: the floored wind carries the prescribed flux. Without the
        // floor u* = 0 and the source is round-off (theta = rho theta / rho is inexact). Under
        // heating zeta reaches -33 and psi_h passes ln(z / z0): a cap of ln(z / z0) itself would
        // flatten the profile and stop the transfer, the 0.9 ln(z / z0) cap keeps it.
        amrex::Gpu::streamSynchronize();
        EXPECT_GT(std::abs(col.src.const_array()(1, 1, 1, RhoTheta_comp)), Real(1.e-6)) << "tflux = " << tflux;
    }
}
