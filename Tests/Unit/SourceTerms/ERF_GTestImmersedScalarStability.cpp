// erf.if_stability_correction and the immersed-forcing temperature kernels (erf-model/ERF#4931).
// The momentum kernels gate the MOST similarity functions on the flag (off by default: the log law,
// #4751); the temperature kernels applied them whatever the flag said, so a deck with a prescribed
// heat flux and the default flag ran the log law for momentum and a stability-corrected profile for
// theta. Now, with the flag off, every temperature branch is the log law: psi_m = psi_h = 0, so
// u* = kappa U / ln(1.5 dz / z0) and the target is theta_above - (theta* / kappa) ln 3, with theta*
// from the prescribed flux (theta* = -q / u*, through the bounded Obukhov length of the kernel) or
// from the prescribed Obukhov length (theta u*^2 / (kappa g L)). With the flag on, the sources differ
// from these neutral ones, so a kernel that ignored the flag fails the neutral comparison; the last
// test pins the flag-on source as the stability-corrected law of before.

#include <cmath>
#include <limits>

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_TableData.H>

#include <gtest/gtest.h>

#include "ERF_Constants.H"
#include "ERF_DataStruct.H"
#include "ERF_ImmersedForcing.H"
#include "ERF_ImmersedWallStability.H"
#include "ERF_IndexDefines.H"
#include "ERF_MOSTUtils.H"

using amrex::Real;

namespace {

constexpr Real dz_   = Real(10.0);
constexpr Real z0_   = Real(0.1);
constexpr Real rho_  = Real(1.2);
constexpr Real Cd_s  = Real(5.0);
constexpr Real u_    = Real(3.0);   // x wind
constexpr Real v_    = Real(2.0);   // y wind: the roof law sees sqrt(u^2 + v^2), a west wall sees v (and w = 0)
const Real roof_wind = std::sqrt(u_ * u_ + v_ * v_);

// theta of cell (i, k): varies in z (the air above differs from the wall cell) and in x (the cell
// away from a west wall differs from the wall cell)
Real theta_of (int i, int k) { return Real(300.0) + Real(2.0) * k + Real(0.5) * i; }

// A 4 x 4 x 6 box of 10 m cells, periodic in x and y, with u = 3 m/s, v = 2 m/s, w = 0 and
// rho = 1.2 everywhere; the blanking comes from the caller.
struct Column {
    amrex::Box bx{amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 5)};
    amrex::Box gbx = amrex::grow(bx, 2);
    amrex::Geometry geom;
    amrex::FArrayBox u, v, w, cell, blank, src;

    template <class F>
    explicit Column (F&& blank_of)
        : u(gbx, 1, amrex::The_Managed_Arena()), v(gbx, 1, amrex::The_Managed_Arena()),
          w(gbx, 1, amrex::The_Managed_Arena()), cell(gbx, NDRY, amrex::The_Managed_Arena()),
          blank(gbx, 1, amrex::The_Managed_Arena()), src(gbx, NDRY, amrex::The_Managed_Arena())
    {
        const amrex::RealBox rb({0.0, 0.0, 0.0}, {40.0, 40.0, 60.0});
        const amrex::Array<int, AMREX_SPACEDIM> periodic{1, 1, 0};
        geom.define(bx, rb, amrex::CoordSys::cartesian, periodic);
        u.setVal<amrex::RunOn::Host>(u_);
        v.setVal<amrex::RunOn::Host>(v_);
        w.setVal<amrex::RunOn::Host>(0.0);
        src.setVal<amrex::RunOn::Host>(0.0);
        const auto c = cell.array();
        const auto b = blank.array();
        amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) {
            c(i, j, k, Rho_comp)      = rho_;
            c(i, j, k, RhoTheta_comp) = rho_ * theta_of(i, k);
            c(i, j, k, 2)             = Real(0.0);
            b(i, j, k) = blank_of(i, j, k);
        });
    }

    Real source_at (int i, int j, int k) const
    {
        amrex::Gpu::streamSynchronize();
        return src.const_array()(i, j, k, RhoTheta_comp);
    }
};

// terrain: solid below k = 1, a half-blanked wall cell at k = 1, fluid above
Real terrain_blank (int, int, int k) { return (k < 1) ? Real(1.0) : ((k == 1) ? Real(0.5) : Real(0.0)); }

// a building: cells i = 1-2, j = 1-2 solid up to k = 2 (roof cells at k = 2), half-solid cells on
// its west side (i = 0, j = 1-2, k <= 2); the west wall law acts on (0, j, k <= 1), and (0, j, 2) is a
// roof cell too
Real building_blank (int i, int j, int k)
{
    if (k > 2 || j < 1 || j > 2) { return Real(0.0); }
    if (i == 1 || i == 2) { return Real(1.0); }
    if (i == 0) { return Real(0.5); }
    return Real(0.0);
}

SolverChoice choice (bool stability, Real tflux, Real Olen)
{
    SolverChoice sc{};
    sc.if_use_most             = true;
    sc.if_stability_correction = stability;
    sc.if_z0                   = z0_;
    sc.if_Cd_momentum          = Real(50.0);
    sc.if_Cd_scalar            = Cd_s;
    sc.if_surf_temp_flux       = tflux;
    sc.if_Olen_in              = Olen;
    return sc;
}

// The rho theta source of the log-law temperature forcing on a wall cell (theta, air theta_nb at
// 1.5 delta, tangential wind U at 1.5 delta): u* = kappa U / ln(1.5 delta / z0), held in [0, 2] m/s;
// theta* from the flux through the kernel's bounded L, or from the given L; and the drag
// Cd_s / cbrt(dx dy dz) with dx = dy = dz here.
Real neutral_source (Real theta, Real theta_nb, Real U, Real delta, Real tflux, Real Olen_in)
{
    const Real kappa = KAPPA;
    const Real ggg   = CONST_GRAV;
    const Real tiny  = std::numeric_limits<Real>::epsilon();
    const Real lnfar = std::log(Real(1.5) * delta / z0_);
    const Real lnnear = std::log(Real(0.5) * delta / z0_);
    Real ustar = U * kappa / lnfar;
    const Real Olen = (Olen_in != Real(1e-8)) ? Olen_in
                    : ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny), Real(1.5) * delta);
    ustar = ib_stability::clamped_ustar(ustar);
    const Real thetastar = theta * ustar * ustar / (kappa * ggg * Olen);
    const Real surf      = theta_nb - thetastar / kappa * lnfar;
    const Real target    = surf + thetastar / kappa * lnnear;
    const Real drag      = Cd_s / dz_;
    return drag * (rho_ * target - rho_ * theta);
}

void run_terrain (Column& col, const SolverChoice& sc)
{
    ImmersedForcingTerrain_Scalar(col.bx, col.u.const_array(), col.v.const_array(),
                                  col.cell.const_array(), col.blank.const_array(),
                                  amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                  amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
}

void run_buildings (Column& col, const SolverChoice& sc)
{
    ImmersedForcingBuildings_Scalar(col.bx, col.u.const_array(), col.v.const_array(), col.w.const_array(),
                                    col.cell.const_array(), col.blank.const_array(),
                                    amrex::Array4<const Real>{}, col.src.array(), col.geom, sc,
                                    amrex::Table1D<Real>{}, amrex::Table1D<Real>{}, Real(0.0));
}

constexpr Real rel_tol = Real(1.e-11);

} // namespace

TEST(ImmersedScalarStability, TerrainHeatFluxIsTheLogLawWithTheFlagOff)
{
    // heating and cooling fluxes; the wall cell (1, 1, 1), air above at k = 2, wind 3.6 m/s at 1.5 dz
    for (Real q : {Real(0.1), Real(-0.02)}) {
        Column off(terrain_blank), on(terrain_blank);
        run_terrain(off, choice(false, q, Real(1e-8)));
        run_terrain(on,  choice(true,  q, Real(1e-8)));
        const Real expect = neutral_source(theta_of(1, 1), theta_of(1, 2), roof_wind, dz_, q, Real(1e-8));
        EXPECT_NE(expect, Real(0.0));
        EXPECT_NEAR(off.source_at(1, 1, 1), expect, rel_tol * std::abs(expect)) << "q = " << q;
        // the flag changes the source, so a kernel ignoring it fails the line above
        EXPECT_GT(std::abs(on.source_at(1, 1, 1) - expect), Real(1.e-6) * std::abs(expect)) << "q = " << q;
        // fluid and solid cells get no wall-law source either way
        EXPECT_EQ(off.source_at(1, 1, 2), Real(0.0));
        EXPECT_EQ(off.source_at(1, 1, 0), Real(0.0));
    }
}

TEST(ImmersedScalarStability, TerrainObukhovLengthIsTheLogLawWithTheFlagOff)
{
    // a prescribed unstable and a prescribed stable length (theta* = theta u*^2 / (kappa g L) still)
    for (Real L : {Real(-20.0), Real(20.0)}) {
        Column off(terrain_blank), on(terrain_blank);
        run_terrain(off, choice(false, Real(1e-8), L));
        run_terrain(on,  choice(true,  Real(1e-8), L));
        const Real expect = neutral_source(theta_of(1, 1), theta_of(1, 2), roof_wind, dz_, Real(1e-8), L);
        EXPECT_NE(expect, Real(0.0));
        EXPECT_NEAR(off.source_at(1, 1, 1), expect, rel_tol * std::abs(expect)) << "L = " << L;
        EXPECT_GT(std::abs(on.source_at(1, 1, 1) - expect), Real(1.e-6) * std::abs(expect)) << "L = " << L;
    }
}

TEST(ImmersedScalarStability, BuildingRoofAndWallAreTheLogLawWithTheFlagOff)
{
    for (Real q : {Real(0.1), Real(-0.02)}) {
        Column off(building_blank), on(building_blank);
        run_buildings(off, choice(false, q, Real(1e-8)));
        run_buildings(on,  choice(true,  q, Real(1e-8)));
        // roof cell (1, 1, 2): air above at k = 3, wind 3.6 m/s at 1.5 dz
        const Real roof = neutral_source(theta_of(1, 2), theta_of(1, 3), roof_wind, dz_, q, Real(1e-8));
        EXPECT_NE(roof, Real(0.0));
        EXPECT_NEAR(off.source_at(1, 1, 2), roof, rel_tol * std::abs(roof)) << "roof, q = " << q;
        EXPECT_GT(std::abs(on.source_at(1, 1, 2) - roof), Real(1.e-6) * std::abs(roof)) << "roof, q = " << q;
        // west wall cell (0, 1, 0): the cell away is i = -1, the tangential wind is v = 2 m/s
        const Real wall = neutral_source(theta_of(0, 0), theta_of(-1, 0), v_, dz_, q, Real(1e-8));
        EXPECT_NE(wall, Real(0.0));
        EXPECT_NEAR(off.source_at(0, 1, 0), wall, rel_tol * std::abs(wall)) << "wall, q = " << q;
        EXPECT_GT(std::abs(on.source_at(0, 1, 0) - wall), Real(1.e-6) * std::abs(wall)) << "wall, q = " << q;
        // a cell with no wall face (interior air) gets nothing
        EXPECT_EQ(off.source_at(3, 3, 4), Real(0.0));
    }
}

TEST(ImmersedScalarStability, FlagOnIsTheStabilityCorrectedLaw)
{
    // with the flag on the terrain flux branch is the law of before: u* from the neutral bracket,
    // L from it (bounded), psi_m and psi_h of the Dyer functions at 0.5 dz and 1.5 dz, u* re-formed
    // with psi_m and clamped, psi_h capped at ln(z / z0), theta* = theta u*^2 / (kappa g L); pinned
    // here against the same formula written out, so a change of the flag-on arithmetic shows
    const Real q = Real(0.1);
    Column on(terrain_blank);
    run_terrain(on, choice(true, q, Real(1e-8)));
    const Real kappa = KAPPA, ggg = CONST_GRAV, tiny = std::numeric_limits<Real>::epsilon();
    const Real theta = theta_of(1, 1), theta_nb = theta_of(1, 2);
    const Real lnfar = std::log(Real(1.5) * dz_ / z0_), lnnear = std::log(Real(0.5) * dz_ / z0_);
    similarity_funs sfuns;
    Real ustar = roof_wind * kappa / lnfar;
    const Real Olen = ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta / (kappa * ggg * q + tiny), Real(1.5) * dz_);
    const Real zeta = ib_stability::bounded_zeta(Real(0.5) * dz_, Olen);
    const Real zeta_nb = ib_stability::bounded_zeta(Real(1.5) * dz_, Olen);
    const Real psi_m = sfuns.calc_psi_m(zeta);
    ustar = ib_stability::clamped_ustar(roof_wind * kappa / (lnfar - psi_m));
    const Real psi_h    = ib_stability::capped_psi(sfuns.calc_psi_h(zeta), Real(0.5) * dz_, z0_, Real(1.0));
    const Real psi_h_nb = ib_stability::capped_psi(sfuns.calc_psi_h(zeta_nb), Real(1.5) * dz_, z0_, Real(1.0));
    const Real thetastar = theta * ustar * ustar / (kappa * ggg * Olen);
    const Real target = theta_nb - thetastar / kappa * (lnfar - psi_h_nb) + thetastar / kappa * (lnnear - psi_h);
    const Real expect = Cd_s / dz_ * (rho_ * target - rho_ * theta);
    EXPECT_NE(psi_m, Real(0.0));
    EXPECT_NEAR(on.source_at(1, 1, 1), expect, rel_tol * std::abs(expect));
    // the heating flux makes the air unstable: psi_h > 0, so the target lies closer to the air above
    // than the log-law one
    EXPECT_LT(std::abs(target - theta_nb), std::abs(theta_nb - thetastar / kappa * std::log(Real(3.0)) - theta_nb));
}
