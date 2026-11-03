#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_DataStruct.H"
#include "ERF_ImmersedForcing.H"
#include "ERF_ImmersedWallCell.H"
#include "ERF_IndexDefines.H"
#include "ERF_Utils.H"

// Motivation: the fraction-stress wall law of the terrain immersed forcing
// (erf.if_wall_form = fraction_stress, kynema-sgf form 7) must put the wall stress in the wall
// cell of each column - the partial cell under the air, or on a face-aligned surface the fluid
// cell on the solid - close that cell's bottom face to diffusion, and measure the RANS wall
// distance from the surface itself. With the old law a face-aligned surface got no wall forcing
// at all and the wind of a neutral ABL came out 40-50 % slow.

using namespace amrex;

namespace {

constexpr int  nx = 4;
constexpr int  ny = 2;
constexpr int  nz = 12;
constexpr Real dz = Real(10.);
constexpr Real dxy = Real(20.);
constexpr Real small = Real(0.005);
constexpr Real z0 = Real(0.1);

Real tol_rel () { return Real(64.) * std::numeric_limits<Real>::epsilon(); }

// Solid fraction of a flat surface at height h in cell k of height dz
Real plane_beta (Real h, int k)
{
    return amrex::min(Real(1.), amrex::max(Real(0.), (h - dz*k) / dz));
}

// ---------------------------------------------------------------- wall cell and the law

TEST(ImmersedFractionStress, WallCellIsThePartialCellOrTheFluidCellOnTheSolid)
{
    // partial cell under fluid
    EXPECT_TRUE (ib_wall::is_wall_cell(Real(1.),   Real(0.25), Real(0.),  small));
    EXPECT_TRUE (ib_wall::is_wall_cell(Real(0.4),  Real(0.25), Real(0.),  small));   // partial on partial
    // face-aligned surface: the fluid cell on the solid cell
    EXPECT_TRUE (ib_wall::is_wall_cell(Real(1.),   Real(0.),   Real(0.),  small));
    EXPECT_TRUE (ib_wall::is_wall_cell(Real(0.996),Real(0.),   Real(0.),  small));   // solid within small of one
    // not wall cells
    EXPECT_FALSE(ib_wall::is_wall_cell(Real(1.),   Real(1.),   Real(0.),  small));   // the solid cell itself
    EXPECT_FALSE(ib_wall::is_wall_cell(Real(0.25), Real(0.),   Real(0.),  small));   // fluid above a partial cell
    EXPECT_FALSE(ib_wall::is_wall_cell(Real(1.),   Real(0.25), Real(0.2), small));   // partial under a partial
    EXPECT_FALSE(ib_wall::is_wall_cell(Real(0.),   Real(0.),   Real(0.),  small));   // open air
    // the thresholds
    EXPECT_TRUE (ib_wall::is_wall_cell(Real(1.),   small,      Real(0.),  small));   // at small: partial
    EXPECT_FALSE(ib_wall::is_wall_cell(Real(0.99), Real(0.004),Real(0.),  small));   // below small: fluid, and 0.99 is not solid
    EXPECT_DOUBLE_EQ(static_cast<double>(ib_wall::wall_beta(Real(0.004), small)), 0.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(ib_wall::wall_beta(Real(0.3),   small)), static_cast<double>(Real(0.3)));
}

TEST(ImmersedFractionStress, FrictionVelocityMixesOwnAndReferenceByTheFraction)
{
    const Real ut = Real(4.), ut_ref = Real(6.);
    // beta = 0: the cell's own velocity at dz/2, the flat-ground surface-layer law
    const Real us0 = ib_wall::friction_velocity(ut, ut_ref, Real(0.), dz, z0);
    EXPECT_NEAR(us0, KAPPA * ut / std::log(Real(5.) / z0), tol_rel() * us0);
    // general beta: u*^2 = (1-b) u*_own^2 + b u*_ref^2 with d1 = (1-b) dz/2, d2 = (3/2-b) dz
    const Real b = Real(0.3);
    const Real own = KAPPA * ut     / std::log(Real(0.5) * (Real(1.) - b) * dz / z0);
    const Real ref = KAPPA * ut_ref / std::log((Real(1.5) - b) * dz / z0);
    const Real usb = ib_wall::friction_velocity(ut, ut_ref, b, dz, z0);
    EXPECT_NEAR(usb, std::sqrt((Real(1.) - b)*own*own + b*ref*ref), tol_rel() * usb);
    // a sliver: finite, and taken almost entirely from the cell above
    const Real us1 = ib_wall::friction_velocity(ut, ut_ref, Real(0.9999), dz, z0);
    EXPECT_TRUE(std::isfinite(us1));
    EXPECT_NEAR(us1, KAPPA * ut_ref / std::log(Real(0.5001) * dz / z0), Real(1e-2) * us1);
}

TEST(ImmersedFractionStress, StressRateIsTheExactIntegralOfTheWallStress)
{
    const Real us = Real(0.4), ut = Real(5.);
    const Real r = us * us / (dz * ut);
    EXPECT_NEAR(ib_wall::stress_rate(us, ut, dz, Real(2.)), (Real(1.) - std::exp(-r * Real(2.))) / Real(2.),
                tol_rel() * r);
    // small dt: the rate itself; dt = 0 stays finite
    EXPECT_NEAR(ib_wall::stress_rate(us, ut, dz, Real(1e-8)), r, Real(1e-6) * r);
    EXPECT_TRUE(std::isfinite(ib_wall::stress_rate(us, ut, dz, Real(0.))));
    // never removes more than the velocity in one step: C dt <= 1
    EXPECT_LE(ib_wall::stress_rate(us, Real(1e-12), dz, Real(1.)) * Real(1.), Real(1.));
}

TEST(ImmersedFractionStress, FaceAboveFactorKeepsFlatGroundAndCorrectsPartialCells)
{
    for (int sp : {0, 1, 2}) {
        SCOPED_TRACE(testing::Message() << "spacing = " << sp);
        EXPECT_NEAR(ib_wall::face_above_factor(Real(0.), dz, z0, sp), Real(1.), tol_rel());
    }
    const Real b = Real(0.5);
    const Real d1 = Real(0.5) * (Real(1.) - b) * dz, d2 = (Real(1.5) - b) * dz;
    EXPECT_EQ(ib_wall::face_above_factor(b, dz, z0, 0), Real(1.));
    EXPECT_NEAR(ib_wall::face_above_factor(b, dz, z0, 1), dz / (d2 - d1), tol_rel() * Real(4.));
    EXPECT_NEAR(ib_wall::face_above_factor(b, dz, z0, 2),
                dz / ((d1 + d2) * std::log(d2 / d1) / (Real(2.) * std::log(Real(3.)))), tol_rel() * Real(8.));
}

// ---------------------------------------------------------------- stratified wall state

TEST(ImmersedFractionStress, NeutralWallStateIsTheNeutralFrictionVelocity)
{
    ib_wall::WallCond c;   // neutral
    for (Real b : {Real(0.), Real(0.3), Real(0.8)}) {
        const ib_wall::WallState w = ib_wall::wall_state(Real(4.), Real(6.), Real(300.), Real(300.), b, dz, z0, c);
        SCOPED_TRACE(testing::Message() << "beta = " << b);
        EXPECT_NEAR(w.ustar, ib_wall::friction_velocity(Real(4.), Real(6.), b, dz, z0), tol_rel() * w.ustar);
        EXPECT_EQ(w.thetastar, Real(0.));
        EXPECT_EQ(w.Linv, Real(0.));
    }
}

// A given heat flux: theta* = -q/u*, and the Obukhov length the iteration settles on is the one
// of that flux and u*, L = -u*^3 theta / (kappa g q); u* then carries the unstable correction
TEST(ImmersedFractionStress, HeatFluxWallStateIsSelfConsistent)
{
    ib_wall::WallCond c;
    c.type = ib_wall::WallCond::heat_flux;
    c.q = Real(0.05);
    const Real th = Real(300.);
    const ib_wall::WallState w = ib_wall::wall_state(Real(8.), Real(10.), th, th, Real(0.25), dz, z0, c);
    EXPECT_NEAR(w.thetastar, -c.q / w.ustar, tol_rel() * std::abs(w.thetastar));
    const Real L = -w.ustar * w.ustar * w.ustar * th / (KAPPA * CONST_GRAV * c.q);
    EXPECT_LT(w.Linv, Real(0.));
    EXPECT_NEAR(Real(1.) / w.Linv, L, Real(1e-3) * std::abs(L));
    EXPECT_GT(w.ustar, ib_wall::friction_velocity(Real(8.), Real(10.), Real(0.25), dz, z0));
}

// A given surface temperature: theta* mixes the wall cell's and the cell above's estimates by the
// fraction, with the stability functions of the settled L; a warm surface gives theta* < 0, L < 0
TEST(ImmersedFractionStress, SurfaceTemperatureWallStateIsSelfConsistent)
{
    ib_wall::WallCond c;
    c.type = ib_wall::WallCond::surface_temp;
    c.theta_s = Real(303.);
    const Real b = Real(0.5), th = Real(300.2), th_ref = Real(300.);
    const ib_wall::WallState w = ib_wall::wall_state(Real(6.), Real(8.), th, th_ref, b, dz, z0, c);
    similarity_funs sf;
    const Real d1 = Real(0.5) * (Real(1.) - b) * dz, d2 = (Real(1.5) - b) * dz;
    const Real ts = (Real(1.) - b) * KAPPA * (th - c.theta_s) / (std::log(d1/z0) - sf.calc_psi_h(d1 * w.Linv))
                  +            b  * KAPPA * (th_ref - c.theta_s) / (std::log(d2/z0) - sf.calc_psi_h(d2 * w.Linv));
    EXPECT_LT(w.thetastar, Real(0.));
    EXPECT_LT(w.Linv, Real(0.));
    EXPECT_NEAR(w.thetastar, ts, Real(1e-3) * std::abs(ts));
    EXPECT_NEAR(w.Linv, KAPPA * CONST_GRAV * w.thetastar / (w.ustar * w.ustar * th), Real(1e-3) * std::abs(w.Linv));
}

// A given Obukhov length is used as it stands: theta* = theta u*^2 / (kappa g L)
TEST(ImmersedFractionStress, ObukhovLengthWallStateUsesTheGivenLength)
{
    ib_wall::WallCond c;
    c.type = ib_wall::WallCond::obukhov;
    c.Linv = Real(1.) / Real(150.);
    const Real th = Real(265.);
    const ib_wall::WallState w = ib_wall::wall_state(Real(5.), Real(6.), th, th, Real(0.), dz, z0, c);
    similarity_funs sf;
    const Real us = KAPPA * Real(5.) / (std::log(Real(0.5) * dz / z0) - sf.calc_psi_m(Real(0.5) * dz * c.Linv));
    EXPECT_EQ(w.Linv, c.Linv);
    EXPECT_NEAR(w.ustar, us, tol_rel() * us);
    EXPECT_NEAR(w.thetastar, th * us * us * c.Linv / (KAPPA * CONST_GRAV), tol_rel() * w.thetastar * Real(4.));
}

// With a surface temperature the wall heat source is a decay toward it, integrated exactly: it can
// never carry the air past the surface temperature, even with a step far above the decay time
TEST(ImmersedFractionStress, SurfaceTemperatureHeatSourceNeverOvershoots)
{
    ib_wall::WallCond c;
    c.type = ib_wall::WallCond::surface_temp;
    c.theta_s = Real(303.);
    const Real th = Real(300.);
    const ib_wall::WallState w = ib_wall::wall_state(Real(6.), Real(8.), th, th, Real(0.9), dz, z0, c);
    for (Real dt : {Real(1.), Real(100.), Real(1e6)}) {
        const Real s = ib_wall::heat_source(w, th, c, dz, dt);
        SCOPED_TRACE(testing::Message() << "dt = " << dt);
        EXPECT_GT(s, Real(0.));
        EXPECT_LE(th + s * dt, c.theta_s * (Real(1.) + tol_rel()));
    }
    // a given flux is applied as it stands
    ib_wall::WallCond cq; cq.type = ib_wall::WallCond::heat_flux; cq.q = Real(0.05);
    const ib_wall::WallState wq = ib_wall::wall_state(Real(6.), Real(8.), th, th, Real(0.), dz, z0, cq);
    EXPECT_NEAR(ib_wall::heat_source(wq, th, cq, dz, Real(1.)), cq.q / dz, tol_rel() * cq.q / dz * Real(4.));
}

// ---------------------------------------------------------------- fields on a mesh

struct Mesh
{
    Box domain;
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
};

Mesh make_mesh ()
{
    Mesh m;
    m.domain = Box(IntVect(0,0,0), IntVect(nx-1,ny-1,nz-1));
    const RealBox rb(Real(0.), Real(0.), Real(0.), dxy*nx, dxy*ny, dz*nz);
    m.geom = Geometry(m.domain, rb, CoordSys::cartesian, {1, 1, 0});
    m.ba = BoxArray(m.domain);
    m.ba.maxSize(IntVect(2,2,nz));
    m.dm = DistributionMapping(m.ba);
    return m;
}

// beta(i,k) from a table (the same in every j), on cells with ghosts
void fill_beta (MultiFab& beta, const std::vector<Real>& table)
{
    Gpu::DeviceVector<Real> bt(table.size());
    Gpu::copy(Gpu::hostToDevice, table.begin(), table.end(), bt.begin());
    const Real* bp = bt.data();
    const int n_x = nx, n_z = nz;
    for (MFIter mfi(beta); mfi.isValid(); ++mfi) {
        auto const b = beta.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ii = ((i % n_x) + n_x) % n_x;                    // periodic in x
            const int kk = amrex::max(0, amrex::min(k, n_z-1));
            b(i,j,k) = bp[ii*n_z + kk];
        });
    }
    Gpu::streamSynchronize();
}

// One valid-region copy of a MultiFab on the host, (i,j,k) -> ((k-klo)*ny + j)*nxx + i
std::vector<Real> to_host (const MultiFab& mf, const Box& region)
{
    MultiFab one(BoxArray(region), DistributionMapping(BoxArray(region)), 1, 0);
    one.ParallelCopy(mf, 0, 0, 1);
    std::vector<Real> out(static_cast<std::size_t>(region.numPts()), Real(-99.));
    for (MFIter mfi(one); mfi.isValid(); ++mfi) {
        Gpu::copy(Gpu::deviceToHost, one[mfi].dataPtr(), one[mfi].dataPtr() + one[mfi].size(), out.begin());
    }
    Gpu::streamSynchronize();
    return out;
}

Real at (const std::vector<Real>& a, const Box& region, int i, int j, int k)
{
    const IntVect lo = region.smallEnd(), len = region.length();
    return a[static_cast<std::size_t>(((k-lo[2])*len[1] + (j-lo[1]))*len[0] + (i-lo[0]))];
}

TEST(ImmersedFractionStress, WallFaceMasksCloseTheBottomFaceOfTheWallCell)
{
    const Mesh m = make_mesh();
    // column i: 0 face-aligned at 40 m, 1 partial 0.25 at k = 4, 2 partial 0.75 at k = 4, 3 no solid
    std::vector<Real> table(static_cast<std::size_t>(nx*nz), Real(0.));
    for (int k = 0; k < nz; ++k) {
        table[0*nz + k] = plane_beta(Real(40.), k);
        table[1*nz + k] = plane_beta(Real(42.5), k);
        table[2*nz + k] = plane_beta(Real(47.5), k);
    }
    MultiFab beta(m.ba, m.dm, 1, 3);
    fill_beta(beta, table);
    MultiFab m13(convert(m.ba, IntVect(1,0,1)), m.dm, 1, 1);
    MultiFab m23(convert(m.ba, IntVect(0,1,1)), m.dm, 1, 1);
    MultiFab m33(convert(m.ba, IntVect(0,0,1)), m.dm, 1, 1);
    make_ib_wall_face_masks(beta, m13, m23, m33, m.domain, small);

    // x-face column i sees the mean of cell columns i-1 and i (periodic: face 0 pairs 3 and 0)
    const Box r13 = convert(m.domain, IntVect(1,0,1));
    const Box r23 = convert(m.domain, IntVect(0,1,1));
    const std::vector<Real> h13 = to_host(m13, r13);
    const std::vector<Real> h23 = to_host(m23, r23);
    auto face_beta = [&] (int ia, int ib, int k) {
        const int kc = std::max(0, std::min(k, nz-1));
        const Real v = Real(0.5) * (table[ia*nz + kc] + table[ib*nz + kc]);
        return (v < small) ? Real(0.) : v;
    };
    for (int i = 0; i <= nx; ++i) {
        const int ia = ((i-1) % nx + nx) % nx, ib = i % nx;
        for (int k = 0; k <= nz; ++k) {
            const bool wall = (k < nz) && ib_wall::is_wall_cell(face_beta(ia,ib,k-1), face_beta(ia,ib,k), face_beta(ia,ib,k+1), small);
            SCOPED_TRACE(testing::Message() << "xz edge i = " << i << ", k = " << k);
            EXPECT_EQ(at(h13, r13, i, 0, k), wall ? Real(0.) : Real(1.));
        }
    }
    // y-face columns and the cell columns (z faces, theta and TKE) sit over a single cell column:
    // the wall face is the bottom of its wall cell
    const Box r33 = convert(m.domain, IntVect(0,0,1));
    const std::vector<Real> h33 = to_host(m33, r33);
    const int expect_kw[nx] = {4, 4, 4, -1};   // 40 m: fluid cell 4 on solid cell 3; partials at k = 4
    for (int i = 0; i < nx; ++i) {
        for (int k = 0; k <= nz; ++k) {
            SCOPED_TRACE(testing::Message() << "yz edge / z face i = " << i << ", k = " << k);
            EXPECT_EQ(at(h23, r23, i, 1, k), (k == expect_kw[i]) ? Real(0.) : Real(1.));
            EXPECT_EQ(at(h33, r33, i, 1, k), (k == expect_kw[i]) ? Real(0.) : Real(1.));
        }
    }
}

// Node and cell-centre heights of the flat mesh, dz apart
void fill_flat_mesh (MultiFab& z_nd, MultiFab& z_cc)
{
    const Real h = dz;
    for (MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
        auto const z = z_nd.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { z(i,j,k) = h * k; });
    }
    for (MFIter mfi(z_cc); mfi.isValid(); ++mfi) {
        auto const z = z_cc.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { z(i,j,k) = h * (k + Real(0.5)); });
    }
    Gpu::streamSynchronize();
}

TEST(ImmersedFractionStress, WallDistanceIsMeasuredFromTheSurface)
{
    const Mesh m = make_mesh();
    const std::vector<Real> hs = {Real(40.), Real(42.5), Real(45.), Real(47.5)};
    std::vector<Real> table(static_cast<std::size_t>(nx*nz), Real(0.));
    for (int i = 0; i < nx; ++i) for (int k = 0; k < nz; ++k) { table[i*nz + k] = plane_beta(hs[i], k); }
    MultiFab beta(m.ba, m.dm, 1, 1);
    fill_beta(beta, table);
    BoxArray ba_nd(m.ba); ba_nd.surroundingNodes();
    MultiFab z_nd(ba_nd, m.dm, 1, 1), z_cc(m.ba, m.dm, 1, 1), wd(m.ba, m.dm, 1, 1);
    fill_flat_mesh(z_nd, z_cc);
    immersed_wall_dist(wd, beta, z_nd, z_cc, m.domain, small);
    const std::vector<Real> d = to_host(wd, m.domain);

    for (int i = 0; i < nx; ++i) {
        const Real h = hs[i];
        const int kw = 4;                         // 40 m: fluid cell 4; others: partial cell 4
        for (int k = 0; k < nz; ++k) {
            const Real zc = dz * (k + Real(0.5));
            const Real expect = (k == kw) ? Real(0.5) * (dz * (kw + 1) - h) : std::abs(zc - h);
            SCOPED_TRACE(testing::Message() << "h = " << h << ", k = " << k);
            EXPECT_NEAR(at(d, m.domain, i, 0, k), expect, tol_rel() * dz * nz);
        }
    }
}

// ---------------------------------------------------------------- the momentum source

struct MomFields
{
    MultiFab cons, u, v, w, beta, src;
};

// u(k) = 1 + 0.3 k, v = 0.5, w = 0, rho = 1.2, surface at height h everywhere
MomFields make_mom_fields (const Mesh& m, Real h)
{
    MomFields f{MultiFab(m.ba, m.dm, NVAR_max, 2),
                MultiFab(convert(m.ba, IntVect(1,0,0)), m.dm, 1, 2),
                MultiFab(convert(m.ba, IntVect(0,1,0)), m.dm, 1, 2),
                MultiFab(convert(m.ba, IntVect(0,0,1)), m.dm, 1, 2),
                MultiFab(m.ba, m.dm, 1, 2),
                MultiFab(convert(m.ba, IntVect(1,0,0)), m.dm, 1, 0)};
    f.cons.setVal(Real(0.));
    f.cons.setVal(Real(1.2), Rho_comp, 1, 2);
    f.v.setVal(Real(0.5));
    f.w.setVal(Real(0.));
    f.src.setVal(Real(0.));
    for (MFIter mfi(f.u); mfi.isValid(); ++mfi) {
        auto const u = f.u.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = Real(1.) + Real(0.3) * k; });
    }
    std::vector<Real> table(static_cast<std::size_t>(nx*nz));
    for (int i = 0; i < nx; ++i) for (int k = 0; k < nz; ++k) { table[i*nz + k] = plane_beta(h, k); }
    fill_beta(f.beta, table);
    return f;
}

std::vector<Real> run_xmom (const Mesh& m, MomFields& f, const SolverChoice& sc, Real dt)
{
    for (MFIter mfi(f.src); mfi.isValid(); ++mfi) {
        ImmersedForcingTerrain_Xmom_FractionStress(mfi.validbox(), f.u.const_array(mfi), f.v.const_array(mfi),
                                                   f.w.const_array(mfi), f.cons.const_array(mfi),
                                                   f.beta.const_array(mfi), Array4<const Real>{},
                                                   f.src.array(mfi), m.geom, sc, dt, Real(0.));
    }
    Gpu::streamSynchronize();
    return to_host(f.src, convert(m.domain, IntVect(1,0,0)));
}

SolverChoice make_choice ()
{
    SolverChoice sc;
    sc.if_Cd_momentum   = Real(50.);
    sc.if_z0            = z0;
    sc.if_implicit_drag = true;
    return sc;
}

// On a face-aligned surface the wall cell is the fluid cell on the solid cell, and the stress it
// feels is the surface-layer stress of flat ground at dz/2: rho u*^2 u/|u_t| / dz for small dt.
TEST(ImmersedFractionStress, FaceAlignedSurfaceGetsTheFlatGroundStress)
{
    const Mesh m = make_mesh();
    MomFields f = make_mom_fields(m, Real(40.));
    const SolverChoice sc = make_choice();
    const Real dt = Real(1e-6);
    const std::vector<Real> s = run_xmom(m, f, sc, dt);
    const Box rx = convert(m.domain, IntVect(1,0,0));

    const Real u4 = Real(1.) + Real(0.3) * 4, ut4 = std::sqrt(u4*u4 + Real(0.25));
    const Real us = KAPPA * ut4 / std::log(Real(0.5) * dz / z0);
    const Real flat = -Real(1.2) * us * us / dz * u4 / ut4;
    for (int i = 0; i <= nx; ++i) {
        SCOPED_TRACE(testing::Message() << "i = " << i);
        EXPECT_NEAR(at(s, rx, i, 0, 4), flat, Real(1e-5) * std::abs(flat));          // wall cell
        EXPECT_EQ  (at(s, rx, i, 0, 5), Real(0.));                                   // air above
        EXPECT_LT  (at(s, rx, i, 0, 3), Real(0.));                                   // solid: drag
    }
}

// A partial cell is the wall cell: the fraction-stress source there, the plain drag in the solid
// below, nothing in the air above.
TEST(ImmersedFractionStress, PartialCellCarriesTheStressAndTheSolidTheDrag)
{
    const Mesh m = make_mesh();
    const Real h = Real(45.), b = Real(0.5);
    MomFields f = make_mom_fields(m, h);
    const SolverChoice sc = make_choice();
    const Real dt = Real(0.5);
    const std::vector<Real> s = run_xmom(m, f, sc, dt);
    const Box rx = convert(m.domain, IntVect(1,0,0));

    auto u_at = [] (int k) { return Real(1.) + Real(0.3) * k; };
    const Real ut4 = std::sqrt(u_at(4)*u_at(4) + Real(0.25));
    const Real ut5 = std::sqrt(u_at(5)*u_at(5) + Real(0.25));
    const Real us  = ib_wall::friction_velocity(ut4, ut5, b, dz, z0);
    const Real wall = -ib_wall::stress_rate(us, ut4, dz, dt) * Real(1.2) * u_at(4);

    // drag in the solid cell 3: lambda = beta Cd/Delta |U| capped, point-implicit
    const Real ut3 = std::sqrt(u_at(3)*u_at(3) + Real(0.25));
    const Real cd  = Real(50.) / std::cbrt(dxy * dxy * dz);
    const Real lam = std::min(cd / ut3, cd) * ut3;
    const Real drag = -lam / (Real(1.) + lam * dt) * Real(1.2) * u_at(3);

    for (int i = 0; i <= nx; ++i) {
        SCOPED_TRACE(testing::Message() << "i = " << i);
        EXPECT_NEAR(at(s, rx, i, 0, 4), wall, tol_rel() * std::abs(wall) * Real(16.));
        EXPECT_NEAR(at(s, rx, i, 0, 3), drag, tol_rel() * std::abs(drag) * Real(16.));
        EXPECT_EQ  (at(s, rx, i, 0, 5), Real(0.));
        EXPECT_EQ  (at(s, rx, i, 0, 8), Real(0.));
    }
}


// The wall heat source of the scalar kernel in the wall cell (cell 4 over a plane at 45 m), with the
// surface temperature erf.if_init_surf_temp changing at erf.if_surf_heating_rate (stored in K/s)
Real wall_heat_source_at (const Mesh& m, MomFields& f, const SolverChoice& sc, Real dt, Real time)
{
    MultiFab src(m.ba, m.dm, NVAR_max, 0);
    src.setVal(Real(0.));
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        ImmersedForcingTerrain_Scalar_FractionStress(mfi.validbox(), f.u.const_array(mfi), f.v.const_array(mfi),
                                                     f.cons.const_array(mfi), f.beta.const_array(mfi),
                                                     src.array(mfi), Array4<Real>{}, Array4<Real>{}, m.geom, sc,
                                                     Table1D<Real>{}, Table1D<Real>{}, dt, time);
    }
    Gpu::streamSynchronize();
    MultiFab one_comp(m.ba, m.dm, 1, 0);
    MultiFab::Copy(one_comp, src, RhoTheta_comp, 0, 1, 0);
    return at(to_host(one_comp, m.domain), m.domain, 1, 0, 4);
}

void set_theta (MultiFab& cons, Real theta)
{
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        auto const c = cons.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            c(i,j,k,RhoTheta_comp) = c(i,j,k,Rho_comp) * theta;
        });
    }
    Gpu::streamSynchronize();
}

TEST(ImmersedFractionStress, SurfaceTemperatureFollowsTheHeatingRate)
{
    const Mesh m = make_mesh();
    MomFields f = make_mom_fields(m, Real(45.));
    set_theta(f.cons, Real(300.));
    SolverChoice sc = make_choice();
    sc.if_init_surf_temp    = Real(303.);
    sc.if_surf_heating_rate = Real(-2.) / Real(3600.);   // -2 K/h, as stored after reading
    const Real dt = Real(1.);

    auto expected = [&] (Real theta_s) {
        ib_wall::WallCond c; c.type = ib_wall::WallCond::surface_temp; c.theta_s = theta_s;
        auto speed = [] (int k) { const Real u = Real(1.) + Real(0.3) * k; return std::sqrt(u*u + Real(0.25)); };
        // the kernel averages the staggered u to the cell centre: u(k) is uniform in x here
        const ib_wall::WallState w = ib_wall::wall_state(speed(4), speed(5), Real(300.), Real(300.),
                                                         Real(0.5), dz, z0, c);
        return Real(1.2) * ib_wall::heat_source(w, Real(300.), c, dz, dt);
    };
    const Real s0 = wall_heat_source_at(m, f, sc, dt, Real(0.));
    const Real s1 = wall_heat_source_at(m, f, sc, dt, Real(3600.));
    EXPECT_NEAR(s0, expected(Real(303.)), Real(1e-10) * std::abs(s0));
    EXPECT_NEAR(s1, expected(Real(301.)), Real(1e-10) * std::abs(s1));
    EXPECT_GT(s0, s1);   // a cooler surface heats the air less
}

TEST(ImmersedFractionStress, WallFaceMasksCarryTheSpacingFactorAboveTheWallCell)
{
    const Mesh m = make_mesh();
    std::vector<Real> table(static_cast<std::size_t>(nx*nz), Real(0.));
    for (int i = 0; i < nx; ++i) for (int k = 0; k < nz; ++k) { table[i*nz + k] = plane_beta(Real(45.), k); }
    MultiFab beta(m.ba, m.dm, 1, 3);
    fill_beta(beta, table);
    MultiFab m13(convert(m.ba, IntVect(1,0,1)), m.dm, 1, 1);
    MultiFab m23(convert(m.ba, IntVect(0,1,1)), m.dm, 1, 1);
    MultiFab m33(convert(m.ba, IntVect(0,0,1)), m.dm, 1, 1);
    make_ib_wall_face_masks(beta, m13, m23, m33, m.domain, small, 1, dz, z0);
    const Box r33 = convert(m.domain, IntVect(0,0,1));
    const std::vector<Real> h33 = to_host(m33, r33);
    const Real g = ib_wall::face_above_factor(Real(0.5), dz, z0, 1);
    for (int k = 0; k <= nz; ++k) {
        SCOPED_TRACE(testing::Message() << "k = " << k);
        EXPECT_NEAR(at(h33, r33, 2, 1, k), (k == 4) ? Real(0.) : ((k == 5) ? g : Real(1.)), tol_rel());
    }
}

} // namespace

// erf.if_wall_tke: the wall value of the k-eqn TKE is Axell & Liungman's Eq. 16, the value
// erf.dirichlet_k holds over flat ground. Neutral and stable: u*^2 / Cmu0^2; unstable adds the
// buoyancy term kappa B d1 with d1 the fluid-centroid height of the wall cell (floored at 2 z0).
TEST(ImmersedFractionStress, WallTkeIsTheFlatGroundWallValue)
{
    const Real Cmu0 = Real(0.5562), theta = Real(300.);
    ib_wall::WallState w{};
    w.ustar = Real(0.4); w.thetastar = Real(0.);
    const Real k_neutral = w.ustar * w.ustar / (Cmu0 * Cmu0);
    for (Real beta : {Real(0.), Real(0.3), Real(0.9)}) {
        EXPECT_NEAR(ib_wall::wall_tke(w, beta, dz, z0, theta, Cmu0), k_neutral, Real(1.e-12));
    }
    w.thetastar = Real(0.05);   // stable: no buoyancy term
    EXPECT_NEAR(ib_wall::wall_tke(w, Real(0.3), dz, z0, theta, Cmu0), k_neutral, Real(1.e-12));
    w.thetastar = Real(-0.05);  // unstable
    const Real B  = CONST_GRAV * w.ustar * Real(0.05) / theta;
    for (Real beta : {Real(0.), Real(0.5)}) {
        const Real d1 = Real(0.5) * (Real(1.) - beta) * dz;
        const Real expect = std::pow(w.ustar*w.ustar*w.ustar + KAPPA * B * d1, Real(2.)/Real(3.)) / (Cmu0*Cmu0);
        EXPECT_NEAR(ib_wall::wall_tke(w, beta, dz, z0, theta, Cmu0), expect, Real(1.e-12));
        EXPECT_GT(ib_wall::wall_tke(w, beta, dz, z0, theta, Cmu0), k_neutral);
    }
    // a sliver of fluid: d1 is floored at 2 z0 as in the wall state
    const Real expect_floor = std::pow(w.ustar*w.ustar*w.ustar + KAPPA * B * Real(2.) * z0, Real(2.)/Real(3.)) / (Cmu0*Cmu0);
    EXPECT_NEAR(ib_wall::wall_tke(w, Real(0.999), dz, z0, theta, Cmu0), expect_floor, Real(1.e-12));
}
