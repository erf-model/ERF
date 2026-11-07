#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_DataStruct.H"
#include "ERF_ImmersedBuildingMasks.H"
#include "ERF_ImmersedForcing.H"
#include "ERF_IndexDefines.H"

// Motivation (erf-model/ERF#4887): with erf.if_implicit_projection the linear drag of the
// building immersed forcing toward zero is applied inside the anelastic projection
// (ERF::make_if_projection_sigma), so the source kernels must drop it: the faces inside the
// building and the partial faces without a wall law get no source, the partial faces carrying
// the MOST wall law keep exactly the source of the explicit path, and a face between two solid
// cells, which the explicit path zeroes after the projection, gets no wall law either. If a kernel
// kept its drag the solid faces would be dragged twice; if it dropped the wall law the walls would
// lose their log-law target. The last test checks that the drag the kernels drop is, face by
// face, the weight the projection reads through ERF_ImmersedBuildingMasks.H, with and without the
// snap of the partial cells and with and without face-centred blanking.

using namespace amrex;

namespace {

constexpr int  nx = 8, ny = 6, nz = 6;
constexpr Real dx = Real(10.), dz = Real(10.);
constexpr Real small = Real(0.005);
constexpr Real rho0 = Real(1.1);

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
    const RealBox rb(Real(0.), Real(0.), Real(0.), dx*nx, dx*ny, dz*nz);
    m.geom = Geometry(m.domain, rb, CoordSys::cartesian, {1, 1, 0});
    m.ba = BoxArray(m.domain);
    m.ba.maxSize(IntVect(4,3,nz));
    m.dm = DistributionMapping(m.ba);
    return m;
}

// Solid fraction of a box building: cells 3-4 in x and 2-3 in y fully solid, half-solid cells on
// its west (i = 2), east (i = 5) and north (j = 4) faces, an exact south wall on the face between
// cells j = 1 and 2 (so the fully solid x-faces of the wall row carry the south wall law), solid up
// to k = 2 and a half-solid roof row at k = 3; a sliver (0.006, so every face it touches is below
// the threshold small and counts as fluid) in cell (7, 0, 0); air elsewhere
AMREX_GPU_HOST_DEVICE inline Real building_beta (int i, int j, int k)
{
    if (i == 7 && j == 0 && k == 0) { return Real(0.006); }
    const Real fx = (i == 3 || i == 4) ? Real(1.) : ((i == 2 || i == 5) ? Real(0.5) : Real(0.));
    const Real fy = (j == 2 || j == 3) ? Real(1.) : ((j == 4) ? Real(0.5) : Real(0.));
    const Real fz = (k <= 2) ? Real(1.) : ((k == 3) ? Real(0.5) : Real(0.));
    return fx * fy * fz;
}

// A face-centred blanking that differs from the mean of the two cells: the exact solid fraction of
// the face of the box (x from 25 to 55 m, y from 20 to 45 m, z up to 35 m), so a face at the edge
// of a half cell is open where the cell mean is a quarter, and the south wall face is solid
AMREX_GPU_HOST_DEVICE inline Real building_face_beta (int i, int j, int k, int dir)
{
    const Real fx = (dir == 0) ? ((i >= 3 && i <= 5) ? Real(1.) : Real(0.))
                               : ((i == 3 || i == 4) ? Real(1.) : ((i == 2 || i == 5) ? Real(0.5) : Real(0.)));
    const Real fy = (dir == 1) ? ((j >= 2 && j <= 4) ? Real(1.) : Real(0.))
                               : ((j == 2 || j == 3) ? Real(1.) : ((j == 4) ? Real(0.5) : Real(0.)));
    const Real fz = (dir == 2) ? ((k <= 3) ? Real(1.) : Real(0.))
                               : ((k <= 2) ? Real(1.) : ((k == 3) ? Real(0.5) : Real(0.)));
    return fx * fy * fz;
}

// The wind of the test: sheared and swirling, every component nonzero everywhere and above 1 m/s
// in the solid, so the kernels' cap min(1, |u|) is off there
AMREX_GPU_HOST_DEVICE inline Real u_of (int, int j, int k) { return Real(2.) + Real(0.4) * k + Real(0.1) * std::sin(Real(0.7) * j); }
AMREX_GPU_HOST_DEVICE inline Real v_of (int i, int, int k) { return Real(0.8) + Real(0.1) * k + Real(0.2) * std::cos(Real(0.5) * i); }
AMREX_GPU_HOST_DEVICE inline Real w_of (int i, int j, int k) { return Real(0.3) + Real(0.05) * std::sin(Real(0.6) * i + Real(0.4) * j + Real(0.2) * k); }

struct Fields
{
    MultiFab cons, u, v, w, beta, bx, by, bz, sx, sy, sz;
};

Fields make_fields (const Mesh& m, bool with_face_blanking)
{
    Fields f{MultiFab(m.ba, m.dm, NVAR_max, 2),
             MultiFab(convert(m.ba, IntVect(1,0,0)), m.dm, 1, 2),
             MultiFab(convert(m.ba, IntVect(0,1,0)), m.dm, 1, 2),
             MultiFab(convert(m.ba, IntVect(0,0,1)), m.dm, 1, 2),
             MultiFab(m.ba, m.dm, 1, 2),
             MultiFab(convert(m.ba, IntVect(1,0,0)), m.dm, 1, 2),
             MultiFab(convert(m.ba, IntVect(0,1,0)), m.dm, 1, 2),
             MultiFab(convert(m.ba, IntVect(0,0,1)), m.dm, 1, 2),
             MultiFab(convert(m.ba, IntVect(1,0,0)), m.dm, 1, 0),
             MultiFab(convert(m.ba, IntVect(0,1,0)), m.dm, 1, 0),
             MultiFab(convert(m.ba, IntVect(0,0,1)), m.dm, 1, 0)};
    f.cons.setVal(Real(0.));
    f.cons.setVal(rho0, Rho_comp, 1, 2);
    f.cons.setVal(rho0 * Real(300.), RhoTheta_comp, 1, 2);
    for (MFIter mfi(f.u); mfi.isValid(); ++mfi) {
        auto const u = f.u.array(mfi);
        auto const v = f.v.array(mfi);
        auto const w = f.w.array(mfi);
        auto const b = f.beta.array(mfi);
        auto const bx = f.bx.array(mfi);
        auto const by = f.by.array(mfi);
        auto const bz = f.bz.array(mfi);
        ParallelFor(f.u[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            u(i,j,k) = u_of(i,j,k);
            bx(i,j,k) = building_face_beta(i,j,k,0);
        });
        ParallelFor(f.v[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            v(i,j,k) = v_of(i,j,k);
            by(i,j,k) = building_face_beta(i,j,k,1);
        });
        ParallelFor(f.w[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            w(i,j,k) = w_of(i,j,k);
            bz(i,j,k) = building_face_beta(i,j,k,2);
        });
        ParallelFor(f.beta[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            b(i,j,k) = building_beta(i, j, k);
        });
    }
    if (!with_face_blanking) { f.bx.clear(); f.by.clear(); f.bz.clear(); }
    f.sx.setVal(Real(0.)); f.sy.setVal(Real(0.)); f.sz.setVal(Real(0.));
    Gpu::streamSynchronize();
    return f;
}

void run_kernels (const Mesh& m, Fields& f, const SolverChoice& sc, Real dt)
{
    const bool fb = f.bx.ok();
    // the kernels limit the drag when they run on the slow step
    const bool slow = !sc.immersed_forcing_substep;
    for (MFIter mfi(f.sx); mfi.isValid(); ++mfi) {
        ImmersedForcingBuildings_Xmom(mfi.nodaltilebox(0), f.u.const_array(mfi), f.v.const_array(mfi), f.w.const_array(mfi),
                                      f.cons.const_array(mfi), f.beta.const_array(mfi),
                                      fb ? f.bx.const_array(mfi) : Array4<const Real>{},
                                      Array4<const Real>{}, f.sx.array(mfi), m.geom, sc, dt, slow);
        ImmersedForcingBuildings_Ymom(mfi.nodaltilebox(1), f.u.const_array(mfi), f.v.const_array(mfi), f.w.const_array(mfi),
                                      f.cons.const_array(mfi), f.beta.const_array(mfi),
                                      fb ? f.by.const_array(mfi) : Array4<const Real>{},
                                      Array4<const Real>{}, f.sy.array(mfi), m.geom, sc, dt, slow);
        ImmersedForcingBuildings_Zmom(mfi.nodaltilebox(2), f.u.const_array(mfi), f.v.const_array(mfi), f.w.const_array(mfi),
                                      f.cons.const_array(mfi), f.beta.const_array(mfi),
                                      fb ? f.bz.const_array(mfi) : Array4<const Real>{},
                                      Array4<const Real>{}, f.sz.array(mfi), m.geom, sc, dt, slow);
    }
    Gpu::streamSynchronize();
}

// One host copy of a MultiFab over region (valid cells, or with the source's ghost cells too)
std::vector<Real> to_host (const MultiFab& mf, const Box& region, bool with_ghosts = false)
{
    MultiFab gathered(BoxArray(region), DistributionMapping(BoxArray(region)), 1, 0);
    gathered.ParallelCopy(mf, 0, 0, 1, with_ghosts ? mf.nGrowVect() : IntVect(0), IntVect(0),
                          IntVect(0), Periodicity::NonPeriodic());
    std::vector<Real> out(static_cast<std::size_t>(region.numPts()), Real(-99.));
    for (MFIter mfi(gathered); mfi.isValid(); ++mfi) {
        Gpu::copy(Gpu::deviceToHost, gathered[mfi].dataPtr(), gathered[mfi].dataPtr() + gathered[mfi].size(), out.begin());
    }
    Gpu::streamSynchronize();
    return out;
}

Real at (const std::vector<Real>& a, const Box& region, int i, int j, int k)
{
    const IntVect lo = region.smallEnd(), len = region.length();
    return a[static_cast<std::size_t>(((k-lo[2])*len[1] + (j-lo[1]))*len[0] + (i-lo[0]))];
}

// A host Array4 over a host copy of region
Array4<const Real> host_view (const std::vector<Real>& a, const Box& region)
{
    const Dim3 lo = amrex::lbound(region);
    const Dim3 hi = amrex::ubound(region);
    return Array4<const Real>(a.data(), lo, Dim3{hi.x+1, hi.y+1, hi.z+1}, 1);
}

// Face blanking as the kernels read it without the snap and without face arrays (the mean of the
// two cells, thresholded): the independent reference of the first test
Real face_beta (int i, int j, int k, int dir)
{
    const int io = (dir == 0), jo = (dir == 1), ko = (dir == 2);
    const Real b = Real(0.5) * (building_beta(i,j,k) + building_beta(i-io,j-jo,k-ko));
    return (b < small) ? Real(0.) : b;
}

SolverChoice make_choice (bool implicit_projection, bool use_most)
{
    SolverChoice sc;
    sc.buildings_type         = BuildingsType::ImmersedForcing;
    sc.if_Cd_momentum         = Real(50.);
    sc.if_z0                  = Real(0.1);
    sc.if_use_most            = use_most;
    sc.if_implicit_drag       = true;
    sc.if_implicit_projection = implicit_projection;
    sc.if_small_volfrac       = small;
    return sc;
}

// The classes of faces of the building above, from the kernels' conditions (dir 0: lo/hi are the
// south/north neighbours). law: the wall law alone (a horizontal wall-law face carries no linear
// drag); law_and_drag: a z-face on the roof row that carries a side wall law and the blanked roof
// drag together; unforced: a partial z-face under MOST with neither (not a side, not under fluid),
// which the kernels leave alone on both paths
enum class Face { air, unforced, law, law_and_drag, drag_only, frozen };

Face classify (int i, int j, int k, int dir)
{
    const int io = (dir == 0), jo = (dir == 1), ko = (dir == 2);
    const Real raw = Real(0.5) * (building_beta(i,j,k) + building_beta(i-io,j-jo,k-ko));
    const Real b = face_beta(i, j, k, dir);
    if (b == Real(0.)) { return Face::air; }
    if (raw == Real(1.)) { return Face::frozen; }
    Real law = Real(0.);
    if (dir < 2) {
        const Real below = (k == 0) ? Real(0.) : face_beta(i,j,k-1,dir);
        const Real above = face_beta(i,j,k+1,dir);
        const Real lo = (dir == 0) ? face_beta(i,j-1,k,dir) : face_beta(i-1,j,k,dir);
        const Real hi = (dir == 0) ? face_beta(i,j+1,k,dir) : face_beta(i+1,j,k,dir);
        // the kernels' roof / lo / hi conditions without the snap
        law += (b < below && above == Real(0.)) ? Real(1.) : Real(0.);
        law += (b <= hi && lo == Real(0.)) ? Real(1.) : Real(0.);
        law += (b <= lo && hi == Real(0.)) ? Real(1.) : Real(0.);
    } else {
        if (k >= 1) {
            const Real s = face_beta(i,j-1,k,dir), n = face_beta(i,j+1,k,dir);
            const Real w = face_beta(i-1,j,k,dir), e = face_beta(i+1,j,k,dir);
            law += (b <= n && s == Real(0.)) ? Real(1.) : Real(0.);
            law += (b <= s && n == Real(0.)) ? Real(1.) : Real(0.);
            law += (b <= e && w == Real(0.)) ? Real(1.) : Real(0.);
            law += (b <= w && e == Real(0.)) ? Real(1.) : Real(0.);
        }
    }
    const bool roof_drag = (dir == 2) && (face_beta(i,j,k+1,dir) == Real(0.));
    if (law == Real(0.)) { return (dir < 2 || roof_drag) ? Face::drag_only : Face::unforced; }
    return roof_drag ? Face::law_and_drag : Face::law;
}

// The kernels' windspeed and normal velocity on a face, from the analytic wind
void face_wind (int i, int j, int k, int dir, Real& un, Real& windspeed)
{
    const Real fourth = Real(0.25);
    Real ux, uy, uz;
    if (dir == 0) {
        ux = u_of(i,j,k);
        uy = fourth * (v_of(i,j,k) + v_of(i-1,j,k) + v_of(i,j+1,k) + v_of(i-1,j+1,k));
        uz = fourth * (w_of(i,j,k) + w_of(i-1,j,k) + w_of(i,j,k+1) + w_of(i-1,j,k+1));
        un = ux;
    } else if (dir == 1) {
        ux = fourth * (u_of(i,j,k) + u_of(i,j-1,k) + u_of(i+1,j,k) + u_of(i+1,j-1,k));
        uy = v_of(i,j,k);
        uz = fourth * (w_of(i,j,k) + w_of(i,j-1,k) + w_of(i,j,k+1) + w_of(i,j-1,k+1));
        un = uy;
    } else {
        ux = fourth * (u_of(i,j,k) + u_of(i+1,j,k) + u_of(i,j,k-1) + u_of(i+1,j,k-1));
        uy = fourth * (v_of(i,j,k) + v_of(i,j+1,k) + v_of(i,j,k-1) + v_of(i,j+1,k-1));
        uz = w_of(i,j,k);
        un = uz;
    }
    windspeed = std::sqrt(ux*ux + uy*uy + uz*uz);
}

} // namespace

TEST(ImmersedBuildingImplicit, KernelsKeepTheWallLawOnlyUnderTheImplicitProjection)
{
    const Mesh m = make_mesh();
    const Real dt = Real(0.5);
    Fields fe = make_fields(m, false);
    Fields fi = make_fields(m, false);
    run_kernels(m, fe, make_choice(false, true), dt);
    run_kernels(m, fi, make_choice(true,  true), dt);

    int n_law = 0, n_law_and_drag = 0, n_drag = 0, n_frozen = 0;
    for (int dir = 0; dir < 3; ++dir) {
        const Box r = convert(m.domain, IntVect::TheDimensionVector(dir));
        const MultiFab& se = (dir == 0) ? fe.sx : (dir == 1) ? fe.sy : fe.sz;
        const MultiFab& si = (dir == 0) ? fi.sx : (dir == 1) ? fi.sy : fi.sz;
        const std::vector<Real> he = to_host(se, r);
        const std::vector<Real> hi = to_host(si, r);
        for (int k = r.smallEnd(2); k <= r.bigEnd(2); ++k) {
        for (int j = r.smallEnd(1); j <= r.bigEnd(1); ++j) {
        for (int i = r.smallEnd(0); i <= r.bigEnd(0); ++i) {
            SCOPED_TRACE(testing::Message() << "dir " << dir << " face " << i << " " << j << " " << k);
            const Real e = at(he, r, i, j, k), im = at(hi, r, i, j, k);
            switch (classify(i, j, k, dir)) {
            case Face::air:
            case Face::unforced:
                EXPECT_EQ(e, Real(0.));
                EXPECT_EQ(im, Real(0.));
                break;
            case Face::law:
                // the wall law of the explicit path, and nothing else: a horizontal wall-law face
                // carries no partial-wall drag (side = 0 there), so the two agree exactly
                EXPECT_NE(e, Real(0.));
                EXPECT_EQ(im, e);
                ++n_law;
                break;
            case Face::law_and_drag:
                // the wall law stays, the roof drag moves to the projection
                EXPECT_NE(e, Real(0.));
                EXPECT_NE(im, Real(0.));
                EXPECT_NE(im, e);
                ++n_law_and_drag;
                break;
            case Face::drag_only:
                EXPECT_NE(e, Real(0.));
                EXPECT_EQ(im, Real(0.));
                ++n_drag;
                break;
            case Face::frozen:
                EXPECT_NE(e, Real(0.));
                EXPECT_EQ(im, Real(0.));
                ++n_frozen;
                break;
            }
        }}}
    }
    // the building has faces of every class
    EXPECT_GT(n_law, 0);
    EXPECT_GT(n_law_and_drag, 0);
    EXPECT_GT(n_drag, 0);
    EXPECT_GT(n_frozen, 0);
}

// Without the wall law every solid face is a drag face, and the kernels then add nothing
TEST(ImmersedBuildingImplicit, KernelsAddNothingWithoutTheWallLaw)
{
    const Mesh m = make_mesh();
    Fields fe = make_fields(m, false);
    Fields fi = make_fields(m, false);
    run_kernels(m, fe, make_choice(false, false), Real(0.5));
    run_kernels(m, fi, make_choice(true,  false), Real(0.5));
    EXPECT_GT(fe.sx.norm0(), Real(0.));
    EXPECT_GT(fe.sy.norm0(), Real(0.));
    EXPECT_GT(fe.sz.norm0(), Real(0.));
    EXPECT_EQ(fi.sx.norm0(), Real(0.));
    EXPECT_EQ(fi.sy.norm0(), Real(0.));
    EXPECT_EQ(fi.sz.norm0(), Real(0.));
}

// The drag the kernels drop under the implicit projection is, on every face, the linear weight
// the projection reads through the shared header times rho CdM |u| u_n (the kernels' explicit
// drag with erf.if_implicit_drag = false and the forcing on the substep, so no limiter and no
// point-implicit rescale), and a face with raw blanking 1 loses the wall law too. Four readings
// of the blanking: with and without the snap, with and without face-centred arrays.
TEST(ImmersedBuildingImplicit, DroppedDragIsTheProjectionWeightOnEveryFace)
{
    const Mesh m = make_mesh();
    const Real dt = Real(0.5);
    const Real drag_coefficient = Real(50.) / std::cbrt(dx*dx*dz);
    const Real tol = Real(64.) * std::numeric_limits<Real>::epsilon();
    for (bool snap : {false, true}) {
    for (bool with_faces : {false, true}) {
        SCOPED_TRACE(testing::Message() << "snap " << snap << " face arrays " << with_faces);
        Fields fe = make_fields(m, with_faces);
        Fields fi = make_fields(m, with_faces);
        SolverChoice sce = make_choice(false, true);
        SolverChoice sci = make_choice(true,  true);
        for (SolverChoice* sc : {&sce, &sci}) {
            sc->if_implicit_drag          = false;
            sc->immersed_forcing_substep  = true;   // the plain branch: src -= drag
            sc->if_snap_partial_cells     = snap;
        }
        run_kernels(m, fe, sce, dt);
        run_kernels(m, fi, sci, dt);

        // the blanking arrays the kernels read, on the host, ghost cells included
        const Box gb = amrex::grow(m.domain, 2);
        const std::vector<Real> hb = to_host(fe.beta, gb, true);
        const Array4<const Real> t = host_view(hb, gb);
        std::vector<Real> hbf[3];
        Box gf[3];
        for (int dir = 0; dir < 3; ++dir) {
            gf[dir] = convert(gb, IntVect::TheDimensionVector(dir));
            if (with_faces) {
                hbf[dir] = to_host((dir == 0) ? fe.bx : (dir == 1) ? fe.by : fe.bz, gf[dir], true);
            }
        }

        int n_weighted = 0, n_raw_solid = 0, n_raw_solid_with_law = 0;
        for (int dir = 0; dir < 3; ++dir) {
            const Array4<const Real> tf = with_faces ? host_view(hbf[dir], gf[dir]) : Array4<const Real>{};
            const Box r = convert(m.domain, IntVect::TheDimensionVector(dir));
            const MultiFab& se = (dir == 0) ? fe.sx : (dir == 1) ? fe.sy : fe.sz;
            const MultiFab& si = (dir == 0) ? fi.sx : (dir == 1) ? fi.sy : fi.sz;
            const std::vector<Real> he = to_host(se, r);
            const std::vector<Real> hi = to_host(si, r);
            for (int k = r.smallEnd(2); k <= r.bigEnd(2); ++k) {
            for (int j = r.smallEnd(1); j <= r.bigEnd(1); ++j) {
            for (int i = r.smallEnd(0); i <= r.bigEnd(0); ++i) {
                SCOPED_TRACE(testing::Message() << "dir " << dir << " face " << i << " " << j << " " << k);
                Real linear, raw, b;
                if (dir == 0) {
                    const if_bld::HorizFace f = if_bld::horiz_face<0>(t, tf, i, j, k, 0, snap, small);
                    linear = if_bld::horiz_masks(f, true, snap).linear(f.t); raw = f.raw; b = f.t;
                } else if (dir == 1) {
                    const if_bld::HorizFace f = if_bld::horiz_face<1>(t, tf, i, j, k, 0, snap, small);
                    linear = if_bld::horiz_masks(f, true, snap).linear(f.t); raw = f.raw; b = f.t;
                } else {
                    const if_bld::VertFace f = if_bld::vert_face(t, tf, i, j, k, snap, small);
                    linear = if_bld::vert_masks(f, k >= 1, true, snap).linear(f.t); raw = f.raw; b = f.t;
                }
                Real un, windspeed;
                face_wind(i, j, k, dir, un, windspeed);
                const Real CdM = std::min(drag_coefficient / (windspeed + std::numeric_limits<Real>::epsilon()), drag_coefficient);
                const Real dropped = linear * rho0 * CdM * un * windspeed;
                const Real e = at(he, r, i, j, k), im = at(hi, r, i, j, k);
                if (linear > Real(0.)) { ++n_weighted; }
                if (raw == Real(1.)) {
                    // the projection holds a raw-solid face by the full rate; the kernels add
                    // nothing there, neither the drag nor a wall law the explicit path carried
                    EXPECT_EQ(if_bld::implicit_weight(linear, raw), Real(1.));
                    EXPECT_EQ(im, Real(0.));
                    ++n_raw_solid;
                    // the exact south wall: solid x-faces whose explicit source holds the wall law
                    // on top of the drag (e = -(law + dropped))
                    if (e + dropped != Real(0.)) { ++n_raw_solid_with_law; }
                } else {
                    // explicit = implicit - dropped drag (the kernels subtract the drag), the wall
                    // law being the same on both paths
                    EXPECT_EQ(if_bld::implicit_weight(linear, raw), linear);
                    EXPECT_NEAR(e - im, -dropped, tol * (std::abs(e) + std::abs(dropped) + Real(1.)));
                }
                // the face blanking the kernels used is the header's: an open face gets no source
                if (b == Real(0.)) { EXPECT_EQ(e, Real(0.)); }
            }}}
        }
        EXPECT_GT(n_weighted, 0);
        EXPECT_GT(n_raw_solid, 0);
        EXPECT_GT(n_raw_solid_with_law, 0);
    }}
}
