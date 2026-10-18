#include <cmath>
#include <limits>
#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <gtest/gtest.h>

#include "ERF_TerrainFileChecks.H"

// Motivation: a terrain or buildings height file that is malformed, mis-ordered or too small
// for the domain used to be read without a word -- a short x list held the edge height over
// the rest of the domain, lists out of order bracketed the wrong points -- and a file written
// with x varying fastest put a north-south ridge east-west.  TableProblem names the first such
// problem, TableSummary prints where the highest point and the steepest slopes sit, and
// FaceBandSlopes / OpenFaceSlopeWarning report steep terrain next to inflow and outflow faces.

using namespace amrex;
using namespace erf_terrain_file;

namespace {

const Real tol = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-5);

GpuArray<Real,AMREX_SPACEDIM> lo_of (Real x, Real y) { return {x, y, Real(0.)}; }
GpuArray<Real,AMREX_SPACEDIM> hi_of (Real x, Real y) { return {x, y, Real(1000.)}; }

// An ERF-order table over [0, 800] x [0, 400] with nodes every 50 m and h = a x + b y
struct Table
{
    int nx = 17, ny = 9;
    Vector<Real> x, y, z;
    Table (Real a, Real b, bool usgs = false)
    {
        for (int i = 0; i < nx; ++i) { x.push_back(Real(50.) * Real(i)); }
        for (int j = 0; j < ny; ++j) { y.push_back(Real(50.) * Real(j)); }
        z.resize(nx*ny);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                const Real h = a * x[i] + b * y[j];
                if (usgs) { z[j*nx + i] = h; } else { z[i*ny + j] = h; }
            }
        }
        if (usgs) {
            Vector<Real> xu;
            for (int j = 0; j < ny; ++j) { for (int i = 0; i < nx; ++i) { xu.push_back(x[i]); } }
            x = xu;
        }
    }
};

std::string
problem (const Table& t, bool usgs = false, Real xlo = 0, Real xhi = 800, Real ylo = 0, Real yhi = 400)
{
    return TableProblem(usgs, t.nx, t.ny, t.x, t.y, t.z, lo_of(xlo, ylo), hi_of(xhi, yhi));
}

// Nodal heights on 50 m cells: a ramp of slope 0.6 in x from x = 400 m to 600 m, flat before
// and after it, and 100 m layers in k (at namespace scope: nvcc rejects a device lambda in a
// TEST body)
void
fill_ramp (MultiFab& z)
{
    for (MFIter mfi(z); mfi.isValid(); ++mfi) {
        auto a = z.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            const Real x = Real(50.) * Real(i);
            const Real h = Real(0.6) * amrex::min(amrex::max(x - Real(400.), Real(0.)), Real(200.));
            a(i,j,k) = h + Real(100.) * Real(k);
        });
    }
}

} // namespace

TEST(TerrainFileChecks, AcceptsANodalAndAPerCellTable)
{
    Table t(Real(0.1), Real(0.));
    EXPECT_EQ(problem(t), "");
    // One height per cell starting at prob_lo (x = 0 .. 750 over a 0 .. 800 domain): each end
    // may stop one table spacing inside the domain, the layout of the buildings files
    EXPECT_EQ(problem(t, false, 0, 850, 0, 450), "");
    // A table larger than the domain is fine
    EXPECT_EQ(problem(t, false, 100, 700, 100, 300), "");
}

TEST(TerrainFileChecks, CoverageLimitIsOneTableSpacing)
{
    Table t(Real(0.), Real(0.));
    // 50 m spacing: a 50 m gap at the top passes, 51 m does not; the same at the bottom
    EXPECT_EQ(problem(t, false, 0, 850), "");
    EXPECT_NE(problem(t, false, 0, 851).find("do not cover the domain"), std::string::npos);
    EXPECT_EQ(problem(t, false, -50, 800), "");
    EXPECT_NE(problem(t, false, -51, 800).find("x values 0 .. 800 do not cover"), std::string::npos);
    EXPECT_NE(problem(t, false, 0, 800, 0, 451).find("y values 0 .. 400 do not cover"), std::string::npos);
}

TEST(TerrainFileChecks, AcceptsPerCellTablesWithRoundedCoordinates)
{
    // 128 cells of 46.875 m over 0 .. 6000 m, one height per cell, coordinates written with
    // six significant digits (Python's %g): the last x is 5953.12, so the gap at the top
    // (46.88 m) exceeds the last spacing (46.87 m) by rounding alone
    Table t(Real(0.), Real(0.));
    t.nx = 128; t.ny = 2;
    t.x.clear();
    for (int i = 0; i < t.nx; ++i) {
        t.x.push_back(std::round(Real(46.875) * Real(i) * Real(100.)) / Real(100.));
    }
    t.x[t.nx-1] = Real(5953.12);
    t.x[t.nx-2] = Real(5906.25);
    t.y = {Real(0.), Real(10.)};
    t.z.assign(t.nx*t.ny, Real(0.));
    EXPECT_EQ(problem(t, false, 0, 6000, 0, 20), "");
    // ... while stopping two cells short still fails
    EXPECT_NE(problem(t, false, 0, 6046.875, 0, 20).find("do not cover the domain"), std::string::npos);

    // Map coordinates: x from 500000 m, 48 cells of 13.3333 m; one cell short at the top
    Table u(Real(0.), Real(0.));
    u.nx = 48; u.ny = 2;
    u.x.clear();
    for (int i = 0; i < u.nx; ++i) { u.x.push_back(Real(500000.) + Real(640.) * Real(i) / Real(48.)); }
    u.y = {Real(0.), Real(10.)};
    u.z.assign(u.nx*u.ny, Real(0.));
    EXPECT_EQ(problem(u, false, 500000, 500640, 0, 20), "");
}

TEST(TerrainFileChecks, RejectsSwappedLists)
{
    // nx/ny and the x/y lists exchanged: on a domain that is not square the x list falls short
    Table t(Real(0.), Real(0.));
    std::swap(t.nx, t.ny);
    std::swap(t.x, t.y);
    EXPECT_NE(problem(t).find("x values 0 .. 400 do not cover the domain 0 .. 800"), std::string::npos);
}

TEST(TerrainFileChecks, RejectsListsThatDoNotIncrease)
{
    Table t(Real(0.), Real(0.));
    std::swap(t.x[5], t.x[6]);
    EXPECT_NE(problem(t).find("x values must increase strictly, but x[6] = 250 follows x[5] = 300"),
              std::string::npos);
    Table u(Real(0.), Real(0.));
    u.y[3] = u.y[2];   // a repeated value is not an increase
    EXPECT_NE(problem(u).find("y values must increase strictly, but y[3]"), std::string::npos);
}

TEST(TerrainFileChecks, RejectsSmallTablesWrongSizesAndNonFiniteValues)
{
    Table t(Real(0.), Real(0.));
    // One row (ny = 1) is a profile constant in y in the ERF format, whatever the y extent;
    // nx = 0 is not a table; the USGS layout needs two points in each direction
    Table one_row = t; one_row.ny = 1; one_row.y.resize(1); one_row.z.resize(t.nx);
    EXPECT_EQ(problem(one_row), "");
    Table none = one_row; none.nx = 0; none.x.clear(); none.z.clear();
    EXPECT_NE(problem(none).find("must both be at least 1"), std::string::npos);
    Table u1(Real(0.), Real(0.), true); u1.ny = 1; u1.y.resize(1); u1.x.resize(u1.nx); u1.z.resize(u1.nx);
    EXPECT_NE(problem(u1, true).find("must both be at least 2 in the USGS format"), std::string::npos);
    Table shrt = t; shrt.z.pop_back();
    EXPECT_NE(problem(shrt).find("expected 17 x values, 9 y values and 153 heights"), std::string::npos);
    Table bad = t; bad.z[40] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_NE(problem(bad).find("height value 40 is not finite"), std::string::npos);
    Table inf = t; inf.x[3] = std::numeric_limits<Real>::infinity();
    EXPECT_NE(problem(inf).find("x value 3 is not finite"), std::string::npos);
}

TEST(TerrainFileChecks, UsgsLayoutChecksEveryRow)
{
    Table t(Real(0.), Real(0.), true);
    EXPECT_EQ(problem(t, true), "");
    // An ERF-order (y fastest) file read as USGS: the x values of a row do not increase
    Table e(Real(0.), Real(0.), true);
    for (int j = 0; j < e.ny; ++j) {
        for (int i = 0; i < e.nx; ++i) { e.x[j*e.nx + i] = Real(50.) * Real((j*e.nx + i) / e.ny); }
    }
    EXPECT_NE(problem(e, true).find("(row 0 of the x-fastest USGS layout)"), std::string::npos);
    // One short row is enough
    Table s(Real(0.), Real(0.), true);
    for (int i = 0; i < s.nx; ++i) { s.x[4*s.nx + i] *= Real(0.9); }
    EXPECT_NE(problem(s, true).find("(row 4 of the x-fastest USGS layout)"), std::string::npos);
}

TEST(TerrainFileChecks, SummaryGivesOrderExtremesAndSlopes)
{
    // h = 0.2 x + 0.05 y: lowest at (0,0), highest at (800,400), slopes 0.2 and 0.05 everywhere
    Table t(Real(0.2), Real(0.05));
    const std::string s = TableSummary("hill.txt", false, t.nx, t.ny, t.x, t.y, t.z);
    EXPECT_NE(s.find("Height file hill.txt: 17 x 9 heights over x 0 .. 800, y 0 .. 400"), std::string::npos) << s;
    EXPECT_NE(s.find("y varies fastest, height n is at (x[n / ny], y[n % ny])"), std::string::npos) << s;
    EXPECT_NE(s.find("lowest  0 m at (x, y) = (0, 0)"), std::string::npos) << s;
    EXPECT_NE(s.find("highest 180 m at (x, y) = (800, 400)"), std::string::npos) << s;
    EXPECT_NE(s.find("largest |dh/dx| 0.2 at"), std::string::npos) << s;
    EXPECT_NE(s.find("largest |dh/dy| 0.05 at"), std::string::npos) << s;

    // The same surface in the USGS layout gives the same extremes and slopes
    Table u(Real(0.2), Real(0.05), true);
    const std::string su = TableSummary("hill.txt", true, u.nx, u.ny, u.x, u.y, u.z);
    EXPECT_NE(su.find("x varies fastest, height n is at (x[n], y[n / nx])"), std::string::npos) << su;
    EXPECT_NE(su.find("highest 180 m at (x, y) = (800, 400)"), std::string::npos) << su;
    EXPECT_NE(su.find("largest |dh/dx| 0.2 at"), std::string::npos) << su;

    // A ridge along y read with x fastest by mistake shows up as slope in y, not in x
    Table r(Real(0.), Real(0.));
    for (int i = 0; i < r.nx; ++i) {
        for (int j = 0; j < r.ny; ++j) {
            r.z[j*r.nx + i] = (i == 8) ? Real(100.) : Real(0.);   // written x fastest
        }
    }
    const std::string sr = TableSummary("ridge.txt", false, r.nx, r.ny, r.x, r.y, r.z);
    EXPECT_EQ(sr.find("largest |dh/dy| 0 at"), std::string::npos) << sr;

    // A single row: no y slope to report
    Vector<Real> z1;
    for (int i = 0; i < t.nx; ++i) { z1.push_back(Real(0.2) * t.x[i]); }
    const std::string s1 = TableSummary("row.txt", false, t.nx, 1, t.x, {Real(0.)}, z1);
    EXPECT_NE(s1.find("largest |dh/dx| 0.2 at"), std::string::npos) << s1;
    EXPECT_NE(s1.find("one y value: constant in y"), std::string::npos) << s1;
}

TEST(TerrainFileChecks, CellTerrainSlopeOfAPlaneWithMapFactors)
{
    // A plane h = a x + b y in physical coordinates on a cell of grid spacing dx, dy with map
    // factors mfx, mfy: the physical spacing is dx/mfx, so the drop across the cell is
    // a dx/mfx and the slope is sqrt(a^2 + b^2) whatever the map factors
    const Box nodes(IntVect(0,0,0), IntVect(1,1,1), IntVect(1,1,1));
    FArrayBox fab(nodes, 1, The_Pinned_Arena());   // written and read on the host
    const Real a = Real(0.3), b = Real(-0.4), dx = Real(20.), dy = Real(10.);
    const Real mfx = Real(1.3), mfy = Real(0.8);
    auto z = fab.array();
    for (int k = 0; k <= 1; ++k) {
        for (int j = 0; j <= 1; ++j) {
            for (int i = 0; i <= 1; ++i) {
                z(i,j,k) = a * (dx/mfx) * Real(i) + b * (dy/mfy) * Real(j) + Real(50.) * Real(k);
            }
        }
    }
    EXPECT_NEAR(CellTerrainSlope(0, 0, 0, fab.const_array(), mfx, mfy, one/dx, one/dy), Real(0.5), tol);
    // Ignoring the map factors would report sqrt((a/mfx)^2 + (b/mfy)^2) = 0.5507 instead
    EXPECT_GT(std::abs(CellTerrainSlope(0, 0, 0, fab.const_array(), one, one, one/dx, one/dy) - Real(0.5)),
              Real(0.03));
}

TEST(TerrainFileChecks, FaceBandSlopesSeeOnlyTheirBands)
{
    // 24 x 12 cells of 50 m, split into 6 x 6 boxes.  A ramp of slope 0.6 in x between
    // x = 400 m and 600 m, flat before and after it: the x faces' bands (3 cells) are flat,
    // the y faces' bands cross the ramp.
    const Box domain(IntVect(0,0,0), IntVect(23,11,3));
    const RealBox rb({Real(0.), Real(0.), Real(0.)}, {Real(1200.), Real(600.), Real(400.)});
    const Geometry geom(domain, rb, CoordSys::cartesian, {0, 0, 0});
    BoxArray ba(domain);
    ba.maxSize(6);
    const DistributionMapping dm(ba);
    MultiFab z(convert(ba, IntVect(1)), dm, 1, 1);
    fill_ramp(z);
    // map factors on the same boxes flattened to one k, as ERF's ba2d
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) { b.setRange(2, 0); }
    const BoxArray ba2d(std::move(bl2d));
    MultiFab mf1(ba2d, dm, 1, 0), mf2(ba2d, dm, 1, 0);
    mf1.setVal(one);
    mf2.setVal(Real(2.));
    const auto s3 = FaceBandSlopes(z, mf1, mf1, geom, 3);
    EXPECT_NEAR(s3[0], Real(0.), tol);    // xlo
    EXPECT_NEAR(s3[1], Real(0.), tol);    // xhi
    EXPECT_NEAR(s3[2], Real(0.6), tol);   // ylo
    EXPECT_NEAR(s3[3], Real(0.6), tol);   // yhi

    // A band of 8 cells (400 m) stops just short of the ramp; 9 cells (450 m) reach its
    // first cell (400 .. 450 m) from xlo
    EXPECT_NEAR(FaceBandSlopes(z, mf1, mf1, geom, 8)[0], Real(0.), tol);
    const auto s9 = FaceBandSlopes(z, mf1, mf1, geom, 9);
    EXPECT_NEAR(s9[0], Real(0.6), tol);
    EXPECT_NEAR(s9[1], Real(0.), tol);    // xhi band 750 .. 1200 m stays flat

    // A map factor of 2 halves the physical cell: the same drops are twice as steep
    EXPECT_NEAR(FaceBandSlopes(z, mf2, mf2, geom, 3)[2], Real(1.2), tol);
}

TEST(TerrainFileChecks, OpenFaceWarningOnlyForSteepOpenFaces)
{
    const Array<std::string,4> names{"Inflow", "Outflow", "SlipWall", "SlipWall"};
    const Array<bool,4> open{true, true, false, false};
    // Below or at the threshold: nothing
    EXPECT_EQ(OpenFaceSlopeWarning(0, 3, {Real(0.2), open_face_slope_warn, Real(-1.), Real(-1.)}, open, names), "");
    // Steep, but only next to walls: nothing
    EXPECT_EQ(OpenFaceSlopeWarning(0, 3, {Real(0.1), Real(0.1), Real(0.9), Real(0.9)}, open, names), "");
    // Steep next to the outflow face
    const std::string w = OpenFaceSlopeWarning(1, 3, {Real(0.1), Real(0.45), Real(0.), Real(0.)}, open, names);
    EXPECT_NE(w.find("WARNING: at level 1 the terrain slope |grad h| reaches 0.45 within 3 cells of the xhi face (Outflow)"),
              std::string::npos) << w;
    EXPECT_EQ(w.find("xlo face"), std::string::npos) << w;
    EXPECT_NE(w.find("Flattening the terrain within 3 cells"), std::string::npos) << w;
    EXPECT_NE(w.find("none of them is on here"), std::string::npos) << w;
}
