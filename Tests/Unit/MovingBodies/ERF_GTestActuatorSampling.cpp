// Contract of the actuator sampler (Core/ERF_ActuatorSampling).
//
// - sample_velocity: each velocity component is interpolated from its own staggered grid,
//   linearly in the physical height, so a field linear in x, y and physical z is reproduced
//   exactly on a uniform-dz mesh, on a stretched mesh (z_phys_nd with equal columns) and on a
//   terrain-following mesh (z_phys_nd varying with x and y), at points inside the cells, on faces,
//   near the ground and near the top; every point is found by exactly one box, however the
//   domain is split, in z too; a point outside the domain aborts, naming it.
//   A point on the domain's upper x or y face belongs to the last cell; in a periodic direction a
//   point beyond the domain is sampled where its image lies inside (sample_velocity, sample_cell_scalar
//   and terrain_heights alike).
// - points_covered_by: the reach counts, and wraps across a periodic seam; on a terrain or stretched mesh
//   (CoverZ::Footprint) the boxes' footprints must hold the reach across, at any height, and for the
//   spreading (CoverZ::Column) the grids must hold it in every cell the heights within reach may lie in
//   (from mesh_z_bounds, exact on a flat uniform mesh), or over the whole height without the bounds.
//   Whether a point is covered does not depend on how the level is cut into boxes.
// - column_cell_ranges: the cells the heights within reach lie in over a point's own columns, never fewer than
//   the cells themselves; points_covered_by with them accepts a point the mesh-wide bounds refuse beside a
//   canyon. Coverage on footprints built once (CoverageGrids) is coverage on a BoxArray.
// - points_covered_by refuses bounds without a positive spacing and one rank's height table.
// - sample_velocity_local with kept buffers gives the answers of a call with its own, one kernel per box.
// - mesh_z_bounds(.., false): the table's cells for a height range hold every cell that reaches into it.
// - wrap_periodic: a point a few ulps below a non-zero prob_lo lands on prob_lo, not on prob_hi.
// - terrain_heights: the bilinear k = 0 node surface.
// - region_touches_grids: a region's cells on a level, through periodic images, far outside the domain, over
//   terrain, and on a level with boxes that do not reach the ground (erf-model/ERF#5056).
// - a point held by two boxes of one rank is counted twice and refused by sample_velocity.

#include <limits>
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <memory>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include "ERF_ActuatorSampling.H"
#include "ERF_GTestThrowOnAbort.H"

namespace {

using amrex::Real;

// The linear field: (u,v,w) = c0 + cx x + cy y + cz z_phys, different coefficients per component
struct LinearField {
    std::array<Real,4> u{{ 3.0,  0.010, -0.020, 0.030}};
    std::array<Real,4> v{{-1.0, -0.005,  0.015, 0.010}};
    std::array<Real,4> w{{ 0.5,  0.002,  0.003, -0.004}};
    static Real eval (const std::array<Real,4>& c, Real x, Real y, Real z) { return c[0] + c[1]*x + c[2]*y + c[3]*z; }
};

// basic terrain-following (BTF) heights: z(x,y,eta) = h(x,y) + (H - h(x,y)) * eta, with eta the
// stretched nominal level fraction; h = 0 gives a flat stretched mesh
struct MeshSpec {
    int nx = 12, ny = 10, nz = 8;
    Real Lx = 1200.0, Ly = 1000.0, H = 400.0;
    Real stretch = 1.15;         // 1.0: uniform nominal levels
    Real hill = 0.0;             // terrain height amplitude (m); 0: flat
    std::array<int,3> max_grid{{1024, 1024, 1024}};

    Real dx () const { return Lx / nx; }
    Real dy () const { return Ly / ny; }
    Real h (Real x, Real y) const
    {
        constexpr Real pi = 3.14159265358979323846;
        return hill * (0.5 + 0.5 * std::cos(2.0 * pi * x / Lx)) * (0.5 + 0.5 * std::sin(2.0 * pi * y / Ly));
    }
    Real eta (int k) const
    {
        if (stretch == 1.0) { return Real(k) / nz; }
        return (std::pow(stretch, k) - 1.0) / (std::pow(stretch, nz) - 1.0);
    }
    Real z_node (int i, int j, int k) const
    {
        const Real x = i * dx(), y = j * dy();
        const Real hs = h(x, y);
        return hs + (H - hs) * eta(k);
    }
};

struct Fields {
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab u, v, w;
    std::unique_ptr<amrex::MultiFab> znd;
    MeshSpec spec;
    LinearField lin;

    // build the mesh and fill the face velocities with the linear field at the face centres;
    // uniform_dz leaves z_phys_nd null so the sampler takes the nominal heights
    Fields (const MeshSpec& m, bool uniform_dz)
        : spec(m)
    {
        const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx-1, m.ny-1, m.nz-1));
        const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(m.Lx, m.Ly, m.H)});
        const std::array<int,3> periodic{{0, 0, 0}};
        geom = amrex::Geometry(domain, &rb, 0, periodic.data());
        ba = amrex::BoxArray(domain);
        ba.maxSize(amrex::IntVect(m.max_grid[0], m.max_grid[1], m.max_grid[2]));
        dm = amrex::DistributionMapping(ba);
        const int ng = 2;
        u.define(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, ng);
        v.define(amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, ng);
        w.define(amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, ng);
        if (!uniform_dz) {
            znd = std::make_unique<amrex::MultiFab>(amrex::convert(ba, amrex::IntVect(1,1,1)), dm, 1, ng);
        }
        const Real dz = m.H / m.nz;
        auto node_z = [&](int i, int j, int k) { return uniform_dz ? k * dz : m.z_node(i, j, k); };
        for (amrex::MFIter mfi(u, false); mfi.isValid(); ++mfi) {
            // fill valid and ghost cells alike: the sampler reads one ghost row across box faces
            auto ua = u.array(mfi); auto va = v.array(mfi); auto wa = w.array(mfi);
            amrex::Array4<Real> za;
            if (znd) { za = znd->array(mfi); }
            const amrex::Box gbx = amrex::grow(mfi.validbox(), ng);   // x-face box grown
            const amrex::Box cbx = amrex::enclosedCells(gbx);
            amrex::LoopOnCpu(amrex::surroundingNodes(cbx), [&](int i, int j, int k) {
                if (za) { za(i,j,k) = node_z(i, j, k); }
            });
            amrex::LoopOnCpu(amrex::convert(cbx, amrex::IntVect(1,0,0)), [&](int i, int j, int k) {
                const Real x = i * m.dx(), y = (j + 0.5) * m.dy();
                const Real z = 0.25 * (node_z(i,j,k) + node_z(i,j+1,k) + node_z(i,j,k+1) + node_z(i,j+1,k+1));
                ua(i,j,k) = LinearField::eval(lin.u, x, y, z);
            });
            amrex::LoopOnCpu(amrex::convert(cbx, amrex::IntVect(0,1,0)), [&](int i, int j, int k) {
                const Real x = (i + 0.5) * m.dx(), y = j * m.dy();
                const Real z = 0.25 * (node_z(i,j,k) + node_z(i+1,j,k) + node_z(i,j,k+1) + node_z(i+1,j,k+1));
                va(i,j,k) = LinearField::eval(lin.v, x, y, z);
            });
            amrex::LoopOnCpu(amrex::convert(cbx, amrex::IntVect(0,0,1)), [&](int i, int j, int k) {
                const Real x = (i + 0.5) * m.dx(), y = (j + 0.5) * m.dy();
                const Real z = 0.25 * (node_z(i,j,k) + node_z(i+1,j,k) + node_z(i,j+1,k) + node_z(i+1,j+1,k));
                wa(i,j,k) = LinearField::eval(lin.w, x, y, z);
            });
        }
    }

    // sample the points and check every component against the analytic field
    void check (const std::vector<Real>& pos, Real tol) const
    {
        std::vector<Real> vel;
        erf_actuator::sample_velocity(u, v, w, znd.get(), geom, pos, vel);
        ASSERT_EQ(vel.size(), pos.size());
        for (std::size_t p = 0; p < pos.size() / 3; ++p) {
            SCOPED_TRACE("point " + std::to_string(p) + " at (" + std::to_string(pos[3*p]) + ", " +
                         std::to_string(pos[3*p+1]) + ", " + std::to_string(pos[3*p+2]) + ")");
            const Real x = pos[3*p], y = pos[3*p+1], z = pos[3*p+2];
            const Real eu = LinearField::eval(lin.u, x, y, z);
            const Real ev = LinearField::eval(lin.v, x, y, z);
            const Real ew = LinearField::eval(lin.w, x, y, z);
            EXPECT_NEAR(vel[3*p],   eu, tol * std::max(Real(1.0), std::abs(eu)));
            EXPECT_NEAR(vel[3*p+1], ev, tol * std::max(Real(1.0), std::abs(ev)));
            EXPECT_NEAR(vel[3*p+2], ew, tol * std::max(Real(1.0), std::abs(ew)));
        }
    }
};

// Points inside cells, on cell faces, on the terrain-following mesh's lowest and highest
// centres, and at the hub-height of a rotor: all inside the domain.
std::vector<Real> test_points (const MeshSpec& m)
{
    std::vector<Real> pos;
    auto add = [&](Real x, Real y, Real z) { pos.push_back(x); pos.push_back(y); pos.push_back(z); };
    add(0.37 * m.Lx, 0.61 * m.Ly, 0.55 * m.H);   // interior, generic
    add(0.50 * m.Lx, 0.50 * m.Ly, 0.50 * m.H);   // on x and y cell faces
    add(3.0 * m.dx(), 2.0 * m.dy(), 0.25 * m.H); // on a cell corner column
    add(0.12 * m.Lx, 0.90 * m.Ly, 0.02 * m.H + m.h(0.12 * m.Lx, 0.90 * m.Ly) * 0.98); // just above the ground
    add(0.81 * m.Lx, 0.23 * m.Ly, 0.97 * m.H);   // near the top
    add(0.66 * m.Lx, 0.44 * m.Ly, 0.375 * m.H);  // rotor hub height
    add(0.66 * m.Lx, 0.44 * m.Ly + 0.3 * m.H, 0.375 * m.H); // rotor tip in y
    return pos;
}

// relative tolerance: each interpolation weight is exact for a linear field, so only the
// roundoff of a few dozen operations on the face values remains: far below 1e-9 in double,
// and a few units of 1e-7 in single precision, where 2e-4 leaves a wide margin
constexpr Real tol = (sizeof(Real) == 8) ? Real(1.0e-9) : Real(2.0e-4);

} // namespace

TEST(ActuatorSampling, ExactForLinearFieldOnUniformMesh)
{
    MeshSpec m;
    Fields f(m, true);
    f.check(test_points(m), tol);
}

TEST(ActuatorSampling, ExactForLinearFieldOnStretchedMesh)
{
    MeshSpec m;
    m.stretch = 1.15;
    Fields f(m, false);
    f.check(test_points(m), tol);
}

TEST(ActuatorSampling, ExactForLinearFieldOverTerrain)
{
    MeshSpec m;
    m.stretch = 1.1;
    m.hill = 60.0;
    Fields f(m, false);
    f.check(test_points(m), tol);
}

// The same points sampled with the domain in one box and split into several boxes agree with
// the analytic field, so a point on or near a box face is found once and read through the
// ghost cells rather than from a neighbour's valid data.
// sample_cell_scalar: bilinear between the four cell centres and linear in the physical height
// between the two centres that bracket the point, so a field linear in x, y and z_phys is exact on
// a terrain-following mesh (cell-centre heights the mean of the cell's eight nodes); below the
// first and above the last centre the nearest pair extrapolates, which a linear field does not notice
TEST(ActuatorSampling, CellScalarIsExactForALinearFieldOverTerrain)
{
    MeshSpec m;
    m.hill = 60.0;
    Fields f(m, false);
    const std::array<Real,4> c{{ 300.0, 0.004, -0.003, 0.0098 }};
    amrex::MultiFab s(f.ba, f.dm, 2, 2);
    for (amrex::MFIter mfi(s, false); mfi.isValid(); ++mfi) {
        auto a = s.array(mfi);
        amrex::LoopOnCpu(amrex::grow(mfi.validbox(), 2), [&](int i, int j, int k) {
            const Real x = (i + 0.5) * m.dx(), y = (j + 0.5) * m.dy();
            Real z = 0.0;
            for (int dk = 0; dk < 2; ++dk) { for (int dj = 0; dj < 2; ++dj) { for (int di = 0; di < 2; ++di) {
                z += 0.125 * m.z_node(i + di, j + dj, k + dk);
            } } }
            a(i,j,k,0) = -7.0;                               // another component, never read
            a(i,j,k,1) = LinearField::eval(c, x, y, z);
        });
    }
    std::vector<Real> pos = test_points(m);
    pos.push_back(0.3 * m.Lx); pos.push_back(0.7 * m.Ly); pos.push_back(0.9995 * m.H);   // above the last centre
    std::vector<Real> val;
    erf_actuator::sample_cell_scalar(s, 1, f.znd.get(), f.geom, pos, val);
    ASSERT_EQ(val.size(), pos.size() / 3);
    for (std::size_t p = 0; p < val.size(); ++p) {
        const Real e = LinearField::eval(c, pos[3*p], pos[3*p+1], pos[3*p+2]);
        EXPECT_NEAR(val[p], e, tol * std::max(Real(1.0), std::abs(e))) << "point " << p;
    }
}

// The point list must be the same on every rank: the reductions are collective. One rank passing a
// shorter list is caught by the collective count check and every rank stops, naming the counts
// (in the parallel executable; with one rank any list agrees with itself and the call succeeds)
TEST(ActuatorSampling, EveryRankMustPassTheSamePoints)
{
    MeshSpec m;
    Fields f(m, true);
    std::vector<Real> pos = test_points(m);
    std::vector<Real> vel;
    erf_gtest::ThrowOnAbort guard;
    if (amrex::ParallelDescriptor::NProcs() == 1) {
        EXPECT_NO_THROW(erf_actuator::sample_velocity(f.u, f.v, f.w, f.znd.get(), f.geom, pos, vel));
        return;
    }
    if (amrex::ParallelDescriptor::MyProc() != 0) { pos.resize(pos.size() - 3); }
    try {
        erf_actuator::sample_velocity(f.u, f.v, f.w, f.znd.get(), f.geom, pos, vel);
        ADD_FAILURE() << "a rank-dependent point list was accepted";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("ranks disagree on the number of points"), std::string::npos) << e.what();
    }
}

TEST(ActuatorSampling, IndependentOfBoxDecomposition)
{
    MeshSpec m;
    m.hill = 60.0;
    m.max_grid = {{4, 5, 1024}};     // 3 x 2 boxes of 4 x 5 cells, never split in z
    Fields f(m, false);
    std::vector<Real> pos = test_points(m);
    // points exactly on the box faces x = 4 dx and y = 5 dy
    auto add = [&](Real x, Real y, Real z) { pos.push_back(x); pos.push_back(y); pos.push_back(z); };
    add(4 * m.dx(), 5 * m.dy(), Real(0.5) * m.H);
    add(4 * m.dx(), Real(0.31) * m.Ly, Real(0.3) * m.H);
    add(Real(0.73) * m.Lx, 5 * m.dy(), Real(0.7) * m.H);
    f.check(pos, tol);
    // split in z as well (boxes of 3, 3 and 2 cells: a fine level stopping below the top is the lowest of
    // them): each box searches only its own column and first ghost layer, and the points on and near the
    // boxes' z faces sample as in one whole column
    MeshSpec mz = m;
    mz.max_grid = {{4, 5, 3}};
    Fields fz(mz, false);
    for (int k : {3, 6}) {
        const Real zf = mz.z_node(5, 5, k);
        add(5.5 * mz.dx(), 5.5 * mz.dy(), zf);
        add(5.5 * mz.dx(), 5.5 * mz.dy(), zf - Real(0.01) * mz.H);
        add(5.5 * mz.dx(), 5.5 * mz.dy(), zf + Real(0.01) * mz.H);
    }
    fz.check(pos, tol);
}

TEST(ActuatorSampling, AnEmptyPointListGivesNoVelocities)
{
    MeshSpec m;
    Fields f(m, true);
    std::vector<Real> vel{1.0, 2.0, 3.0};
    erf_actuator::sample_velocity(f.u, f.v, f.w, nullptr, f.geom, {}, vel);
    EXPECT_TRUE(vel.empty());
}

// On a refined level the grids are not the whole domain: a point is usable only if the cells it
// and its reach need are all on the level.
TEST(ActuatorSampling, CoverageByALevelIncludesTheReach)
{
    // a 3000 x 1200 x 600 m domain on 50 m cells, periodic in x and y; the level covers only
    // x in [500, 1500), y in [300, 900), z in [0, 450) as two boxes
    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(59, 23, 11));
    amrex::RealBox rb({0.0, 0.0, 0.0}, {3000.0, 1200.0, 600.0});
    amrex::Geometry geom(domain, rb, 0, {1, 1, 0});
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(10, 6, 0), amrex::IntVect(19, 17, 8)));
    bl.push_back(amrex::Box(amrex::IntVect(20, 6, 0), amrex::IntVect(29, 17, 8)));
    amrex::BoxArray ba(bl);
    std::string outside;
    // a point well inside with a 150 m reach: covered
    EXPECT_TRUE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 150.0}, 150.0, outside));
    EXPECT_TRUE(outside.empty());
    // the same point with a reach that crosses the level's x edge (500 m): not covered, named
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 150.0}, 300.0, outside));
    EXPECT_NE(outside.find("750"), std::string::npos);
    // a point across the two boxes' shared face: the union covers it
    EXPECT_TRUE(erf_actuator::points_covered_by(ba, geom, {1000.0, 600.0, 150.0}, 100.0, outside));
    // a point outside the level, and one above its top
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {2000.0, 600.0, 150.0}, 0.0, outside));
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 500.0}, 0.0, outside));
    // several points: the first outside is the one named
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 150.0, 2500.0, 600.0, 150.0, 2600.0, 600.0, 150.0}, 0.0, outside));
    EXPECT_NE(outside.find("2500"), std::string::npos);
    // a level covering the whole domain covers everything, reach or not
    amrex::BoxArray whole(domain);
    EXPECT_TRUE(erf_actuator::points_covered_by(whole, geom, {10.0, 10.0, 10.0, 2990.0, 1190.0, 590.0}, 400.0, outside));
    // an empty point list is covered
    EXPECT_TRUE(erf_actuator::points_covered_by(ba, geom, {}, 100.0, outside));
}

// Whether the grids cover a point does not depend on how the level is cut into boxes: the same region as one
// box and as unequal, non-mirror boxes (at most 7 x 5 x 3 cells) gives the same answer for every point of a sweep
// across its edges and seams, every reach and every z rule
TEST(ActuatorSampling, CoverageIsIndependentOfTheBoxDecomposition)
{
    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(59, 23, 11));
    amrex::RealBox rb({0.0, 0.0, 0.0}, {3000.0, 1200.0, 600.0});
    amrex::Geometry geom(domain, rb, 0, {1, 1, 0});
    const amrex::Box region(amrex::IntVect(10, 6, 0), amrex::IntVect(29, 17, 8));
    const amrex::BoxArray one(region);
    amrex::BoxArray cut(region);
    cut.maxSize(amrex::IntVect(7, 5, 3));
    ASSERT_GT(cut.size(), 20);
    int points = 0, covered = 0;
    for (const auto z : {erf_actuator::CoverZ::Reach, erf_actuator::CoverZ::Footprint, erf_actuator::CoverZ::Column}) {
        for (const amrex::Real reach : {amrex::Real(0.0), amrex::Real(60.0), amrex::Real(130.0)}) {
            for (amrex::Real x = 380.0; x < 1650.0; x += 37.0) {
                for (amrex::Real y = 190.0; y < 1010.0; y += 41.0) {
                    for (amrex::Real h = 5.0; h < 560.0; h += 53.0) {
                        std::string a, b;
                        const bool ca = erf_actuator::points_covered_by(one, geom, {x, y, h}, reach, a, z);
                        const bool cb = erf_actuator::points_covered_by(cut, geom, {x, y, h}, reach, b, z);
                        ASSERT_EQ(ca, cb) << "(" << x << ", " << y << ", " << h << ") reach " << reach << " rule " << static_cast<int>(z);
                        EXPECT_EQ(a, b);
                        ++points;
                        if (ca) { ++covered; }
                    }
                }
            }
        }
    }
    // the sweep holds points on both sides of the answer
    EXPECT_GT(covered, points / 10);
    EXPECT_LT(covered, points - points / 10);
}

// In a periodic direction the cells a point needs wrap across the seam: a level that covers part
// of the periodic extent must also hold the cells on the far side of the seam
TEST(ActuatorSampling, CoverageWrapsAcrossAPeriodicSeam)
{
    // 3000 m in x on 50 m cells (60 cells), periodic in x and y
    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(59, 23, 11));
    amrex::RealBox rb({0.0, 0.0, 0.0}, {3000.0, 1200.0, 600.0});
    amrex::Geometry geom(domain, rb, 0, {1, 1, 0});
    std::string outside;
    // the level holds the first half in x: a point 25 m from the seam with a 100 m reach needs the
    // cells 58 and 59 across it, which the level does not hold
    const amrex::BoxArray half(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(29, 23, 11)));
    EXPECT_FALSE(erf_actuator::points_covered_by(half, geom, {25.0, 600.0, 150.0}, 100.0, outside));
    EXPECT_NE(outside.find("25"), std::string::npos) << outside;
    // the same point well inside the half: covered
    EXPECT_TRUE(erf_actuator::points_covered_by(half, geom, {750.0, 600.0, 150.0}, 100.0, outside));
    // a level holding both sides of the seam covers the point, from either side
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(9, 23, 11)));
    bl.push_back(amrex::Box(amrex::IntVect(50, 0, 0), amrex::IntVect(59, 23, 11)));
    const amrex::BoxArray seam(bl);
    EXPECT_TRUE(erf_actuator::points_covered_by(seam, geom, {25.0, 600.0, 150.0}, 100.0, outside));
    EXPECT_TRUE(erf_actuator::points_covered_by(seam, geom, {2975.0, 600.0, 150.0}, 100.0, outside));
    // a reach wider than the period needs the whole periodic extent
    EXPECT_FALSE(erf_actuator::points_covered_by(seam, geom, {25.0, 600.0, 150.0}, 4000.0, outside));
    // the wrap is in y too, and a non-periodic z is clamped to the domain
    const amrex::BoxArray low_y(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(59, 11, 11)));
    EXPECT_FALSE(erf_actuator::points_covered_by(low_y, geom, {750.0, 25.0, 10.0}, 100.0, outside));
    EXPECT_TRUE(erf_actuator::points_covered_by(amrex::BoxArray(domain), geom, {750.0, 600.0, 10.0}, 100.0, outside));
}

// On a terrain-following or stretched mesh a point's k is not its height over dz: for the sampler a level
// whose grids stand over the point's column, at any height, covers it, and one beside it does not; for the
// spreading the grids must cover the whole height
TEST(ActuatorSampling, CoverageWithoutTheZCheckNeedsTheFootprints)
{
    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(59, 23, 11));
    amrex::RealBox rb({0.0, 0.0, 0.0}, {3000.0, 1200.0, 600.0});
    amrex::Geometry geom(domain, rb, 0, {1, 1, 0});
    // two stacked boxes over x in [500, 1500), the upper one not reaching the ground, and one beside them
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(10, 6, 0), amrex::IntVect(29, 17, 3)));
    bl.push_back(amrex::Box(amrex::IntVect(10, 6, 4), amrex::IntVect(29, 17, 8)));
    const amrex::BoxArray ba(bl);
    std::string outside;
    // 500 m up, above the level's top (450 m) by its nominal cells: refused with the z check, not without it
    using erf_actuator::CoverZ;
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 500.0}, 100.0, outside));
    EXPECT_TRUE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 500.0}, 100.0, outside, CoverZ::Footprint));
    // across, the reach still counts
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 500.0}, 300.0, outside, CoverZ::Footprint));
    // a box raised off the ground over the column covers it too
    const amrex::BoxArray raised(amrex::Box(amrex::IntVect(10, 6, 5), amrex::IntVect(29, 17, 8)));
    EXPECT_TRUE(erf_actuator::points_covered_by(raised, geom, {750.0, 600.0, 50.0}, 100.0, outside, CoverZ::Footprint));
    EXPECT_FALSE(erf_actuator::points_covered_by(raised, geom, {1750.0, 600.0, 50.0}, 100.0, outside, CoverZ::Footprint));
    // for the spreading neither covers it: the stacked boxes stop at 450 m, the raised one starts at 250 m
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 100.0}, 100.0, outside, CoverZ::Column));
    EXPECT_FALSE(erf_actuator::points_covered_by(raised, geom, {750.0, 600.0, 300.0}, 100.0, outside, CoverZ::Column));
    // two boxes stacked over the whole height do
    amrex::BoxList full;
    full.push_back(amrex::Box(amrex::IntVect(10, 6, 0), amrex::IntVect(29, 17, 5)));
    full.push_back(amrex::Box(amrex::IntVect(10, 6, 6), amrex::IntVect(29, 17, 11)));
    EXPECT_TRUE(erf_actuator::points_covered_by(amrex::BoxArray(full), geom, {750.0, 600.0, 100.0}, 100.0, outside, CoverZ::Column));
    // with the mesh's bounds: on this flat 50 m mesh the cells within reach of 300 m (k 4..8) are on the stacked
    // boxes (to k 8), as with the reach in z, and those of 400 m (k 6..10) are not; with ground from 0 to 200 m and
    // levels 30 to 50 m apart the heights within reach of 300 m may lie in any of k 0..11, above the boxes' top
    erf_actuator::ZBounds flat;
    flat.ground_lo = 0.0; flat.ground_hi = 0.0; flat.dz_lo = 50.0; flat.dz_hi = 50.0;
    EXPECT_TRUE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 300.0}, 100.0, outside, CoverZ::Column, &flat));
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 400.0}, 100.0, outside, CoverZ::Column, &flat));
    erf_actuator::ZBounds hilly;
    hilly.ground_lo = 0.0; hilly.ground_hi = 200.0; hilly.dz_lo = 30.0; hilly.dz_hi = 50.0;
    EXPECT_FALSE(erf_actuator::points_covered_by(ba, geom, {750.0, 600.0, 300.0}, 100.0, outside, CoverZ::Column, &hilly));
}

namespace {

// A terrain-following level z = h(i) + (H - h(i)) k / nz on the boxes of bl (8 x 4 x 10 cells of 100 m across,
// H = 500 m, periodic in y): its node heights, its bounds over all ranks and its table over all ranks
struct FittedLevel {
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::MultiFab znd;
    erf_actuator::ZBounds bounds, table;
    template <class Ground>
    FittedLevel (const amrex::BoxList& bl, const Ground& h)
    {
        const int nz = 10;
        const Real H = 500.0;
        const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 3, nz - 1));
        const amrex::RealBox rb({0.0, 0.0, 0.0}, {800.0, 400.0, H});
        geom = amrex::Geometry(domain, rb, 0, {0, 1, 0});
        ba = amrex::BoxArray(bl);
        const amrex::DistributionMapping dm(ba);
        znd.define(amrex::convert(ba, amrex::IntVect(1, 1, 1)), dm, 1, 0);
        for (amrex::MFIter mfi(znd, false); mfi.isValid(); ++mfi) {
            auto z = znd.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { z(i,j,k) = h(i) + (H - h(i)) * Real(k) / Real(nz); });
        }
        bounds = erf_actuator::mesh_z_bounds(znd, geom);
        table = erf_actuator::mesh_z_levels(znd, geom);
    }
    bool touches (std::array<Real,2> x, std::array<Real,2> z, const erf_actuator::ZBounds& zb, std::array<Real,2> y = {{100.0, 200.0}}) const
    {
        return erf_actuator::region_touches_grids(ba, geom, amrex::RealBox({x[0], y[0], z[0]}, {x[1], y[1], z[1]}), &zb);
    }
};

} // namespace

// region_touches_grids on a uniform mesh: the region's own cells, its periodic images however far, and a region
// that reaches to infinity or lies far outside the domain (its cell indices are bounded before they become
// integers).
TEST(ActuatorSampling, ARegionTouchesALevelThroughItsOwnCellsAndItsPeriodicImages)
{
    // 3000 x 1200 x 600 m on 50 m cells, periodic in x; the level holds x in [500, 1500), y in [300, 900),
    // z in [200, 450)
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(59, 23, 11));
    const amrex::RealBox prob({0.0, 0.0, 0.0}, {3000.0, 1200.0, 600.0});
    const amrex::Geometry geom(domain, prob, 0, {1, 0, 0});
    const amrex::BoxArray ba(amrex::Box(amrex::IntVect(10, 6, 4), amrex::IntVect(29, 17, 8)));
    auto touches = [&](std::array<Real,3> lo, std::array<Real,3> hi) {
        return erf_actuator::region_touches_grids(ba, geom, amrex::RealBox(lo.data(), hi.data()));
    };
    std::feclearexcept(FE_ALL_EXCEPT);
    EXPECT_TRUE(touches({700.0, 500.0, 300.0}, {800.0, 600.0, 350.0}));       // inside
    EXPECT_TRUE(touches({400.0, 500.0, 300.0}, {510.0, 600.0, 350.0}));       // across the x edge
    EXPECT_FALSE(touches({400.0, 500.0, 300.0}, {499.0, 600.0, 350.0}));      // beside it in x
    EXPECT_FALSE(touches({700.0, 100.0, 300.0}, {800.0, 299.0, 350.0}));      // beside it in y
    EXPECT_FALSE(touches({700.0, 500.0, 0.0}, {800.0, 600.0, 199.0}));        // under it
    EXPECT_TRUE(touches({700.0, 500.0, 0.0}, {800.0, 600.0, 201.0}));         // reaching its lowest cell
    EXPECT_FALSE(touches({700.0, 500.0, 450.0}, {800.0, 600.0, 600.0}));      // above it
    // x is periodic: the region's images one and several periods on (3000 m each) lie in the level; an image
    // beside the level stays beside it, and a region across the seam reaches the level from the other side
    EXPECT_TRUE(touches({-2300.0, 500.0, 300.0}, {-2200.0, 600.0, 350.0}));
    EXPECT_TRUE(touches({3700.0, 500.0, 300.0}, {3800.0, 600.0, 350.0}));
    EXPECT_TRUE(touches({6700.0, 500.0, 300.0}, {6800.0, 600.0, 350.0}));
    EXPECT_TRUE(touches({-8300.0, 500.0, 300.0}, {-8200.0, 600.0, 350.0}));
    EXPECT_FALSE(touches({7700.0, 500.0, 300.0}, {7800.0, 600.0, 350.0}));    // x 1700 to 1800: beside
    EXPECT_FALSE(touches({-7300.0, 500.0, 300.0}, {-7200.0, 600.0, 350.0}));  // x 1700 to 1800: beside
    EXPECT_TRUE(touches({5900.0, 500.0, 300.0}, {6600.0, 600.0, 350.0}));     // x 2900 to 3600: wraps to 0 to 600
    EXPECT_FALSE(touches({5900.0, 500.0, 300.0}, {6400.0, 600.0, 350.0}));    // wraps to 0 to 400: beside
    EXPECT_TRUE(touches({-9000.0, 500.0, 300.0}, {-5500.0, 600.0, 350.0}));   // longer than a period: every cell
    // y is not periodic: an image there does not count
    EXPECT_FALSE(touches({700.0, -700.0, 300.0}, {800.0, -600.0, 350.0}));
    // far outside and infinite extents: an answer, not an out-of-range conversion
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real big = std::numeric_limits<Real>::max();
    EXPECT_TRUE(touches({-inf, -inf, -inf}, {inf, inf, inf}));
    EXPECT_TRUE(touches({-big, -big, -big}, {big, big, big}));
    EXPECT_FALSE(touches({700.0, 500.0, Real(1.0e30)}, {800.0, 600.0, inf}));
    EXPECT_FALSE(touches({700.0, 500.0, -inf}, {800.0, 600.0, Real(-1.0e30)}));
    EXPECT_FALSE(touches({700.0, Real(1.0e30), 300.0}, {800.0, big, 350.0}));
    EXPECT_TRUE(touches({Real(1.0e30), 500.0, 300.0}, {inf, 600.0, 350.0}));  // periodic and endless: every cell
    EXPECT_NO_THROW(touches({Real(1.0e30), 500.0, 300.0}, {Real(1.0e30), 600.0, 350.0}));
    EXPECT_TRUE(touches({inf, 500.0, 300.0}, {inf, 600.0, 350.0}));           // both ends endless: every cell too
    EXPECT_TRUE(touches({-inf, 500.0, 300.0}, {-inf, 600.0, 350.0}));
    // a level with no box touches nothing
    EXPECT_FALSE(erf_actuator::region_touches_grids(amrex::BoxArray(), geom, amrex::RealBox({0.0, 0.0, 0.0}, {3000.0, 1200.0, 600.0})));
    // no invalid operation so far, the infinite regions included (it would trap with amrex.fpe_trap_invalid = 1)
    EXPECT_EQ(std::fetestexcept(FE_INVALID), 0);
    // an inverted or NaN region is refused
    EXPECT_NE(erf_gtest::abort_message([&] { touches({800.0, 500.0, 300.0}, {700.0, 600.0, 350.0}); }).find("lo <= hi"), std::string::npos);
    const Real qnan = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_NE(erf_gtest::abort_message([&] { touches({qnan, 500.0, 300.0}, {700.0, 600.0, 350.0}); }).find("lo <= hi"), std::string::npos);
}

// A height range with an infinite end, against the level's table (ZBounds::k_range): the cells are those of
// the finite end alone, a range at one infinity holds no cell, and no invalid operation is raised (the margin
// must not take inf - inf, which traps with amrex.fpe_trap_invalid = 1). Boxes on cells 0 to 2: ground from 0
// to 280 m, the top node level between 150 and 346 m.
TEST(ActuatorSampling, AnInfiniteHeightRangeTakesTheCellsOfItsFiniteEnd)
{
    const FittedLevel f(amrex::BoxList(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 3, 2))),
                        [](int i) { return Real(40.0) * Real(std::min(i, 7)); });
    ASSERT_TRUE(f.table.has_table());
    const Real inf = std::numeric_limits<Real>::infinity();
    std::feclearexcept(FE_ALL_EXCEPT);
    EXPECT_FALSE(f.touches({100.0, 200.0}, {inf, inf}, f.table));
    EXPECT_FALSE(f.touches({100.0, 200.0}, {-inf, -inf}, f.table));
    EXPECT_TRUE(f.touches({100.0, 200.0}, {-inf, inf}, f.table));
    EXPECT_FALSE(f.touches({100.0, 200.0}, {-inf, -10.0}, f.table));   // below the lowest ground
    EXPECT_TRUE(f.touches({100.0, 200.0}, {-inf, 5.0}, f.table));
    EXPECT_FALSE(f.touches({100.0, 200.0}, {400.0, inf}, f.table));    // above the highest node (346 m)
    EXPECT_TRUE(f.touches({100.0, 200.0}, {300.0, inf}, f.table));
    EXPECT_EQ(std::fetestexcept(FE_INVALID), 0);
    int kb = 0, kt = 0;
    f.table.k_range(-inf, Real(5.0), 0, 9, kb, kt);
    EXPECT_EQ(kb, 0);
    EXPECT_EQ(kt, 0);   // 5 m is below node level 1 everywhere (50 m at its lowest)
    f.table.k_range(Real(330.0), inf, 0, 9, kb, kt);
    EXPECT_EQ(kb, 2);   // node level 2 reaches 280 + 44 = 324 m, level 3 reaches 346 m: cell 2 only
    EXPECT_EQ(kt, 2);
    EXPECT_EQ(std::fetestexcept(FE_INVALID), 0);
}

// LevelHeights keeps the table of mesh_z_levels: built once for a set of boxes, again for other boxes or
// another domain, and on the same boxes only when asked (heights that change in time).
TEST(ActuatorSampling, TheLevelsTableIsKeptUntilTheBoxesChange)
{
    auto ground = [](int i) { return Real(40.0) * Real(i); };
    FittedLevel a(amrex::BoxList(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 3, 2))), ground);
    FittedLevel b(amrex::BoxList(amrex::Box(amrex::IntVect(4, 0, 5), amrex::IntVect(7, 3, 9))), ground);
    erf_actuator::LevelHeights h;
    EXPECT_EQ(h.builds(), 0);
    EXPECT_EQ(h.get(a.znd, a.geom).node_lo, a.table.node_lo);
    EXPECT_EQ(h.get(a.znd, a.geom).has_cell, a.table.has_cell);
    EXPECT_EQ(h.builds(), 1);
    // other boxes
    EXPECT_EQ(h.get(b.znd, b.geom).has_cell, b.table.has_cell);
    EXPECT_EQ(h.get(b.znd, b.geom).node_hi, b.table.node_hi);
    EXPECT_EQ(h.builds(), 2);
    // the heights 10 m higher on the same boxes: the kept table, until the caller asks
    b.znd.plus(Real(10.0), 0, 1, 0);
    EXPECT_EQ(h.get(b.znd, b.geom).node_lo, b.table.node_lo);
    EXPECT_EQ(h.builds(), 2);
    EXPECT_NEAR(h.get(b.znd, b.geom, true).node_lo[5], b.table.node_lo[5] + Real(10.0), 1.0e-3);
    EXPECT_EQ(h.builds(), 3);
    // the same boxes in a taller domain
    const amrex::Box tall(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 3, 11));
    const amrex::Geometry gtall(tall, amrex::RealBox({0.0, 0.0, 0.0}, {800.0, 400.0, 600.0}), 0, {0, 1, 0});
    EXPECT_EQ(h.get(b.znd, gtall).node_lo.size(), std::size_t(13));
    EXPECT_EQ(h.builds(), 4);
}

// A terrain-following level whose boxes all stand on the ground and hold cells 0 to 2. Ground from 0 (i = 0) to
// 280 m (i = 7), node levels 22 to 50 m apart. With the bounds alone the heights' cells come from the ground and
// spacing extremes, in whatever column; with the table from the node heights of the level itself, which is tighter.
TEST(ActuatorSampling, ARegionOverTerrainTouchesTheCellsItsHeightsMayLieIn)
{
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 2)));
    bl.push_back(amrex::Box(amrex::IntVect(4, 0, 0), amrex::IntVect(7, 3, 2)));
    const FittedLevel f(bl, [](int i) { return Real(40.0) * Real(std::min(i, 7)); });
    ASSERT_TRUE(f.bounds.all_ground);
    ASSERT_TRUE(f.bounds.has_ground());
    ASSERT_FALSE(f.bounds.has_table());
    ASSERT_TRUE(f.table.has_table());
    // node level 3, the boxes' top, lies between 150 m (i = 0) and 346 m (i = 7)
    EXPECT_NEAR(f.table.node_lo[3], 150.0, 1.0e-3);
    EXPECT_NEAR(f.table.node_hi[3], 346.0, 1.0e-3);
    for (const erf_actuator::ZBounds* zb : {&f.bounds, &f.table}) {
        SCOPED_TRACE(zb->has_table() ? "table" : "bounds");
        // heights 100 to 150 m: bounds, cells floor((100 - 280) / 50) = -4 to floor(150 / 22) = 6; table, cells 0 to 2
        EXPECT_TRUE(f.touches({100.0, 200.0}, {100.0, 150.0}, *zb));
        // heights 450 to 480 m, inside the domain: bounds, no cell below floor((450 - 280) / 50) = 3; table, above
        // the highest node of the boxes (346 m)
        EXPECT_FALSE(f.touches({100.0, 200.0}, {450.0, 480.0}, *zb));
        // below every ground height
        EXPECT_FALSE(f.touches({100.0, 200.0}, {-500.0, -400.0}, *zb));
        // the footprint decides first
        EXPECT_TRUE(f.touches({-50.0, 10.0}, {100.0, 150.0}, *zb));
        EXPECT_FALSE(f.touches({-500.0, -10.0}, {100.0, 150.0}, *zb));
    }
    // 350 to 400 m is above the boxes' highest node (346 m), which only the table knows: the bounds allow
    // cell floor((350 - 280) / 50) = 1
    EXPECT_TRUE(f.touches({100.0, 200.0}, {350.0, 400.0}, f.bounds));
    EXPECT_FALSE(f.touches({100.0, 200.0}, {350.0, 400.0}, f.table));
}

// erf-model/ERF#5056. A refined level aloft (no box at the ground) has no ground heights: its bounds hold the
// largest and the lowest Real, and (z - ground) / dz is no cell index. The table of its node heights still tells
// which heights its cells hold: a region in the box touches, one under it or beside it does not. With the bounds
// alone the whole column is taken, and a region under the box touches too.
TEST(ActuatorSampling, ALevelAloftTouchesARegionAtItsHeights)
{
    // cells 5 to 9 over x in [400, 800) m; flat ground, so the box holds z in [250, 500) m
    const FittedLevel f(amrex::BoxList(amrex::Box(amrex::IntVect(4, 0, 5), amrex::IntVect(7, 3, 9))), [](int) { return Real(0.0); });
    ASSERT_FALSE(f.bounds.all_ground);
    ASSERT_FALSE(f.bounds.has_ground());
    ASSERT_TRUE(f.table.has_table());
    for (const erf_actuator::ZBounds* zb : {&f.bounds, &f.table}) {
        SCOPED_TRACE(zb->has_table() ? "table" : "bounds");
        EXPECT_TRUE(f.touches({500.0, 600.0}, {300.0, 350.0}, *zb));    // in the box
        EXPECT_TRUE(f.touches({500.0, 600.0}, {200.0, 260.0}, *zb));    // reaching its lowest cell
        EXPECT_FALSE(f.touches({100.0, 300.0}, {300.0, 350.0}, *zb));   // beside it
    }
    // under the box, in its columns
    EXPECT_FALSE(f.touches({500.0, 600.0}, {10.0, 60.0}, f.table));
    EXPECT_FALSE(f.touches({500.0, 600.0}, {10.0, 240.0}, f.table));
    EXPECT_TRUE(f.touches({500.0, 600.0}, {10.0, 60.0}, f.bounds));
}

// Boxes on the ground and a box aloft together. The ground heights of the bounds are those of the box on the
// ground (0 m); the spacings are over all boxes (20 to 50 m). They say nothing of the cells under the box aloft:
// over its 300 m plateau a height of 425 m is in cell 6, which the box aloft holds (cells 5 to 7), while the
// bounds give cells floor(420 / 50) = 8 and up. With the bounds alone the whole column must be taken; the table
// gives the cells themselves.
TEST(ActuatorSampling, ABoxAloftBesideABoxOnTheGroundIsFoundAtItsOwnHeights)
{
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 9)));
    bl.push_back(amrex::Box(amrex::IntVect(4, 0, 5), amrex::IntVect(7, 3, 7)));
    const FittedLevel f(bl, [](int i) { return (i <= 4) ? Real(0.0) : ((i == 5) ? Real(150.0) : Real(300.0)); });
    ASSERT_FALSE(f.bounds.all_ground);
    ASSERT_TRUE(f.bounds.has_ground());
    EXPECT_NEAR(f.bounds.ground_hi, 0.0, 1.0e-6);
    EXPECT_NEAR(f.bounds.dz_lo, 20.0, 1.0e-4);
    EXPECT_NEAR(f.bounds.dz_hi, 50.0, 1.0e-4);
    // x in [650, 750): columns i = 6 and 7, on the plateau; heights 420 to 430 m are in cell 6 there
    int kb = 0, kt = 0;
    f.bounds.k_range(420.0, 430.0, 0, 9, kb, kt);
    EXPECT_GT(kb, 7) << "the ground and spacing extremes miss the box aloft (cells 5 to 7)";
    f.table.k_range(420.0, 430.0, 0, 9, kb, kt);
    EXPECT_LE(kb, 6);
    EXPECT_GE(kt, 6);
    for (const erf_actuator::ZBounds* zb : {&f.bounds, &f.table}) {
        SCOPED_TRACE(zb->has_table() ? "table" : "bounds");
        EXPECT_TRUE(f.touches({650.0, 750.0}, {420.0, 430.0}, *zb));
        // one period on in y (y is periodic on this level) the region is the same one
        EXPECT_TRUE(f.touches({650.0, 750.0}, {420.0, 430.0}, *zb, {{500.0, 600.0}}));
    }
    // under the box aloft, over the plateau (300 to 390 m: below its lowest node there, 400 m, but at the heights
    // of the ground box's cells): the level holds cells at these heights only elsewhere. The table answers for the
    // level as a whole, so it says yes here; a region below every node of the level's cells over x >= 400 m does
    // not exist on this mesh, as the ground box holds every height
    EXPECT_TRUE(f.touches({650.0, 750.0}, {300.0, 390.0}, f.table));
    // above the box aloft and beside the ground box: cells 8 and 9 are the ground box's, over x < 400 m only
    EXPECT_FALSE(f.touches({650.0, 750.0}, {470.0, 490.0}, f.table));
}

// The table over all ranks is that of the boxes whichever rank holds them: the same from every rank, and equal
// to a one-box table of the same nodes.
TEST(ActuatorSampling, TheLevelsTableIsTheSameWhateverTheBoxes)
{
    auto ground = [](int i) { return Real(40.0) * Real(std::min(i, 7)); };
    amrex::BoxList split;
    split.push_back(amrex::Box(amrex::IntVect(0, 0, 2), amrex::IntVect(3, 3, 5)));
    split.push_back(amrex::Box(amrex::IntVect(4, 0, 2), amrex::IntVect(7, 1, 5)));
    split.push_back(amrex::Box(amrex::IntVect(4, 2, 2), amrex::IntVect(7, 3, 5)));
    split.push_back(amrex::Box(amrex::IntVect(0, 0, 6), amrex::IntVect(7, 3, 7)));
    const FittedLevel a(split, ground);
    const FittedLevel b(amrex::BoxList(amrex::Box(amrex::IntVect(0, 0, 2), amrex::IntVect(7, 3, 7))), ground);
    ASSERT_EQ(a.table.node_lo.size(), b.table.node_lo.size());
    EXPECT_EQ(a.table.k0, b.table.k0);
    for (std::size_t n = 0; n < a.table.has_cell.size(); ++n) {
        EXPECT_EQ(a.table.has_cell[n], b.table.has_cell[n]) << "cell " << n;
        EXPECT_EQ(a.table.has_cell[n] != 0, n >= 2 && n <= 7) << "cell " << n;
    }
    for (std::size_t n = 2; n <= 8; ++n) {
        EXPECT_EQ(a.table.node_lo[n], b.table.node_lo[n]) << "node level " << n;
        EXPECT_EQ(a.table.node_hi[n], b.table.node_hi[n]) << "node level " << n;
        // z = h + (500 - h) n / 10 with h from 0 to 280 m (node i = 7; i = 8 has the same h)
        EXPECT_NEAR(a.table.node_lo[n], Real(50.0) * Real(n), 1.0e-3);
        EXPECT_NEAR(a.table.node_hi[n], Real(280.0) + Real(22.0) * Real(n), 1.0e-3);
    }
}

// mesh_z_bounds: the lowest and highest ground and the smallest and largest level spacing of a
// terrain-following mesh z = h + (H - h) k / nz over a ramp h = 4 i (m), found across boxes
TEST(ActuatorSampling, TheMeshsZBoundsAreItsGroundAndSpacingExtremes)
{
    const int nx = 8, nz = 10;
    const Real H = 500.0;
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, 3, nz - 1));
    const amrex::RealBox rb({0.0, 0.0, 0.0}, {800.0, 400.0, 500.0});
    const amrex::Geometry geom(domain, rb, 0, {0, 0, 0});
    amrex::BoxArray ba(domain);
    ba.maxSize(4);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab znd(amrex::convert(ba, amrex::IntVect(1, 1, 1)), dm, 1, 0);
    for (amrex::MFIter mfi(znd, false); mfi.isValid(); ++mfi) {
        auto z = znd.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const Real h = Real(4.0) * Real(i);
            z(i,j,k) = h + (H - h) * Real(k) / Real(nz);
            (void) j;
        });
    }
    const erf_actuator::ZBounds b = erf_actuator::mesh_z_bounds(znd, geom);
    EXPECT_NEAR(b.ground_lo, 0.0, 1.0e-4);
    EXPECT_NEAR(b.ground_hi, 4.0 * nx, 1.0e-4);
    EXPECT_NEAR(b.dz_lo, (H - Real(4.0 * nx)) / nz, 1.0e-3);
    EXPECT_NEAR(b.dz_hi, H / nz, 1.0e-3);
}

// A point a few ulps below a non-zero prob_lo wraps to lo + (len - tiny), which rounds onto prob_hi, owned by
// no box: it lands on prob_lo instead
TEST(ActuatorSampling, AWrappedPointNeverLandsOnTheUpperFace)
{
    const Real lo = Real(-177.507), hi = Real(8056.149);
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(63, 7, 7));
    const amrex::RealBox rb({lo, 0.0, 0.0}, {hi, 800.0, 800.0});
    const amrex::Geometry geom(domain, rb, 0, {1, 0, 0});
    const Real below = std::nextafter(static_cast<Real>(lo), static_cast<Real>(-1.0e30));
    const std::vector<Real> p = erf_actuator::wrap_periodic({below, 400.0, 400.0}, geom);
    EXPECT_LT(p[0], static_cast<Real>(geom.ProbHi(0))) << "wrapped onto the upper face";
    EXPECT_GE(p[0], static_cast<Real>(geom.ProbLo(0)));
    // every point a few ulps either side of the seam stays inside [lo, hi)
    Real x = static_cast<Real>(lo);
    for (int n = 0; n < 8; ++n) { x = std::nextafter(x, static_cast<Real>(-1.0e30)); }
    for (int n = 0; n < 16; ++n, x = std::nextafter(x, static_cast<Real>(1.0e30))) {
        const Real w = erf_actuator::wrap_periodic({x, 400.0, 400.0}, geom)[0];
        EXPECT_TRUE(w >= static_cast<Real>(geom.ProbLo(0)) && w < static_cast<Real>(geom.ProbHi(0))) << n << " " << w;
    }
}

// The domain's upper x, y and z faces belong to the last cells, on every mesh type; in a periodic
// direction a point past the seam takes its image's value
TEST(ActuatorSampling, TheUpperFaceIsTheLastCellAndPeriodicPointsTheirImage)
{
    MeshSpec m;
    m.stretch = 1.0;
    m.max_grid = {{5, 4, 1024}};
    Fields f(m, true);
    // on the domain's upper x, y and z faces (an end of a line placed on the boundary, a point at the top)
    f.check({m.Lx, Real(0.37) * m.Ly, Real(0.5) * m.H, Real(0.41) * m.Lx, m.Ly, Real(0.5) * m.H, m.Lx, m.Ly, Real(0.25) * m.H,
             Real(0.29) * m.Lx, Real(0.53) * m.Ly, m.H, m.Lx, m.Ly, m.H}, tol);
    {
        // the top on a stretched and on a terrain-following mesh (z_phys_nd present), the box split in z
        MeshSpec ms;
        ms.max_grid = {{1024, 1024, 3}};
        Fields fs(ms, false);
        fs.check({Real(0.29) * ms.Lx, Real(0.53) * ms.Ly, ms.H}, tol);
        MeshSpec mt;
        mt.hill = 60.0;
        Fields ft(mt, false);
        ft.check({Real(0.29) * mt.Lx, Real(0.53) * mt.Ly, mt.H, mt.Lx, mt.Ly, mt.H}, tol);
    }
    {
        // 1000 m on 88 cells: prob_hi / dx rounds one ulp above 88 in double, as 5000 m on 96 cells does in float
        MeshSpec m88 = m;
        m88.Lx = (sizeof(Real) == 8) ? 1000.0 : 5000.0;
        m88.nx = (sizeof(Real) == 8) ? 88 : 96;
        Fields f88(m88, true);
        f88.check({m88.Lx, Real(0.37) * m88.Ly, Real(0.5) * m88.H}, tol);
    }
    // periodic in x: a point a third of a cell past either seam takes its image's value
    const std::array<int,3> periodic{{1, 0, 0}};
    const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(m.Lx, m.Ly, m.H)});
    f.geom = amrex::Geometry(f.geom.Domain(), &rb, 0, periodic.data());
    const Real past = m.dx() / 3.0;
    const std::vector<Real> outside{m.Lx + past, Real(0.37) * m.Ly, Real(0.5) * m.H, -past, Real(0.61) * m.Ly, Real(0.25) * m.H};
    const std::vector<Real> image{past, Real(0.37) * m.Ly, Real(0.5) * m.H, m.Lx - past, Real(0.61) * m.Ly, Real(0.25) * m.H};
    std::vector<Real> vel, ref;
    erf_actuator::sample_velocity(f.u, f.v, f.w, nullptr, f.geom, outside, vel);
    erf_actuator::sample_velocity(f.u, f.v, f.w, nullptr, f.geom, image, ref);
    ASSERT_EQ(vel.size(), ref.size());
    for (std::size_t i = 0; i < vel.size(); ++i) { EXPECT_NEAR(vel[i], ref[i], tol * std::max(Real(1.0), std::abs(ref[i]))) << i; }
    const auto wrapped = erf_actuator::wrap_periodic(outside, f.geom);
    for (std::size_t i = 0; i < image.size(); ++i) { EXPECT_NEAR(wrapped[i], image[i], 1.0e-4 * m.Lx) << i; }
}

// A point outside the domain is found by no box: the sampler aborts and names it
TEST(ActuatorSampling, APointOutsideTheDomainIsRefusedNamingIt)
{
    MeshSpec m;
    Fields f(m, true);
    std::vector<Real> vel;
    const std::string msg = erf_gtest::abort_message([&] {
        erf_actuator::sample_velocity(f.u, f.v, f.w, nullptr, f.geom, {Real(5000.0), Real(500.0), Real(100.0)}, vel);
    });
    EXPECT_NE(msg.find("5000"), std::string::npos) << msg;
    EXPECT_NE(msg.find("sampled by 0 boxes"), std::string::npos) << msg;
}

// Valid boxes never overlap, so a point lies in one box. If a caller passes fields whose boxes do overlap, a point
// in the overlap is held by two boxes of the same rank: the count must say two, and sample_velocity must stop.
TEST(ActuatorSampling, APointHeldByTwoBoxesOfOneRankIsRefused)
{
    MeshSpec m;
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx-1, m.ny-1, m.nz-1));
    const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(m.Lx, m.Ly, m.H)});
    const std::array<int,3> periodic{{0, 0, 0}};
    const amrex::Geometry geom(domain, &rb, 0, periodic.data());
    // two boxes that share the cells i = 5..6, both on rank 0
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(6, m.ny-1, m.nz-1)));
    bl.push_back(amrex::Box(amrex::IntVect(5, 0, 0), amrex::IntVect(m.nx-1, m.ny-1, m.nz-1)));
    const amrex::BoxArray ba(bl);
    const amrex::DistributionMapping dm(amrex::Vector<int>{0, 0});
    amrex::MultiFab u(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 2);
    amrex::MultiFab v(amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 2);
    amrex::MultiFab w(amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 2);
    u.setVal(1.0); v.setVal(2.0); w.setVal(3.0);

    // one point in the shared cells (x = 5.5 dx) and one in the first box only
    const std::vector<Real> pos{Real(5.5) * m.dx(), Real(500.0), Real(100.0), Real(1.5) * m.dx(), Real(500.0), Real(100.0)};
    std::vector<Real> vel;
    std::vector<int> cnt;
    erf_actuator::sample_velocity_local(u, v, w, nullptr, geom, pos, vel, cnt);
    amrex::ParallelDescriptor::ReduceIntSum(cnt.data(), static_cast<int>(cnt.size()));
    ASSERT_EQ(cnt.size(), 2u);
    EXPECT_EQ(cnt[0], 2);
    EXPECT_EQ(cnt[1], 1);

    const std::string msg = erf_gtest::abort_message([&] {
        erf_actuator::sample_velocity(u, v, w, nullptr, geom, pos, vel);
    });
    EXPECT_NE(msg.find("sampled by 2 boxes"), std::string::npos) << msg;
}

// The terrain surface under a point: the k = 0 node plane of z_phys_nd, bilinear between the four
// nodes around (x, y); prob_lo z on a uniform-dz mesh
TEST(ActuatorSampling, TerrainHeightIsTheBilinearNodeSurface)
{
    MeshSpec m; m.hill = 80.0; m.max_grid = {{4, 5, 1024}};   // several boxes: each point has one owner
    Fields f(m, false);
    std::vector<Real> pos, expect;
    // at nodes: exact
    for (int i : {0, 3, 7, m.nx}) { for (int j : {0, 4, m.ny}) {
        pos.insert(pos.end(), {i * m.dx(), j * m.dy(), 999.0}); expect.push_back(m.z_node(i, j, 0));
    } }
    // at a cell centre and at general points: bilinear between the four nodes
    auto bilin = [&](Real x, Real y) {
        const int i = static_cast<int>(std::floor(x / m.dx())), j = static_cast<int>(std::floor(y / m.dy()));
        const Real wx = x / m.dx() - i, wy = y / m.dy() - j;
        return (1 - wy) * ((1 - wx) * m.z_node(i, j, 0) + wx * m.z_node(i+1, j, 0)) + wy * ((1 - wx) * m.z_node(i, j+1, 0) + wx * m.z_node(i+1, j+1, 0));
    };
    for (auto xy : {std::array<Real,2>{{250.0, 350.0}}, std::array<Real,2>{{612.5, 137.25}}, std::array<Real,2>{{1199.0, 999.0}}}) {
        pos.insert(pos.end(), {xy[0], xy[1], 0.0}); expect.push_back(bilin(xy[0], xy[1]));
    }
    std::vector<Real> h;
    erf_actuator::terrain_heights(f.znd.get(), f.geom, pos, h);
    ASSERT_EQ(h.size(), expect.size());
    // heights up to 80 m: 1e-9 m in double, a few float spacings (8e-6 m at 80 m) in single precision
    const Real htol = (sizeof(Real) == 8) ? Real(1.0e-9) : Real(5.0e-5);
    for (std::size_t p = 0; p < h.size(); ++p) { EXPECT_NEAR(h[p], expect[p], htol) << "point " << p; }
    EXPECT_GT(*std::max_element(h.begin(), h.end()), 40.0);   // the hill is really there
    // a uniform-dz mesh: the domain floor everywhere
    Fields flat(m, true);
    erf_actuator::terrain_heights(flat.znd.get(), flat.geom, pos, h);
    for (const Real v : h) { EXPECT_DOUBLE_EQ(v, 0.0); }
    // no points: nothing to do
    erf_actuator::terrain_heights(f.znd.get(), f.geom, {}, h);
    EXPECT_TRUE(h.empty());
}

// On a slope the ground under a point (bilinear between the four nodes of its cell) lies below the
// cell's bottom face (their mean) on one side: a point just above the ground there is sampled from
// the bottom cell, by the column's extrapolation, still exactly for a linear field, on any split
TEST(ActuatorSampling, APointJustAboveSlopingGroundIsSampled)
{
    MeshSpec m; m.stretch = 1.1; m.hill = 120.0; m.max_grid = {{4, 5, 1024}};
    Fields f(m, false);
    std::vector<Real> xy;
    for (int i = 0; i < 40; ++i) { for (int j = 0; j < 30; ++j) { xy.insert(xy.end(), {Real((i + 0.37) * m.Lx / 40), Real((j + 0.61) * m.Ly / 30), Real(0.0)}); } }
    std::vector<Real> h;
    erf_actuator::terrain_heights(f.znd.get(), f.geom, xy, h);
    std::vector<Real> pos;
    int below_face = 0;
    for (std::size_t p = 0; p < h.size(); ++p) {
        const Real x = xy[3*p], y = xy[3*p+1];
        pos.insert(pos.end(), {x, y, h[p] + Real(0.5)});
        const int i = static_cast<int>(std::floor(x / m.dx())), j = static_cast<int>(std::floor(y / m.dy()));
        const Real face = 0.25 * (m.z_node(i, j, 0) + m.z_node(i+1, j, 0) + m.z_node(i, j+1, 0) + m.z_node(i+1, j+1, 0));
        if (h[p] + 0.5 < face) { ++below_face; }
    }
    EXPECT_GT(below_face, 50) << "the slope must put points below their cell's bottom face";
    f.check(pos, tol);
}

namespace {

template <class F>
amrex::MultiFab make_znd (const amrex::BoxArray& level_cells, const F& z)
{
    const amrex::DistributionMapping dm(level_cells);
    amrex::MultiFab znd(amrex::convert(level_cells, amrex::IntVect(1, 1, 1)), dm, 1, 0);
    for (amrex::MFIter mfi(znd, false); mfi.isValid(); ++mfi) {
        auto za = znd.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { za(i,j,k) = z(i, j, k); });
    }
    return znd;
}

} // namespace

// column_cell_ranges: the cells a kernel's heights may lie in, over the point's own columns. A plateau
// beside a deep canyon: a level of 32 of the domain's 64 cells over a plateau 300 m
// up and a canyon floor 40 m up, on a terrain-following mesh z = h + (H - h) k / nz. A point on the
// plateau at 700 m with a reach of 60 m needs the cells 22 to 30 of its column, which the level holds; the
// mesh's bounds (the canyon's ground with the plateau's spacing) put the same heights in cells up to 47 and
// refuse it. A point above its own column's top (790 m) is refused by both, as is one over the canyon,
// whose column's top is 660 m; the per-column answer is sharp to the top node's height.
TEST(ActuatorSampling, ColumnCellsAcceptAColumnBesideADeepCanyonAndRefuseOneAboveItsTop)
{
    using erf_actuator::CoverZ;
    const int nx = 40, ny = 4, nz = 64;
    const Real Lx = 4000.0, Ly = 400.0, H = 1280.0;
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    const amrex::RealBox rb({0.0, 0.0, 0.0}, {Lx, Ly, H});
    const amrex::Geometry geom(domain, rb, 0, {0, 0, 0});
    amrex::BoxArray level(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, 31)));
    level.maxSize(amrex::IntVect(8, 8, 1024));   // columns whole, as ERF's levels keep them (amr.no_box_split_dir)
    auto h = [&](Real x) { return (x <= 1500.0) ? Real(300.0) : (x >= 2500.0) ? Real(40.0) : Real(300.0) - Real(260.0) * (x - Real(1500.0)) / Real(1000.0); };
    const Real dx = Lx / nx;
    const amrex::MultiFab znd = make_znd(level, [&](int i, int, int k) { const Real hs = h(Real(i) * dx); return hs + (H - hs) * Real(k) / Real(nz); });
    // the mesh's bounds over the level: the canyon's ground, the plateau's spacing
    const erf_actuator::ZBounds old = erf_actuator::mesh_z_bounds(znd, geom);
    ASSERT_TRUE(old.all_ground);
    EXPECT_NEAR(old.ground_lo, 40.0, 1.0e-6);
    EXPECT_NEAR(old.dz_lo, (H - 300.0) / nz, 1.0e-6);
    const std::vector<Real> pos{500.0, 200.0, 700.0,    // on the plateau, within its column's top of 790 m
                                500.0, 200.0, 800.0,    // on the plateau, above it
                                3500.0, 200.0, 700.0};  // over the canyon, whose column tops at 660 m
    const std::vector<Real> reach(3, 60.0);
    erf_actuator::ColumnCells cells;
    erf_actuator::column_cell_ranges(znd, geom, pos, reach, cells);
    ASSERT_EQ(cells.kb.size(), 3u);
    // 640 to 760 m in a column of 15.3125 m cells from 300 m: cells 22 to 30
    EXPECT_EQ(cells.kb[0], 22);
    EXPECT_EQ(cells.kt[0], 30);
    // 740 to 860 m: cells 28 to 31, and the cell above the level's top node (790 m), which it crosses
    EXPECT_EQ(cells.kb[1], 28);
    EXPECT_EQ(cells.kt[1], 32);
    // over the canyon (19.375 m cells from 40 m): 640 to 760 m is cells 30 and 31 and above the top
    EXPECT_EQ(cells.kb[2], 30);
    EXPECT_EQ(cells.kt[2], 32);
    std::string outside;
    for (int p = 0; p < 3; ++p) {
        const std::vector<Real> one(pos.begin() + 3*p, pos.begin() + 3*p + 3);
        const erf_actuator::ColumnCells c = cells.slice(p, 1);
        const bool per_column = erf_actuator::points_covered_by(level, geom, one, 60.0, outside, CoverZ::Column, nullptr, &c);
        const bool bounds = erf_actuator::points_covered_by(level, geom, one, 60.0, outside, CoverZ::Column, &old);
        EXPECT_EQ(per_column, p == 0) << "point " << p;
        EXPECT_FALSE(bounds) << "point " << p;
    }
    // sharp: a reach of 89 m from 700 m stays under the top node, 91 m crosses it
    for (const Real r : {Real(89.0), Real(91.0)}) {
        erf_actuator::ColumnCells c;
        erf_actuator::column_cell_ranges(znd, geom, {500.0, 200.0, 700.0}, {r}, c);
        EXPECT_EQ(erf_actuator::points_covered_by(level, geom, {500.0, 200.0, 700.0}, r, outside, CoverZ::Column, nullptr, &c), r < Real(90.0)) << r;
    }
    // the ranges do not depend on how the level is cut into boxes, in z too
    for (const int ms : {4, 16, 1024}) {
        amrex::BoxArray cut(level.minimalBox());
        cut.maxSize(ms);
        const amrex::MultiFab z2 = make_znd(cut, [&](int i, int, int k) { const Real hs = h(Real(i) * dx); return hs + (H - hs) * Real(k) / Real(nz); });
        erf_actuator::ColumnCells c;
        erf_actuator::column_cell_ranges(z2, geom, pos, reach, c);
        EXPECT_EQ(c.kb, cells.kb) << "max box size " << ms;
        EXPECT_EQ(c.kt, cells.kt) << "max box size " << ms;
    }
}

// column_cell_ranges against the cells themselves: a point is covered when, at every node column its reach
// across touches (a periodic x wrapping) and at every height within its reach between that column's floor
// and top, the cell holding the height is on the level in every cell column next to the node column. On a
// level whose boxes span the same cells in every column (at the ground, aloft, part of the domain across),
// over this gentle hill, the per-column check gives exactly that; on one whose boxes stop at different
// heights in neighbouring columns it never covers a point the cells do not. The answers do not depend on
// how the level is cut into boxes.
TEST(ActuatorSampling, ColumnCellsAreTheCellsThemselvesAndNeverFewer)
{
    using erf_actuator::CoverZ;
    MeshSpec m;
    m.stretch = 1.12;
    m.hill = 120.0;
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx - 1, m.ny - 1, m.nz - 1));
    const amrex::RealBox rb({0.0, 0.0, 0.0}, {m.Lx, m.Ly, m.H});
    const amrex::Geometry geom(domain, rb, 0, {1, 0, 0});
    auto znode = [&](int i, int j, int k) { return m.z_node(((i % m.nx) + m.nx) % m.nx, j, k); };
    // the levels: the whole domain across to cell 5; aloft from cell 2 to 6; part of the domain across;
    // boxes to cell 3 everywhere and to cell 6 over x >= 600 m only (the uneven one)
    struct Level { const char* name; amrex::BoxList boxes; bool uniform; };
    std::vector<Level> levels;
    levels.push_back({"to cell 5", amrex::BoxList(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx - 1, m.ny - 1, 5))), true});
    levels.push_back({"aloft", amrex::BoxList(amrex::Box(amrex::IntVect(0, 0, 2), amrex::IntVect(m.nx - 1, m.ny - 1, 6))), true});
    levels.push_back({"part across", amrex::BoxList(amrex::Box(amrex::IntVect(2, 1, 0), amrex::IntVect(8, m.ny - 2, 6))), true});
    {
        amrex::BoxList bl;
        bl.push_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx - 1, m.ny - 1, 3)));
        bl.push_back(amrex::Box(amrex::IntVect(6, 0, 4), amrex::IntVect(m.nx - 1, m.ny - 1, 6)));
        levels.push_back({"uneven", bl, false});
    }
    // the points: a fixed pseudo-random sweep across the domain (beyond it in x, which wraps), the whole
    // height and below and above it, with reaches up to 1.5 cells across
    std::vector<Real> pos, reach;
    unsigned long long seed = 12345;
    auto next = [&seed]() { seed = seed * 6364136223846793005ULL + 1442695040888963407ULL; return Real((seed >> 11) % 1000000) / Real(1.0e6); };
    const int npts = 400;
    for (int p = 0; p < npts; ++p) {
        pos.push_back(-200.0 + (m.Lx + 400.0) * next());
        pos.push_back(m.Ly * next());
        pos.push_back(-30.0 + (m.H + 60.0) * next());
        reach.push_back(1.5 * m.dx() * next());
    }
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    for (const auto& lv : levels) {
        SCOPED_TRACE(lv.name);
        const amrex::BoxArray level(lv.boxes);
        // the reference, from the heights of every node column within reach
        std::vector<char> ref(npts, 1);
        for (int p = 0; p < npts; ++p) {
            const Real zl = pos[3*p+2] - reach[p], zh = pos[3*p+2] + reach[p];
            int clo[2], chi[2];
            for (int d = 0; d < 2; ++d) {
                clo[d] = static_cast<int>(std::floor((pos[3*p+d] - reach[p] - plo[d]) * dxi[d]));
                chi[d] = static_cast<int>(std::floor((pos[3*p+d] + reach[p] - plo[d]) * dxi[d]));
            }
            clo[1] = std::max(clo[1], 0); chi[1] = std::min(chi[1], m.ny - 1);   // y is not periodic: clamped, as the check clamps it
            bool covered = true, touched = false;
            for (int ic = clo[0]; ic <= chi[0] && covered; ++ic) {
                for (int jc = clo[1]; jc <= chi[1] && covered; ++jc) {
                    const int i = ((ic % m.nx) + m.nx) % m.nx;   // the cell column, wrapped
                    for (int ni = i; ni <= i + 1 && covered; ++ni) {
                        for (int nj = jc; nj <= jc + 1 && covered; ++nj) {
                            // every height within reach between this node column's floor and top lies in one of its cells
                            for (int k = 0; k < m.nz && covered; ++k) {
                                const Real z0 = znode(ni, nj, k), z1 = znode(ni, nj, k + 1);
                                if (z0 <= zh && z1 >= zl) {
                                    touched = true;
                                    if (!level.contains(amrex::Box(amrex::IntVect(i, jc, k), amrex::IntVect(i, jc, k)))) { covered = false; }
                                }
                            }
                        }
                    }
                }
            }
            // a reach wholly above the domain's top or below its floor lies in no cell: not covered
            ref[p] = (covered && touched) ? 1 : 0;
        }
        int covered_ref = 0;
        for (const char c : ref) { covered_ref += c; }
        EXPECT_GT(covered_ref, 20);
        EXPECT_LT(covered_ref, npts - 20);
        // the check, on the level cut into boxes three ways
        std::vector<char> answer_first;
        for (const std::array<int,3> ms : {std::array<int,3>{{1024, 1024, 1024}}, std::array<int,3>{{4, 5, 3}}, std::array<int,3>{{3, 2, 2}}}) {
            amrex::BoxArray cut(lv.boxes);
            cut.maxSize(amrex::IntVect(ms[0], ms[1], ms[2]));
            const amrex::MultiFab znd = make_znd(cut, znode);
            erf_actuator::ColumnCells cells;
            erf_actuator::column_cell_ranges(znd, geom, pos, reach, cells);
            std::vector<char> answer(npts, 0);
            std::string outside;
            int covered_new = 0;
            for (int p = 0; p < npts; ++p) {
                const std::vector<Real> one(pos.begin() + 3*p, pos.begin() + 3*p + 3);
                const erf_actuator::ColumnCells c = cells.slice(p, 1);
                const bool covered = erf_actuator::points_covered_by(level, geom, one, reach[p], outside, CoverZ::Column, nullptr, &c);
                covered_new += covered ? 1 : 0;
                answer[p] = covered ? 1 : 0;
                if (lv.uniform) {
                    EXPECT_EQ(covered, ref[p] != 0) << "point " << p << " at (" << one[0] << ", " << one[1] << ", " << one[2]
                                                    << ") reach " << reach[p] << " cells " << c.kb[0] << ".." << c.kt[0];
                } else if (covered) {
                    EXPECT_TRUE(ref[p] != 0) << "point " << p << " at (" << one[0] << ", " << one[1] << ", " << one[2]
                                             << ") reach " << reach[p] << " covered, but not by the cells";
                }
            }
            // never everything refused: the points low over the taller boxes are covered on every level
            EXPECT_GT(covered_new, 20);
            // the same answers however the level is cut
            if (answer_first.empty()) { answer_first = answer; }
            EXPECT_EQ(answer, answer_first) << "boxes of " << ms[0] << " " << ms[1] << " " << ms[2];
        }
    }
}

// The level's boxes built once (CoverageGrids): its cells and footprints are the level's (checked directly),
// and a sweep across the level's edges and a periodic seam, with every reach and every z rule, on boxes
// cut unevenly, runs both entry points (the BoxArray one builds the same grids, so the sweep alone
// would not tell them apart)
TEST(ActuatorSampling, CoverageOnGridsBuiltOnceIsCoverageOnGridsBuiltEachTime)
{
    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(59, 23, 11));
    amrex::RealBox rb({0.0, 0.0, 0.0}, {3000.0, 1200.0, 600.0});
    amrex::Geometry geom(domain, rb, 0, {1, 1, 0});
    // a level across the x seam (cells 50..59 and 0..9) and an island, cut unevenly
    amrex::BoxList bl;
    bl.push_back(amrex::Box(amrex::IntVect(50, 4, 0), amrex::IntVect(59, 19, 8)));
    bl.push_back(amrex::Box(amrex::IntVect(0, 4, 0), amrex::IntVect(9, 19, 8)));
    bl.push_back(amrex::Box(amrex::IntVect(25, 8, 0), amrex::IntVect(34, 15, 11)));
    amrex::BoxArray ba(bl);
    ba.maxSize(amrex::IntVect(7, 5, 3));
    const erf_actuator::CoverageGrids grids(ba, geom);
    // the grids themselves: the cell-centred boxes, and each flattened onto the bottom layer, built here apart
    ASSERT_EQ(grids.cells.size(), ba.size());
    ASSERT_EQ(grids.footprints.size(), ba.size());
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
        EXPECT_EQ(grids.cells[b], ba[b]);
        amrex::Box flat = ba[b];
        flat.setSmall(2, 0);
        flat.setBig(2, 0);
        EXPECT_EQ(grids.footprints[b], flat);
    }
    erf_actuator::ZBounds zb;
    zb.ground_lo = 0.0; zb.ground_hi = 0.0; zb.dz_lo = 50.0; zb.dz_hi = 50.0;
    int points = 0, covered = 0;
    for (const auto z : {erf_actuator::CoverZ::Reach, erf_actuator::CoverZ::Footprint, erf_actuator::CoverZ::Column}) {
        for (const amrex::Real reach : {amrex::Real(0.0), amrex::Real(60.0), amrex::Real(130.0)}) {
            for (amrex::Real x = -90.0; x < 3090.0; x += 71.0) {
                for (amrex::Real y = 110.0; y < 1110.0; y += 83.0) {
                    for (amrex::Real h = 5.0; h < 590.0; h += 97.0) {
                        const std::vector<amrex::Real> pos{x, y, h, 1500.0, 600.0, 300.0};
                        std::string a, b;
                        const bool once = erf_actuator::points_covered_by(grids, geom, pos, reach, a, z, &zb);
                        const bool fresh = erf_actuator::points_covered_by(ba, geom, pos, reach, b, z, &zb);
                        ASSERT_EQ(once, fresh) << "(" << x << ", " << y << ", " << h << ") reach " << reach << " rule " << static_cast<int>(z);
                        EXPECT_EQ(a, b);
                        ++points;
                        if (once) { ++covered; }
                    }
                }
            }
        }
    }
    // the sweep is not vacuous: some points are covered and some are not
    EXPECT_GT(covered, 0);
    EXPECT_LT(covered, points);
}

// The local sampler with buffers kept from call to call: one kernel per box, the same answers as with
// its own buffers, also for a call with fewer points after one with more (the buffers keep their size):
// it gets its own answers, not the previous call's
TEST(ActuatorSampling, KeptBuffersGiveTheSameAnswersAndOneKernelPerBox)
{
    MeshSpec m;
    m.hill = 60.0;
    m.max_grid = {{4, 5, 1024}};
    Fields f(m, false);
    const std::vector<Real> many = test_points(m);
    // the first two points: the second call reads only the start of the buffers the first one filled
    const std::vector<Real> few(many.begin(), many.begin() + 6);
    erf_actuator::ActuatorDeviceBuffers buf;
    int nlocal = 0;
    for (amrex::MFIter mfi(f.u, false); mfi.isValid(); ++mfi) { ++nlocal; }
    for (const std::vector<Real>* pts : {&many, &few, &many}) {
        std::vector<Real> vel_kept, vel_own;
        std::vector<int> cnt_kept, cnt_own;
        const long before = buf.kernel_launches;
        erf_actuator::sample_velocity_local(f.u, f.v, f.w, f.znd.get(), f.geom, *pts, vel_kept, cnt_kept, &buf);
        EXPECT_EQ(buf.kernel_launches - before, nlocal);
        erf_actuator::sample_velocity_local(f.u, f.v, f.w, f.znd.get(), f.geom, *pts, vel_own, cnt_own);
        EXPECT_EQ(vel_kept, vel_own);
        EXPECT_EQ(cnt_kept, cnt_own);
        EXPECT_GE(buf.sample_vel.device.size(), many.size());
        EXPECT_GE(buf.sample_cnt.device.size(), many.size() / 3);
    }
}

// The table of one rank's node heights (mesh_z_bounds over this rank's boxes): the cells it gives for a
// range of heights hold every cell of those boxes whose corners reach into the range, over terrain on a
// stretched mesh with each column split in z, and none outside the domain's cells; above the mesh, none.
// (That the cells of no box of the rank are left out is the spreader's test of a rank holding only boxes aloft.)
TEST(ActuatorSampling, TheHeightTableHoldsEveryCellAHeightMayLieIn)
{
    MeshSpec m;
    m.stretch = 1.15;
    m.hill = 100.0;
    m.max_grid = {{4, 5, 3}};
    Fields f(m, false);
    const erf_actuator::ZBounds zb = erf_actuator::mesh_z_bounds(*f.znd, f.geom, false);
    if (amrex::ParallelDescriptor::NProcs() == 1) { ASSERT_TRUE(zb.has_table()); }
    if (!zb.has_table()) { return; }   // a rank whose boxes all reach the ground (or that holds none) has no table
    const int klo = 0, khi = m.nz - 1;
    int checked = 0;
    for (int a = 0; a < 40; ++a) {
        const Real zl = m.H * Real(a) / Real(41);
        for (const Real span : {Real(1.0), Real(20.0), Real(90.0)}) {
            const Real zh = zl + span;
            int kb = 0, kt = 0;
            zb.k_range(zl, zh, klo, khi, kb, kt);
            // every cell of this rank's boxes with a corner height in [zl, zh] or spanning it
            for (amrex::MFIter mfi(*f.znd, false); mfi.isValid(); ++mfi) {
                const amrex::Box cells = amrex::enclosedCells(mfi.validbox());
                const auto z = f.znd->const_array(mfi);
                amrex::LoopOnCpu(cells, [&](int i, int j, int k) {
                    Real lo = std::numeric_limits<Real>::max(), hi = std::numeric_limits<Real>::lowest();
                    for (int c = 0; c < 8; ++c) {
                        const Real zc = z(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));
                        lo = std::min(lo, zc); hi = std::max(hi, zc);
                    }
                    if (lo <= zh && hi >= zl) {
                        ++checked;
                        EXPECT_TRUE(k >= kb && k <= kt) << "cell k " << k << " for heights " << zl << " .. " << zh
                                                        << " outside [" << kb << ", " << kt << "]";
                    }
                });
            }
            if (kb <= kt) { EXPECT_TRUE(kb >= klo && kt <= khi); }
        }
    }
    EXPECT_GT(checked, 100);
    // above the mesh's top: no cell
    int kb = 0, kt = 0;
    zb.k_range(m.H + Real(10.0), m.H + Real(50.0), klo, khi, kb, kt);
    EXPECT_GT(kb, kt);
}

// points_covered_by takes the bounds over all ranks: bounds whose spacing is not positive (a z_phys_nd
// not yet filled) and one rank's table of node heights are refused, the first as before
TEST(ActuatorSampling, CoverageRefusesBoundsWithoutASpacingAndARanksTable)
{
    MeshSpec m;
    m.hill = 60.0;
    m.max_grid = {{1024, 1024, 3}};   // columns split in z: a rank with a box aloft gets a table
    Fields f(m, false);
    const amrex::BoxArray ba = amrex::convert(f.u.boxArray(), amrex::IntVect(0, 0, 0));
    const std::vector<Real> pos{Real(0.5 * m.Lx), Real(0.5 * m.Ly), Real(0.5 * m.H)};
    std::string outside;
    erf_actuator::ZBounds flat;   // a ground, and no spacing
    flat.ground_lo = 0.0;
    flat.ground_hi = 0.0;
    flat.dz_lo = 0.0;
    flat.dz_hi = 0.0;
    std::string msg = erf_gtest::abort_message([&] {
        erf_actuator::points_covered_by(ba, f.geom, pos, Real(50.0), outside, erf_actuator::CoverZ::Column, &flat);
    });
    EXPECT_NE(msg.find("the levels' spacing must be positive"), std::string::npos) << msg;
    const erf_actuator::ZBounds table = erf_actuator::mesh_z_bounds(*f.znd, f.geom, false);
    if (amrex::ParallelDescriptor::NProcs() == 1) { ASSERT_TRUE(table.has_table()); }
    if (table.has_table()) {   // a rank whose boxes all reach the ground has none
        msg = erf_gtest::abort_message([&] {
            erf_actuator::points_covered_by(ba, f.geom, pos, Real(50.0), outside, erf_actuator::CoverZ::Column, &table);
        });
        EXPECT_NE(msg.find("not one rank's table"), std::string::npos) << msg;
    }
    // the bounds over all ranks are taken
    const erf_actuator::ZBounds all = erf_actuator::mesh_z_bounds(*f.znd, f.geom);
    EXPECT_FALSE(all.has_table());
    EXPECT_TRUE(erf_actuator::points_covered_by(ba, f.geom, pos, Real(50.0), outside, erf_actuator::CoverZ::Column, &all)) << outside;
}
