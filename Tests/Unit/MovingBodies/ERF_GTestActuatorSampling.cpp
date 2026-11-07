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
// - wrap_periodic: a point a few ulps below a non-zero prob_lo lands on prob_lo, not on prob_hi.
// - terrain_heights: the bilinear k = 0 node surface.

#include <algorithm>
#include <array>
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
