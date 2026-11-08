// Contract of erf_actuator::spread_forces: the source spread from a point integrates back to
// the point's force exactly on every component, on a uniform mesh, on a terrain-following mesh
// with its cell volumes, when the kernel is cut off by the ground, across a periodic boundary,
// for several points at once, and however the domain is split; the source is zero beyond the
// kernel's reach; over a hill taller than the kernel's reach the source is centred on the point
// (the faces are found by their physical heights, not their index times dz); each point is spread
// with its own kernel width; and a non-finite position or force, or a width that is not positive,
// aborts, naming the point.

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include "ERF_ActuatorSpreading.H"
#include "ERF_GTestThrowOnAbort.H"

namespace {

using amrex::Real;

struct Mesh {
    int nx = 24, ny = 12, nz = 10;
    Real Lx = 2400.0, Ly = 1200.0, H = 500.0;
    Real hill = 0.0;                       // basic terrain-following (BTF) terrain amplitude (m); 0: flat
    std::array<int,3> max_grid{{1024, 1024, 1024}};
    std::array<int,3> periodic{{1, 1, 0}};

    Real dx () const { return Lx / nx; }
    Real dy () const { return Ly / ny; }
    Real h (Real x, Real y) const
    {
        constexpr Real pi = 3.14159265358979323846;
        return hill * (0.5 + 0.5 * std::cos(2.0 * pi * x / Lx)) * (0.5 + 0.5 * std::sin(2.0 * pi * y / Ly));
    }
    // the horizontal spacings, divided once for the loops that fill or read a field. Called inside those loops,
    // dx() and dy() let clang -O2 divide two-lane vectors with an unused lane: {Lx, 0} / {nx, 0} and
    // {Ly, H} / {ny, 0}. The unused lanes' 0 / 0 and H / 0 are discarded, but they raise the flags that
    // amrex.fpe_trap_invalid and fpe_trap_zero stop on
    struct Spacing { Real dx, dy; };
    Spacing spacing () const { return {dx(), dy()}; }
    // basic terrain-following (BTF) with uniform nominal levels: the Jacobian of a column is (H - h) / H
    Real z_node (const Spacing& s, int i, int j, int k) const
    {
        const Real hs = h(i * s.dx, j * s.dy);
        return hs + (H - hs) * k / nz;
    }
    Real detj (const Spacing& s, int i, int j) const { return (H - h((i + 0.5) * s.dx, (j + 0.5) * s.dy)) / H; }
};

struct Sources {
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab sx, sy, sz;
    std::unique_ptr<amrex::MultiFab> znd, detj;
    Mesh mesh;

    Sources (const Mesh& m, bool terrain)
        : mesh(m)
    {
        const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx-1, m.ny-1, m.nz-1));
        const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(m.Lx, m.Ly, m.H)});
        geom = amrex::Geometry(domain, &rb, 0, m.periodic.data());
        ba = amrex::BoxArray(domain);
        ba.maxSize(amrex::IntVect(m.max_grid[0], m.max_grid[1], m.max_grid[2]));
        dm = amrex::DistributionMapping(ba);
        sx.define(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 0);
        sy.define(amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 0);
        sz.define(amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 0);
        if (terrain) {
            const int ng = 1;
            znd = std::make_unique<amrex::MultiFab>(amrex::convert(ba, amrex::IntVect(1,1,1)), dm, 1, ng);
            detj = std::make_unique<amrex::MultiFab>(ba, dm, 1, ng);
            const Mesh::Spacing sp = m.spacing();
            for (amrex::MFIter mfi(*detj, false); mfi.isValid(); ++mfi) {
                auto za = znd->array(mfi);
                auto ja = detj->array(mfi);
                const amrex::Box gbx = amrex::grow(mfi.validbox(), ng);
                amrex::LoopOnCpu(amrex::surroundingNodes(gbx),
                                 [&](int i, int j, int k) { za(i,j,k) = m.z_node(sp, i, j, k); });
                amrex::LoopOnCpu(gbx, [&](int i, int j, int k) { ja(i,j,k) = m.detj(sp, i, j); (void) k; });
            }
        }
    }

    // spread and return the integrated source per component
    std::array<Real,3> spread (const std::vector<Real>& pos, const std::vector<Real>& force, Real eps)
    {
        erf_actuator::spread_forces(pos, force, eps, znd.get(), detj.get(), geom, sx, sy, sz);
        return totals();
    }
    // the same with a kernel width per point
    std::array<Real,3> spread (const std::vector<Real>& pos, const std::vector<Real>& force, const std::vector<Real>& eps)
    {
        erf_actuator::spread_forces(pos, force, eps, znd.get(), detj.get(), geom, sx, sy, sz);
        return totals();
    }
    std::array<Real,3> totals () const
    {
        return {{erf_actuator::integrate_source(0, sx, detj.get(), geom),
                 erf_actuator::integrate_source(1, sy, detj.get(), geom),
                 erf_actuator::integrate_source(2, sz, detj.get(), geom)}};
    }
};

// exact up to the roundoff of sums of a few thousand terms (a few units of 1e-7 in single
// precision); the scale of a component is the sum of the magnitudes of the forces spread,
// since the totals may cancel to zero
constexpr Real rel_tol = (sizeof(Real) == 8) ? Real(1.0e-11) : Real(1.0e-5);

std::vector<Real> xyz (double x, double y, double z)
{
    return {static_cast<Real>(x), static_cast<Real>(y), static_cast<Real>(z)};
}

void expect_total (const std::array<Real,3>& got, const std::array<Real,3>& want,
                   const std::array<Real,3>& scale)
{
    for (int d = 0; d < 3; ++d) {
        SCOPED_TRACE("component " + std::to_string(d));
        EXPECT_NEAR(got[d], want[d], rel_tol * std::max(Real(1.0), scale[d]));
    }
}

std::array<Real,3> abs_sum (const std::vector<Real>& force)
{
    std::array<Real,3> s{{0.0, 0.0, 0.0}};
    for (std::size_t p = 0; p < force.size() / 3; ++p) { for (int d = 0; d < 3; ++d) { s[d] += std::abs(force[3*p+d]); } }
    return s;
}

} // namespace

TEST(ActuatorSpreading, PointForceIntegratesBackOnUniformMesh)
{
    Mesh m;
    Sources s(m, false);
    const Real eps = 2.0 * m.dx();
    const std::vector<Real> pos = xyz(0.4 * m.Lx, 0.55 * m.Ly, 0.5 * m.H);
    const std::vector<Real> force = xyz(-2.0e6, 3.0e5, -1.0e4);
    expect_total(s.spread(pos, force, eps), {{force[0], force[1], force[2]}}, abs_sum(force));
}

TEST(ActuatorSpreading, SourceIsZeroBeyondThreeEpsilon)
{
    Mesh m;
    Sources s(m, false);
    const Real eps = 2.0 * m.dx();
    const std::vector<Real> pos = xyz(0.4 * m.Lx, 0.55 * m.Ly, 0.5 * m.H);
    const std::vector<Real> force = xyz(-2.0e6, 0.0, 0.0);
    s.spread(pos, force, eps);
    int nonzero_far = 0, nonzero_near = 0;
    std::string first_far;
    const Mesh::Spacing sp = m.spacing();
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) {
        auto a = s.sx.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const Real x = i * sp.dx, y = (j + 0.5) * sp.dy, z = (k + 0.5) * m.H / m.nz;
            // minimum-image distance: the domain is periodic in x and y
            Real ddx = x - pos[0], ddy = y - pos[1];
            ddx -= m.Lx * std::round(ddx / m.Lx);
            ddy -= m.Ly * std::round(ddy / m.Ly);
            const Real r = std::sqrt(ddx*ddx + ddy*ddy + (z-pos[2])*(z-pos[2]));
            if (a(i,j,k) != 0.0) {
                if (r > 3.0 * eps) {
                    ++nonzero_far;
                    if (first_far.empty()) {
                        first_far = "face (" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) +
                                    ") at r = " + std::to_string(r) + " has " + std::to_string(a(i,j,k));
                    }
                } else {
                    ++nonzero_near;
                }
            }
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(nonzero_far);
    amrex::ParallelDescriptor::ReduceIntSum(nonzero_near);
    EXPECT_EQ(nonzero_far, 0) << first_far;
    EXPECT_GT(nonzero_near, 0);
}

TEST(ActuatorSpreading, PointForceIntegratesBackOverTerrain)
{
    Mesh m;
    m.hill = 100.0;
    Sources s(m, true);
    const Real eps = 2.0 * m.dx();
    // over the hill's flank, where neighbouring columns have different volumes
    const std::vector<Real> pos = xyz(0.2 * m.Lx, 0.3 * m.Ly, 0.4 * m.H);
    const std::vector<Real> force = xyz(-1.5e6, 2.0e5, 4.0e4);
    expect_total(s.spread(pos, force, eps), {{force[0], force[1], force[2]}}, abs_sum(force));
}

// The face volumes are half of each neighbouring cell's dx dy dz detJ: with a detJ that varies in
// k as well as in x and y, the source summed with the test's own face volumes equals the force.
// integrate_source is not used here, so a wrong neighbour in face_volume (the cell above instead
// of the cell below a z face, a single cell instead of the average) would not cancel out.
TEST(ActuatorSpreading, TerrainVolumesAreTheFaceAveragesOfDetJ)
{
    Mesh m;
    m.hill = 100.0;
    Sources s(m, true);
    const Mesh::Spacing sp = m.spacing();
    auto dj = [&](int i, int j, int k) { return m.detj(sp, i, j) * (Real(1.0) + Real(0.5) * k / m.nz); };
    for (amrex::MFIter mfi(*s.detj, false); mfi.isValid(); ++mfi) {
        auto ja = s.detj->array(mfi);
        amrex::LoopOnCpu(amrex::grow(mfi.validbox(), 1), [&](int i, int j, int k) { ja(i,j,k) = dj(i, j, k); });
    }
    const Real eps = 2.0 * m.dx();
    const std::vector<Real> pos = xyz(0.2 * m.Lx, 0.3 * m.Ly, 0.4 * m.H);
    const std::vector<Real> force = xyz(-1.5e6, 2.0e5, 4.0e4);
    erf_actuator::spread_forces(pos, force, eps, s.znd.get(), s.detj.get(), s.geom, s.sx, s.sy, s.sz);
    const double dv = static_cast<double>(sp.dx) * static_cast<double>(sp.dy) * static_cast<double>(m.H / m.nz);
    std::array<double,3> total{{0.0, 0.0, 0.0}};
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) {
        const auto a = s.sx.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            if (i == m.nx) { return; }                                   // the periodic image of face 0
            total[0] += static_cast<double>(a(i,j,k)) * 0.5 * (dj(i-1,j,k) + dj(i,j,k)) * dv;
        });
    }
    for (amrex::MFIter mfi(s.sy, false); mfi.isValid(); ++mfi) {
        const auto a = s.sy.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            if (j == m.ny) { return; }
            total[1] += static_cast<double>(a(i,j,k)) * 0.5 * (dj(i,j-1,k) + dj(i,j,k)) * dv;
        });
    }
    for (amrex::MFIter mfi(s.sz, false); mfi.isValid(); ++mfi) {
        const auto a = s.sz.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            total[2] += static_cast<double>(a(i,j,k)) * 0.5 * (dj(i,j,k-1) + dj(i,j,k)) * dv;   // boundary faces hold no source
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(total.data(), 3);
    const std::array<Real,3> scale = abs_sum(force);
    for (int d = 0; d < 3; ++d) {
        EXPECT_NEAR(static_cast<Real>(total[d]), force[d], rel_tol * std::max(Real(1.0), scale[d])) << "component " << d;
    }
}

// Over a 400 m hill on 40 levels (15 m apart at its crest, 25 m nominal) a point 120 m above the crest
// sits at nominal index 20 but between the crest's levels 7 and 8: the kernel's faces are those within
// 3 eps of it in physical height, so the source's centroid is the point's height
TEST(ActuatorSpreading, OverATallHillTheSourceIsCentredOnThePoint)
{
    Mesh m;
    m.nz = 40;
    m.H = 1000.0;
    m.hill = 400.0;
    Sources s(m, true);
    const Real eps = 0.5 * m.dx();
    // the crest: x = 0, y = Ly / 4
    const Real zp = m.hill + 120.0;
    const std::vector<Real> pos = xyz(0.0, 0.25 * m.Ly, zp);
    const std::vector<Real> force = xyz(-1.0e5, 0.0, 0.0);
    expect_total(s.spread(pos, force, eps), {{force[0], force[1], force[2]}}, abs_sum(force));
    // the x source's centroid in physical height, each face at the mean of its four nodes, its weight scaled by
    // the column's Jacobian
    double num = 0.0, den = 0.0;
    const Mesh::Spacing sp = m.spacing();
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) {
        const auto a = s.sx.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            if (i == m.nx) { return; }   // the periodic image of face 0
            const double zc = 0.25 * (m.z_node(sp, i, j, k) + m.z_node(sp, i, j + 1, k)
                                      + m.z_node(sp, i, j, k + 1) + m.z_node(sp, i, j + 1, k + 1));
            const double wv = static_cast<double>(a(i,j,k)) * static_cast<double>(m.detj(sp, i, j));
            num += wv * zc;
            den += wv;
        });
    }
    // every rank's boxes: a rank may hold none near the point
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(den);
    ASSERT_NE(den, 0.0);
    EXPECT_NEAR(num / den, static_cast<double>(zp), 0.1 * static_cast<double>(eps));
}

TEST(ActuatorSpreading, GroundCutKernelStillIntegratesToTheForce)
{
    Mesh m;
    Sources s(m, false);
    const Real eps = 2.0 * m.dx();
    // a point 30 m above the ground: most of the 3 epsilon = 300 m sphere is underground
    const std::vector<Real> pos = xyz(0.5 * m.Lx, 0.5 * m.Ly, 30.0);
    const std::vector<Real> force = xyz(-8.0e5, 0.0, -2.0e5);
    expect_total(s.spread(pos, force, eps), {{force[0], force[1], force[2]}}, abs_sum(force));
}

// A point half a cell from the periodic x boundary: the kernel wraps, the faces near x = 0
// carry source, and the total is still the force.
TEST(ActuatorSpreading, KernelWrapsAcrossAPeriodicBoundary)
{
    Mesh m;
    Sources s(m, false);
    const Real eps = 2.0 * m.dx();
    const std::vector<Real> pos = xyz(m.Lx - 0.5 * m.dx(), 0.5 * m.Ly, 0.5 * m.H);
    const std::vector<Real> force = xyz(-1.0e6, 2.0e5, -3.0e4);
    expect_total(s.spread(pos, force, eps), {{force[0], force[1], force[2]}}, abs_sum(force));
    // the x face at i = 1 (x = dx, 1.5 dx past the periodic boundary) carries source, and the
    // periodic image faces i = 0 and i = nx hold the same value
    Real v1 = 0.0, v0 = 0.0, vn = 0.0;
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) {
        auto a = s.sx.const_array(mfi);
        const amrex::IntVect f1(1, m.ny / 2, m.nz / 2), f0(0, m.ny / 2, m.nz / 2), fn(m.nx, m.ny / 2, m.nz / 2);
        if (mfi.validbox().contains(f1)) { v1 = a(f1); }
        if (mfi.validbox().contains(f0)) { v0 = a(f0); }
        if (mfi.validbox().contains(fn)) { vn = a(fn); }
    }
    amrex::ParallelDescriptor::ReduceRealSum(v1);
    amrex::ParallelDescriptor::ReduceRealSum(v0);
    amrex::ParallelDescriptor::ReduceRealSum(vn);
    EXPECT_LT(v1, 0.0);
    EXPECT_NEAR(v0, vn, rel_tol * std::abs(v0));
    EXPECT_LT(v0, v1);   // the face at x = 0 is closer to the point than the one at x = dx
}

TEST(ActuatorSpreading, SeveralPointsSuperpose)
{
    Mesh m;
    Sources s(m, false);
    const Real eps = 2.0 * m.dx();
    std::vector<Real> pos, force;
    std::array<Real,3> total{{0.0, 0.0, 0.0}};
    for (int p = 0; p < 12; ++p) {
        // a ring of points, the rotor of radius 120 m at hub height 150 m
        constexpr Real pi = 3.14159265358979323846;
        const Real th = 2.0 * pi * p / 12;
        const auto q = xyz(0.3 * m.Lx, 0.5 * m.Ly + 120.0 * std::cos(th), 150.0 + 120.0 * std::sin(th));
        pos.insert(pos.end(), q.begin(), q.end());
        const std::array<Real,3> f{{static_cast<Real>(-1.0e5 * (1 + p)), static_cast<Real>(2.0e4 * std::sin(th)), static_cast<Real>(-2.0e4 * std::cos(th))}};
        force.insert(force.end(), f.begin(), f.end());
        for (int d = 0; d < 3; ++d) { total[d] += f[d]; }
    }
    expect_total(s.spread(pos, force, eps), total, abs_sum(force));
}

// One box and 3 x 3 boxes of 8 x 4 cells (never split in z) give the same integrated force and
// the same source at an x face that two boxes share.
TEST(ActuatorSpreading, IndependentOfBoxDecomposition)
{
    Mesh m;
    m.hill = 60.0;
    const Real eps = 2.0 * m.dx();
    const std::vector<Real> pos = xyz(8.0 * m.dx(), 5.0 * m.dy(), 0.5 * m.H);   // on the x box face i = 8 of the split
    const std::vector<Real> force = xyz(-1.2e6, 1.0e5, -5.0e4);
    Sources one(m, true);
    const auto t1 = one.spread(pos, force, eps);
    Mesh ms = m;
    ms.max_grid = {{8, 5, 1024}};
    Sources split(ms, true);
    const auto t2 = split.spread(pos, force, eps);
    expect_total(t1, {{force[0], force[1], force[2]}}, abs_sum(force));
    expect_total(t2, {{force[0], force[1], force[2]}}, abs_sum(force));
    // the value at the x face (8, 5, k) on the x box face i = 8, which two boxes hold in the split
    // layout: both copies must equal the single-box value
    const amrex::IntVect fc(8, 5, m.nz / 2);
    Real v_one = 0.0;
    for (amrex::MFIter mfi(one.sx, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(fc)) { v_one = one.sx.const_array(mfi)(fc); }
    }
    amrex::ParallelDescriptor::ReduceRealSum(v_one);
    EXPECT_NE(v_one, 0.0);
    int copies = 0;
    for (amrex::MFIter mfi(split.sx, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(fc)) {
            ++copies;
            EXPECT_NEAR(split.sx.const_array(mfi)(fc), v_one, rel_tol * std::abs(v_one));
        }
    }
    amrex::ParallelDescriptor::ReduceIntSum(copies);
    EXPECT_EQ(copies, 2);
}

TEST(ActuatorSpreading, ANonFinitePositionOrForceIsRefusedNamingThePoint)
{
    Mesh m;
    Sources s(m, false);
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    std::vector<Real> pos = xyz(0.4 * m.Lx, 0.55 * m.Ly, 0.5 * m.H), force = xyz(-2.0e6, 0.0, 0.0);
    pos.insert(pos.end(), {Real(0.5 * m.Lx), Real(0.5 * m.Ly), Real(0.5 * m.H)});
    force.insert(force.end(), {Real(1.0e5), nan, Real(0.0)});
    std::string msg = erf_gtest::abort_message([&] { s.spread(pos, force, 2.0 * m.dx()); });
    EXPECT_NE(msg.find("point 1 (0-based)"), std::string::npos) << msg;
    force[4] = 0.0;
    pos[0] = std::numeric_limits<Real>::infinity();
    msg = erf_gtest::abort_message([&] { s.spread(pos, force, 2.0 * m.dx()); });
    EXPECT_NE(msg.find("point 0 (0-based)"), std::string::npos) << msg;
}

TEST(ActuatorSpreading, EmptyPointListGivesZeroSources)
{
    Mesh m;
    Sources s(m, false);
    s.sx.setVal(7.0);
    const auto t = s.spread({}, {}, 2.0 * m.dx());
    expect_total(t, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}});
    EXPECT_EQ(s.sx.max(0), 0.0);
}

// Two points, 480 m apart on 20 m cells, spread with widths of 40 m and 80 m: the x-source of each,
// summed over its own half of the domain, integrates to its own force and has the second moment of
// its own kernel about it, in y and in z. For exp(-r^2/eps^2) cut off at 3 eps that is
// eps^2/2 * (1 - 0.0025077) (the cut-off's share, from the radial integrals); the face sums lie
// within 0.04 % of it at two and four cells per width. One width for both points (the first's)
// would put the second's moment four times too low.
TEST(ActuatorSpreading, EachPointIsSpreadWithItsOwnWidth)
{
    Mesh m;
    m.nx = 48; m.ny = 48; m.nz = 48;
    m.Lx = 960.0; m.Ly = 960.0; m.H = 960.0;
    m.periodic = {{0, 0, 0}};
    m.max_grid = {{16, 24, 48}};   // six boxes: the kernels cross box faces, and two ranks both own some
    Sources s(m, false);
    const std::array<Real,2> px{{247.0, 727.0}}, eps{{40.0, 80.0}}, fx{{-3.0e5, 2.0e5}};
    const Real py = 473.0, pz = 483.0;   // off the cell centres and faces
    std::vector<Real> pos, force;
    for (int p = 0; p < 2; ++p) {
        pos.insert(pos.end(), {px[p], py, pz});
        force.insert(force.end(), {fx[p], Real(0.0), Real(0.0)});
    }
    s.spread(pos, force, std::vector<Real>{eps[0], eps[1]});
    // per point (the half x < 480 m or above): sum of src dV, of src dV (y - py)^2 and of src dV (z - pz)^2
    // in double: some 15 000 faces carry the wider kernel, which a float sum would blur to 1e-5
    std::array<double,6> sums{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    const Real dv = m.dx() * m.dy() * (m.H / m.nz);
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) {
        // each x face once: a box's cells' low faces
        amrex::Box fb = amrex::convert(amrex::enclosedCells(mfi.validbox()), amrex::IntVect(1,0,0));
        fb.growHi(0, -1);
        const auto a = s.sx.const_array(mfi);
        amrex::LoopOnCpu(fb, [&](int i, int j, int k) {
            const int p = (i * m.dx() < 0.5 * m.Lx) ? 0 : 1;
            const double w = static_cast<double>(a(i,j,k)) * dv;
            const double y = (j + 0.5) * m.dy() - py, z = (k + 0.5) * (m.H / m.nz) - pz;
            sums[3*p] += w; sums[3*p+1] += w * y * y; sums[3*p+2] += w * z * z;
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(sums.data(), 6);
    constexpr double cut = 0.0025077;
    for (int p = 0; p < 2; ++p) {
        SCOPED_TRACE("point " + std::to_string(p) + ", width " + std::to_string(eps[p]) + " m");
        // the spreader's own normalisation is a sum over some 7 000 faces in Real: 1e-4 in single precision
        const double tol_total = (sizeof(Real) == 8) ? double(rel_tol) : 1.0e-4;
        EXPECT_NEAR(sums[3*p], fx[p], tol_total * std::abs(fx[p]));
        const double want = 0.5 * eps[p] * eps[p] * (1.0 - cut);
        EXPECT_NEAR(sums[3*p+1] / sums[3*p], want, 1.0e-3 * want);
        EXPECT_NEAR(sums[3*p+2] / sums[3*p], want, 1.0e-3 * want);
    }
}

TEST(ActuatorSpreading, AWidthThatIsMissingOrNotPositiveIsRefused)
{
    Mesh m;
    Sources s(m, false);
    std::vector<Real> pos = xyz(0.4 * m.Lx, 0.55 * m.Ly, 0.5 * m.H), force = xyz(-2.0e6, 0.0, 0.0);
    pos.insert(pos.end(), {Real(0.5 * m.Lx), Real(0.5 * m.Ly), Real(0.5 * m.H)});
    force.insert(force.end(), {Real(1.0e5), Real(0.0), Real(0.0)});
    std::string msg = erf_gtest::abort_message([&] { s.spread(pos, force, std::vector<Real>{Real(2.0 * m.dx())}); });
    EXPECT_NE(msg.find("one kernel width per point"), std::string::npos) << msg;
    for (const Real bad : {Real(0.0), Real(-1.0), std::numeric_limits<Real>::quiet_NaN()}) {
        msg = erf_gtest::abort_message([&] { s.spread(pos, force, std::vector<Real>{Real(2.0 * m.dx()), bad}); });
        EXPECT_NE(msg.find("point 1 (0-based) has the kernel width"), std::string::npos) << msg;
    }
    // one width for every point: refused when it is not positive, without points too
    for (const Real bad : {Real(0.0), Real(-1.0)}) {
        msg = erf_gtest::abort_message([&] { s.spread({}, {}, bad); });
        EXPECT_NE(msg.find("the kernel width"), std::string::npos) << msg;
    }
}
