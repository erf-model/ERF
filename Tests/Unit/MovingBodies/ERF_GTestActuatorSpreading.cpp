// Contract of erf_actuator::spread_forces: the source spread from a point integrates back to
// the point's force exactly on every component, on a uniform mesh, on a terrain-following mesh
// with its cell volumes, when the kernel is cut off by the ground, across a periodic boundary,
// for several points at once, and however the domain is split; the source is zero beyond the
// kernel's reach; over a hill taller than the kernel's reach the source is centred on the point
// (the faces are found by their physical heights, not their index times dz); each point is spread
// with its own kernel width; and a non-finite position or force, or a width that is not positive,
// aborts, naming the point. The two passes on their own (spread_partial_norms, add_scaled_kernels) give
// spread_forces' source term by term, against a point-by-point reference, with the columns split in z and
// on a rank holding only boxes aloft;
// add_spread_forces adds to what the sources hold.

#include <AMReX_ParallelReduce.H>
#include "ERF_ActuatorGeometry.H"
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

    // dm_in: the distribution of another Sources on the same mesh, to compare the two box by box
    Sources (const Mesh& m, bool terrain, const amrex::DistributionMapping* dm_in = nullptr)
        : mesh(m)
    {
        const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx-1, m.ny-1, m.nz-1));
        const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(m.Lx, m.Ly, m.H)});
        geom = amrex::Geometry(domain, &rb, 0, m.periodic.data());
        ba = amrex::BoxArray(domain);
        ba.maxSize(amrex::IntVect(m.max_grid[0], m.max_grid[1], m.max_grid[2]));
        dm = (dm_in != nullptr) ? *dm_in : amrex::DistributionMapping(ba);
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

namespace {

// A point's periodic images, in the spreader's order
std::vector<std::array<Real,3>> reference_images (const amrex::Geometry& geom, Real px, Real py, Real pz, Real eps)
{
    int n[2] = {0, 0};
    Real len[2] = {0.0, 0.0};
    for (int d = 0; d < 2; ++d) {
        if (geom.isPeriodic(d)) {
            len[d] = geom.ProbLength(d);
            n[d] = static_cast<int>(std::ceil(Real(3.0) * eps / len[d]));
        }
    }
    std::vector<std::array<Real,3>> q;
    for (int sx = -n[0]; sx <= n[0]; ++sx) {
        for (int sy = -n[1]; sy <= n[1]; ++sy) { q.push_back({{px + sx * len[0], py + sy * len[1], pz}}); }
    }
    return q;
}

// The kernel weight of face (i,j,k) of grid dir for one image of a point, as the spreader defines it
Real reference_weight (int dir, int i, int j, int k, const std::array<Real,3>& q, const amrex::Geometry& geom,
                       amrex::Array4<Real const> const& znd, bool has_znd, Real eps)
{
    const auto plo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const Real inv_eps2 = Real(1.0) / (eps * eps);
    const Real reach2 = Real(9.0) * eps * eps;
    const Real ddx = plo[0] + (i + ((dir == 0) ? Real(0.0) : Real(0.5))) * dx[0] - q[0];
    const Real ddy = plo[1] + (j + ((dir == 1) ? Real(0.0) : Real(0.5))) * dx[1] - q[1];
    const Real ddz = erf_actuator::face_height(dir, i, j, k, znd, has_znd, plo[2], dx[2]) - q[2];
    const Real r2 = ddx*ddx + ddy*ddy + ddz*ddz;
    return (r2 > reach2) ? Real(0.0) : std::exp(-r2 * inv_eps2);
}

bool reference_boundary_face (int dir, int i, int j, int k, const amrex::Geometry& geom)
{
    if (geom.isPeriodic(dir)) { return false; }
    const int f = (dir == 0) ? i : ((dir == 1) ? j : k);
    return (f == geom.Domain().smallEnd(dir)) || (f == geom.Domain().bigEnd(dir) + 1);
}

// spread_forces point by point, as before the kernels were batched: every face of every box visited for
// every image of every point (not only those within reach), one point after the other, and the terms
// added in that order (box, point, image, then the faces k, j, i). Its sums take the batched kernels'
// order, so on the host the two agree to the last bit.
void reference_spread (const std::vector<Real>& pos, const std::vector<Real>& force, const std::vector<Real>& eps,
                       Sources& s)
{
    const int npts = static_cast<int>(eps.size());
    const bool has_znd = (s.znd != nullptr);
    const auto dx = s.geom.CellSizeArray();
    const Real dxdydz = dx[0] * dx[1] * dx[2];
    amrex::MultiFab* src[3] = {&s.sx, &s.sy, &s.sz};
    std::vector<Real> norm(3 * npts, 0.0);
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) {
        const amrex::Box vbx = amrex::enclosedCells(mfi.validbox());
        amrex::Array4<Real const> znd = has_znd ? s.znd->const_array(mfi) : amrex::Array4<Real const>{};
        amrex::Array4<Real const> detj = has_znd ? s.detj->const_array(mfi) : amrex::Array4<Real const>{};
        for (int dir = 0; dir < 3; ++dir) {
            amrex::Box fbx = amrex::convert(vbx, amrex::IntVect::TheDimensionVector(dir));
            fbx.growHi(dir, -1);
            for (int p = 0; p < npts; ++p) {
                for (const auto& q : reference_images(s.geom, pos[3*p], pos[3*p+1], pos[3*p+2], eps[p])) {
                    amrex::LoopOnCpu(fbx, [&](int i, int j, int k) {
                        if (reference_boundary_face(dir, i, j, k, s.geom)) { return; }
                        const Real w = reference_weight(dir, i, j, k, q, s.geom, znd, has_znd, eps[p]);
                        if (w > Real(0.0)) { const Real term = w * erf_actuator::face_volume(dir, i, j, k, detj, has_znd, dxdydz); norm[3*p+dir] += term; }
                    });
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(norm.data(), 3 * npts);
    for (int dir = 0; dir < 3; ++dir) { src[dir]->setVal(0.0); }
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) {
        amrex::Array4<Real const> znd = has_znd ? s.znd->const_array(mfi) : amrex::Array4<Real const>{};
        for (int dir = 0; dir < 3; ++dir) {
            const amrex::Box fbx = mfi.validbox().convert(amrex::IntVect::TheDimensionVector(dir));
            const auto a = src[dir]->array(mfi);
            for (int p = 0; p < npts; ++p) {
                const Real scale = force[3*p+dir] / norm[3*p+dir];
                if (scale == Real(0.0)) { continue; }
                for (const auto& q : reference_images(s.geom, pos[3*p], pos[3*p+1], pos[3*p+2], eps[p])) {
                    amrex::LoopOnCpu(fbx, [&](int i, int j, int k) {
                        if (reference_boundary_face(dir, i, j, k, s.geom)) { return; }
                        const Real w = reference_weight(dir, i, j, k, q, s.geom, znd, has_znd, eps[p]);
                        if (w > Real(0.0)) { const Real term = scale * w; a(i,j,k) += term; }
                    });
                }
            }
        }
    }
}

// Forty points of two widths over terrain, some near the periodic corner (their kernels wrap in x and
// y at once), some on box faces, some with zero force components, on boxes of 8 x 5 cells
std::vector<Real> many_points (const Mesh& m, std::vector<Real>& force, std::vector<Real>& eps)
{
    std::vector<Real> pos;
    force.clear();
    eps.clear();
    for (int p = 0; p < 40; ++p) {
        const Real fx = static_cast<Real>(p % 7) / 7, fy = static_cast<Real>((3 * p) % 11) / 11, fz = static_cast<Real>((5 * p) % 13) / 13;
        if (p < 6) {
            pos.insert(pos.end(), {Real((0.98 + 0.003 * p) * m.Lx), Real((0.01 + 0.002 * p) * m.Ly), Real((0.3 + 0.05 * p) * m.H)});
        } else if (p < 10) {
            pos.insert(pos.end(), {Real(8.0 * m.dx()), Real((5.0 + p - 6) * m.dy()), Real(0.45 * m.H)});
        } else {
            pos.insert(pos.end(), {Real((0.1 + 0.8 * fx) * m.Lx), Real((0.05 + 0.9 * fy) * m.Ly), Real((0.25 + 0.5 * fz) * m.H)});
        }
        force.insert(force.end(), {Real(-1.0e5 * (1 + p % 5)), (p % 4 == 0) ? Real(0.0) : Real(2.0e4 * (fy - 0.5)),
                                   (p % 3 == 0) ? Real(0.0) : Real(-1.0e4 * fz)});
        eps.push_back((p % 2 == 0) ? Real(1.5 * m.dx()) : Real(2.5 * m.dx()));
    }
    return pos;
}

} // namespace

// Columns split in z (boxes stacked over the same columns, some aloft): the spreading bounds each kernel's
// cells from the table of this rank's node heights, where without it the whole column was swept. Only
// faces of exactly zero weight are skipped, so the norms and every face of the source must be the
// full-column answer bit for bit, on the same layout.
TEST(ActuatorSpreading, AColumnSplitInZGivesTheFullColumnsAnswerBitForBit)
{
    Mesh m;
    m.hill = 60.0;
    m.nz = 40;
    m.H = 2000.0;                 // 50 m cells: a kernel reaches a few hundred metres of the column
    m.max_grid = {{24, 12, 10}};  // four boxes up each column, three of them aloft
    const Real eps = 1.5 * m.dx();
    // near the ground, mid-height and near the top, each with a kernel a few cells tall
    std::vector<Real> pos, force, epsv;
    for (const Real zf : {Real(0.18), Real(0.5), Real(0.8)}) {
        const auto p = xyz(9.3 * m.dx(), 5.6 * m.dy(), zf * m.H);
        pos.insert(pos.end(), p.begin(), p.end());
        const auto f = xyz(-1.0e6 * zf, 2.0e5, -3.0e4);
        force.insert(force.end(), f.begin(), f.end());
        epsv.push_back(eps);
    }
    Sources table(m, true);
    Sources column(m, true, &table.dm);
    // each rank's bounds: the table where one of its boxes is aloft, the ground-and-spacing bounds where all of
    // its boxes reach the ground (on several ranks a rank may hold only the bottom boxes)
    const erf_actuator::ZBounds zb = erf_actuator::mesh_z_bounds(*table.znd, table.geom, false);
    EXPECT_EQ(zb.has_table(), !zb.all_ground);
    int any_table = zb.has_table() ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntMax(any_table);
    EXPECT_EQ(any_table, 1);
    // the table bounds each kernel to part of the column (without it: all of it), so faces are skipped
    if (zb.has_table()) {
        for (std::size_t p = 0; p < epsv.size(); ++p) {
            int kb = 0, kt = 0;
            zb.k_range(pos[3*p+2] - Real(3.0) * eps, pos[3*p+2] + Real(3.0) * eps, 0, m.nz - 1, kb, kt);
            if (kb <= kt) { EXPECT_LT(kt - kb + 1, 3 * m.nz / 4) << "point " << p; }
        }
    }
    const erf_actuator::ZBounds none{};
    auto run = [&] (Sources& s, const erf_actuator::ZBounds* bounds, std::vector<Real>& norm) {
        erf_actuator::spread_partial_norms(pos, epsv, s.znd.get(), s.detj.get(), s.geom, s.sx, norm, bounds);
        std::vector<Real> f_over_s(force.size());
        for (std::size_t q = 0; q < force.size(); ++q) { f_over_s[q] = force[q] / norm[q]; }
        s.sx.setVal(0.0); s.sy.setVal(0.0); s.sz.setVal(0.0);
        erf_actuator::add_scaled_kernels(pos, f_over_s, epsv, s.znd.get(), s.detj.get(), s.geom, s.sx, s.sy, s.sz, bounds);
    };
    std::vector<Real> n_table, n_column;
    run(table, &zb, n_table);
    run(column, &none, n_column);
    ASSERT_EQ(n_table.size(), n_column.size());
    for (std::size_t q = 0; q < n_table.size(); ++q) { EXPECT_EQ(n_table[q], n_column[q]) << "norm " << q; }
    long nonzero = 0, faces = 0;
    for (int dir = 0; dir < 3; ++dir) {
        const amrex::MultiFab& a = (dir == 0) ? table.sx : (dir == 1) ? table.sy : table.sz;
        const amrex::MultiFab& b = (dir == 0) ? column.sx : (dir == 1) ? column.sy : column.sz;
        for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
            const auto x = a.const_array(mfi);
            const auto y = b.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                ++faces;
                if (x(i,j,k) != Real(0.0)) { ++nonzero; }
                if (x(i,j,k) != y(i,j,k)) { ADD_FAILURE() << "dir " << dir << " face (" << i << "," << j << "," << k << "): "
                                                          << x(i,j,k) << " vs " << y(i,j,k); }
            });
        }
    }
    EXPECT_GT(nonzero, 0);
    EXPECT_LT(nonzero, faces);
}

// A rank that holds only boxes aloft (the upper half of every column here): the table has none of the lower
// cells, so a kernel near the ground reaches none of its faces (an empty range, a zero partial norm) and one
// aloft gives the partial norm of the full layout's upper cells, as without the table
TEST(ActuatorSpreading, ARankHoldingOnlyBoxesAloftSkipsTheKernelsBelowThem)
{
    Mesh m;
    m.hill = 60.0;
    m.nz = 40;
    m.H = 2000.0;
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(m.nx-1, m.ny-1, m.nz-1));
    const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(m.Lx, m.Ly, m.H)});
    const amrex::Geometry geom(domain, &rb, 0, m.periodic.data());
    const amrex::BoxArray ba(amrex::Box(amrex::IntVect(0, 0, m.nz / 2), amrex::IntVect(m.nx-1, m.ny-1, m.nz-1)));
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab znd(amrex::convert(ba, amrex::IntVect(1,1,1)), dm, 1, 1);
    amrex::MultiFab detj(ba, dm, 1, 1);
    amrex::MultiFab sx(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 0);
    const Mesh::Spacing sp = m.spacing();
    for (amrex::MFIter mfi(detj, false); mfi.isValid(); ++mfi) {
        auto za = znd.array(mfi);
        auto ja = detj.array(mfi);
        const amrex::Box gbx = amrex::grow(mfi.validbox(), 1);
        amrex::LoopOnCpu(amrex::surroundingNodes(gbx), [&](int i, int j, int k) { za(i,j,k) = m.z_node(sp, i, j, k); });
        amrex::LoopOnCpu(gbx, [&](int i, int j, int k) { ja(i,j,k) = m.detj(sp, i, j); (void) k; });
    }
    // the rank that holds the box (the others hold none: no table, and nothing to spread); every call
    // below is this rank's alone (no collective)
    if (dm[0] != amrex::ParallelDescriptor::MyProc()) { return; }
    const erf_actuator::ZBounds zb = erf_actuator::mesh_z_bounds(znd, geom, false);
    ASSERT_TRUE(zb.has_table());
    for (int k = 0; k < m.nz; ++k) { EXPECT_EQ(zb.has_cell[k] != 0, k >= m.nz / 2) << "cell " << k; }
    const Real eps = 1.5 * m.dx();
    const Real x = 9.3 * m.dx(), y = 5.6 * m.dy();
    int kb = 0, kt = 0;
    zb.k_range(Real(150.0) - Real(3.0) * eps, Real(150.0) + Real(3.0) * eps, 0, m.nz - 1, kb, kt);
    EXPECT_GT(kb, kt);
    // the partial norms on this rank: none for the low point; for the high one, as without the table
    const std::vector<Real> pos{x, y, Real(150.0), x, y, Real(1500.0)};
    const std::vector<Real> epsv{eps, eps};
    std::vector<Real> with, without;
    erf_actuator::spread_partial_norms(pos, epsv, &znd, &detj, geom, sx, with, &zb);
    const erf_actuator::ZBounds none{};
    erf_actuator::spread_partial_norms(pos, epsv, &znd, &detj, geom, sx, without, &none);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(with[d], Real(0.0)) << "low point, dir " << d;
        EXPECT_EQ(without[d], Real(0.0)) << "low point, dir " << d;
        EXPECT_GT(with[3+d], Real(0.0)) << "high point, dir " << d;
        EXPECT_EQ(with[3+d], without[3+d]) << "high point, dir " << d;
    }
}

// The batched kernels (one per box and face grid in each pass, no atomics) give the sources of the
// point-by-point reference: to the last bit on the host, where the reference adds its terms in the
// batched kernels' order, and integrate back to the forces. With buffers kept from call to call, each
// pass launches at most one kernel per box and face grid, however many points and images, and a second
// call with fewer points after the first (the buffers keep their size) is unaffected by the first's.
TEST(ActuatorSpreading, BatchedKernelsMatchThePointByPointReference)
{
    Mesh m;
    m.hill = 60.0;
    m.max_grid = {{8, 5, 1024}};
    std::vector<Real> force, eps;
    const std::vector<Real> pos = many_points(m, force, eps);
    Sources s(m, true);
    Sources ref(m, true, &s.dm);
    int nlocal = 0;
    for (amrex::MFIter mfi(s.sx, false); mfi.isValid(); ++mfi) { ++nlocal; }

    erf_actuator::ActuatorDeviceBuffers buf;
    for (const int n : {40, 9, 40}) {
        SCOPED_TRACE(std::to_string(n) + " points");
        const std::vector<Real> p_pos(pos.begin(), pos.begin() + 3 * n), p_force(force.begin(), force.begin() + 3 * n);
        const std::vector<Real> p_eps(eps.begin(), eps.begin() + n);
        // the two passes of add_spread_forces, through the kept buffers
        s.sx.setVal(0.0);
        s.sy.setVal(0.0);
        s.sz.setVal(0.0);
        std::vector<Real> norm;
        long before = buf.kernel_launches;
        erf_actuator::spread_partial_norms(p_pos, p_eps, s.znd.get(), s.detj.get(), s.geom, s.sx, norm, nullptr, &buf);
        EXPECT_LE(buf.kernel_launches - before, 3 * nlocal);
        amrex::ParallelDescriptor::ReduceRealSum(norm.data(), static_cast<int>(norm.size()));
        std::vector<Real> f_over_s(norm.size());
        for (std::size_t q = 0; q < norm.size(); ++q) { f_over_s[q] = p_force[q] / norm[q]; }
        before = buf.kernel_launches;
        // the face-by-face pass, whose launches this counts (on the host the images scatter without kernels)
        buf.force_face_major = true;
        erf_actuator::add_scaled_kernels(p_pos, f_over_s, p_eps, s.znd.get(), s.detj.get(), s.geom, s.sx, s.sy, s.sz, nullptr, &buf);
        buf.force_face_major = false;
        EXPECT_GT(buf.kernel_launches - before, 0);
        EXPECT_LE(buf.kernel_launches - before, 3 * nlocal);
        EXPECT_GE(buf.image_pos.device.size(), 3 * p_eps.size());

        reference_spread(p_pos, p_force, p_eps, ref);
        std::array<Real,3> p_total{{0.0, 0.0, 0.0}};
        for (std::size_t q = 0; q < p_force.size(); ++q) { p_total[q % 3] += p_force[q]; }
        expect_total(s.totals(), p_total, abs_sum(p_force));
        expect_total(ref.totals(), p_total, abs_sum(p_force));
        int differ = 0, nonzero = 0;
        Real worst = 0.0;
        for (int dir = 0; dir < 3; ++dir) {
            const amrex::MultiFab& got = (dir == 0) ? s.sx : ((dir == 1) ? s.sy : s.sz);
            const amrex::MultiFab& want = (dir == 0) ? ref.sx : ((dir == 1) ? ref.sy : ref.sz);
            const Real scale = want.norminf(0);
            for (amrex::MFIter mfi(got, false); mfi.isValid(); ++mfi) {
                const auto a = got.const_array(mfi);
                const auto b = want.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    if (b(i,j,k) != Real(0.0)) { ++nonzero; }
                    if (a(i,j,k) != b(i,j,k)) {
                        ++differ;
                        worst = std::max(worst, std::abs(a(i,j,k) - b(i,j,k)) / scale);
                    }
                });
            }
        }
        amrex::ParallelDescriptor::ReduceIntSum(nonzero);
        amrex::ParallelDescriptor::ReduceIntSum(differ);
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        EXPECT_GT(nonzero, 1000);
#if !defined(AMREX_USE_GPU) && (!defined(__FP_FAST_FMA) || (defined(__clang__) && !defined(__INTEL_LLVM_COMPILER)))
        // The same operations in the same order: bit for bit on any host build whose compiler keeps a multiply
        // and an add of two statements apart. That is every target without a fused multiply-add, and clang on
        // any target (it fuses within one expression only).
        EXPECT_EQ(differ, 0) << "worst relative difference " << worst;
#else
        // A device, or a compiler that fuses a multiply and an add across statements where the target has the
        // instruction (GCC and the Intel compiler by default), may round otherwise.
        EXPECT_LE(worst, rel_tol);
#endif
    }
}

// The second pass on one host thread adds each image into the faces of its own reach box; on a device (or with
// several threads) one thread per face adds the images binned to its tile. Both add every face's terms in image
// order: every face of every box must be bit-identical, over terrain with several boxes, several points per face
// and periodic images.
TEST(ActuatorSpreading, TheHostScatterIsTheFaceByFacePassBitForBit)
{
#ifdef AMREX_USE_GPU
    GTEST_SKIP() << "a device has only the face-by-face pass";
#else
    // over terrain and on a flat mesh, several boxes, several points per face, periodic images in x and y, zero
    // scales, kernels reaching the ground and the top, and sources that already hold values (several bodies add
    // into the same sources)
    for (const Real hill : {Real(60.0), Real(0.0)}) {
        SCOPED_TRACE("hill " + std::to_string(hill));
        Mesh m;
        m.hill = hill;
        m.max_grid = {{8, 6, 1024}};
        std::vector<Real> pos, scale, epsv;
        for (int p = 0; p < 12; ++p) {
            const auto q = xyz((0.7 + 1.9 * p) * m.dx(), (0.4 + 0.83 * p) * m.dy(), (0.2 + 0.05 * p) * m.H);
            pos.insert(pos.end(), q.begin(), q.end());
            const auto f = xyz(1.0e3 * (p + 1), -2.0e2 * p, (p % 3 == 0) ? 0.0 : 5.0e1);   // a zero scale skips a component
            scale.insert(scale.end(), f.begin(), f.end());
            epsv.push_back((1.2 + 0.1 * (p % 4)) * m.dx());
        }
        Sources scatter(m, true);
        Sources by_face(m, true, &scatter.dm);
        for (Sources* s : {&scatter, &by_face}) {
            // the same non-zero start in both
            for (int dir = 0; dir < 3; ++dir) {
                amrex::MultiFab& f = (dir == 0) ? s->sx : (dir == 1) ? s->sy : s->sz;
                for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
                    auto a = f.array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i,j,k) = Real(0.37) * (i - 2 * j + 3 * k + dir) / Real(7.0); });
                }
            }
            erf_actuator::ActuatorDeviceBuffers buf;
            buf.force_face_major = (s == &by_face);
            erf_actuator::add_scaled_kernels(pos, scale, epsv, s->znd.get(), s->detj.get(), s->geom, s->sx, s->sy, s->sz, nullptr, &buf);
            // the face-by-face pass launches a kernel per box and face grid it reaches; the host scatter none
            long launches = buf.kernel_launches;
            amrex::ParallelDescriptor::ReduceLongSum(launches);
            if (s == &scatter) { EXPECT_EQ(launches, 0); } else { EXPECT_GT(launches, 0); }
        }
        long changed = 0;
        for (int dir = 0; dir < 3; ++dir) {
            const amrex::MultiFab& a = (dir == 0) ? scatter.sx : (dir == 1) ? scatter.sy : scatter.sz;
            const amrex::MultiFab& b = (dir == 0) ? by_face.sx : (dir == 1) ? by_face.sy : by_face.sz;
            for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
                const auto x = a.const_array(mfi);
                const auto y = b.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    if (x(i,j,k) != Real(0.37) * (i - 2 * j + 3 * k + dir) / Real(7.0)) { ++changed; }
                    if (x(i,j,k) != y(i,j,k)) { ADD_FAILURE() << "dir " << dir << " face (" << i << "," << j << "," << k << "): "
                                                              << x(i,j,k) << " vs " << y(i,j,k); }
                });
            }
        }
        amrex::ParallelDescriptor::ReduceLongSum(changed);
        EXPECT_GT(changed, 1000);
    }
#endif
}


// add_spread_forces adds to what the sources hold: two sets spread one after the other give the sum of the
// two sets spread alone, and the integral of the sum is the sum of the forces
TEST(ActuatorSpreading, AddingASecondSetLeavesTheFirstInPlace)
{
    Mesh m;
    m.hill = 80.0;
    m.max_grid = {{8, 6, 1024}};
    Sources a(m, true), b(m, true, &a.dm), both(m, true, &a.dm);
    const std::vector<Real> p1 = xyz(0.3 * m.Lx, 0.4 * m.Ly, 0.5 * m.H), f1 = xyz(-2.0e6, 3.0e5, 1.0e5);
    const std::vector<Real> p2 = xyz(0.36 * m.Lx, 0.45 * m.Ly, 0.55 * m.H), f2 = xyz(4.0e5, -1.0e5, 2.0e5);
    const std::vector<Real> e1{Real(2.0 * m.dx())}, e2{Real(1.5 * m.dx())};
    erf_actuator::spread_forces(p1, f1, e1, a.znd.get(), a.detj.get(), a.geom, a.sx, a.sy, a.sz);
    erf_actuator::spread_forces(p2, f2, e2, b.znd.get(), b.detj.get(), b.geom, b.sx, b.sy, b.sz);
    erf_actuator::spread_forces(p1, f1, e1, both.znd.get(), both.detj.get(), both.geom, both.sx, both.sy, both.sz);
    erf_actuator::add_spread_forces(p2, f2, e2, both.znd.get(), both.detj.get(), both.geom, both.sx, both.sy, both.sz);
    amrex::MultiFab* sa[3] = {&a.sx, &a.sy, &a.sz};
    amrex::MultiFab* sb[3] = {&b.sx, &b.sy, &b.sz};
    amrex::MultiFab* sc[3] = {&both.sx, &both.sy, &both.sz};
    for (int d = 0; d < 3; ++d) {
        ASSERT_GT(sb[d]->norm0(), Real(0.0));
        // the second set alone is what was added: the same additions onto a zero and onto the first set
        amrex::MultiFab::Subtract(*sc[d], *sa[d], 0, 0, 1, 0);
        amrex::MultiFab::Subtract(*sc[d], *sb[d], 0, 0, 1, 0);
        // to rounding: a face of both kernels takes (first + second) - first - second
        EXPECT_LE(sc[d]->norm0(), Real(16.0) * std::numeric_limits<Real>::epsilon() * (sa[d]->norm0() + sb[d]->norm0()))
            << "component " << d;
    }
}
