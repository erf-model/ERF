#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBParams.H"
#include "ERF_IBSEBSolar.H"

// The ray cast of the building faces reads a per-rank patch of the column
// map instead of a map of the whole domain replicated on every rank
// (erf-model/ERF#4020): the halo of the patch is the reach of the shallowest
// ray the map must follow, a ray that leaves the patch is in open country,
// and a patch that spans a periodic domain wraps round. The building ids
// come from a label exchange across the boxes and the periodic seams, in
// the same scan-order numbering whatever the decomposition. The blankings
// here are set by hand (1 solid, 0 fluid), no embedded boundary.

namespace {

// A ray map over a hand-made patch of 10 m columns: tops[] is indexed
// (i - lo_i) + nx_patch * (j - lo_j) over the patch box, as an Array4 is.
struct PatchMap
{
    std::vector<int> tops;
    ibseb::RayMap m;
    PatchMap (const amrex::Box& patch, int nx, int ny, bool wrap_x, amrex::Real z_max)
        : tops(static_cast<size_t>(patch.numPts()), -1)
    {
        m.top = amrex::Array4<const int>(tops.data(), amrex::begin(patch), amrex::end(patch), 1);
        m.ilo = 0; m.jlo = 0; m.nx = nx; m.ny = ny;
        m.wrap_x = wrap_x; m.wrap_y = false;
        m.x_lo = 0.0; m.y_lo = 0.0; m.dx = 10.0; m.dy = 10.0; m.dz = 10.0;
        m.z_ground = 0.0; m.z_max = z_max;
        m.max_steps = 1000;
    }
    void set (int i, int j, int kt)
    {
        const amrex::Box b(amrex::IntVect(m.top.begin[0], m.top.begin[1], 0), amrex::IntVect(m.top.end[0] - 1, m.top.end[1] - 1, 0));
        tops[static_cast<size_t>(b.index(amrex::IntVect(i, j, 0)))] = kt;
    }
};

// One level's blanking on the given boxes over a 160 m x 160 m x 40 m
// domain of nx x ny columns, periodic in x and y, one ghost layer filled.
amrex::MultiFab blanking (const amrex::Geometry& geom, const amrex::BoxArray& ba,
                          const std::vector<std::array<int, 5>>& solid)    // ilo, ihi, jlo, jhi, khi
{
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab b(ba, dm, 1, 1);
    b.setVal(0.0);
    for (amrex::MFIter mfi(b); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& a = b.array(mfi);
        for (const auto& s : solid) {
            const amrex::Box sb(amrex::IntVect(s[0], s[2], 0), amrex::IntVect(s[1], s[3], s[4]));
            const amrex::Box ov = sb & bx;
            if (ov.isEmpty()) { continue; }
            amrex::ParallelFor(ov, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = 1.0; });
        }
    }
    b.FillBoundary(geom.periodicity());
    return b;
}

amrex::Geometry level_geom (int nx, int ny)
{
    const amrex::RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 40.0});
    const amrex::Box dom(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, 3));
    return amrex::Geometry(dom, rb, 0, {1, 1, 0});
}

template <class T>
std::vector<T> host (const amrex::Gpu::DeviceVector<T>& d)
{
    std::vector<T> h(d.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d.begin(), d.end(), h.begin());
    amrex::Gpu::streamSynchronize();
    return h;
}

} // namespace

// The halo is the horizontal reach of the shallowest ray over the tallest
// column, height / tan(elevation), rounded up plus the starting column; and
// the shallowest ray of the hemisphere sampling is the wall ray nearest the
// horizontal, sin(theta_max) sin(pi / n_az), unless a roof ray is lower.
TEST(IBSEBRayMap, HaloIsTheReachOfTheShallowestRay)
{
    // 140 m over 5 m cells for a sun 30 deg up: 140 / tan(30) = 242.5 m = 48.5 cells,
    // ceil(48.5) + 1 = 50 with the starting column.
    EXPECT_EQ(ibseb::ray_halo_cells(140.0, 5.0, std::sin(30.0 * PI / 180.0)), 50);
    // Nothing built: one cell, for the solid neighbour of a face.
    EXPECT_EQ(ibseb::ray_halo_cells(0.0, 5.0, 0.5), 1);
    // A level ray reaches for ever; the caller caps it at the domain.
    EXPECT_GT(ibseb::ray_halo_cells(10.0, 5.0, 0.0), 1000000);
    // 16 x 8 rays: the wall ray nearest the normal, theta_min = asin(sqrt(1/16)),
    // at the azimuth nearest the horizontal, phi_min = pi/16: 2.8 deg.
    const double expect = std::sqrt(1.0 / 16.0) * std::sin(PI / 16.0);
    EXPECT_NEAR(ibseb::hemisphere_min_sz(16, 8), expect, 1.0e-12);
    // 4 azimuths: the wall ray is 45 deg round, sqrt(1/16) sin(pi/4), still below the
    // roof ray's cos(theta_max) = sqrt(1/16) = 0.25.
    EXPECT_NEAR(ibseb::hemisphere_min_sz(4, 8), 0.25 * std::sin(PI / 4.0), 1.0e-12);
    // 2 x 1 rays: the wall rays go straight up and down the tangent plane, the roof
    // ray 45 deg from the normal: both sqrt(1/2).
    EXPECT_NEAR(ibseb::hemisphere_min_sz(2, 1), std::sqrt(0.5), 1.0e-12);
    // Brute force over the samples of hemisphere_direction() for a few sizes.
    for (int n_az : {2, 3, 4, 8, 16}) {
        for (int n_el : {1, 2, 8}) {
            double lowest = 1.0;
            for (int dir : {0, 1, 2}) {
                for (int ia = 0; ia < n_az; ++ia) {
                    for (int ie = 0; ie < n_el; ++ie) {
                        amrex::Real sx, sy, sz;
                        ibseb::hemisphere_direction(dir, 1, ia, ie, n_az, n_el, sx, sy, sz);
                        if (std::abs(sz) > 1.0e-9) { lowest = std::min(lowest, double(std::abs(sz))); }
                    }
                }
            }
            EXPECT_NEAR(ibseb::hemisphere_min_sz(n_az, n_el), lowest, 1.0e-12) << n_az << " x " << n_el;
        }
    }
}

// A tall column inside the patch blocks a low ray; cut the patch short of
// it and the ray has left the map, where the country is open: not blocked,
// and a falling ray reaches the ground.
TEST(IBSEBRayMap, LeavingThePatchIsOpenGround)
{
    const double el = 10.0 * PI / 180.0;
    const double sx = std::cos(el), sz = std::sin(el);
    // Column (8, 5) is 40 m tall (top index 3); the ray from (15, 55, 5) toward +x
    // enters it at x = 80 m at 5 + 65 tan(10) = 16.5 m: blocked.
    PatchMap whole(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(9, 9, 0)), 10, 10, false, 40.0);
    whole.set(8, 5, 3);
    EXPECT_TRUE(ibseb::ray_blocked(15.0, 55.0, 5.0, sx, 0.0, sz, whole.m));
    EXPECT_EQ(ibseb::ray_hit(15.0, 55.0, 35.0, sx, 0.0, -sz, whole.m), ibseb::RAY_BUILDING);
    EXPECT_EQ(ibseb::ray_hit(15.0, 55.0, 5.0, sx, 0.0, sz, whole.m), ibseb::RAY_BUILDING);
    // The same ray on a patch that ends at column 6: open ground beyond. The
    // patch's storage is solid everywhere but on the ray's own columns, so a
    // walk that read past the patch's bounds would find a building there.
    PatchMap cut(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(6, 9, 0)), 10, 10, false, 40.0);
    std::fill(cut.tops.begin(), cut.tops.end(), 3);
    for (int i = 0; i <= 6; ++i) { cut.set(i, 5, -1); }
    EXPECT_FALSE(ibseb::ray_blocked(15.0, 55.0, 5.0, sx, 0.0, sz, cut.m));
    EXPECT_EQ(ibseb::ray_hit(15.0, 55.0, 35.0, sx, 0.0, -sz, cut.m), ibseb::RAY_GROUND);
    EXPECT_EQ(ibseb::ray_hit(15.0, 55.0, 5.0, sx, 0.0, sz, cut.m), ibseb::RAY_SKY);
    // A ray above the tallest column is never blocked.
    EXPECT_FALSE(ibseb::ray_blocked(15.0, 55.0, 40.0, sx, 0.0, sz, whole.m));
}

// A patch that spans the periodic domain wraps the ray round: from (95, 55, 5)
// in the fluid column 9 toward +x the ray leaves at x = 100 m, re-enters at
// x = 0 and meets column 8 at x = 180 m, 85 m on, at 5 + 85 tan(10) = 20.0 m:
// blocked. Without the wrap the same patch ends at the domain's edge and the
// ray is free.
TEST(IBSEBRayMap, APatchSpanningThePeriodicDomainWrapsRound)
{
    const double el = 10.0 * PI / 180.0;
    const double sx = std::cos(el), sz = std::sin(el);
    PatchMap wrap(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(9, 9, 0)), 10, 10, true, 40.0);
    wrap.set(8, 5, 3);
    EXPECT_TRUE(ibseb::ray_blocked(95.0, 55.0, 5.0, sx, 0.0, sz, wrap.m));
    PatchMap open(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(9, 9, 0)), 10, 10, false, 40.0);
    open.set(8, 5, 3);
    EXPECT_FALSE(ibseb::ray_blocked(95.0, 55.0, 5.0, sx, 0.0, sz, open.m));
    // The wrapped walk still ends: a ray that clears the tallest column is free
    // (85 m on at a slope of 0.8 / 0.98 it is 74 m up).
    EXPECT_FALSE(ibseb::ray_blocked(95.0, 55.0, 5.0, sx, 0.0, 0.8, wrap.m));
}

// Three blocks on a 16 x 16 periodic level of four 8 x 8 boxes: A runs across
// the box edge at i = 8, B across the periodic seam in y (rows 15 and 0), C
// sits alone. The labelling must give three buildings, numbered in scan
// order of their first column (C at (2, 12), A at (6, 2), B at (12, 0): C 1,
// A 2, B 3), with every face carrying its building's id, and the same on one
// box as on four.
TEST(IBSEBRayMap, BuildingsAreLabelledAcrossBoxesAndSeams)
{
    IBSEBParams params;
    const amrex::Geometry g = level_geom(16, 16);
    const std::vector<std::array<int, 5>> blocks = {{6, 9, 2, 3, 1},      // A
                                                    {12, 13, 15, 15, 1},  // B, row 15
                                                    {12, 13, 0, 0, 1},    // B, row 0
                                                    {2, 2, 12, 12, 0}};   // C
    amrex::BoxArray ba_four(g.Domain());
    ba_four.maxSize(8);
    ASSERT_EQ(ba_four.size(), 4);
    const amrex::BoxArray ba_one(g.Domain());
    IBFaceSet f4(params, 0), f1(params, 0);
    f4.build(blanking(g, ba_four, blocks), g);
    f1.build(blanking(g, ba_one, blocks), g);
    EXPECT_EQ(f4.n_buildings(), 3);
    EXPECT_EQ(f1.n_buildings(), 3);
    for (const IBFaceSet* f : {&f4, &f1}) {
        const auto i = host(f->d_i), j = host(f->d_j), bid = host(f->d_bid), dir = host(f->d_dir), side = host(f->d_side);
        int nA = 0, nB = 0, nC = 0;
        for (size_t n = 0; n < i.size(); ++n) {
            // The solid neighbour of the face.
            int si = i[n], sj = j[n];
            if (dir[n] == 0) { si += side[n]; } else if (dir[n] == 1) { sj += side[n]; }
            sj = (sj + 16) % 16;
            int expect = 0;
            if (si >= 6 && si <= 9 && sj >= 2 && sj <= 3) { expect = 2; ++nA; }
            else if (si >= 12 && si <= 13 && (sj == 15 || sj == 0)) { expect = 3; ++nB; }
            else if (si == 2 && sj == 12) { expect = 1; ++nC; }
            EXPECT_EQ(bid[n], expect) << "face " << n << " at (" << i[n] << "," << j[n] << ") dir " << dir[n] << " side " << side[n];
        }
        // A: 2 x 4 columns of 2 cells, roof 8, walls 2 (2 + 4) 2 = 24; B: 2 x 2 of 2, roof 4,
        // walls 16; C: 1 column of 1, roof 1, walls 4.
        EXPECT_EQ(nA, 32);
        EXPECT_EQ(nB, 20);
        EXPECT_EQ(nC, 5);
    }
    EXPECT_EQ(f4.n_faces(), f1.n_faces());
}

// A building outside a refined level is seen through the ray map, from the
// level below: level 0 (20 m cells) holds a 30 m block at x = 120-160 m,
// y = 60-100 m; level 1 (10 m cells) covers x < 80 m with its own 20 m block
// at x = 20-40 m, y = 70-90 m. Its east wall sees only sky and ground until
// the level below is merged into its ray map, then the block too. (The
// level's map wraps, so the wall's shallowest rays also come round to the
// back of its own block, 20 m to the side of their path: a miss by 20 m.)
TEST(IBSEBRayMap, TheLevelBelowOccludesThroughTheRayMap)
{
    IBSEBParams params;
    const amrex::Geometry g0 = level_geom(8, 8);
    const amrex::Geometry g1 = level_geom(16, 16);
    const amrex::BoxArray ba0(g0.Domain());
    const amrex::BoxArray ba1(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 15, 3)));
    IBFaceSet f0(params, 0), f1(params, 1);
    f0.build(blanking(g0, ba0, {{6, 7, 3, 4, 2}}), g0);
    f1.build(blanking(g1, ba1, {{2, 3, 7, 8, 1}}), g1);
    ASSERT_EQ(f0.n_buildings(), 1);
    ASSERT_EQ(f1.n_buildings(), 1);
    auto east_wall_sees_building = [&] () {
        f1.compute_view_fractions();
        const auto i = host(f1.d_i), dir = host(f1.d_dir), side = host(f1.d_side);
        const auto fb = host(f1.d_f_bldg);
        double sum = 0.0; int n_east = 0;
        for (size_t n = 0; n < i.size(); ++n) {
            if (dir[n] == 0 && side[n] == -1 && i[n] == 4) { sum += fb[n]; ++n_east; }
        }
        EXPECT_EQ(n_east, 4);
        return sum;
    };
    EXPECT_EQ(east_wall_sees_building(), 0.0);
    f1.add_outside_occluders({&f0});
    EXPECT_GT(east_wall_sees_building(), 0.0);
    // The ray map of level 1 spans the periodic domain, so the block is in reach.
    EXPECT_TRUE(f1.ray_map_spans_domain());
}
