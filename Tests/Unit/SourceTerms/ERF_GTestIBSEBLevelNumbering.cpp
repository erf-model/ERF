#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBParams.H"

// A refined level numbers its buildings from its own cells; inputs that go by
// building (erf.ibseb.material_by_building) and the report's building_level0
// use level 0's numbers, which IBFaceSet::map_buildings_to_level0() supplies.
// The blankings here are set by hand (1 solid, 0 fluid), no embedded boundary.

namespace {

struct Cells { int ilo, ihi, jlo, jhi, khi; };

// One level's blanking on a single box over a 160 m x 160 m x 40 m domain,
// periodic in x and y, with one ghost layer filled.
amrex::MultiFab blanking (const amrex::Geometry& geom, const std::vector<Cells>& solid)
{
    const amrex::BoxArray ba(geom.Domain());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab b(ba, dm, 1, 1);
    b.setVal(0.0);
    for (amrex::MFIter mfi(b); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& a = b.array(mfi);
        for (const Cells& s : solid) {
            const amrex::Box sb(amrex::IntVect(s.ilo, s.jlo, 0), amrex::IntVect(s.ihi, s.jhi, s.khi));
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

} // namespace

// Level 0 (20 m cells) has a block A (40 m x 40 m, 30 m) and a block B
// (20 m x 20 m, 20 m). Level 1 (10 m cells, ratio 2 x 2 x 1) has the same two
// blocks, A with a 40 m spire one 10 m column wide on its east side, alone in
// a 20 m cell that level 0 leaves fluid, and a 10 m block C only level 1
// resolves. A building's tallest column is the spire, so looking a building up
// by that column alone maps A to nothing; the vote of its columns maps it to
// level 0's A. B maps to level 0's B, and C, which no labelled coarse column
// lies under, to 0.
TEST(IBSEBLevelNumbering, ColumnsVoteForTheLevelZeroBuilding)
{
    IBSEBParams params;
    const amrex::Geometry g0 = level_geom(8, 8);
    const amrex::Geometry g1 = level_geom(16, 16);
    const amrex::MultiFab b0 = blanking(g0, {{2, 3, 2, 3, 2}, {5, 5, 5, 5, 1}});
    const amrex::MultiFab b1 = blanking(g1, {{4, 7, 4, 7, 2}, {8, 8, 4, 4, 3},     // A and its spire
                                             {10, 11, 10, 11, 1},                  // B
                                             {14, 14, 0, 0, 0}});                  // C, level 1 only
    IBFaceSet f0(params, 0), f1(params, 1);
    f0.build(b0, g0);
    f1.build(b1, g1);
    ASSERT_EQ(f0.n_buildings(), 2);
    ASSERT_EQ(f1.n_buildings(), 3);
    f1.map_buildings_to_level0(f0, amrex::IntVect(2, 2, 1));
    // Scan order (i outer): level 1 numbers A 1, B 2, C 3; level 0 A 1, B 2.
    EXPECT_EQ(f1.building_level0(1), 1);
    EXPECT_EQ(f1.building_level0(2), 2);
    EXPECT_EQ(f1.building_level0(3), 0);
    // Level 0 is its own numbering.
    EXPECT_EQ(f0.building_level0(1), 1);
    EXPECT_EQ(f0.building_level0(2), 2);
}

// The numbering follows the buildings, not the scan order of each level: with
// a level 1 that holds only B, level 1's building 1 is level 0's building 2.
TEST(IBSEBLevelNumbering, ALevelHoldingSomeBuildingsKeepsLevelZeroNumbers)
{
    IBSEBParams params;
    const amrex::Geometry g0 = level_geom(8, 8);
    const amrex::Geometry g1 = level_geom(16, 16);
    const amrex::MultiFab b0 = blanking(g0, {{2, 3, 2, 3, 2}, {5, 5, 5, 5, 1}});
    const amrex::MultiFab b1 = blanking(g1, {{10, 11, 10, 11, 1}});
    IBFaceSet f0(params, 0), f1(params, 1);
    f0.build(b0, g0);
    f1.build(b1, g1);
    ASSERT_EQ(f1.n_buildings(), 1);
    f1.map_buildings_to_level0(f0, amrex::IntVect(2, 2, 1));
    EXPECT_EQ(f1.building_level0(1), 2);
}
