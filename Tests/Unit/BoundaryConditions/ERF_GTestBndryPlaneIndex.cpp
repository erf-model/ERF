// Contract of bndry_plane_index, which picks the cell of a boundary plane (read from a precursor's
// boundary registers) that fills each ghost cell beside the face: inside the domain the cell
// itself; along a periodic direction the periodic image, as across the seam inside the domain, so
// that the ghost cells beside an inflow face at y = -1 hold the plane at y = ny - 1, not at y = 0;
// along a non-periodic direction the nearest cell of the plane. Also plane_index, which places a
// fine ghost cell or node between two plane cells for refined levels, with the same periodic image
// and clamp: across a periodic seam it blends the cells either side of the seam.

#include <gtest/gtest.h>

#include "ERF_BndryPlaneIndex.H"

TEST(BndryPlaneIndex, InsideTheDomainACellIsItself)
{
    for (int i = 0; i < 96; ++i) {
        EXPECT_EQ(bndry_plane_index(i, 0, 95, true), i);
        EXPECT_EQ(bndry_plane_index(i, 0, 95, false), i);
    }
}

TEST(BndryPlaneIndex, APeriodicDirectionTakesThePeriodicImage)
{
    // the ghost cells either side of a periodic seam, as many as a high-order stencil reaches
    EXPECT_EQ(bndry_plane_index(-1, 0, 95, true), 95);
    EXPECT_EQ(bndry_plane_index(-3, 0, 95, true), 93);
    EXPECT_EQ(bndry_plane_index(96, 0, 95, true), 0);
    EXPECT_EQ(bndry_plane_index(98, 0, 95, true), 2);
    // a domain that does not start at zero, and a ghost further out than the domain is wide
    EXPECT_EQ(bndry_plane_index(3, 4, 7, true), 7);
    EXPECT_EQ(bndry_plane_index(-6, 4, 7, true), 6);
    EXPECT_EQ(bndry_plane_index(12, 4, 7, true), 4);
}

TEST(BndryPlaneIndex, ANonPeriodicDirectionTakesTheNearestCell)
{
    EXPECT_EQ(bndry_plane_index(-1, 0, 47, false), 0);
    EXPECT_EQ(bndry_plane_index(-3, 0, 47, false), 0);
    EXPECT_EQ(bndry_plane_index(48, 0, 47, false), 47);
    EXPECT_EQ(bndry_plane_index(50, 0, 47, false), 47);
}

// plane_index places a fine index among the level-0 plane indices (a plane of 4 cells, 0..3, at
// ratio 2 below): the two neighbours to blend and the weight on the second, with the neighbours
// wrapped across a periodic seam and clamped otherwise
namespace {
void expect_plane_index (const PlaneIndex& p, int lo, int hi, double w)
{
    EXPECT_EQ(p.lo, lo);
    EXPECT_EQ(p.hi, hi);
    EXPECT_NEAR(static_cast<double>(p.w), w, 1.0e-6);
}
} // namespace

TEST(BndryPlaneIndex, AtLevelZeroAPlaneIndexIsTheCellItselfWithNoWeight)
{
    for (int j = -3; j < 7; ++j) {
        for (bool nodal : {false, true}) {
            for (bool periodic : {false, true}) {
                const PlaneIndex p = plane_index(j, 1, nodal, 0, 3, periodic);
                EXPECT_EQ(p.lo, bndry_plane_index(j, 0, 3, periodic)) << j;
                EXPECT_EQ(p.w, 0.0) << j;
            }
        }
    }
}

TEST(BndryPlaneIndex, AFineCellAcrossAPeriodicSeamBlendsTheCellsEitherSideOfIt)
{
    // fine cell 0 sits a quarter of a coarse cell inside the domain: between coarse cell 3
    // (its periodic image at -1) and coarse cell 0
    expect_plane_index(plane_index(0, 2, false, 0, 3, true), 3, 0, 0.75);
    expect_plane_index(plane_index(7, 2, false, 0, 3, true), 3, 0, 0.25);
    expect_plane_index(plane_index(-1, 2, false, 0, 3, true), 3, 0, 0.25);
    // away from the seam, the same as without the periodic image
    expect_plane_index(plane_index(3, 2, false, 0, 3, true), 1, 2, 0.25);
    expect_plane_index(plane_index(3, 2, false, 0, 3, false), 1, 2, 0.25);
}

TEST(BndryPlaneIndex, AFineNodeAcrossAPeriodicSeamBlendsTheNodesEitherSideOfIt)
{
    // node 8 of the fine level is node 4 of the coarse one, the periodic image of node 0
    expect_plane_index(plane_index(8, 2, true, 0, 3, true), 0, 1, 0.0);
    expect_plane_index(plane_index(7, 2, true, 0, 3, true), 3, 0, 0.5);
    expect_plane_index(plane_index(-1, 2, true, 0, 3, true), 3, 0, 0.5);
}

TEST(BndryPlaneIndex, ANonPeriodicDirectionClampsBothNeighboursToThePlane)
{
    expect_plane_index(plane_index(0, 2, false, 0, 3, false), 0, 0, 0.75);
    expect_plane_index(plane_index(7, 2, false, 0, 3, false), 3, 3, 0.25);
    expect_plane_index(plane_index(8, 2, true, 0, 3, false), 3, 3, 0.0);
    expect_plane_index(plane_index(-1, 2, true, 0, 3, false), 0, 0, 0.5);
}
