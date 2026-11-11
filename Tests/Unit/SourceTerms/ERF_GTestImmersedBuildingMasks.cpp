#include <gtest/gtest.h>

#include <AMReX_REAL.H>

#include "ERF_ImmersedBuildingMasks.H"

// Motivation (erf-model/ERF#4887): the building source kernels and the implicit projection
// (erf.if_implicit_projection) must agree on which momentum faces carry the MOST wall law and
// which the linear drag toward zero, so the masks are shared (ERF_ImmersedBuildingMasks.H). The
// tests pin the masks the kernels computed inline before the refactor, case by case, and the two
// invariants the projection relies on: a face carries the wall law or the drag, never both when
// the snap is on, and the drag weight of a fully solid face is one.

using namespace amrex;

namespace {

// zero and one come from ERF_NumericalConstants.H through the header under test
constexpr Real half = Real(0.5), quarter = Real(0.25);

} // namespace

// Without the snap, as the kernels did: a fully solid horizontal face is an interior face (the roof
// law needs t < below, so the top solid row of an exact box gets the interior drag only), a partial
// face carries the roof or side wall law alone, and a partial face without a law gets the blanked
// drag
TEST(ImmersedBuildingMasks, WithoutSnapSolidFacesAreInteriorAndPartialFacesCarryOneForcing)
{
    // in the body of the building: solid on every side
    auto body = if_bld::horiz_masks(one, one, one, one, one, false, true, false);
    EXPECT_EQ(body.most(), zero);
    EXPECT_EQ(body.interior, one);
    EXPECT_EQ(body.linear(one), one);

    // the top solid row: fluid above, solid below. Without the snap the roof law needs t < below,
    // so a fully solid face in the top row keeps the interior drag only
    auto top = if_bld::horiz_masks(one, one, zero, one, one, false, true, false);
    EXPECT_EQ(top.roof, zero);
    EXPECT_EQ(top.interior, one);

    // a partial roof face (t = 0.5) over a solid row, fluid above: the roof law alone
    auto roof = if_bld::horiz_masks(half, one, zero, half, half, false, true, false);
    EXPECT_EQ(roof.roof, one);
    EXPECT_EQ(roof.side, zero);
    EXPECT_EQ(roof.interior, zero);
    EXPECT_EQ(roof.linear(half), zero);

    // a partial face on the south wall (fluid to the south, as solid or more to the north)
    auto south = if_bld::horiz_masks(half, half, half, zero, one, false, true, false);
    EXPECT_EQ(south.lo, one);
    EXPECT_EQ(south.hi, zero);
    EXPECT_EQ(south.side, zero);
    EXPECT_EQ(south.linear(half), zero);

    // a partial face with solid on both neighbouring sides and no fluid above: not a wall-law
    // face, so it gets the partial-wall drag weighted by its blanking (the east/west walls of
    // the x-momentum)
    auto side = if_bld::horiz_masks(half, half, half, half, half, false, true, false);
    EXPECT_EQ(side.most(), zero);
    EXPECT_EQ(side.side, one);
    EXPECT_EQ(side.linear(half), half);
    EXPECT_EQ(side.linear(quarter), quarter);

    // without MOST every partial face is a wall face with the blanked drag
    auto nomost = if_bld::horiz_masks(half, one, zero, zero, one, false, false, false);
    EXPECT_EQ(nomost.most(), zero);
    EXPECT_EQ(nomost.wall, one);
    EXPECT_EQ(nomost.linear(half), half);

    // open air: nothing
    auto air = if_bld::horiz_masks(zero, zero, zero, zero, zero, false, true, false);
    EXPECT_EQ(air.most(), zero);
    EXPECT_EQ(air.linear(zero), zero);
}

// With the snap every solid face has t = 1: the roof face lies in the top solid row (t <= below),
// a wall-law face is not also an interior face, and a wall-normal face gets the interior drag
TEST(ImmersedBuildingMasks, SnapGivesEachFaceOneForcing)
{
    // the top solid row: roof law, no interior drag
    auto top = if_bld::horiz_masks(one, one, zero, one, one, false, true, true);
    EXPECT_EQ(top.roof, one);
    EXPECT_EQ(top.interior, zero);
    EXPECT_EQ(top.linear(one), zero);

    // the south wall column of a solid row: wall law (lo), no interior drag
    auto south = if_bld::horiz_masks(one, one, one, zero, one, false, true, true);
    EXPECT_EQ(south.lo, one);
    EXPECT_EQ(south.interior, zero);

    // a wall-normal face (joining a solid and a fluid cell): the interior drag, no law
    auto normal = if_bld::horiz_masks(one, one, zero, zero, one, true, true, true);
    EXPECT_EQ(normal.most(), zero);
    EXPECT_EQ(normal.interior, one);
    EXPECT_EQ(normal.linear(one), one);

    // the body: interior drag
    auto body = if_bld::horiz_masks(one, one, one, one, one, false, true, true);
    EXPECT_EQ(body.interior, one);
    EXPECT_EQ(body.most(), zero);

    // the invariant the projection relies on: with the snap a face carries the law or the drag
    for (Real below : {zero, one}) {
    for (Real above : {zero, one}) {
    for (Real lo : {zero, one}) {
    for (Real hi : {zero, one}) {
    for (bool normal_face : {false, true}) {
        auto m = if_bld::horiz_masks(one, below, above, lo, hi, normal_face, true, true);
        EXPECT_TRUE(m.most() == zero || m.linear(one) == zero)
            << "below " << below << " above " << above << " lo " << lo << " hi " << hi << " normal " << normal_face;
        // a fully solid face always has one of the two
        EXPECT_TRUE(m.most() > zero || m.linear(one) == one)
            << "below " << below << " above " << above << " lo " << lo << " hi " << hi << " normal " << normal_face;
    }}}}}
}

// z-faces: the four sides carry the wall law, the roof face the blanked drag, the domain bottom
// no wall law
TEST(ImmersedBuildingMasks, VerticalFacesSidesCarryTheLawAndTheRoofTheDrag)
{
    // a partial z-face on the west wall (fluid to the west), above the floor
    auto west = if_bld::vert_masks(half, half, half, half, zero, one, true, false, true, false);
    EXPECT_EQ(west.west, one);
    EXPECT_EQ(west.east, zero);
    EXPECT_EQ(west.most(), one);
    EXPECT_EQ(west.linear(half), zero);

    // the same face at the domain bottom: no wall law, and no drag either (not interior)
    auto floor = if_bld::vert_masks(half, half, half, half, zero, one, false, false, true, false);
    EXPECT_EQ(floor.most(), zero);
    EXPECT_EQ(floor.linear(half), zero);

    // a partial roof z-face (fluid above, solid on every side): the blanked drag, no law
    auto roof = if_bld::vert_masks(half, zero, one, one, one, one, true, false, true, false);
    EXPECT_EQ(roof.roof, one);
    EXPECT_EQ(roof.most(), zero);
    EXPECT_EQ(roof.linear(half), half);

    // a fully solid roof face without the snap: roof drag and interior drag, weight 1 + t
    auto solid_roof = if_bld::vert_masks(one, zero, one, one, one, one, true, false, true, false);
    EXPECT_EQ(solid_roof.roof, one);
    EXPECT_EQ(solid_roof.interior, one);
    // with the snap: the roof drag only
    auto snap_roof = if_bld::vert_masks(one, zero, one, one, one, one, true, false, true, true);
    EXPECT_EQ(snap_roof.roof, one);
    EXPECT_EQ(snap_roof.interior, zero);
    EXPECT_EQ(snap_roof.linear(one), one);

    // without MOST: the blanked wall drag on a partial face
    auto nomost = if_bld::vert_masks(half, half, half, half, zero, one, true, false, false, false);
    EXPECT_EQ(nomost.wall, one);
    EXPECT_EQ(nomost.most(), zero);
    EXPECT_EQ(nomost.linear(half), half);

    // the body
    auto body = if_bld::vert_masks(one, one, one, one, one, one, true, false, true, false);
    EXPECT_EQ(body.interior, one);
    EXPECT_EQ(body.most(), zero);
    EXPECT_EQ(body.linear(one), one);
}
