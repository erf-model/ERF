#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>

#include <ERF_Advection.H>

#include <gtest/gtest.h>

#include <set>
#include <tuple>
#include <vector>

using amrex::Box;
using amrex::BoxArray;
using amrex::IntVect;
using amrex::surroundingNodes;

namespace {

using Face    = std::tuple<int,int,int>;
using FaceSet = std::set<Face>;

Box test_domain () { return Box(IntVect(0,0,0), IntVect(7,5,15)); }

// Splits of the domain across ranks, including several that cut in z. Whether a grid
// is split in z must not change which w-faces get a source term.
std::vector<BoxArray> decompositions (const Box& domain)
{
    std::vector<BoxArray> out;
    for (const IntVect& max_size : {IntVect(64,64,64), IntVect(64,64, 8),
                                    IntVect(64,64, 4), IntVect( 4, 3, 4)})
    {
        BoxArray ba(domain);
        ba.maxSize(max_size);
        out.push_back(ba);
    }
    return out;
}

void insert (FaceSet& s, const Box& b)
{
    for (int k(b.smallEnd(2)); k <= b.bigEnd(2); ++k) {
    for (int j(b.smallEnd(1)); j <= b.bigEnd(1); ++j) {
    for (int i(b.smallEnd(0)); i <= b.bigEnd(0); ++i) {
        s.insert(Face(i,j,k));
    }}}
}

// Every w-face that should carry a source term: all of them but the two on the
// bottom and top of the domain.
FaceSet expected_faces (const Box& domain)
{
    Box all = surroundingNodes(domain,2);
    all.growLo(2,-1);
    all.growHi(2,-1);

    FaceSet s;
    insert(s, all);
    return s;
}

// ---------------------------------------------------------------------------------
// Motivation: the shrink exists to skip the bottom and top of the *domain*. The boxes
// it is handed are per-grid, and on a grid decomposed in z their ends are usually in
// the interior, where nothing should be skipped.
// ---------------------------------------------------------------------------------
TEST(ZmomBoxDecomposition, ShrinkActsOnlyWhereTheBoxReachesTheDomainEnd)
{
    const Box domain = test_domain();

    const Box interior = surroundingNodes(Box(IntVect(0,0,4), IntVect(7,5,11)),2);
    EXPECT_EQ(ShrinkZmomBoxAtDomainEnds(interior, domain), interior);

    const Box at_bottom = surroundingNodes(Box(IntVect(0,0,0), IntVect(7,5,7)),2);
    const Box cut_bottom = ShrinkZmomBoxAtDomainEnds(at_bottom, domain);
    EXPECT_EQ(cut_bottom.smallEnd(2), at_bottom.smallEnd(2)+1);
    EXPECT_EQ(cut_bottom.bigEnd(2),   at_bottom.bigEnd(2));

    const Box at_top = surroundingNodes(Box(IntVect(0,0,8), IntVect(7,5,15)),2);
    const Box cut_top = ShrinkZmomBoxAtDomainEnds(at_top, domain);
    EXPECT_EQ(cut_top.smallEnd(2), at_top.smallEnd(2));
    EXPECT_EQ(cut_top.bigEnd(2),   at_top.bigEnd(2)-1);
}

// ---------------------------------------------------------------------------------
// Motivation: the bug. An unconditional shrink removes the shared w-face at every
// internal z split -- and it removes it from the grid below *and* the grid above, so
// no grid computes it. Anything run on such a box, in particular the open-boundary
// advection kernels, silently skips those faces and leaves whatever the interior
// stencil put there.
// ---------------------------------------------------------------------------------
TEST(ZmomBoxDecomposition, UnconditionalShrinkWouldLeaveAGapAtEveryInternalSplit)
{
    const Box domain = test_domain();

    BoxArray ba(domain);
    ba.maxSize(IntVect(64,64,8));
    ASSERT_GT(ba.size(), 1);

    FaceSet guarded, unconditional;
    for (int ibox(0); ibox < ba.size(); ++ibox) {
        Box nodal = surroundingNodes(ba[ibox],2);
        insert(guarded, ShrinkZmomBoxAtDomainEnds(nodal, domain));

        Box bad(nodal); bad.growLo(2,-1); bad.growHi(2,-1);
        insert(unconditional, bad);
    }

    EXPECT_EQ(guarded, expected_faces(domain));

    // The split is at k=8, and the unconditional shrink loses that whole plane
    EXPECT_LT(unconditional.size(), guarded.size());
    for (int j(domain.smallEnd(1)); j <= domain.bigEnd(1); ++j) {
    for (int i(domain.smallEnd(0)); i <= domain.bigEnd(0); ++i) {
        EXPECT_EQ(unconditional.count(Face(i,j,8)), 0u);
        EXPECT_EQ(guarded.count(Face(i,j,8)),       1u);
    }}
}

// ---------------------------------------------------------------------------------
// Motivation: the property that has to hold however the grid is split. Note the faces
// are covered, not owned: a w-face shared by two grids in z is in both grids' boxes
// and both compute the same value into their own copy, which is how nodal data works
// here. What must not happen is a face covered by neither.
// ---------------------------------------------------------------------------------
TEST(ZmomBoxDecomposition, EveryInteriorWFaceIsCoveredForEveryDecomposition)
{
    const Box domain = test_domain();
    const FaceSet expected = expected_faces(domain);

    for (const BoxArray& ba : decompositions(domain)) {
        FaceSet covered;
        for (int ibox(0); ibox < ba.size(); ++ibox) {
            insert(covered, ShrinkZmomBoxAtDomainEnds(surroundingNodes(ba[ibox],2), domain));
        }
        EXPECT_EQ(covered, expected);
    }
}

} // namespace
