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
using amrex::makeSlab;
using amrex::surroundingNodes;

namespace {

using Cell    = std::tuple<int,int,int>;
using CellSet = std::set<Cell>;

Box test_domain () { return Box(IntVect(0,0,0), IntVect(11,9,5)); }

// The ways a domain can be handed to the advection routines: as one box, or split
// across ranks in any combination of directions. Every open-boundary answer below
// must come out the same for all of them.
std::vector<BoxArray> decompositions (const Box& domain)
{
    std::vector<BoxArray> out;
    for (const IntVect& max_size : {IntVect(64,64,64), IntVect( 4,64,64),
                                    IntVect(64, 5,64), IntVect(64,64, 3),
                                    IntVect( 3, 4, 2)})
    {
        BoxArray ba(domain);
        ba.maxSize(max_size);
        out.push_back(ba);
    }
    return out;
}

void insert (CellSet& s, const Box& b)
{
    for (int k(b.smallEnd(2)); k <= b.bigEnd(2); ++k) {
    for (int j(b.smallEnd(1)); j <= b.bigEnd(1); ++j) {
    for (int i(b.smallEnd(0)); i <= b.bigEnd(0); ++i) {
        s.insert(Cell(i,j,k));
    }}}
}

// Returns false if any cell of b is already in s -- i.e. two open-boundary kernels
// would both write it, and since they assign rather than accumulate, one would win.
bool insert_disjoint (CellSet& s, const Box& b)
{
    bool disjoint = true;
    for (int k(b.smallEnd(2)); k <= b.bigEnd(2); ++k) {
    for (int j(b.smallEnd(1)); j <= b.bigEnd(1); ++j) {
    for (int i(b.smallEnd(0)); i <= b.bigEnd(0); ++i) {
        if (!s.insert(Cell(i,j,k)).second) { disjoint = false; }
    }}}
    return disjoint;
}

// Every cell of `b` that lies on an open lateral boundary; what the patches must cover.
CellSet expected_band (const Box& b, const Box& domain,
                       bool xlo_open, bool xhi_open, bool ylo_open, bool yhi_open)
{
    const int dom_xhi = domain.bigEnd(0) + (b.ixType().nodeCentered(0) ? 1 : 0);
    const int dom_yhi = domain.bigEnd(1) + (b.ixType().nodeCentered(1) ? 1 : 0);

    CellSet band;
    for (int k(b.smallEnd(2)); k <= b.bigEnd(2); ++k) {
    for (int j(b.smallEnd(1)); j <= b.bigEnd(1); ++j) {
    for (int i(b.smallEnd(0)); i <= b.bigEnd(0); ++i) {
        const bool on_open_bndry = (xlo_open && i == domain.smallEnd(0)) ||
                                   (xhi_open && i == dom_xhi)            ||
                                   (ylo_open && j == domain.smallEnd(1)) ||
                                   (yhi_open && j == dom_yhi);
        if (on_open_bndry) { band.insert(Cell(i,j,k)); }
    }}}
    return band;
}

// ---------------------------------------------------------------------------------
// Motivation: this is the guard the fix turns on. The open flags are properties of
// the domain, but the boxes handed to the shrink are per-rank slabs whose ends are
// usually in the interior. Trimming those would drop interior cells from the
// open-boundary treatment and make the answer depend on the decomposition.
// ---------------------------------------------------------------------------------
TEST(OpenBCCorners, TrimOnlyActsWhereTheBoxReachesTheDomainBoundary)
{
    const Box domain = test_domain();

    const Box interior(IntVect(0,3,0), IntVect(0,6,5));
    EXPECT_EQ(TrimOpenBCCorner(interior, 1, domain, true, true), interior);

    const Box touches_lo(IntVect(0,domain.smallEnd(1),0), IntVect(0,6,5));
    const Box trimmed_lo = TrimOpenBCCorner(touches_lo, 1, domain, true, true);
    EXPECT_EQ(trimmed_lo.smallEnd(1), touches_lo.smallEnd(1)+1);
    EXPECT_EQ(trimmed_lo.bigEnd(1),   touches_lo.bigEnd(1));

    const Box touches_hi(IntVect(0,3,0), IntVect(0,domain.bigEnd(1),5));
    const Box trimmed_hi = TrimOpenBCCorner(touches_hi, 1, domain, true, true);
    EXPECT_EQ(trimmed_hi.smallEnd(1), touches_hi.smallEnd(1));
    EXPECT_EQ(trimmed_hi.bigEnd(1),   touches_hi.bigEnd(1)-1);

    // A closed boundary owns nothing, so nothing is given up to it
    EXPECT_EQ(TrimOpenBCCorner(touches_lo, 1, domain, false, false), touches_lo);
}

// ---------------------------------------------------------------------------------
// Motivation: a node-centered box's high boundary sits one index past the domain's,
// so inferring the index space from the box is what keeps the trim on the boundary
// face rather than one face inside it.
// ---------------------------------------------------------------------------------
TEST(OpenBCCorners, TrimUsesTheIndexSpaceOfTheBoxItIsGiven)
{
    const Box domain = test_domain();

    const Box cc(IntVect(0,0,0), IntVect(domain.bigEnd(0),0,0));
    EXPECT_EQ(TrimOpenBCCorner(cc, 0, domain, false, true).bigEnd(0), domain.bigEnd(0)-1);

    const Box nd = surroundingNodes(cc,0);
    ASSERT_EQ(nd.bigEnd(0), domain.bigEnd(0)+1);
    EXPECT_EQ(TrimOpenBCCorner(nd, 0, domain, false, true).bigEnd(0), domain.bigEnd(0));

    // A node-centered box that stops one face short of the boundary is left alone
    Box nd_short(nd); nd_short.growHi(0,-1);
    EXPECT_EQ(TrimOpenBCCorner(nd_short, 0, domain, false, true), nd_short);
}

// ---------------------------------------------------------------------------------
// Motivation: the bug itself, for w and the state. Those are tangent to both lateral
// boundaries, so nothing else claims the corner cell; before the fix both directions'
// tangential kernels wrote it and the later launch won, using a stencil that
// differenced across the open boundary it was not treating.
// ---------------------------------------------------------------------------------
TEST(OpenBCCorners, TangentPatchesCoverEveryOpenBoundaryCellExactlyOnce)
{
    const Box domain = test_domain();

    for (const BoxArray& ba : decompositions(domain)) {
        CellSet covered;
        CellSet expected;
        bool disjoint = true;

        for (int ibox(0); ibox < ba.size(); ++ibox) {
            const Box b = ba[ibox];
            for (const OpenBCPatch& patch : OpenBCTangentPatches(b, domain, true, true, true, true)) {
                disjoint = insert_disjoint(covered, patch.box) && disjoint;
            }
            const CellSet band = expected_band(b, domain, true, true, true, true);
            expected.insert(band.begin(), band.end());
        }

        EXPECT_TRUE(disjoint) << "an open-boundary cell is written by two kernels";
        EXPECT_EQ(covered, expected);
    }
}

// ---------------------------------------------------------------------------------
// Motivation: a corner must be tagged open on both sides. That tag is what makes the
// kernel use the open-boundary form in x and in y, so it reaches across neither.
// ---------------------------------------------------------------------------------
TEST(OpenBCCorners, CornersAreTaggedOpenInBothDirections)
{
    const Box domain = test_domain();
    const auto patches = OpenBCTangentPatches(domain, domain, true, true, true, true);

    int n_corner = 0;
    for (const OpenBCPatch& patch : patches) {
        const bool is_corner = (patch.x_side != OpenSide::none) && (patch.y_side != OpenSide::none);
        if (!is_corner) { continue; }
        ++n_corner;
        EXPECT_EQ(patch.box.length(0), 1);
        EXPECT_EQ(patch.box.length(1), 1);
        const int i = patch.box.smallEnd(0);
        const int j = patch.box.smallEnd(1);
        EXPECT_EQ(i, (patch.x_side == OpenSide::lo) ? domain.smallEnd(0) : domain.bigEnd(0));
        EXPECT_EQ(j, (patch.y_side == OpenSide::lo) ? domain.smallEnd(1) : domain.bigEnd(1));
    }
    EXPECT_EQ(n_corner, 4);

    // With only one direction open there is no corner, and the single edge is untrimmed
    const auto one_side = OpenBCTangentPatches(domain, domain, false, true, false, false);
    ASSERT_EQ(one_side.size(), 1u);
    EXPECT_EQ(one_side[0].x_side, OpenSide::hi);
    EXPECT_EQ(one_side[0].y_side, OpenSide::none);
    EXPECT_EQ(one_side[0].box.length(1), domain.length(1));
}

// ---------------------------------------------------------------------------------
// Motivation: the x-momentum case, mirroring AdvectionSrcForMom. A face on an open x
// boundary belongs to the boundary-normal (radiation) kernel even where it also lies
// on an open y boundary; before the fix the y-direction tangential kernel ran last
// and replaced the radiation condition with a centered difference reaching outside
// the domain. The two sets of faces must not overlap, and must not depend on how the
// grid is split.
// ---------------------------------------------------------------------------------
TEST(OpenBCCorners, NormalAndTangentMomentumFacesAreDisjointForEveryDecomposition)
{
    const Box domain = test_domain();
    const bool xlo_open = true, xhi_open = true, ylo_open = true, yhi_open = true;

    CellSet reference_tangent;
    bool first = true;

    for (const BoxArray& ba : decompositions(domain)) {
        CellSet normal, tangent;

        for (int ibox(0); ibox < ba.size(); ++ibox) {
            const Box tbx = surroundingNodes(ba[ibox],0);

            if (xlo_open && tbx.smallEnd(0) == domain.smallEnd(0)) {
                insert(normal, makeSlab(tbx,0,domain.smallEnd(0)));
            }
            if (xhi_open && tbx.bigEnd(0) == domain.bigEnd(0)+1) {
                insert(normal, makeSlab(tbx,0,domain.bigEnd(0)+1));
            }
            if (ylo_open && tbx.smallEnd(1) == domain.smallEnd(1)) {
                insert(tangent, TrimOpenBCCorner(makeSlab(tbx,1,domain.smallEnd(1)),
                                                 0, domain, xlo_open, xhi_open));
            }
            if (yhi_open && tbx.bigEnd(1) == domain.bigEnd(1)) {
                insert(tangent, TrimOpenBCCorner(makeSlab(tbx,1,domain.bigEnd(1)),
                                                 0, domain, xlo_open, xhi_open));
            }
        }

        for (const Cell& c : tangent) {
            EXPECT_EQ(normal.count(c), 0u)
                << "face (" << std::get<0>(c) << "," << std::get<1>(c) << "," << std::get<2>(c)
                << ") is written by both the normal and the tangential open-BC kernel";
        }

        if (first) { reference_tangent = tangent; first = false; }
        else       { EXPECT_EQ(tangent, reference_tangent); }
    }

    // The corner faces really are in play: they belong to the normal kernel, not
    // the tangential one
    for (int k(domain.smallEnd(2)); k <= domain.bigEnd(2); ++k) {
        EXPECT_EQ(reference_tangent.count(Cell(domain.bigEnd(0)+1, domain.bigEnd(1), k)), 0u);
        EXPECT_EQ(reference_tangent.count(Cell(domain.smallEnd(0), domain.smallEnd(1), k)), 0u);
    }
}

} // namespace
