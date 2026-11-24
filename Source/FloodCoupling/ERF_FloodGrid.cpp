#include <ERF_FloodGrid.H>
#include <AMReX_BoxList.H>

using namespace amrex;

FloodGrid
create_flood_grid (const BoxArray& ba_atm,
                   const DistributionMapping& dm_atm,
                   const Geometry& geom_atm,
                   int C)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(C >= 1, "[FLOOD] create_flood_grid: erf.flood.grid_ratio must be >= 1");
    FloodGrid fg;
    fg.C = C;

    // Step 1: the k = 0 slab of each atmospheric box, in the order of ba_atm, so
    // that flood box n refines atmospheric box n. The k index here is the flood
    // grid's own single layer; the atmospheric slab read in FloodCoupler is at
    // the bottom of the level's domain, which is k = 0 in ERF.
    Vector<Box> box_list_2d;
    for (int i = 0; i < ba_atm.size(); ++i) {
        Box b = ba_atm[i];
        b.setSmall(2, 0);
        b.setBig(2, 0);
        box_list_2d.push_back(b);
    }
    BoxList bl_2d(std::move(box_list_2d));
    BoxArray ba_2d(bl_2d);

    // Step 2: the region the flood grid covers, the bounding rectangle of the
    // level's boxes. The flood grid treats the region's edges as its domain
    // edges, so the boxes must cover it without holes
    // (verify_flood_prerequisites says so first with the inputs to change).
    const Box region = ba_2d.minimalBox();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(region.numPts() == ba_2d.numPts(),
        "[FLOOD] create_flood_grid: the boxes of the atmospheric level do not form one rectangle");
    const IntVect atm_lo(region.smallEnd(0), region.smallEnd(1), 0);

    // Step 3: shift to the region's corner and refine by {C, C, 1}. refine()
    // and shift() keep the number and order of the boxes, so dm_atm maps
    // identically (box n of the flood grid is owned by the rank of atmospheric
    // box n).
    BoxArray ba_flood = ba_2d;
    ba_flood.shift(-atm_lo);
    ba_flood.refine(IntVect(C, C, 1));

    // Step 4: the 2D Geometry over the region. Where the region spans the
    // domain the domain's own bounds are used; elsewhere the bounds are the
    // atmospheric cell faces at the region's edges. Periodic only in a
    // direction the region spans and the atmosphere is periodic in.
    const Box& dom_atm = geom_atm.Domain();
    Box flood_domain(IntVect(0, 0, 0),
                     IntVect(region.length(0) * C - 1, region.length(1) * C - 1, 0));

    Array<Real, AMREX_SPACEDIM> lo {geom_atm.ProbLo(0), geom_atm.ProbLo(1), geom_atm.ProbLo(2)};
    Array<Real, AMREX_SPACEDIM> hi {geom_atm.ProbHi(0), geom_atm.ProbHi(1), geom_atm.ProbLo(2) + Real(1.0)};  // z extent = 1 m (dummy)
    Array<int, AMREX_SPACEDIM> is_per {0, 0, 0};
    for (int d = 0; d < 2; ++d) {
        const bool spans = (region.smallEnd(d) == dom_atm.smallEnd(d)) && (region.bigEnd(d) == dom_atm.bigEnd(d));
        if (!spans) {
            lo[d] = geom_atm.ProbLo(d) + Real(region.smallEnd(d) - dom_atm.smallEnd(d))     * geom_atm.CellSize(d);
            hi[d] = geom_atm.ProbLo(d) + Real(region.bigEnd(d)   - dom_atm.smallEnd(d) + 1) * geom_atm.CellSize(d);
        }
        is_per[d] = (spans && geom_atm.isPeriodic(d)) ? 1 : 0;
    }
    RealBox prob_domain_2d(lo, hi);

    fg.ba     = ba_flood;
    fg.dm     = dm_atm;
    fg.geom   = Geometry(flood_domain, prob_domain_2d, CoordSys::cartesian, is_per);
    fg.atm_lo = atm_lo;
    return fg;
}
