#include <ERF_FireGrid.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxList.H>

#include <string>

using namespace amrex;

FireGrid
create_fire_grid(const BoxArray& ba_atm,
                 const DistributionMapping& dm_atm,
                 const Geometry& geom_atm,
                 int C)
{
    FireGrid fg;
    fg.C = C;

    // Step 1: the k = 0 slab of each atmospheric box, in the order of ba_atm, so
    // that fire box n refines atmospheric box n.
    Vector<Box> box_list_2d;
    for (int i = 0; i < ba_atm.size(); ++i) {
        Box b = ba_atm[i];
        b.setSmall(2, 0);
        b.setBig(2, 0);
        box_list_2d.push_back(b);
    }
    BoxList bl_2d(std::move(box_list_2d));   // rvalue — matches Vector<Box>&&
    BoxArray ba_2d(bl_2d);

    // Step 2: the region the fire grid covers, the bounding rectangle of the
    // level's boxes: the domain on level 0, the refined region on a finer level.
    // The fire grid treats the region's edges as its domain edges, so the boxes
    // must cover it without holes (verify_fire_prerequisites says so first with
    // the inputs to change).
    const Box region = ba_2d.minimalBox();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(region.numPts() == ba_2d.numPts(),
        "[FIRE] create_fire_grid: the boxes of the atmospheric level do not form one rectangle");
    const IntVect atm_lo(region.smallEnd(0), region.smallEnd(1), 0);

    // Step 3: shift to the region's corner and refine by {C, C, 1}. The fire
    // grid's index space starts at zero there, as every AMReX Geometry does, so
    // positions from ProbLo + (i + 1/2) dx, the fuel map rows and the plotfile
    // header are those of the region as they stand; only the map back to the
    // atmosphere adds atm_lo. refine() and shift() keep the number and order of
    // the boxes, so dm_atm maps identically (box n of the fire grid is owned by
    // the rank of atmospheric box n).
    BoxArray ba_fire = ba_2d;
    ba_fire.shift(-atm_lo);
    ba_fire.refine(IntVect(C, C, 1));
    DistributionMapping dm_fire = dm_atm;

    // Step 4: the 2D Geometry over the region. Where the region spans the
    // domain the domain's own bounds are used, so a level-0 fire grid is the
    // same to the last bit as before the region was introduced; elsewhere the
    // bounds are the atmospheric cell faces at the region's edges.
    const Box& dom_atm = geom_atm.Domain();
    Box fire_domain(IntVect(0, 0, 0),
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
        // Inherit the atmospheric periodicity in x and y where the region spans
        // the domain. Hard-coding this non-periodic silently disabled every
        // FillBoundary on the fire grid: the fire grid is one box spanning the
        // domain, so all of its ghost cells are domain-boundary ghosts, and a
        // non-periodic Periodicity() leaves them untouched. The level-set stage
        // fields carry three ghost cells past the box edge and were reading
        // whatever the allocator happened to supply. A region narrower than the
        // domain has edges of its own and is not periodic there.
        is_per[d] = (spans && geom_atm.isPeriodic(d)) ? 1 : 0;
    }
    RealBox prob_domain_2d(lo, hi);

    Geometry geom_fire_2d(fire_domain, prob_domain_2d, CoordSys::cartesian, is_per);

    fg.ba = ba_fire;
    fg.dm = dm_fire;
    fg.geom = geom_fire_2d;
    fg.atm_lo = atm_lo;

    return fg;
}
