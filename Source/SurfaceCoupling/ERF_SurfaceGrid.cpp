#include <ERF_SurfaceGrid.H>
#include <AMReX_BoxList.H>

#include <string>

using namespace amrex;

SurfaceGrid
create_surface_grid (const BoxArray& ba_atm,
                     const DistributionMapping& dm_atm,
                     const Geometry& geom_atm,
                     int C,
                     const SurfaceCouplingNames& names)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(C >= 1,
        (names.tag() + " Building the " + names.grid() + ": " + names.key("grid_ratio")
         + " must be >= 1").c_str());
    SurfaceGrid sg;
    sg.C = C;

    // Step 1: the k = 0 slab of each atmospheric box, in the order of ba_atm, so
    // that surface box n refines atmospheric box n. The k index here is the
    // surface grid's own single layer.
    Vector<Box> box_list_2d;
    for (int i = 0; i < ba_atm.size(); ++i) {
        Box b = ba_atm[i];
        b.setSmall(2, 0);
        b.setBig(2, 0);
        box_list_2d.push_back(b);
    }
    BoxList bl_2d(std::move(box_list_2d));
    BoxArray ba_2d(bl_2d);

    // Step 2: the region the grid covers, the bounding rectangle of the level's
    // boxes: the domain on level 0, the refined region on a finer level. The
    // grid treats the region's edges as its domain edges, so the boxes must
    // cover it without holes (verify_surface_grid_level says so first with the
    // inputs to change).
    const Box region = ba_2d.minimalBox();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(region.numPts() == ba_2d.numPts(),
        (names.tag() + " Building the " + names.grid()
         + ": the boxes of the atmospheric level do not form one rectangle").c_str());
    const IntVect atm_lo(region.smallEnd(0), region.smallEnd(1), 0);

    // Step 3: shift to the region's corner and refine by {C, C, 1}. The grid's
    // index space starts at zero there, as every AMReX Geometry does, so
    // positions from ProbLo + (i + 1/2) dx and the plotfile header are those of
    // the region as they stand; only the map back to the atmosphere adds
    // atm_lo. refine() and shift() keep the number and order of the boxes, so
    // dm_atm maps identically (box n of the grid is owned by the rank of
    // atmospheric box n).
    BoxArray ba_surf = ba_2d;
    ba_surf.shift(-atm_lo);
    ba_surf.refine(IntVect(C, C, 1));

    // Step 4: the 2D Geometry over the region. Where the region spans the
    // domain the domain's own bounds are used; elsewhere the bounds are the
    // atmospheric cell faces at the region's edges. The grid inherits the
    // atmospheric periodicity in x and y where the region spans the domain, so
    // FillBoundary exchanges across a periodic seam; a region narrower than the
    // domain has edges of its own and is not periodic there.
    const Box& dom_atm = geom_atm.Domain();
    Box surf_domain(IntVect(0, 0, 0),
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

    sg.ba     = ba_surf;
    sg.dm     = dm_atm;
    sg.geom   = Geometry(surf_domain, prob_domain_2d, CoordSys::cartesian, is_per);
    sg.atm_lo = atm_lo;
    return sg;
}
