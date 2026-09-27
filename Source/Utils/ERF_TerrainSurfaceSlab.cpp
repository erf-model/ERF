/**
 * \file ERF_TerrainSurfaceSlab.cpp
 */
#include "ERF_TerrainSurfaceSlab.H"
#include "ERF_ProbCommon.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>

using namespace amrex;

BoxArray
bottom_node_slab (const BoxArray& ba, const Geometry& geom)
{
    const int klo = geom.Domain().smallEnd(2);
    BoxList bl = ba.boxList();
    bl.convert(IndexType::TheNodeType());   // the list's index type as well as each box's
    for (auto& b : bl) {
        b.setRange(2, klo);
    }
    return BoxArray(std::move(bl));
}

void
fill_terrain_surface_slab (MultiFab& zs, const Geometry& geom, ProblemBase& prob, double time)
{
    // The whole domain at once, as ERF::fill_terrain_surface and the immersed
    // boundary build it, so that the surface is the same whichever box asks
    const int klo = geom.Domain().smallEnd(2);
    const Box bx(surroundingNodes(geom.Domain()));
    FArrayBox terrain_fab(makeSlab(bx, 2, klo), 1);
    prob.init_terrain_surface(geom, terrain_fab, time);
    for (MFIter mfi(zs); mfi.isValid(); ++mfi) {
        const Box isect = terrain_fab.box() & zs[mfi].box();
        if (!isect.isEmpty()) {
            zs[mfi].template copy<RunOn::Device>(terrain_fab, isect, 0, isect, 0, 1);
        }
    }
}
