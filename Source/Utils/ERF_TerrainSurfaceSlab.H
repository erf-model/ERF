/**
 * \file ERF_TerrainSurfaceSlab.H
 */
#ifndef ERF_TERRAINSURFACESLAB_H_
#define ERF_TERRAINSURFACESLAB_H_

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

class ProblemBase;

//
// The terrain height under the grids of a level, on the nodes of the bottom
// slab (k = the bottom of the domain) of each box.
//
// With a terrain-fitted mesh the terrain is the bottom of the mesh, z_phys_nd at
// the bottom nodes.  With immersed-forcing terrain the mesh is flat and the
// terrain is the surface the immersed boundary is built from; these functions
// give that surface at a level's resolution for code that measures heights
// above it.
//

// The nodal boxes of the bottom slab under each cell-centred box of ba
amrex::BoxArray bottom_node_slab (const amrex::BoxArray& ba, const amrex::Geometry& geom);

// Fill zs (defined on bottom_node_slab) with the terrain surface that prob gives
// at the resolution of geom, as the immersed boundary is built
void fill_terrain_surface_slab (amrex::MultiFab& zs, const amrex::Geometry& geom,
                                ProblemBase& prob, double time);

#endif
