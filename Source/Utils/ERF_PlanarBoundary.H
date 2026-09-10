#ifndef ERF_PLANAR_BOUNDARY_H_
#define ERF_PLANAR_BOUNDARY_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Vector.H>

/*
 * Planar (2D) MultiFabs, such as the surface-layer fields and the MOST averages, are
 * built on the z-collapse of a 3D BoxArray, one 2D box per 3D box.  When the 3D
 * BoxArray is split in z, several 3D boxes collapse onto the same 2D box, so the
 * planar BoxArray holds duplicate boxes with overlapping valid regions, and only the
 * copy that belongs to the 3D box touching the surface is ever computed.
 *
 * MakeSurfaceBoxes selects the 2D boxes of the 3D boxes touching the surface, on the
 * ranks that own them, and records the index of each in the planar BoxArray.
 */
void MakeSurfaceBoxes (const amrex::BoxArray& ba3d,
                       const amrex::DistributionMapping& dm3d,
                       int klo,
                       amrex::BoxArray& ba_sfc,
                       amrex::DistributionMapping& dm_sfc,
                       amrex::Vector<int>& src_index);

/*
 * FillBoundary for a planar MultiFab built on the z-collapse of a 3D BoxArray.
 * A FillBoundary on a BoxArray with duplicate boxes may fill a ghost cell from an
 * uncomputed copy, so the computed surface copies are gathered on the BoxArray from
 * MakeSurfaceBoxes and every copy, valid and ghost, is filled from that.
 */
void FillPlanarBoundary (amrex::MultiFab& mf,
                         const amrex::BoxArray& ba_sfc,
                         const amrex::DistributionMapping& dm_sfc,
                         const amrex::Vector<int>& src_index,
                         const amrex::Periodicity& period);

#endif
