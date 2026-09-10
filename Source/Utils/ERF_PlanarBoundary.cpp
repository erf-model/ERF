/**
 * \file ERF_PlanarBoundary.cpp
 */
#include <ERF_PlanarBoundary.H>

using namespace amrex;

/**
 * Select the surface copies of the planar boxes.
 *
 * @param[in]  ba3d      3D BoxArray the planar BoxArray was collapsed from
 * @param[in]  dm3d      DistributionMapping of the 3D BoxArray (shared by the planar MultiFabs)
 * @param[in]  klo       k index of the surface
 * @param[out] ba_sfc    2D boxes of the 3D boxes that touch the surface, without duplicates
 * @param[out] dm_sfc    ranks owning those boxes
 * @param[out] src_index index of each of those boxes in the planar BoxArray
 */
void MakeSurfaceBoxes (const BoxArray& ba3d,
                       const DistributionMapping& dm3d,
                       int klo,
                       BoxArray& ba_sfc,
                       DistributionMapping& dm_sfc,
                       Vector<int>& src_index)
{
    BoxList bl_sfc(ba3d.ixType());
    Vector<int> pmap;
    src_index.clear();
    for (int ib = 0; ib < static_cast<int>(ba3d.size()); ++ib) {
        if (ba3d[ib].smallEnd(2) == klo) {
            Box b = ba3d[ib]; b.setRange(2,0);
            bl_sfc.push_back(b);
            pmap.push_back(dm3d[ib]);
            src_index.push_back(ib);
        }
    }
    ba_sfc = BoxArray(std::move(bl_sfc));
    dm_sfc = DistributionMapping(std::move(pmap));

    // Exactly one computed copy per column: the surface boxes must not overlap
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ba_sfc.isDisjoint(),
        "MakeSurfaceBoxes: the boxes touching the surface overlap in the plane, so a planar "
        "field would have more than one computed copy per column");
}

/**
 * FillBoundary for a planar MultiFab built on the z-collapse of a 3D BoxArray.
 *
 * @param[in,out] mf        planar MultiFab whose ghost cells are to be filled
 * @param[in]     ba_sfc    surface boxes from MakeSurfaceBoxes
 * @param[in]     dm_sfc    ranks owning the surface boxes
 * @param[in]     src_index index of each surface box in the planar BoxArray
 * @param[in]     period    periodicity of the level
 */
void FillPlanarBoundary (MultiFab& mf,
                         const BoxArray& ba_sfc,
                         const DistributionMapping& dm_sfc,
                         const Vector<int>& src_index,
                         const Periodicity& period)
{
    const int nsfc = static_cast<int>(src_index.size());

    // Every planar box is a surface box: no duplicates, FillBoundary is well defined
    if (nsfc == static_cast<int>(mf.size())) {
        mf.FillBoundary(period);
        return;
    }

    // No box on this level touches the surface: nothing to fill from
    if (nsfc == 0) { return; }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nsfc < static_cast<int>(mf.size()),
        "FillPlanarBoundary: more surface boxes than planar boxes");

    const int ncomp = mf.nComp();
    MultiFab sfc(ba_sfc, dm_sfc, ncomp, 0);
    for (MFIter mfi(sfc); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        // The surface copy must be the planar box with the same footprint
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bx == mf.boxArray()[src_index[mfi.index()]],
            "FillPlanarBoundary: surface box does not match its planar box");
        sfc[mfi].template copy<RunOn::Device>(mf[src_index[mfi.index()]], bx, 0, bx, 0, ncomp);
    }

    // Fill every copy, valid and ghost, from the computed surface copies
    mf.ParallelCopy(sfc, 0, 0, ncomp, IntVect(0), mf.nGrowVect(), period);
}
