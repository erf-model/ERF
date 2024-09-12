#ifndef ERF_TILE_NO_Z_H_
#define ERF_TILE_NO_Z_H_

#include <AMReX.H>

/**
 * Function returns explicit tile size that prevents tiling in z
 */

AMREX_FORCE_INLINE
amrex::IntVect TileNoZ ()
{
    if (amrex::TilingIfNotGPU()) {
        return amrex::IntVect{amrex::FabArrayBase::mfiter_tile_size[0],
                              amrex::FabArrayBase::mfiter_tile_size[1],
                              1024000};
    }
    else {
        return amrex::IntVect::TheZeroVector();
    }
}

#endif
