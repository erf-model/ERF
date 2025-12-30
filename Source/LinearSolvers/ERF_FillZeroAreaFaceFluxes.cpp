#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

/**
 * Compute phi gradients where the area of the face is zero
 */
void ERF::FillZeroAreaFaceFluxes (int lev, Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes)
{
    BL_PROFILE("ERF::FillZeroAreaFaceFluxes()");

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
// #ifdef _OPENMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//     {
//         Box const& bx    = mfi.tilebox();


//     } // mfi

    // Now overwrite with periodic fill outside domain and fine-fine fill inside
    // phi.FillBoundary(geom[lev].periodicity());
}
