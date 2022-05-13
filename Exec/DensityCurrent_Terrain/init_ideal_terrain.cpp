#include "AMReX_ParmParse.H"
#include "ERF.H"
#include "ERF_Constants.H"
#include "IndexDefines.H"

using namespace amrex;

#ifdef ERF_USE_TERRAIN
void
ERF::init_ideal_terrain(int lev)
{
    auto dx = geom[lev].CellSizeArray();
    auto ProbLoArr = geom[lev].ProbLoArray();
    auto ProbHiArr = geom[lev].ProbHiArray();

    for ( MFIter mfi(z_phys_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& gbx = mfi.growntilebox(1);
        Array4<Real> z_arr = z_phys_nd[lev].array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real z = k * dx[2];

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k) = 2.0 * z;
        });
    }
    z_phys_nd[lev].FillBoundary(geom[lev].periodicity());
}
#endif
