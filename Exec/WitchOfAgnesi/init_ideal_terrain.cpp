#include "AMReX_ParmParse.H"
#include "ERF.H"
#include "ERF_Constants.H"
#include "IndexDefines.H"

using namespace amrex;

void
ERF::init_ideal_terrain(int lev)
{
    auto dx = geom[lev].CellSize();
    Real dz = dx[2];
    for ( MFIter mfi(z_phys_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> z_arr = z_phys_nd[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            z_arr(i,j,k) = k * dz;
        });
    }
}
