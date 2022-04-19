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

    // 2 a is the high point of the hill
    Real a = 0.5;

    Real num = 8 * a * a * a;
    Real ztop = ProbHiArr[2];

    Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
    Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

    for ( MFIter mfi(z_phys_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& gbx = mfi.growntilebox(1);
        Array4<Real> z_arr = z_phys_nd[lev].array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            // Location of nodes
            Real x = (i * dx[0] - xcen);
            Real y = j * dx[1] - ycen;
            Real z = k * dx[2];

            // WoA Hill in x-direction
            Real woa = num / (x*x + 4 * a * a);

            // WoA Hill in y-direction
            // Real woa = num / (y*y + 4 * a * a);

            // This is the BTF model from p2163 of Klemp2011
            z_arr(i,j,k) = z + (1. - (z/ztop)) * woa;

            // Flat terrain with z = 0 at k = 0
            //z_arr(i,j,k) = z;
        });
    }
}
#endif
