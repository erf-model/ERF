#include "ERF.H"
#include "IndexDefines.H"

using namespace amrex;

#ifdef ERF_USE_TERRAIN
void
ERF::init_ideal_terrain(int lev)
{
    auto dx = geom[lev].CellSizeArray();
    auto ProbLoArr = geom[lev].ProbLoArray();
    auto ProbHiArr = geom[lev].ProbHiArray();

    // Number of ghost cells
    int ngrow = z_phys_nd[lev].nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<amrex::Real> const& z_arr = z_phys_nd[lev].array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}
#endif
