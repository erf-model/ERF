#include <AMReX_ParmParse.H>
#include <ERF.H>
#include <ERF_Constants.H>
#include <IndexDefines.H>

#ifdef ERF_USE_TERRAIN
void
ERF::init_ideal_terrain(int lev)
{
    auto dx = geom[lev].CellSizeArray();
    auto ProbLoArr = geom[lev].ProbLoArray();
    auto ProbHiArr = geom[lev].ProbHiArray();

    // User function parameters
    amrex::Real a = 0.5;
    amrex::Real num = 8 * a * a * a;
    amrex::Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
    amrex::Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

    // Populate bottom plane
    int k0 = 0;

    for ( amrex::MFIter mfi(z_phys_nd[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const amrex::Box& vbx = mfi.validbox();
        amrex::Box xybx = vbx;
        xybx.setRange(2,0);
        
        amrex::Array4<amrex::Real> const& z_arr = z_phys_nd[lev].array(mfi);
        
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Location of nodes
            amrex::Real x = (i  * dx[0] - xcen);
            amrex::Real y = (j  * dx[1] - ycen);
            amrex::Real z =  k0 * dx[2];

            // WoA Hill in x-direction
            amrex::Real height = num / (x*x + 4 * a * a);

            // Populate terrain height
            z_arr(i,j,k0) = height;

            // Flat terrain with z = 0 at k = 0
            //z_arr(i,j,k0) = z;
        });
    }
}
#endif
