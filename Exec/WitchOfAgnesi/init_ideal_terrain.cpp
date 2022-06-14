#include <AMReX_ParmParse.H>
#include <ERF.H>
#include <ERF_Constants.H>
#include <IndexDefines.H>

#ifdef ERF_USE_TERRAIN
void
ERF::init_ideal_terrain(int lev)
{
    // Domain cell size and real bounds
    auto dx = geom[lev].CellSizeArray();
    auto ProbLoArr = geom[lev].ProbLoArray();
    auto ProbHiArr = geom[lev].ProbHiArray();

    // Domain valid box (z_nd is nodal)
    const amrex::Box& domain = geom[lev].Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2); int domhi_z = domain.bigEnd(2) + 1;

    // User function parameters
    amrex::Real a    = 0.5;
    amrex::Real num  = 8 * a * a * a;
    amrex::Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
    amrex::Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

    // Number of ghost cells
    int ngrow = z_phys_nd[lev].nGrow();

    // Populate bottom plane
    int k0 = domlo_z;

    for ( amrex::MFIter mfi(z_phys_nd[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        amrex::Array4<amrex::Real> const& z_arr = z_phys_nd[lev].array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Clip indices for ghost-cells
            int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
            int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

            // Location of nodes
            amrex::Real x = (ii  * dx[0] - xcen);
            amrex::Real y = (jj  * dx[1] - ycen);
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
