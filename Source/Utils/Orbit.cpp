#include <Orbit.H>

using namespace amrex;

void
zenith (int& calday,
        amrex::MultiFab* clat,
        amrex::MultiFab* clon,
        real1d& coszrs,
        int& ncol,
        const Real& eccen,
        const Real& mvelpp,
        const Real& lambm0,
        const Real& obliqr,
        amrex::Real uniform_angle)
{
    Real delta;    // Solar declination angle  in radians
    Real eccf;     // Earth orbit eccentricity factor

    // Populate delta & eccf
    shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf);

    // If we have a valid pointer, go through the whole machinery
    if (clat) {
        for (MFIter mfi(*clat); mfi.isValid(); ++mfi) {
            const auto& tbx = mfi.tilebox();
            auto nx = tbx.length(0);

            auto lat_array = clat->array(mfi);
            auto lon_array = clon->array(mfi);

            // NOTE: lat/lon are 2D multifabs!
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
            {
                auto icol = j*nx+i+1;
                coszrs(icol) = shr_orb_cosz(calday, lat_array(i,j,0), lon_array(i,j,0), delta, uniform_angle);
            });
       }
    }
    // Use constant value near center of USA or the uniform angle
    else {
        Real cons_lat = 40.0;
        Real cons_lon = 100.0;
        Real val = shr_orb_cosz(calday, cons_lat, cons_lon, delta, uniform_angle);
        yakl::memset(coszrs, val);
    }
}
