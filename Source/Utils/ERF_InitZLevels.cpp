#include <ERF_TerrainMetrics.H>
#include <ERF_Utils.H>
#include <AMReX_ParmParse.H>
#include <ERF_Constants.H>
#include <ERF_Interpolation_1D.H>
#include <cmath>

using namespace amrex;

void
init_zlevels (Vector<Vector<Real>>& zlevels_stag,
              Vector<Vector<Real>>& stretched_dz_h,
              Vector<Gpu::DeviceVector<Real>>& stretched_dz_d,
              Vector<Geometry> const& geom,
              Vector<IntVect> const& ref_ratio,
              const Real grid_stretching_ratio,
              const Real zsurf,
              const Real dz0)
{
    int max_level = zlevels_stag.size()-1;

    for (int lev = 0; lev <= max_level; lev++)
    {
        auto dx = geom[lev].CellSizeArray();
        const Box& domain = geom[lev].Domain();
        int nz = domain.length(2);

        zlevels_stag[lev].resize(nz+1);

        stretched_dz_h[lev].resize(domain.length(2));

        if (grid_stretching_ratio == 0) {
            // This is the default for z_levels
            for (int k = 0; k < nz+1; k++)
            {
                zlevels_stag[lev][k] = k * dx[2];
            }
            for (int k = 0; k < nz; k++)
            {
                stretched_dz_h[lev][k] = dx[2];
            }
        } else if (lev == 0) {
            // Create stretched grid based on initial dz and stretching ratio
            zlevels_stag[lev][0] = zsurf;
            Real dz = dz0;
            Print() << "Stretched grid levels at level : " << lev << " is " <<  zsurf;
            for (int k = 1; k < nz+1; k++)
            {
                stretched_dz_h[lev][k-1] = dz;
                zlevels_stag[lev][k] = zlevels_stag[lev][k-1] + dz;
                Print() << " " << zlevels_stag[lev][k];
                dz *= grid_stretching_ratio;
            }
            Print() << std::endl;
        } else if (lev > 0) {
            int rr = ref_ratio[lev-1][2];
            expand_and_interpolate_1d(zlevels_stag[lev], zlevels_stag[lev-1], rr, false);
            for (int k = 0; k < nz; k++)
            {
                stretched_dz_h[lev][k] = (zlevels_stag[lev][k+1] - zlevels_stag[lev][k]);
            }
        }
    }

    // Try reading in terrain_z_levels, which allows arbitrarily spaced grid
    // levels to be specified and will take precedence over the
    // grid_stretching_ratio parameter
    ParmParse pp("erf");
    int n_zlevels = pp.countval("terrain_z_levels");
    if (n_zlevels > 0)
    {
        int nz = geom[0].Domain().length(2);
        if (n_zlevels != nz+1) {
            Print() << "You supplied " << n_zlevels << " staggered terrain_z_levels " << std::endl;
            Print() << "but n_cell+1 in the z-direction is " << nz+1 << std::endl;
            Abort("You must specify a z_level for every value of k");
        }

        if (grid_stretching_ratio > 0) {
            Print() << "Note: Found terrain_z_levels, ignoring grid_stretching_ratio" << std::endl;
        }

        pp.getarr("terrain_z_levels", zlevels_stag[0], 0, nz+1);

        // These levels should range from 0 at the surface to the height of the
        // top of model domain (see the coordinate surface height, zeta, in
        // Klemp 2011)
        AMREX_ALWAYS_ASSERT(zlevels_stag[0][0] == 0);

        for (int lev = 1; lev <= max_level; lev++) {
            int rr = ref_ratio[lev-1][2];
            expand_and_interpolate_1d(zlevels_stag[lev], zlevels_stag[lev-1], rr, false);
        }

        for (int lev = 0; lev <= max_level; lev++) {
            int nz_zlevs = zlevels_stag[lev].size();
            for (int k = 0; k < nz_zlevs-1; k++)
            {
                stretched_dz_h[lev][k] = (zlevels_stag[lev][k+1] - zlevels_stag[lev][k]);
            }
        }
    }

    for (int lev = 0; lev <= max_level; lev++) {
        stretched_dz_d[lev].resize(stretched_dz_h[lev].size());
        Gpu::copy(Gpu::hostToDevice, stretched_dz_h[lev].begin(), stretched_dz_h[lev].end(), stretched_dz_d[lev].begin());
    }
}
