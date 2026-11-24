#include <ERF_FloodGrid.H>

using namespace amrex;

FloodGrid
create_flood_grid (const BoxArray& ba_atm,
                   const DistributionMapping& dm_atm,
                   const Geometry& geom_atm,
                   int C)
{
    FloodGrid fg;
    static_cast<SurfaceGrid&>(fg) = create_surface_grid(ba_atm, dm_atm, geom_atm, C, flood_coupling_names());
    return fg;
}
