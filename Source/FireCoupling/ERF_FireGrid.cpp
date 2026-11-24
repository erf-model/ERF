#include <ERF_FireGrid.H>

using namespace amrex;

FireGrid
create_fire_grid(const BoxArray& ba_atm,
                 const DistributionMapping& dm_atm,
                 const Geometry& geom_atm,
                 int C)
{
    FireGrid fg;
    static_cast<SurfaceGrid&>(fg) = create_surface_grid(ba_atm, dm_atm, geom_atm, C, fire_coupling_names());
    return fg;
}
