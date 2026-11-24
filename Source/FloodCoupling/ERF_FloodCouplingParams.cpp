#include <ERF_FloodCouplingParams.H>
#include <ERF_FloodEngineRegistry.H>
#include <ERF_FloodGrid.H>
#include <ERF_SurfaceCouplingParams.H>
#include <AMReX_ParmParse.H>

#include <string>
#include <vector>

FloodCouplingParams
FloodCouplingParams::read ()
{
    FloodCouplingParams c;
    amrex::ParmParse pp("erf.flood");
    // enable, grid_ratio, anchor_level and engine: the inputs every surface
    // coupling reads (ERF_SurfaceCouplingParams.H)
    read_surface_coupling_common(pp, flood_coupling_names(), flood_engine_names(),
                                 c.enable, c.grid_ratio, c.anchor_level, c.engine);
    return c;
}
