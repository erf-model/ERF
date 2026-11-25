#include <ERF_FloodPrerequisites.H>
#include <ERF_FloodGrid.H>
#include <ERF_SurfacePrerequisites.H>

#include <string>

int
flood_anchor_level (int anchor_level, int finest_level)
{
    return surface_anchor_level(anchor_level, finest_level, flood_coupling_names());
}

void
verify_flood_precip_source (const std::string& engine_name,
                            const SurfacePrecipAccumulationSources& sources)
{
    verify_surface_precip_source(engine_name, sources, "precipitation", flood_coupling_names());
}

void
verify_flood_prerequisites (const ERF& erf, int lev)
{
    // The checks every surface coupling makes on the level its grid refines
    // (ERF_SurfacePrerequisites.H): the level exists, keeps its boxes above
    // level 0, its boxes span the full height and cover one rectangle, and its
    // domain starts at index 0.
    verify_surface_grid_level(erf, lev, flood_coupling_names());
}
