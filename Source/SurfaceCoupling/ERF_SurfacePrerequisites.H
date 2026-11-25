#ifndef ERF_SURFACE_PREREQUISITES_H
#define ERF_SURFACE_PREREQUISITES_H

#include <ERF_SurfaceCouplingNames.H>
#include <ERF_SurfacePrecipitation.H>

#include <string>

class ERF;

/**
 * @file ERF_SurfacePrerequisites.H
 * @brief Start-up checks every surface coupling makes on the level its grid
 *        refines
 */

/**
 * @brief The level the surface grid refines: anchor_level, or the finest level
 *        when it is -1. Aborts, naming erf.<coupling>.anchor_level, when it is
 *        above the finest level.
 */
int surface_anchor_level (int anchor_level, int finest_level, const SurfaceCouplingNames& names);

/**
 * @brief Aborts, naming erf.moisture_model, when no microphysics surface
 *        precipitation accumulator exists for an engine that reads it.
 *
 * @param[in] what  What the engine takes from the accumulators, for the
 *                  message ("rain", "precipitation")
 */
void verify_surface_precip_source (const std::string& engine_name,
                                   const SurfacePrecipAccumulationSources& sources,
                                   const std::string& what,
                                   const SurfaceCouplingNames& names);

/**
 * @brief The checks on the level the surface grid refines, each naming the
 *        inputs to change
 *
 * 1. The level exists (0 to the finest level).
 * 2. Above level 0, the level keeps its boxes (regrid_int <= 0, read with
 *    ERF's own ParmParse prefix): the grid is built once
 *    (verify_surface_static_refinement).
 * 3. Every box spans the level's full height: the grid pairs one surface slab
 *    with each box, and reads whole columns.
 * 4. The boxes cover one rectangle: the grid treats its edges as domain edges.
 * 5. The level's domain starts at index 0 (the map i / C + atm_lo assumes it).
 */
void verify_surface_grid_level (const ERF& erf, int lev, const SurfaceCouplingNames& names);

/**
 * @brief Check 2 of verify_surface_grid_level: above level 0 the level must
 *        keep its boxes
 *
 * @param[in] lev        The level the surface grid refines
 * @param[in] pp_prefix  The ParmParse prefix ERF reads regrid_int with
 *                       (ERF::pp_prefix: "erf", or a MultiBlock block's own)
 * @param[in] names      The coupling's names, for the message
 */
void verify_surface_static_refinement (int lev, const std::string& pp_prefix,
                                       const SurfaceCouplingNames& names);

#endif
