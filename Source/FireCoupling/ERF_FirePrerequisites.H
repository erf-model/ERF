#ifndef ERF_FIRE_PREREQUISITES_H
#define ERF_FIRE_PREREQUISITES_H

#include <AMReX_MultiFab.H>
#include <string>
#include <AMReX_Geometry.H>
#include <ERF_SurfacePrecipitation.H>

// Forward declaration
class SurfaceLayer;
class ERF;

/**
 * @file ERF_FirePrerequisites.H
 * @brief Start-up checks of the atmosphere side of the fire coupling
 *
 * Verifies that ERF is configured so that a fire grid can be coupled to it,
 * whichever fire engine runs on that grid:
 * - MOST surface layer is active with proper boundary conditions
 * - The level the fire grid refines exists and keeps its boxes
 * - No z-direction MPI decomposition (all z-levels on same rank)
 * - The level's boxes cover one rectangle
 * - Grid divisibility by refinement factor
 *
 * Each engine checks its own inputs after these (FireEngine::verify_prerequisites).
 *
 * All checks use AMREX_ALWAYS_ASSERT_WITH_MESSAGE for actionable diagnostics.
 */

/**
 * @brief The AMR level the fire grid refines: erf.fire.anchor_level, or the
 *        finest level when it is -1 (the default).
 *
 * Aborts, naming the input, when the level is above the finest level.
 *
 * @param[in] anchor_level  erf.fire.anchor_level (-1 = the finest level)
 * @param[in] finest_level  Finest level of the run at start-up
 */
int fire_anchor_level(int anchor_level, int finest_level);

/**
 * @brief Start-up check for an engine that reads the rain of every column
 *        (FireEngine::wants_surface_precip): the microphysics scheme must
 *        provide surface precipitation accumulators.
 *
 * Aborts, naming the engine and erf.moisture_model, when the sources are
 * empty (no moisture model, Kessler_NoRain, saturation adjustment only, ...).
 * Called after the microphysics of the fire grid's level is initialised.
 *
 * @param[in] engine_name  The engine's name (FireEngine::name)
 * @param[in] sources      Microphysics::Get_Surface_Precip_Accumulation_Ptrs(lev)
 */
void verify_fire_precip_source(const std::string& engine_name,
                               const SurfacePrecipAccumulationSources& sources);

/**
 * @brief Verify the ERF configuration the fire coupling needs
 *
 * Called once at startup from FireCoupler::initialize() after the SurfaceLayer
 * is ready, before the engine's own checks. Aborts with clear diagnostic
 * messages if any requirement is not met.
 *
 * @param[in] erf              ERF instance (for grids, dmap, Geom access)
 * @param[in] lev              AMR level the fire grid refines (fire_anchor_level)
 * @param[in] surface_layer    Pointer to initialized SurfaceLayer
 * @param[in] grid_ratio       erf.fire.grid_ratio, fire cells per atmospheric cell in x and y
 *
 * Checks performed (in order):
 * 1. SurfaceLayer pointer is not null (zlo.type = surface_layer), unless
 *    needs_surface_layer is false (FireEngine::needs_surface_layer)
 * 2. 0 <= lev <= finest level
 * 3. On a level above 0: static refinement (erf.regrid_int <= 0)
 * 7. Every box of the level spans its full height (no z-decomposition)
 * 7b. The level's boxes cover one rectangle in x and y
 * 8. grid_ratio >= 1
 * 9. All boxes divisible by grid_ratio in x and y
 * 10. Distribution map size matches grid size
 * 11. Domain z, x and y indices start at 0
 */
void verify_fire_prerequisites(const ERF& erf,
                               int lev,
                               const SurfaceLayer* surface_layer,
                               int grid_ratio,
                               bool needs_surface_layer = true);

#endif
