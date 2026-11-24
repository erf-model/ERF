#ifndef ERF_FLOOD_PREREQUISITES_H
#define ERF_FLOOD_PREREQUISITES_H

#include <string>
#include <ERF_SurfacePrecipitation.H>

class ERF;

/**
 * @file ERF_FloodPrerequisites.H
 * @brief Start-up checks of the atmosphere side of the flood coupling
 *
 * Verifies that ERF is configured so that a flood grid can be coupled to it,
 * whichever flood engine runs on that grid. Each engine checks its own inputs
 * after these (FloodEngine::verify_prerequisites). Every check aborts with a
 * message that names the input to change.
 */

/**
 * @brief The AMR level the flood grid refines: erf.flood.anchor_level, or the
 *        finest level when it is -1 (the default).
 *
 * Aborts, naming the input, when the level is above the finest level.
 */
int flood_anchor_level (int anchor_level, int finest_level);

/**
 * @brief Start-up check for an engine that reads the precipitation
 *        (FloodEngine::wants_surface_precip): the microphysics scheme must
 *        provide surface precipitation accumulators.
 *
 * Aborts, naming the engine and erf.moisture_model, when the sources are empty.
 *
 * @param[in] engine_name  The engine's name (FloodEngine::name)
 * @param[in] sources      Microphysics::Get_Surface_Precip_Accumulation_Ptrs(lev)
 */
void verify_flood_precip_source (const std::string& engine_name,
                                 const SurfacePrecipAccumulationSources& sources);

/**
 * @brief Verify the ERF configuration the flood coupling needs
 *
 * Checks, in order:
 * 1. 0 <= lev <= finest level
 * 2. On a level above 0: static refinement (erf.regrid_int <= 0)
 * 3. Every box of the level spans its full height (no z decomposition), so
 *    the level's k = 0 slabs are one box per atmospheric box
 * 4. The level's boxes cover one rectangle in x and y
 * 5. The domain starts at index 0 in x, y and z
 */
void verify_flood_prerequisites (const ERF& erf, int lev);

#endif
