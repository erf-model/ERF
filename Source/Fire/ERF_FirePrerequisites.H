#ifndef ERF_FIRE_PREREQUISITES_H
#define ERF_FIRE_PREREQUISITES_H

#include <ERF_FireParams.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

// Forward declaration
class SurfaceLayer;
class ERF;

/**
 * @file ERF_FirePrerequisites.H
 * @brief Fire module prerequisite checks at startup
 *
 * Verifies that ERF is configured correctly for fire simulation:
 * - MOST surface layer is active with proper boundary conditions
 * - MOST diagnostics (u*, z0, L, wind) are computed
 * - No z-direction MPI decomposition (all z-levels on same rank)
 * - Grid divisibility by refinement factor
 * - Domain height is sufficient for wind extraction
 *
 * All checks use AMREX_ALWAYS_ASSERT_WITH_MESSAGE for actionable diagnostics.
 */

/**
 * @brief Verify fire module prerequisites in ERF configuration
 *
 * Called once at startup from ERF::init_fire_module() after SurfaceLayer is ready.
 * Aborts with clear diagnostic messages if any requirement is not met.
 *
 * @param[in] erf              ERF instance (for grids, dmap, Geom access)
 * @param[in] surface_layer    Pointer to initialized SurfaceLayer
 * @param[in] fire_params      Fire configuration parameters
 *
 * Checks performed (in order):
 * 1. Surface layer BC type is "surface_layer" on z-low domain
 * 2. SurfaceLayer pointer is not null
 * 3. u_star is computed at level 0
 * 4. z0 is set (from most.z0 parameter)
 * 5. Obukhov length L is computed
 * 6. Wind averages (u, v) are computed (from fluxes)
 * 7. All boxes in level 0 have full z-extent (no z-decomposition)
 * 8. grid_ratio >= 1
 * 9. All boxes divisible by grid_ratio in x and y
 * 10. Distribution map size matches grid size
 * 11. Domain z-index starts at 0
 * 12. Domain height > wind_ref_ht
 */
void verify_fire_prerequisites(const ERF& erf,
                               const SurfaceLayer* surface_layer,
                               const FireParams& fire_params);

#endif
