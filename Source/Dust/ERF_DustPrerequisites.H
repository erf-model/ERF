/**
 * @file ERF_DustPrerequisites.H
 * @brief Prerequisite verification for the ERF-Dust module.
 *
 * Declares a function that verifies the ERF configuration meets requirements
 * for dust module initialization.
 *
 * References:
 *   MRF surface layer reference height:
 *     Hong, S.-Y., and H.-L. Pan (1996). Mon. Wea. Rev., 124, 2322-2339.
 *     https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
 */

#ifndef ERF_DUST_PREREQUISITES_H
#define ERF_DUST_PREREQUISITES_H

#include <ERF_DustParams.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

// Forward declarations
class SurfaceLayer;
class ERF;

/**
 * @brief Verify ERF configuration requirements for the dust module.
 *
 * Called once at startup before the DustLayer is allocated. Aborts with
 * actionable diagnostic messages if any requirement is not met.
 *
 * Checks performed (in order):
 *   1. SurfaceLayer pointer is not null (zlo.type = "surface_layer" required).
 *   2. erf.most.z0 is set (roughness length required for u* computation).
 *   3. No z-direction MPI decomposition: all boxes in ba_atm must span the
 *      full z-extent of the domain (set amrex.max_grid_size_z = domain_nz).
 *   4. grid_ratio >= 1.
 *   5. All atmospheric boxes have x,y sizes divisible by grid_ratio.
 *   6. DistributionMapping size matches BoxArray size.
 *   7. Domain z-index starts at 0.
 *   8. Domain physical height exceeds the MRF reference height z_ref
 *      (from erf.most.zref, default 10 m).
 *
 * References:
 *   MRF surface layer reference height:
 *     Hong, S.-Y., and H.-L. Pan (1996). Mon. Wea. Rev., 124, 2322-2339.
 *     https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
 *
 * @param[in] erf           ERF instance (for public AmrCore accessors).
 * @param[in] surface_layer Pointer to SurfaceLayer (may be nullptr, checked first).
 * @param[in] dust_params   Dust configuration parameters.
 */
void verify_dust_prerequisites(const ERF&          erf,
                               const SurfaceLayer* surface_layer,
                               const DustParams&   dust_params);

#endif // ERF_DUST_PREREQUISITES_H
