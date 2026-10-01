#ifndef ERF_SEB_TURBULENT_FLUX_H_
#define ERF_SEB_TURBULENT_FLUX_H_

/**
 * @file ERF_SEBTurbulentFlux.H
 * @brief Where the simplified surface energy balance takes its sensible (H) and
 *        latent (LE) heat fluxes from.
 *
 * The force-restore update removes H + LE from the ground. The air above receives
 * the flux the lower surface layer applies, which it computes from the same skin
 * temperature. Unless the two are the same flux the ground never loses what the
 * air gains, so by default the balance takes the surface layer's flux.
 *
 * Precedence, per field and per level:
 *   1. a land-surface model field of that name ("hfx" or "lh");
 *   2. the surface layer's applied flux, with seb_turbulent_flux_source =
 *      surface_layer and a flux present (a zlo surface layer that applies one,
 *      not on EB terrain; the latent flux exists only with a moisture model);
 *   3. the scalar default (seb_hfx_default or seb_lh_default).
 *
 * The surface-layer flux is the conservative flux on the lowest z face (rho
 * w'theta' and rho w'qv'), converted to W/m^2 by the same helpers as the
 * sensible_heat_flux and latent_heat_flux 2D outputs, so the two cannot
 * disagree. Positive is away from the surface, as in the balance.
 */

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <ERF_RadStruct.H>

/** @brief Which of the three sources a field was filled from. */
enum class SEBFluxOrigin { LandSurface, SurfaceLayer, Default };

/** @brief The turbulent flux being filled: selects the unit conversion. */
enum class SEBTurbulentFlux { Sensible, Latent };

/**
 * @brief The precedence rule: land-surface field, then surface layer, then default.
 *
 * @param[in] land_surface_field_present a land-surface model exposes the field
 * @param[in] source                     erf.radiation.seb_turbulent_flux_source
 * @param[in] surface_layer_flux_present the surface layer applies this flux
 */
[[nodiscard]] inline SEBFluxOrigin
select_seb_flux_origin (bool land_surface_field_present,
                        SEBTurbulentFluxSource source,
                        bool surface_layer_flux_present) noexcept
{
    if (land_surface_field_present) { return SEBFluxOrigin::LandSurface; }
    if (source == SEBTurbulentFluxSource::SurfaceLayer && surface_layer_flux_present) {
        return SEBFluxOrigin::SurfaceLayer;
    }
    return SEBFluxOrigin::Default;
}

/**
 * @brief Fill the SEB's H or LE [W/m^2] on the valid cells of seb_flux.
 *
 * @param[out] seb_flux           2D SEB field (k = 0), laid out like the level's state
 * @param[in]  kind               sensible or latent
 * @param[in]  land_surface_field land-surface field of that name, or nullptr
 * @param[in]  surface_layer_flux conservative z-face flux the surface layer applies
 *                                (rho w'theta' or rho w'qv'), or nullptr
 * @param[in]  surface_k          z index of the surface face in surface_layer_flux
 * @param[in]  source             erf.radiation.seb_turbulent_flux_source
 * @param[in]  fallback           the scalar default [W/m^2]
 * @return the source the field was filled from
 */
SEBFluxOrigin
fill_seb_turbulent_flux (amrex::MultiFab& seb_flux,
                         SEBTurbulentFlux kind,
                         const amrex::MultiFab* land_surface_field,
                         const amrex::MultiFab* surface_layer_flux,
                         int surface_k,
                         SEBTurbulentFluxSource source,
                         amrex::Real fallback);

#endif
