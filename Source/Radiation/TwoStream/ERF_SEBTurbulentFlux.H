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
 * Precedence, per field, per level and per cell:
 *   1. a land-surface model field of that name ("hfx" or "lh") where it holds a
 *      valid value (is_valid_lsm_value: finite and below the lsm_undefined
 *      placeholder a land model leaves in cells it does not process);
 *   2. the surface layer's applied flux, with seb_turbulent_flux_source =
 *      surface_layer and the flux field allocated: a zlo surface layer with any
 *      diffusion or closure, not on EB terrain; the latent flux only with a
 *      moisture model. An adiabatic surface layer has the field and a zero flux,
 *      so H (or LE) is then 0, and ERF warns at start-up if the scalar default it
 *      replaces is nonzero;
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
 * The selected values go to the valid cells; the halo always holds the fallback.
 * The force-restore update, the residual diagnostic and the seb_hfx/seb_lh outputs
 * read valid cells only. A consumer that needs neighbouring values in the halo must
 * fill it first (FillBoundary).
 *
 * @param[out] seb_flux           2D SEB field (k = 0), laid out like the level's state
 * @param[in]  kind               sensible or latent
 * @param[in]  land_surface_field land-surface field of that name, or nullptr; cells
 *                                holding no valid value take the next source
 * @param[in]  surface_layer_flux conservative z-face flux the surface layer applies
 *                                (rho w'theta' or rho w'qv'), or nullptr
 * @param[in]  surface_k          z index of the surface face in surface_layer_flux
 * @param[in]  source             erf.radiation.seb_turbulent_flux_source
 * @param[in]  fallback           the scalar default [W/m^2]
 * @return the highest-precedence source the field was filled from (cells of a
 *         land-surface field without a valid value take the next one)
 */
SEBFluxOrigin
fill_seb_turbulent_flux (amrex::MultiFab& seb_flux,
                         SEBTurbulentFlux kind,
                         const amrex::MultiFab* land_surface_field,
                         const amrex::MultiFab* surface_layer_flux,
                         int surface_k,
                         SEBTurbulentFluxSource source,
                         amrex::Real fallback);


/**
 * @brief Fill one 2D surface-energy-balance input from land-surface model fields,
 *        cell by cell.
 *
 * Each valid cell of seb_field is set to scale * field + add_field where both
 * values are valid land-surface values (is_valid_lsm_value), and to fallback
 * elsewhere. A land model leaves the lsm_undefined placeholder in cells it does
 * not process (Noah-MP: open water and sea ice) and everywhere before its first
 * step, which runs after the first radiation call; those cells take the scalar
 * default rather than the placeholder. The halo always holds the fallback; only
 * valid cells take land-surface values.
 *
 * @param[out] seb_field 2D SEB field (k = 0), laid out like the land-surface fields
 * @param[in]  field     land-surface field, or nullptr
 * @param[in]  fallback  the scalar default
 * @param[in]  scale     factor on field (Noah-MP's fira is positive upward; the
 *                       balance wants the absorbed flux, so -1)
 * @param[in]  add_field second land-surface field added on top, or nullptr
 *                       (Noah-MP splits the absorbed shortwave into sav and sag)
 */
void
fill_seb_field_from_land_surface (amrex::MultiFab& seb_field,
                                  const amrex::MultiFab* field,
                                  amrex::Real fallback,
                                  amrex::Real scale = amrex::Real(1.0),
                                  const amrex::MultiFab* add_field = nullptr);

#endif
