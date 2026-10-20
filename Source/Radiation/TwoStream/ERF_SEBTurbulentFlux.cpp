#include <ERF_SEBTurbulentFlux.H>
#include <ERF_Plotfile2DFill.H>

using namespace amrex;

SEBFluxOrigin
fill_seb_turbulent_flux (MultiFab& seb_flux,
                         SEBTurbulentFlux kind,
                         const MultiFab* land_surface_field,
                         const MultiFab* surface_layer_flux,
                         int surface_k,
                         SEBTurbulentFluxSource source,
                         Real fallback)
{
    const SEBFluxOrigin origin = select_seb_flux_origin(land_surface_field != nullptr,
                                                        source,
                                                        surface_layer_flux != nullptr);
    switch (origin) {
    case SEBFluxOrigin::LandSurface:
        MultiFab::Copy(seb_flux, *land_surface_field, 0, 0, 1, 0);
        break;
    case SEBFluxOrigin::SurfaceLayer:
        // The 2D-output conversions, so the balance's H and LE are exactly the
        // sensible_heat_flux and latent_heat_flux a plotfile reports.
        if (kind == SEBTurbulentFlux::Sensible) {
            plotfile2d::fill_sensible_heat_flux_from_klevel_or_missing(
                seb_flux, 0, surface_layer_flux, surface_k, fallback);
        } else {
            plotfile2d::fill_latent_heat_flux_from_klevel_or_missing(
                seb_flux, 0, surface_layer_flux, surface_k, fallback);
        }
        break;
    case SEBFluxOrigin::Default:
        seb_flux.setVal(fallback);
        break;
    }
    return origin;
}
