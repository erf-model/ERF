#include <ERF_SEBTurbulentFlux.H>
#include <ERF_Plotfile2DFill.H>
#include <ERF_Constants.H>

#include <AMReX_Gpu.H>

using namespace amrex;

namespace {

// Overwrite the valid cells of dst where src (plus add, when given) holds a valid
// land-surface value with scale * src + add; leave the others as they are.
void
overwrite_valid_land_surface_cells (MultiFab& dst,
                                    const MultiFab& src,
                                    Real scale,
                                    const MultiFab* add)
{
    // The kernel reads src and add through dst's MFIter, so the boxes themselves (and
    // their index type), not only their number, must agree.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src.boxArray() == dst.boxArray() &&
        src.DistributionMap() == dst.DistributionMap() &&
        (add == nullptr || (add->boxArray() == dst.boxArray() &&
                            add->DistributionMap() == dst.DistributionMap())),
        "SEB land-surface fill: the land-surface fields and the SEB field are laid out differently");
    const bool has_add = (add != nullptr);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto dst_arr = dst.array(mfi);
        const auto src_arr = src.const_array(mfi);
        const auto add_arr = has_add ? add->const_array(mfi) : Array4<const Real>{};
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real value = src_arr(i, j, k);
            const Real extra = has_add ? add_arr(i, j, k) : Real(0.0);
            if (is_valid_lsm_value(value) && is_valid_lsm_value(extra)) {
                dst_arr(i, j, k) = scale * value + extra;
            }
        });
    }
}

} // namespace

void
fill_seb_field_from_land_surface (MultiFab& seb_field,
                                  const MultiFab* field,
                                  Real fallback,
                                  Real scale,
                                  const MultiFab* add_field)
{
    seb_field.setVal(fallback);
    if (field != nullptr) {
        overwrite_valid_land_surface_cells(seb_field, *field, scale, add_field);
    }
}

SEBFluxOrigin
fill_seb_turbulent_flux (MultiFab& seb_flux,
                         SEBTurbulentFlux kind,
                         const MultiFab* land_surface_field,
                         const MultiFab* surface_layer_flux,
                         int surface_k,
                         SEBTurbulentFluxSource source,
                         Real fallback)
{
    // The fallback fills every cell, halo included; the surface layer's flux then
    // overwrites the valid cells when it is the source, and the land model's valid
    // values overwrite those in turn.
    seb_flux.setVal(fallback);
    const SEBFluxOrigin below = select_seb_flux_origin(false, source,
                                                       surface_layer_flux != nullptr);
    if (below == SEBFluxOrigin::SurfaceLayer) {
        // The 2D-output conversions, so the balance's H and LE are exactly the
        // sensible_heat_flux and latent_heat_flux a plotfile reports.
        if (kind == SEBTurbulentFlux::Sensible) {
            plotfile2d::fill_sensible_heat_flux_from_klevel_or_missing(
                seb_flux, 0, surface_layer_flux, surface_k, fallback);
        } else {
            plotfile2d::fill_latent_heat_flux_from_klevel_or_missing(
                seb_flux, 0, surface_layer_flux, surface_k, fallback);
        }
    }
    if (land_surface_field != nullptr) {
        overwrite_valid_land_surface_cells(seb_flux, *land_surface_field, Real(1.0), nullptr);
    }
    return select_seb_flux_origin(land_surface_field != nullptr, source,
                                  surface_layer_flux != nullptr);
}
