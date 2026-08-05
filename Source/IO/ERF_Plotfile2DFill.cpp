/**
 * \file ERF_Plotfile2DFill.cpp
 */
#include "ERF_Plotfile2DFill.H"

#include <AMReX_Gpu.H>

#include "ERF_Constants.H"
#include "Diagnostics/ERF_SurfaceFluxDiagnostics.H"

using namespace amrex;

namespace plotfile2d
{

namespace
{

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
is_valid_land_surface_value (Real value) noexcept
{
    // Noah-MP's -9999 fill convention and ERF's large internal sentinel are
    // both invalid. This preserves valid zero and negative flux values.
    return std::isfinite(value) && value > Real(-9990.0) &&
           value < Real(0.5) * lsm_undefined;
}

} // namespace

void
fill_component_with_value (MultiFab& dst, int dst_comp, Real value)
{
    dst.setVal(value, dst_comp, 1, 0);
}

void
fill_component_from_klevel (MultiFab& dst,
                            int dst_comp,
                            const MultiFab& src,
                            int src_k,
                            int src_comp)
{
    // Iterate over dst because it defines the 2D output component layout. The
    // source must be box-compatible with dst on the horizontal tile covered by
    // each MFIter.
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& dst_arr = dst.array(mfi);
        const auto& src_arr = src.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            dst_arr(i, j, k, dst_comp) = src_arr(i, j, src_k, src_comp);
        });
    }
}

void
fill_component_from_klevel_or_value (MultiFab& dst,
                                     int dst_comp,
                                     const MultiFab* src,
                                     int src_k,
                                     Real missing_value,
                                     int src_comp)
{
    if (src) {
        fill_component_from_klevel(dst, dst_comp, *src, src_k, src_comp);
    } else {
        fill_component_with_value(dst, dst_comp, missing_value);
    }
}

void
fill_land_surface_component_from_klevel_or_missing (MultiFab& dst,
                                                    int dst_comp,
                                                    const MultiFab* src,
                                                    int src_k,
                                                    Real missing_value)
{
    if (!src) {
        fill_component_with_value(dst, dst_comp, missing_value);
        return;
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& dst_arr = dst.array(mfi);
        const auto& src_arr = src->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real value = src_arr(i, j, src_k, 0);
            dst_arr(i, j, k, dst_comp) = is_valid_land_surface_value(value)
                ? value : missing_value;
        });
    }
}

void
fill_sensible_heat_flux_from_klevel_or_missing (MultiFab& dst,
                                                int dst_comp,
                                                const MultiFab* src,
                                                int src_k,
                                                Real missing_value)
{
    if (!src) {
        fill_component_with_value(dst, dst_comp, missing_value);
        return;
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Keep unit conversion in the surface-flux diagnostic helper so this
        // mechanical layer does not duplicate physical constants or units.
        const Box& bx = mfi.tilebox();
        const auto& dst_arr = dst.array(mfi);
        const auto& src_arr = src->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            dst_arr(i, j, k, dst_comp) =
                surface_flux_diagnostics::sensible_heat_flux_wm2_from_rhotheta_flux(
                    src_arr(i, j, src_k, 0));
        });
    }
}

void
fill_latent_heat_flux_from_klevel_or_missing (MultiFab& dst,
                                              int dst_comp,
                                              const MultiFab* src,
                                              int src_k,
                                              Real missing_value)
{
    if (!src) {
        fill_component_with_value(dst, dst_comp, missing_value);
        return;
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Keep unit conversion in the surface-flux diagnostic helper so this
        // mechanical layer does not duplicate physical constants or units.
        const Box& bx = mfi.tilebox();
        const auto& dst_arr = dst.array(mfi);
        const auto& src_arr = src->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            dst_arr(i, j, k, dst_comp) =
                surface_flux_diagnostics::latent_heat_flux_wm2_from_rhoqv_flux(
                    src_arr(i, j, src_k, 0));
        });
    }
}

} // namespace plotfile2d
