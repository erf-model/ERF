/**
 * \file ERF_Plotfile2DPrecip.cpp
 */
#include "ERF_Plotfile2DPrecip.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>

using namespace amrex;

namespace plotfile2d
{

namespace
{

bool has_rain_component (const MoistureComponentIndices& moisture_indices) noexcept
{
    return moisture_indices.qr >= 0;
}

bool has_snow_component (const MoistureComponentIndices& moisture_indices) noexcept
{
    return moisture_indices.qs >= 0;
}

bool has_graupel_component (const MoistureComponentIndices& moisture_indices) noexcept
{
    return moisture_indices.qg >= 0;
}

bool has_any_precip_component (const MoistureComponentIndices& moisture_indices) noexcept
{
    return has_rain_component(moisture_indices)
        || has_snow_component(moisture_indices)
        || has_graupel_component(moisture_indices);
}

const SurfacePrecipAccumulationSource* first_available_source (const SurfacePrecipAccumulationSources& sources) noexcept
{
    if (surface_precip_has_source(sources.total))   return &sources.total;
    if (surface_precip_has_source(sources.rain))    return &sources.rain;
    if (surface_precip_has_source(sources.snow))    return &sources.snow;
    if (surface_precip_has_source(sources.graupel)) return &sources.graupel;
    if (surface_precip_has_source(sources.hail))    return &sources.hail;
    return nullptr;
}

} // namespace

bool
is_precipitation_accumulation (DiagnosticID id) noexcept
{
    switch (id) {
    case DiagnosticID::PrecipTotalAccum:
    case DiagnosticID::PrecipRainAccum:
    case DiagnosticID::PrecipSnowAccum:
    case DiagnosticID::PrecipGraupelAccum:
    case DiagnosticID::PrecipHailAccum:
    case DiagnosticID::PrecipFrozenAccum:
        return true;
    default:
        return false;
    }
}

bool
is_precipitation_accumulation_name (const std::string& name)
{
    const auto* descriptor = find_diagnostic(name);
    return descriptor && is_precipitation_accumulation(descriptor->id);
}

bool
precipitation_diagnostic_available (DiagnosticID id,
                                    const MoistureComponentIndices& moisture_indices) noexcept
{
    const bool has_any_precip = has_any_precip_component(moisture_indices);

    switch (id) {
    case DiagnosticID::PrecipTotalAccum:
    case DiagnosticID::PrecipFrozenAccum:
        return has_any_precip;
    case DiagnosticID::PrecipRainAccum:
        return has_rain_component(moisture_indices);
    case DiagnosticID::PrecipSnowAccum:
        return has_snow_component(moisture_indices);
    case DiagnosticID::PrecipGraupelAccum:
        return has_graupel_component(moisture_indices);
    case DiagnosticID::PrecipHailAccum:
        return false;
    default:
        return false;
    }
}

SelectedSurfacePrecipAccumulationComponents
selected_precipitation_accumulation_components (const amrex::Vector<std::string>& plot_var_names,
                                                const SurfacePrecipAccumulationSources& sources)
{
    SelectedSurfacePrecipAccumulationComponents selected;

    for (int dst_comp = 0; dst_comp < static_cast<int>(plot_var_names.size()); ++dst_comp) {
        const auto* descriptor = find_diagnostic(plot_var_names[dst_comp]);
        if (!descriptor || !is_precipitation_accumulation(descriptor->id)) {
            continue;
        }

        const bool has_total = surface_precip_has_total_source(sources);
        const bool has_rain = surface_precip_has_source(sources.rain);

        const bool source_available =
            (descriptor->id == DiagnosticID::PrecipTotalAccum)
                ? surface_precip_has_any_source(sources)
                : (descriptor->id == DiagnosticID::PrecipRainAccum)
                    ? (has_rain || has_total)
                    : (descriptor->id == DiagnosticID::PrecipSnowAccum)
                        ? surface_precip_has_source(sources.snow)
                        : (descriptor->id == DiagnosticID::PrecipGraupelAccum)
                            ? surface_precip_has_source(sources.graupel)
                            : (descriptor->id == DiagnosticID::PrecipHailAccum)
                                ? surface_precip_has_source(sources.hail)
                                : (descriptor->id == DiagnosticID::PrecipFrozenAccum)
                                    ? surface_precip_has_any_frozen_source(sources)
                                    : false;

        if (!source_available) {
            continue;
        }

        selected.dst_comp[selected.n] = dst_comp;
        selected.id[selected.n] = descriptor->id;
        ++selected.n;
    }

    return selected;
}

void
fill_precipitation_accumulations (MultiFab& dst,
                                  const SurfacePrecipAccumulationSources& sources,
                                  const SelectedSurfacePrecipAccumulationComponents& selected,
                                  const int klo)
{
    if (selected.n == 0) {
        return;
    }

    AMREX_ALWAYS_ASSERT(selected.n >= 0);
    AMREX_ALWAYS_ASSERT(selected.n <= MaxSurfacePrecipAccumulationComponents);

    for (int n = 0; n < selected.n; ++n) {
        AMREX_ALWAYS_ASSERT(selected.dst_comp[n] >= 0);
        AMREX_ALWAYS_ASSERT(selected.dst_comp[n] < dst.nComp());
        dst.setVal(0.0, selected.dst_comp[n], 1, 0);
    }

    if (!surface_precip_has_any_source(sources)) {
        return;
    }

    // Normalize scheme-native precipitation accumulators to kg/m^2 at the
    // output boundary, then build the derived total and frozen fields from the
    // normalized species values.
    const SurfacePrecipAccumulationSource* anchor = first_available_source(sources);
    AMREX_ALWAYS_ASSERT(anchor != nullptr);
    const MultiFab& anchor_mf = *anchor->accum;
    const bool has_total_source = surface_precip_has_source(sources.total);
    const bool has_rain_source = surface_precip_has_source(sources.rain);
    const bool has_snow_source = surface_precip_has_source(sources.snow);
    const bool has_graupel_source = surface_precip_has_source(sources.graupel);
    const bool has_hail_source = surface_precip_has_source(sources.hail);
    const amrex::Real total_factor = sources.total.native_to_kg_m2;
    const amrex::Real rain_factor = sources.rain.native_to_kg_m2;
    const amrex::Real snow_factor = sources.snow.native_to_kg_m2;
    const amrex::Real graupel_factor = sources.graupel.native_to_kg_m2;
    const amrex::Real hail_factor = sources.hail.native_to_kg_m2;

    for (MFIter mfi(anchor_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto dst_arr = dst.array(mfi);
        const auto total_arr = has_total_source ? sources.total.accum->const_array(mfi) : Array4<const Real>{};
        const auto rain_arr = has_rain_source ? sources.rain.accum->const_array(mfi) : Array4<const Real>{};
        const auto snow_arr = has_snow_source ? sources.snow.accum->const_array(mfi) : Array4<const Real>{};
        const auto graupel_arr = has_graupel_source ? sources.graupel.accum->const_array(mfi) : Array4<const Real>{};
        const auto hail_arr = has_hail_source ? sources.hail.accum->const_array(mfi) : Array4<const Real>{};

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (k != klo) {
                return;
            }

            const amrex::Real total_source =
                has_total_source ? total_arr(i, j, k, 0) * total_factor : amrex::Real(0.0);
            const amrex::Real rain_source =
                has_rain_source ? rain_arr(i, j, k, 0) * rain_factor : amrex::Real(0.0);
            const amrex::Real snow =
                has_snow_source ? snow_arr(i, j, k, 0) * snow_factor : amrex::Real(0.0);
            const amrex::Real graupel =
                has_graupel_source ? graupel_arr(i, j, k, 0) * graupel_factor : amrex::Real(0.0);
            const amrex::Real hail =
                has_hail_source ? hail_arr(i, j, k, 0) * hail_factor : amrex::Real(0.0);
            const amrex::Real frozen = snow + graupel + hail;
            const amrex::Real total = has_total_source ? total_source : rain_source + frozen;
            // Total-plus-frozen-subset schemes derive rain as total - frozen.
            // Clamp to zero to avoid negative output from roundoff or small
            // scheme inconsistencies.
            const amrex::Real rain = has_rain_source ? rain_source
                                                     : amrex::max(amrex::Real(0.0), total - frozen);

            for (int n = 0; n < selected.n; ++n) {
                switch (selected.id[n]) {
                case DiagnosticID::PrecipTotalAccum:
                    dst_arr(i, j, 0, selected.dst_comp[n]) = total;
                    break;
                case DiagnosticID::PrecipRainAccum:
                    dst_arr(i, j, 0, selected.dst_comp[n]) = rain;
                    break;
                case DiagnosticID::PrecipSnowAccum:
                    dst_arr(i, j, 0, selected.dst_comp[n]) = snow;
                    break;
                case DiagnosticID::PrecipGraupelAccum:
                    dst_arr(i, j, 0, selected.dst_comp[n]) = graupel;
                    break;
                case DiagnosticID::PrecipHailAccum:
                    dst_arr(i, j, 0, selected.dst_comp[n]) = hail;
                    break;
                case DiagnosticID::PrecipFrozenAccum:
                    dst_arr(i, j, 0, selected.dst_comp[n]) = frozen;
                    break;
                default:
                    break;
                }
            }
        });
    }
}

} // namespace plotfile2d
