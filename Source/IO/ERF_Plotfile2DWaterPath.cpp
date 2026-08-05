/**
 * \file ERF_Plotfile2DWaterPath.cpp
 */
#include "ERF_Plotfile2DWaterPath.H"
#include "ERF_Plotfile2DPrecip.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <algorithm>

#include "ERF_IndexDefines.H"

using namespace amrex;

namespace plotfile2d
{

namespace
{

struct CondensedWaterPathSpec
{
    DiagnosticID id;
    const char* name;
};

constexpr amrex::GpuArray<CondensedWaterPathSpec, MaxCondensedWaterPathComponents> condensed_specs{{
    {DiagnosticID::IntegratedQc, "integrated_qc"},
    {DiagnosticID::IntegratedQi, "integrated_qi"},
    {DiagnosticID::IntegratedQr, "integrated_qr"},
    {DiagnosticID::IntegratedQs, "integrated_qs"},
    {DiagnosticID::IntegratedQg, "integrated_qg"},
}};

int source_component_for (DiagnosticID id, const MoistureComponentIndices& moisture_indices) noexcept
{
    switch (id) {
    case DiagnosticID::IntegratedQc: return moisture_indices.qc;
    case DiagnosticID::IntegratedQi: return moisture_indices.qi;
    case DiagnosticID::IntegratedQr: return moisture_indices.qr;
    case DiagnosticID::IntegratedQs: return moisture_indices.qs;
    case DiagnosticID::IntegratedQg: return moisture_indices.qg;
    default:                         return -1;
    }
}

} // namespace

bool
is_condensed_water_path (DiagnosticID id) noexcept
{
    switch (id) {
    case DiagnosticID::IntegratedQc:
    case DiagnosticID::IntegratedQi:
    case DiagnosticID::IntegratedQr:
    case DiagnosticID::IntegratedQs:
    case DiagnosticID::IntegratedQg:
        return true;
    default:
        return false;
    }
}

amrex::Vector<WaterPathDescriptor>
active_condensed_water_path_descriptors (const SolverChoice& solver_choice)
{
    amrex::Vector<WaterPathDescriptor> descriptors;
    descriptors.reserve(MaxCondensedWaterPathComponents);

    for (const auto& spec : condensed_specs) {
        const int source_component = source_component_for(spec.id, solver_choice.moisture_indices);
        if (source_component >= 0) {
            descriptors.push_back({spec.id, spec.name, source_component});
        }
    }

    return descriptors;
}

bool
is_condensed_water_path_name (const std::string& name)
{
    const auto* descriptor = find_diagnostic(name);
    return descriptor && is_condensed_water_path(descriptor->id);
}

bool is_noahmp_active (const SolverChoice& solver_choice) noexcept
{
#ifdef ERF_USE_NOAHMP
    return solver_choice.lsm_type == LandSurfaceType::NOAHMP;
#else
    amrex::ignore_unused(solver_choice);
    return false;
#endif
}

bool is_land_surface_provider_field (DiagnosticID id) noexcept
{
    switch (id) {
    case DiagnosticID::LandSurfaceTsfC:
    case DiagnosticID::LandSurfaceEmissivity:
    case DiagnosticID::LandSurfaceAlbDirVis:
    case DiagnosticID::LandSurfaceAlbDirNir:
    case DiagnosticID::LandSurfaceAlbDifVis:
    case DiagnosticID::LandSurfaceAlbDifNir:
    case DiagnosticID::LandSurfaceCosZenith:
    case DiagnosticID::LandSurfaceSwFluxDn:
    case DiagnosticID::LandSurfaceSwFluxDnDirVis:
    case DiagnosticID::LandSurfaceSwFluxDnDirNir:
    case DiagnosticID::LandSurfaceSwFluxDnDifVis:
    case DiagnosticID::LandSurfaceSwFluxDnDifNir:
    case DiagnosticID::LandSurfaceLwFluxDn:
    case DiagnosticID::LandSurfaceGrdflx:
    case DiagnosticID::LandSurfaceFira:
    case DiagnosticID::LandSurfaceSav:
    case DiagnosticID::LandSurfaceSag:
    case DiagnosticID::LandSurfaceAlbedo:
    case DiagnosticID::LandSurfaceSfcrunoff:
    case DiagnosticID::LandSurfaceUdrunoff:
    case DiagnosticID::NoahmpTemperature2mVegetated:
    case DiagnosticID::NoahmpTemperature2mBare:
    case DiagnosticID::NoahmpWaterVaporMixingRatio2mVegetated:
    case DiagnosticID::NoahmpWaterVaporMixingRatio2mBare:
    case DiagnosticID::NoahmpVegetationFraction:
        return true;
    default:
        return false;
    }
}

bool active_lsm_contains (const amrex::Vector<std::string>& active_lsm_names,
                          const char* name)
{
    return std::find(active_lsm_names.begin(), active_lsm_names.end(), name) !=
           active_lsm_names.end();
}

amrex::Vector<std::string>
available_diagnostic_names (const SolverChoice& solver_choice)
{
    return available_diagnostic_names(solver_choice, true);
}

amrex::Vector<std::string>
available_diagnostic_names (const SolverChoice& solver_choice,
                            bool has_surface_layer)
{
    return available_diagnostic_names(solver_choice, has_surface_layer,
                                      amrex::Vector<std::string>{});
}

amrex::Vector<std::string>
available_diagnostic_names (const SolverChoice& solver_choice,
                            bool has_surface_layer,
                            const amrex::Vector<std::string>& active_lsm_names)
{
    amrex::Vector<std::string> names;
    names.reserve(diagnostic_catalog().size());
    const bool has_noahmp = is_noahmp_active(solver_choice);
    const bool has_moisture = solver_choice.moisture_type != MoistureType::None;

    for (const auto& descriptor : diagnostic_catalog()) {
        if (is_land_surface_provider_field(descriptor.id) &&
            ((!active_lsm_names.empty() &&
              !active_lsm_contains(active_lsm_names, descriptor.name)) ||
             (active_lsm_names.empty() && !has_noahmp))) {
            continue;
        }
        if (descriptor.id == DiagnosticID::Temperature2m && !(has_noahmp || has_surface_layer)) {
            continue;
        }
        if ((descriptor.id == DiagnosticID::WaterVaporMixingRatio2m ||
             descriptor.id == DiagnosticID::NearSurfaceDiagnosticSource) &&
            !(has_moisture && (has_noahmp || has_surface_layer))) {
            if (descriptor.id == DiagnosticID::NearSurfaceDiagnosticSource &&
                (has_noahmp || has_surface_layer)) {
                names.push_back(descriptor.name);
            }
            continue;
        }
        if (is_condensed_water_path(descriptor.id)) {
            if (source_component_for(descriptor.id, solver_choice.moisture_indices) >= 0) {
                names.push_back(descriptor.name);
            }
            continue;
        }

        if (is_precipitation_accumulation(descriptor.id)) {
            if (precipitation_diagnostic_available(descriptor.id, solver_choice.moisture_indices)) {
                names.push_back(descriptor.name);
            }
            continue;
        }

        names.push_back(descriptor.name);
    }

    if (has_noahmp || !active_lsm_names.empty()) {
        amrex::ParmParse pp("erf");
        int nsoil = 4;
        const auto soil_names = active_lsm_names.empty()
            ? (pp.query("lsm_nsoil", nsoil), dynamic_soil_diagnostic_names(nsoil))
            : dynamic_soil_diagnostic_names(active_lsm_names);
        names.insert(names.end(), soil_names.begin(), soil_names.end());
    }

    return names;
}

SelectedWaterPathComponents
selected_condensed_water_path_components (const amrex::Vector<std::string>& plot_var_names,
                                          const SolverChoice& solver_choice)
{
    SelectedWaterPathComponents selected;

    for (int dst_comp = 0; dst_comp < static_cast<int>(plot_var_names.size()); ++dst_comp) {
        const auto* descriptor = find_diagnostic(plot_var_names[dst_comp]);
        if (!descriptor || !is_condensed_water_path(descriptor->id)) {
            continue;
        }

        const int src_comp = source_component_for(descriptor->id, solver_choice.moisture_indices);
        if (src_comp < 0) {
            continue;
        }

        selected.dst_comp[selected.n] = dst_comp;
        selected.src_comp[selected.n] = src_comp;
        ++selected.n;
    }

    return selected;
}

void
fill_condensed_water_paths (MultiFab& dst,
                            const MultiFab& cons,
                            const SelectedWaterPathComponents& selected,
                            const Geometry& geom,
                            const MultiFab& detJ)
{
    if (selected.n == 0) {
        return;
    }

    // Validate selected components on the host before device kernels use them.
    AMREX_ALWAYS_ASSERT(selected.n >= 0);
    AMREX_ALWAYS_ASSERT(selected.n <= MaxCondensedWaterPathComponents);
    for (int n = 0; n < selected.n; ++n) {
        AMREX_ALWAYS_ASSERT(selected.dst_comp[n] >= 0);
        AMREX_ALWAYS_ASSERT(selected.dst_comp[n] < dst.nComp());
        AMREX_ALWAYS_ASSERT(selected.src_comp[n] >= 0);
        AMREX_ALWAYS_ASSERT(selected.src_comp[n] < cons.nComp());
    }

    for (int n = 0; n < selected.n; ++n) {
        dst.setVal(0., selected.dst_comp[n], 1, 0);
    }

    const auto& dx = geom.CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto dst_arr = dst.array(mfi);
        const auto src_arr = cons.const_array(mfi);

        if (SolverChoice::mesh_type == MeshType::ConstantDz) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int n = 0; n < selected.n; ++n) {
                    amrex::HostDevice::Atomic::Add(
                        &dst_arr(i, j, 0, selected.dst_comp[n]),
                        src_arr(i, j, k, selected.src_comp[n]));
                }
            });
        } else {
            const auto& detJ_arr = detJ.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real metric = detJ_arr(i, j, k);
                for (int n = 0; n < selected.n; ++n) {
                    amrex::HostDevice::Atomic::Add(
                        &dst_arr(i, j, 0, selected.dst_comp[n]),
                        src_arr(i, j, k, selected.src_comp[n]) * metric);
                }
            });
        }
    }

    for (int n = 0; n < selected.n; ++n) {
        dst.mult(dx[2], selected.dst_comp[n], 1, 0);
    }
}

} // namespace plotfile2d
