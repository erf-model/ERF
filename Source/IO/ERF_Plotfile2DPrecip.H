/**
 * \file ERF_Plotfile2DPrecip.H
 */
#ifndef ERF_PLOTFILE2DPRECIP_H_
#define ERF_PLOTFILE2DPRECIP_H_

#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "ERF_DataStruct.H"
#include "ERF_Plotfile2DCatalog.H"
#include "ERF_SurfacePrecipitation.H"

// Scheme-aware helpers for 2D surface precipitation accumulation diagnostics.
// The output layer normalizes scheme-native precipitation accumulators to the
// public liquid-water-equivalent kg/m^2 basis before writing plotfile slabs.
// Some schemes provide explicit rain/snow/graupel sources, while others
// provide a total source plus frozen subsets. The helper layer preserves both
// contracts and derives the public fields without double-counting.

namespace plotfile2d
{

static constexpr int MaxSurfacePrecipAccumulationComponents = 6;

struct SelectedSurfacePrecipAccumulationComponents
{
    int n = 0;
    amrex::GpuArray<int, MaxSurfacePrecipAccumulationComponents> dst_comp{};
    amrex::GpuArray<DiagnosticID, MaxSurfacePrecipAccumulationComponents> id{};
};

bool is_precipitation_accumulation (DiagnosticID id) noexcept;

bool is_precipitation_accumulation_name (const std::string& name);

bool precipitation_diagnostic_available (DiagnosticID id,
                                         const MoistureComponentIndices& moisture_indices) noexcept;

SelectedSurfacePrecipAccumulationComponents
selected_precipitation_accumulation_components (const amrex::Vector<std::string>& plot_var_names,
                                                const SurfacePrecipAccumulationSources& sources);

void fill_precipitation_accumulations (amrex::MultiFab& dst,
                                       const SurfacePrecipAccumulationSources& sources,
                                       const SelectedSurfacePrecipAccumulationComponents& selected,
                                       const int klo);

} // namespace plotfile2d

#endif
