/**
 * \file ERF_Plotfile2DWaterPath.H
 */
#ifndef ERF_PLOTFILE2DWATERPATH_H_
#define ERF_PLOTFILE2DWATERPATH_H_

#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "ERF_DataStruct.H"
#include "ERF_Plotfile2DCatalog.H"

// Scheme-aware helpers for 2D column water-path diagnostics. This module maps
// active moisture species to conserved ERF state components and reduces those
// components over k. It does not define microphysics species. It uses
// solverChoice.moisture_indices as the source of truth.

namespace plotfile2d
{

static constexpr int MaxCondensedWaterPathComponents = 5;

struct WaterPathDescriptor
{
    DiagnosticID id;
    const char* name;
    int source_component;
};

struct SelectedWaterPathComponents
{
    int n = 0;
    amrex::GpuArray<int, MaxCondensedWaterPathComponents> dst_comp{};
    amrex::GpuArray<int, MaxCondensedWaterPathComponents> src_comp{};
};

// A species is available when the active moisture model exposes a conserved
// mass component for it. Number concentrations are not water paths.
amrex::Vector<WaterPathDescriptor>
active_condensed_water_path_descriptors (const SolverChoice& solver_choice);

bool is_condensed_water_path (DiagnosticID id) noexcept;

bool is_condensed_water_path_name (const std::string& name);

amrex::Vector<std::string>
available_diagnostic_names (const SolverChoice& solver_choice);

// Include the fields supplied by the active land-surface provider and the
// unified near-surface fields that have a possible source in this run.
amrex::Vector<std::string>
available_diagnostic_names (const SolverChoice& solver_choice,
                            bool has_surface_layer);

amrex::Vector<std::string>
available_diagnostic_names (const SolverChoice& solver_choice,
                            bool has_surface_layer,
                            const amrex::Vector<std::string>& active_lsm_names);

SelectedWaterPathComponents
selected_condensed_water_path_components (const amrex::Vector<std::string>& plot_var_names,
                                          const SolverChoice& solver_choice);

// This mirrors volWgtColumnSum() for several selected condensed species in one
// pass. The source components are already rho*q, so the reduction gives kg/m^2
// after multiplication by dz.
void fill_condensed_water_paths (amrex::MultiFab& dst,
                                 const amrex::MultiFab& cons,
                                 const SelectedWaterPathComponents& selected,
                                 const amrex::Geometry& geom,
                                 const amrex::MultiFab& detJ);

} // namespace plotfile2d

#endif
