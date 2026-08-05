/**
 * \file ERF_Plotfile2DUtils.H
 */
#ifndef ERF_PLOTFILE2DUTILS_H_
#define ERF_PLOTFILE2DUTILS_H_

#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>

// Narrow helper utilities for 2D plotfile request handling and message
// formatting. Keep this layer small: it should support validation and tests,
// not become the long-term home for diagnostic science.

namespace plotfile2d
{

// Native SHOC state_update preserves consumed surface fluxes for diagnostics,
// but each flux component remains available only if the corresponding host SFS
// field existed before SHOC consumed and cleared it.
AMREX_FORCE_INLINE
bool use_native_shoc_consumed_flux_source (bool native_shoc_owns_scalar_fluxes,
                                           bool native_shoc_has_consumed_flux_diagnostics,
                                           bool host_flux_field_available) noexcept
{
    return native_shoc_owns_scalar_fluxes &&
           native_shoc_has_consumed_flux_diagnostics &&
           host_flux_field_available;
}

struct PlotVariableSelection
{
    amrex::Vector<std::string> accepted;
    amrex::Vector<std::string> unavailable;
};

PlotVariableSelection select_requested_plot_variables (const amrex::Vector<std::string>& requested,
                                                       const amrex::Vector<std::string>& available);

std::string format_unavailable_2d_plot_var_warning (const std::string& parameter_name,
                                                    const std::string& unavailable_name,
                                                    const amrex::Vector<std::string>& available_names);

std::string format_plot2d_parameter_name (const std::string& pp_prefix,
                                          const std::string& parameter_name);

std::string format_2d_component_count_error (int lev, int filled, int expected);

std::string format_invalid_2d_stream_error (int which);

} // namespace plotfile2d

#endif
