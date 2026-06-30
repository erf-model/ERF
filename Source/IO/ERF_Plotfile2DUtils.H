#ifndef ERF_PLOTFILE2DUTILS_H_
#define ERF_PLOTFILE2DUTILS_H_

#include <string>

#include <AMReX_Vector.H>

// Narrow helper utilities for 2D plotfile request handling and message
// formatting. Keep this layer small: it should support validation and tests,
// not become the long-term home for diagnostic science.

namespace plotfile2d
{

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
