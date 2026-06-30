#include "ERF_Plotfile2DUtils.H"

#include <sstream>
#include <unordered_set>

namespace plotfile2d
{

namespace
{

std::string join_names (const amrex::Vector<std::string>& names)
{
    if (names.empty()) {
        return "(none)";
    }

    std::ostringstream os;
    for (int i = 0; i < static_cast<int>(names.size()); ++i) {
        if (i > 0) {
            os << ", ";
        }
        os << names[i];
    }
    return os.str();
}

} // namespace

PlotVariableSelection
select_requested_plot_variables (const amrex::Vector<std::string>& requested,
                                 const amrex::Vector<std::string>& available)
{
    PlotVariableSelection selection;

    std::unordered_set<std::string> requested_set(requested.begin(), requested.end());
    std::unordered_set<std::string> available_set(available.begin(), available.end());
    std::unordered_set<std::string> seen_unavailable;

    // Preserve the canonical built-in ordering so component indices remain
    // stable even when the user lists variables in a different order.
    for (const auto& name : available) {
        if (requested_set.count(name) != 0) {
            selection.accepted.push_back(name);
        }
    }

    // Report each missing name once, in the order it was requested. The caller
    // uses these names to emit warnings without duplicating messages for repeats.
    for (const auto& name : requested) {
        if (available_set.count(name) == 0) {
            if (seen_unavailable.insert(name).second) {
                selection.unavailable.push_back(name);
            }
        }
    }

    return selection;
}

std::string
format_unavailable_2d_plot_var_warning (const std::string& parameter_name,
                                        const std::string& unavailable_name,
                                        const amrex::Vector<std::string>& available_names)
{
    std::ostringstream os;
    os << "WARNING: Requested 2D plot variable '" << unavailable_name
       << "' from '" << parameter_name
       << "' is not available and will be skipped. Available built-in 2D plot variables are: "
       << join_names(available_names) << ".";
    return os.str();
}

std::string
format_plot2d_parameter_name (const std::string& pp_prefix,
                              const std::string& parameter_name)
{
    if (pp_prefix.empty()) {
        return parameter_name;
    }

    return pp_prefix + "." + parameter_name;
}

std::string
format_2d_component_count_error (int lev, int filled, int expected)
{
    std::ostringstream os;
    os << "Write2DPlotFile internal error at level " << lev
       << ": filled " << filled
       << " components but expected " << expected
       << ". The 2D plot variable list and fill blocks are inconsistent.";
    return os.str();
}

std::string
format_invalid_2d_stream_error (int which)
{
    std::ostringstream os;
    os << "Write2DPlotFile received invalid stream index " << which
       << "; expected 1 or 2.";
    return os.str();
}

} // namespace plotfile2d
