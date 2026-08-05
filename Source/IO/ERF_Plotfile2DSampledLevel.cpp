/**
 * \file ERF_Plotfile2DSampledLevel.cpp
 */
#include "ERF_Plotfile2DSampledLevel.H"

#include <cctype>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>
#include <unordered_set>

#include <AMReX_ParallelDescriptor.H>

namespace plotfile2d
{

namespace
{

const char*
coordinate_units_default (SampledCoordinate coordinate) noexcept
{
    switch (coordinate) {
    case SampledCoordinate::ModelIndex: return "1";
    case SampledCoordinate::HeightMSL:   return "m";
    case SampledCoordinate::HeightAGL:   return "m";
    case SampledCoordinate::Pressure:    return "Pa";
    }
    return "";
}

const char*
coordinate_interpolation_default (SampledCoordinate coordinate) noexcept
{
    switch (coordinate) {
    case SampledCoordinate::ModelIndex: return "none";
    case SampledCoordinate::HeightMSL:   return "linear";
    case SampledCoordinate::HeightAGL:   return "linear";
    case SampledCoordinate::Pressure:    return "linear";
    }
    return "";
}

bool string_iequals (const std::string& lhs, const char* rhs) noexcept
{
    if (lhs.size() != std::char_traits<char>::length(rhs)) {
        return false;
    }

    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(lhs[i])) !=
            std::tolower(static_cast<unsigned char>(rhs[i]))) {
            return false;
        }
    }

    return true;
}

std::string trim_trailing_zeros (std::string value)
{
    auto pos = value.find('.');
    if (pos == std::string::npos) {
        return value;
    }

    while (!value.empty() && value.back() == '0') {
        value.pop_back();
    }
    if (!value.empty() && value.back() == '.') {
        value.pop_back();
    }
    if (value.empty() || value == "-") {
        return "0";
    }
    return value;
}

std::string format_numeric_tag (amrex::Real value)
{
    std::ostringstream os;
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(8);
    os << value;

    std::string tag = trim_trailing_zeros(os.str());
    for (char& c : tag) {
        if (c == '.') {
            c = 'p';
        }
    }
    if (!tag.empty() && tag.front() == '-') {
        tag.front() = 'm';
    }
    return tag;
}

std::string build_duplicate_output_error (const std::string& name)
{
    std::ostringstream os;
    os << "Duplicate sampled-level output variable name '" << name << "'";
    return os.str();
}

std::string build_missing_field_error (const std::string& level_set_name,
                                       const std::string& parameter_name)
{
    std::ostringstream os;
    os << "Sampled-level definition '" << level_set_name << "' is missing required parameter '"
       << parameter_name << "'";
    return os.str();
}

} // namespace

const char*
sampled_coordinate_to_string (SampledCoordinate coordinate) noexcept
{
    switch (coordinate) {
    case SampledCoordinate::ModelIndex: return "model_index";
    case SampledCoordinate::HeightMSL:   return "height_msl";
    case SampledCoordinate::HeightAGL:   return "height_agl";
    case SampledCoordinate::Pressure:    return "pressure";
    }
    return "unknown";
}

const char*
sampled_coordinate_tag (SampledCoordinate coordinate) noexcept
{
    switch (coordinate) {
    case SampledCoordinate::ModelIndex: return "k";
    case SampledCoordinate::HeightMSL:   return "z_msl";
    case SampledCoordinate::HeightAGL:   return "z_agl";
    case SampledCoordinate::Pressure:    return "p";
    }
    return "unknown";
}

const char*
sampled_coordinate_default_units (SampledCoordinate coordinate) noexcept
{
    return coordinate_units_default(coordinate);
}

const char*
sampled_coordinate_default_interpolation (SampledCoordinate coordinate) noexcept
{
    return coordinate_interpolation_default(coordinate);
}

const char*
sampled_interpolation_to_string (SampledInterpolation interpolation) noexcept
{
    switch (interpolation) {
    case SampledInterpolation::None:   return "none";
    case SampledInterpolation::Linear: return "linear";
    }
    return "unknown";
}

SampledCoordinate
sampled_coordinate_from_string (const std::string& value)
{
    if (string_iequals(value, "model_index")) {
        return SampledCoordinate::ModelIndex;
    }
    if (string_iequals(value, "height_msl")) {
        return SampledCoordinate::HeightMSL;
    }
    if (string_iequals(value, "height_agl")) {
        return SampledCoordinate::HeightAGL;
    }
    if (string_iequals(value, "pressure")) {
        return SampledCoordinate::Pressure;
    }

    amrex::Abort("Unknown sampled-level coordinate '" + value + "'");
    return SampledCoordinate::ModelIndex;
}

SampledInterpolation
sampled_interpolation_from_string (const std::string& value)
{
    if (string_iequals(value, "none")) {
        return SampledInterpolation::None;
    }
    if (string_iequals(value, "linear")) {
        return SampledInterpolation::Linear;
    }

    amrex::Abort("Unknown sampled-level interpolation '" + value + "'");
    return SampledInterpolation::None;
}

std::string
validate_sampled_coordinate_string (const std::string& level_set_name,
                                    const std::string& value)
{
    if (string_iequals(value, "isentropic")) {
        return sampled_level_error_prefix(level_set_name, "coordinate") +
               "isentropic output needs a crossing policy";
    }

    if (string_iequals(value, "model_index") ||
        string_iequals(value, "height_msl") ||
        string_iequals(value, "height_agl") ||
        string_iequals(value, "pressure")) {
        return {};
    }

    return "Unknown sampled-level coordinate '" + value + "'";
}

std::string
sampled_level_error_prefix (const std::string& level_set_name,
                            const std::string& parameter_name)
{
    std::ostringstream os;
    os << "Sampled-level definition '" << level_set_name << "' parameter '"
       << parameter_name << "': ";
    return os.str();
}

std::string
sampled_level_value_tag (amrex::Real value, const std::string& units)
{
    std::string tag = format_numeric_tag(value);
    if (!units.empty() && units != "1") {
        tag += units;
    }
    return tag;
}

std::string
sampled_output_name (const std::string& field_name,
                     SampledCoordinate coordinate,
                     amrex::Real value,
                     const std::string& units)
{
    std::ostringstream os;
    os << field_name << "_" << sampled_coordinate_tag(coordinate) << "_"
       << sampled_level_value_tag(value, units);
    return os.str();
}

std::string
validate_sampled_level_definition (const SampledLevelDefinition& level_set)
{
    if (level_set.name.empty()) {
        return "Sampled-level definition name is empty";
    }

    if (level_set.fields.empty()) {
        return sampled_level_error_prefix(level_set.name, "fields") + "missing fields";
    }

    if (level_set.values.empty()) {
        return sampled_level_error_prefix(level_set.name, "values") + "missing values";
    }

    const char* default_units = sampled_coordinate_default_units(level_set.coordinate);
    if (level_set.units.empty() || level_set.units == default_units) {
        // valid
    } else if (level_set.coordinate == SampledCoordinate::Pressure &&
               (level_set.units == "Pa" || level_set.units == "hPa")) {
        // valid
    } else {
        return sampled_level_error_prefix(level_set.name, "units") +
               "unsupported units '" + level_set.units + "'";
    }

    const char* default_interp = sampled_coordinate_default_interpolation(level_set.coordinate);
    if (level_set.interpolation == SampledInterpolation::None) {
        if (std::string(default_interp) != "none") {
            return sampled_level_error_prefix(level_set.name, "interpolation") +
                   "unsupported interpolation 'none' for coordinate '" +
                   sampled_coordinate_to_string(level_set.coordinate) + "'";
        }
    } else if (level_set.interpolation == SampledInterpolation::Linear) {
        if (std::string(default_interp) == "none" && level_set.coordinate == SampledCoordinate::ModelIndex) {
            return sampled_level_error_prefix(level_set.name, "interpolation") +
                   "unsupported interpolation 'linear' for coordinate 'model_index'";
        }
    }

    if (level_set.coordinate == SampledCoordinate::ModelIndex) {
        for (const auto value : level_set.values) {
            const auto rounded = std::nearbyint(value);
            if (std::abs(value - rounded) > amrex::Real(1.0e-12)) {
                return sampled_level_error_prefix(level_set.name, "values") +
                       "non-integer model_index value '" + std::to_string(value) + "'";
            }
        }
    }

    std::unordered_set<std::string> unique_fields(level_set.fields.begin(), level_set.fields.end());
    const auto n_unique_fields    = static_cast<int>(unique_fields.size());
    const auto n_requested_fields = static_cast<int>(level_set.fields.size());
    if (n_unique_fields != n_requested_fields) {
        return sampled_level_error_prefix(level_set.name, "fields") +
               build_duplicate_output_error(level_set.name);
    }

    return {};
}

SampledLevelDefinition
parse_sampled_level_definition (const std::string& level_set_name,
                                const std::string& pp_prefix)
{
    SampledLevelDefinition level_set;
    level_set.name = level_set_name;

    amrex::ParmParse pp(pp_prefix);
    const std::string base = "plot2d.level_set." + level_set_name + ".";

    std::string coordinate_value;
    if (!pp.query((base + "coordinate").c_str(), coordinate_value)) {
        amrex::Abort(build_missing_field_error(level_set_name, "coordinate"));
    }
    const std::string coordinate_error =
        validate_sampled_coordinate_string(level_set_name, coordinate_value);
    if (!coordinate_error.empty()) {
        amrex::Abort(coordinate_error);
    }
    level_set.coordinate = sampled_coordinate_from_string(coordinate_value);

    if (pp.query((base + "units").c_str(), level_set.units) == 0) {
        level_set.units = sampled_coordinate_default_units(level_set.coordinate);
    }

    std::string interpolation_value;
    if (pp.query((base + "interpolation").c_str(), interpolation_value) == 0) {
        level_set.interpolation =
            sampled_interpolation_from_string(sampled_coordinate_default_interpolation(level_set.coordinate));
    } else {
        level_set.interpolation = sampled_interpolation_from_string(interpolation_value);
    }

    if (pp.query((base + "missing_value").c_str(), level_set.missing_value) == 0) {
        level_set.missing_value = amrex::Real(-999.0);
    }

    int n_values = pp.countval((base + "values").c_str());
    if (n_values <= 0) {
        amrex::Abort(build_missing_field_error(level_set_name, "values"));
    }
    level_set.values.resize(n_values);
    pp.getarr((base + "values").c_str(), level_set.values, 0, n_values);

    int n_fields = pp.countval((base + "fields").c_str());
    if (n_fields <= 0) {
        amrex::Abort(build_missing_field_error(level_set_name, "fields"));
    }
    level_set.fields.resize(n_fields);
    pp.getarr((base + "fields").c_str(), level_set.fields, 0, n_fields);

    const std::string error = validate_sampled_level_definition(level_set);
    if (!error.empty()) {
        amrex::Abort(error);
    }

    return level_set;
}

amrex::Vector<std::string>
parse_requested_sampled_level_sets (const std::string& pp_prefix,
                                    int which)
{
    amrex::Vector<std::string> names;
    amrex::ParmParse pp(pp_prefix);

    const std::string param_name = "plot2d_level_sets_" + std::to_string(which);
    if (!pp.contains(param_name.c_str())) {
        return names;
    }

    const int n_sets = pp.countval(param_name.c_str());
    names.resize(n_sets);
    pp.getarr(param_name.c_str(), names, 0, n_sets);
    return names;
}

amrex::Vector<Plotfile2DOutputDescriptor>
build_sampled_level_output_descriptors_from_definitions (
    const amrex::Vector<SampledLevelDefinition>& level_sets,
    const amrex::Vector<std::string>& static_plot_vars,
    const SolverChoice& solver_choice)
{
    amrex::Vector<Plotfile2DOutputDescriptor> descriptors;
    descriptors.reserve(static_plot_vars.size());

    std::unordered_set<std::string> output_names;

    for (const auto& name : static_plot_vars) {
        const auto* static_descriptor = find_diagnostic(name);
        if (static_descriptor == nullptr) {
            static_descriptor = find_dynamic_soil_diagnostic(name);
        }
        if (static_descriptor == nullptr) {
            amrex::Abort("Unknown built-in 2D plotfile diagnostic '" + name + "'");
        }

        Plotfile2DOutputDescriptor descriptor;
        descriptor.name = static_descriptor->name;
        descriptor.long_name = static_descriptor->long_name;
        descriptor.units = static_descriptor->units;
        descriptor.category = static_descriptor->category;
        descriptor.missing_policy = static_descriptor->missing_policy;
        descriptor.missing_value =
            (descriptor.missing_policy == MissingPolicy::FillZeroWhenUnavailable) ? amrex::Real(0.0)
                                                                                  : amrex::Real(-999.0);
        descriptor.static_diagnostic = static_descriptor;

        if (!output_names.insert(descriptor.name).second) {
            amrex::Abort(build_duplicate_output_error(descriptor.name));
        }

        descriptors.push_back(std::move(descriptor));
    }

    for (const auto& level_set : level_sets) {
        const auto field_selection = select_requested_sampled_fields(level_set.fields, solver_choice);
        if (!field_selection.unavailable.empty()) {
            if (amrex::ParallelDescriptor::IOProcessor()) {
                for (const auto& field_name : field_selection.unavailable) {
                    amrex::Warning(sampled_level_error_prefix(level_set.name, "fields") +
                                   "skipping unavailable field '" + field_name + "'");
                }
            }
        }

        if (field_selection.accepted.empty()) {
            amrex::Abort(sampled_level_error_prefix(level_set.name, "fields") +
                         "all requested fields unavailable");
        }

        const std::string coordinate_units = level_set.units.empty()
            ? sampled_coordinate_default_units(level_set.coordinate)
            : level_set.units;
        const std::string interpolation_name = sampled_interpolation_to_string(level_set.interpolation);
        const std::string canonical_units =
            (level_set.coordinate == SampledCoordinate::Pressure && coordinate_units == "hPa") ? "Pa"
                                                                                               : coordinate_units;

        for (const auto& value : level_set.values) {
            amrex::Real canonical_value = value;
            if (level_set.coordinate == SampledCoordinate::Pressure &&
                coordinate_units == "hPa") {
                canonical_value = value * amrex::Real(100.0);
            }

            for (const auto& field : field_selection.accepted) {
                Plotfile2DOutputDescriptor descriptor;
                descriptor.name = sampled_output_name(field.name, level_set.coordinate, value, coordinate_units);
                descriptor.long_name = std::string(field.long_name) + " sampled on " +
                                       sampled_coordinate_to_string(level_set.coordinate) + " levels";
                descriptor.units = field.units;
                descriptor.category = DiagnosticCategory::SampledLevel;
                descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
                descriptor.missing_value = level_set.missing_value;
                descriptor.sampled_level = SampledLevelMetadata{
                    level_set.name,
                    field.name,
                    SampledVerticalCoordinateMetadata{
                        sampled_coordinate_to_string(level_set.coordinate),
                        value,
                        coordinate_units,
                        canonical_value,
                        canonical_units,
                        interpolation_name
                    }
                };

                if (!output_names.insert(descriptor.name).second) {
                    amrex::Abort(build_duplicate_output_error(descriptor.name));
                }

                descriptors.push_back(std::move(descriptor));
            }
        }
    }

    return descriptors;
}

amrex::Vector<Plotfile2DOutputDescriptor>
build_sampled_level_output_descriptors (const std::string& pp_prefix,
                                        int which,
                                        const amrex::Vector<std::string>& static_plot_vars,
                                        const SolverChoice& solver_choice)
{
    amrex::Vector<SampledLevelDefinition> level_sets;
    for (const auto& level_set_name : parse_requested_sampled_level_sets(pp_prefix, which)) {
        level_sets.push_back(parse_sampled_level_definition(level_set_name, pp_prefix));
    }

    return build_sampled_level_output_descriptors_from_definitions(level_sets,
                                                                   static_plot_vars,
                                                                   solver_choice);
}

} // namespace plotfile2d
