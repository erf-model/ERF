/**
 * \file ERF_Plotfile2DMetadata.cpp
 */
#include "ERF_Plotfile2DMetadata.H"

#include <fstream>
#include <sstream>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

namespace plotfile2d
{

namespace
{

void append_escaped_codepoint (std::string& out, unsigned char c)
{
    static constexpr char hex_digits[] = "0123456789abcdef";

    out += "\\u00";
    out += hex_digits[(c >> 4) & 0x0f];
    out += hex_digits[c & 0x0f];
}

void append_json_string (std::ostringstream& os, const std::string& value)
{
    os << '\"' << escape_json_string(value) << '\"';
}

void append_sampled_level_metadata (std::ostringstream& os,
                                    const SampledLevelMetadata& metadata)
{
    os << "      \"source_field\": ";
    append_json_string(os, metadata.source_field);
    os << ",\n";
    os << "      \"vertical_coordinate\": {\n";
    os << "        \"type\": ";
    append_json_string(os, metadata.vertical_coordinate.type);
    os << ",\n";
    os << "        \"value\": " << metadata.vertical_coordinate.value << ",\n";
    os << "        \"units\": ";
    append_json_string(os, metadata.vertical_coordinate.units);
    os << ",\n";
    os << "        \"canonical_value\": " << metadata.vertical_coordinate.canonical_value << ",\n";
    os << "        \"canonical_units\": ";
    append_json_string(os, metadata.vertical_coordinate.canonical_units);
    os << ",\n";
    os << "        \"interpolation\": ";
    append_json_string(os, metadata.vertical_coordinate.interpolation);
    os << "\n";
    os << "      }\n";
}

void append_variable_record (std::ostringstream& os,
                             const Plotfile2DOutputDescriptor& descriptor,
                             int component_index,
                             bool is_last)
{
    os << "    {\n";
    os << "      \"component_index\": " << component_index << ",\n";
    os << "      \"name\": ";
    append_json_string(os, descriptor.name);
    os << ",\n";
    os << "      \"long_name\": ";
    append_json_string(os, descriptor.long_name);
    os << ",\n";
    os << "      \"units\": ";
    append_json_string(os, descriptor.units);
    os << ",\n";
    os << "      \"category\": ";
    append_json_string(os, diagnostic_category_to_string(descriptor.category));
    os << ",\n";
    os << "      \"missing_policy\": ";
    append_json_string(os, missing_policy_to_string(descriptor.missing_policy));
    os << ",\n";
    os << "      \"missing_value\": ";
    if (descriptor.sampled_level) {
        os << descriptor.missing_value;
    } else {
        os << missing_value_json(descriptor.missing_policy);
    }
    if (descriptor.sampled_level) {
        os << ",\n";
        append_sampled_level_metadata(os, *descriptor.sampled_level);
    } else {
        os << "\n";
    }
    os << "    }";
    if (!is_last) {
        os << ",";
    }
    os << "\n";
}

} // namespace

const char*
diagnostic_category_to_string (DiagnosticCategory category) noexcept
{
    switch (category) {
    case DiagnosticCategory::Geometry:        return "Geometry";
    case DiagnosticCategory::SurfaceLayer:     return "SurfaceLayer";
    case DiagnosticCategory::Radiation:        return "Radiation";
    case DiagnosticCategory::SurfaceFlux:      return "SurfaceFlux";
    case DiagnosticCategory::PBL:              return "PBL";
    case DiagnosticCategory::SurfaceState:     return "SurfaceState";
    case DiagnosticCategory::Precipitation:    return "Precipitation";
    case DiagnosticCategory::ColumnIntegral:   return "ColumnIntegral";
    case DiagnosticCategory::LandSurface:      return "LandSurface";
    case DiagnosticCategory::SampledLevel:     return "SampledLevel";
    }

    amrex::Abort("Unhandled DiagnosticCategory in 2D metadata writer");
    return "";
}

const char*
missing_policy_to_string (MissingPolicy policy) noexcept
{
    switch (policy) {
    case MissingPolicy::AlwaysAvailable:            return "AlwaysAvailable";
    case MissingPolicy::FillZeroWhenUnavailable:     return "FillZeroWhenUnavailable";
    case MissingPolicy::FillMinus999WhenUnavailable: return "FillMinus999WhenUnavailable";
    }

    amrex::Abort("Unhandled MissingPolicy in 2D metadata writer");
    return "";
}

std::string
missing_value_json (MissingPolicy policy)
{
    switch (policy) {
    case MissingPolicy::AlwaysAvailable:            return "null";
    case MissingPolicy::FillZeroWhenUnavailable:     return "0";
    case MissingPolicy::FillMinus999WhenUnavailable: return "-999";
    }

    amrex::Abort("Unhandled MissingPolicy in 2D metadata writer");
    return {};
}

std::string
escape_json_string (const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size() + 8);

    for (unsigned char c : value) {
        switch (c) {
        case '\"': escaped += "\\\""; break;
        case '\\': escaped += "\\\\"; break;
        case '\n': escaped += "\\n"; break;
        case '\t': escaped += "\\t"; break;
        case '\r': escaped += "\\r"; break;
        default:
            if (c < 0x20) {
                append_escaped_codepoint(escaped, c);
            } else {
                escaped.push_back(static_cast<char>(c));
            }
            break;
        }
    }

    return escaped;
}

std::string
metadata_json_filename (const std::string& plotfilename)
{
    return plotfilename + "/2DMetadata.json";
}

std::string
format_2d_metadata_json (const amrex::Vector<std::string>& varnames)
{
    amrex::Vector<Plotfile2DOutputDescriptor> descriptors;
    descriptors.reserve(varnames.size());

    for (const auto& name : varnames) {
        const auto* static_descriptor = find_diagnostic(name);
        if (static_descriptor == nullptr) {
            static_descriptor = find_dynamic_soil_diagnostic(name);
        }
        if (static_descriptor == nullptr) {
            amrex::Abort("2D metadata requested for unknown diagnostic '" + name + "'");
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
        descriptors.push_back(std::move(descriptor));
    }

    return format_2d_metadata_json(descriptors);
}

std::string
format_2d_metadata_json (const amrex::Vector<Plotfile2DOutputDescriptor>& descriptors)
{
    // Native AMReX 2D plotfiles get a metadata sidecar for the selected
    // output variables only. The writer formats catalog metadata and sampled-
    // level metadata; it does not compute diagnostics or encode runtime
    // source selection.
    std::ostringstream os;
    os << "{\n";
    os << "  \"format_version\": 2,\n";
    os << "  \"kind\": \"ERF 2D plotfile metadata\",\n";
    os << "  \"n_variables\": " << static_cast<int>(descriptors.size()) << ",\n";
    os << "  \"variables\": [\n";

    for (int i = 0; i < static_cast<int>(descriptors.size()); ++i) {
        append_variable_record(os, descriptors[i], i, i == static_cast<int>(descriptors.size()) - 1);
    }

    os << "  ]\n";
    os << "}\n";
    return os.str();
}

void
write_2d_metadata_json (const std::string& plotfilename,
                        const amrex::Vector<std::string>& varnames)
{
    // Native AMReX 2D plotfiles write this sidecar next to the plotfile
    // directory on the I/O processor only.
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string filename = metadata_json_filename(plotfilename);
    std::ofstream outfile(filename, std::ios::out | std::ios::trunc);
    if (!outfile.good()) {
        amrex::FileOpenFailed(filename);
    }

    outfile << format_2d_metadata_json(varnames);
    if (!outfile.good()) {
        amrex::FileOpenFailed(filename);
    }
}

void
write_2d_metadata_json (const std::string& plotfilename,
                        const amrex::Vector<Plotfile2DOutputDescriptor>& descriptors)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string filename = metadata_json_filename(plotfilename);
    std::ofstream outfile(filename, std::ios::out | std::ios::trunc);
    if (!outfile.good()) {
        amrex::FileOpenFailed(filename);
    }

    outfile << format_2d_metadata_json(descriptors);
    if (!outfile.good()) {
        amrex::FileOpenFailed(filename);
    }
}

} // namespace plotfile2d
