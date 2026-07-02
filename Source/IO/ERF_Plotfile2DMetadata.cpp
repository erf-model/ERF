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

void append_variable_record (std::ostringstream& os,
                             const DiagnosticDescriptor& descriptor,
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
    os << "      \"missing_value\": " << missing_value_json(descriptor.missing_policy) << "\n";
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
    case DiagnosticCategory::SurfaceState:     return "SurfaceState";
    case DiagnosticCategory::ColumnIntegral:   return "ColumnIntegral";
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
    // Native AMReX 2D plotfiles get a metadata sidecar for the selected
    // output variables only. The writer formats catalog metadata; it does not
    // compute diagnostics or encode runtime source selection.
    std::ostringstream os;
    os << "{\n";
    os << "  \"format_version\": 1,\n";
    os << "  \"kind\": \"ERF 2D plotfile metadata\",\n";
    os << "  \"n_variables\": " << static_cast<int>(varnames.size()) << ",\n";
    os << "  \"variables\": [\n";

    for (int i = 0; i < static_cast<int>(varnames.size()); ++i) {
        const auto* descriptor = find_diagnostic(varnames[i]);
        if (descriptor == nullptr) {
            amrex::Abort("2D metadata requested for unknown diagnostic '" + varnames[i] + "'");
        }

        append_variable_record(os, *descriptor, i, i == static_cast<int>(varnames.size()) - 1);
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

} // namespace plotfile2d
