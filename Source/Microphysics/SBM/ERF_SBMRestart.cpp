#include "ERF_SBMRestart.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

SBMCheckpointSchema make_checkpoint_schema(const SBMLayout& layout,
                                           const std::string& constraint_policy,
                                           const std::string& transport_identity,
                                           const std::string& numerical_policy)
{
    SBMCheckpointSchema result;
    result.schema_version = "ERF-SBM-P2-1";
    std::ostringstream ids, moments, grids;
    for (const auto& population : layout.populations()) {
        ids << population.population_id << ':' << population.semantic_id << ';';
        moments << population.population_id << ':' << static_cast<int>(population.moment_mode) << ';';
        grids << population.population_id << ':' << population.grid.identity() << ';';
    }
    result.population_ids = ids.str();
    result.moment_modes = moments.str();
    result.grid_identity = grids.str();
    result.property_identity = layout.schema_identity();
    result.constraint_policy = constraint_policy;
    std::ostringstream projection;
    for (const auto& rule : layout.bulk_projection().rules()) {
        projection << rule.target << ':' << rule.source_begin << ':' << rule.source_count
                   << ':' << static_cast<int>(rule.target_kind) << ':'
                   << static_cast<int>(rule.source_kind) << ';';
    }
    result.projection_identity = projection.str();
    result.transport_identity = transport_identity;
    result.numerical_policy = numerical_policy;
    return result;
}

std::string compare_checkpoint_schema(const SBMCheckpointSchema& expected,
                                      const SBMCheckpointSchema& actual)
{
    const std::pair<const char*, const std::string&> expected_fields[] = {
        {"schema_version", expected.schema_version}, {"population_ids", expected.population_ids},
        {"moment_modes", expected.moment_modes}, {"grid_identity", expected.grid_identity},
        {"property_identity", expected.property_identity}, {"constraint_policy", expected.constraint_policy},
        {"projection_identity", expected.projection_identity}, {"transport_identity", expected.transport_identity},
        {"numerical_policy", expected.numerical_policy}};
    const std::string* actual_fields[] = {&actual.schema_version, &actual.population_ids, &actual.moment_modes,
        &actual.grid_identity, &actual.property_identity, &actual.constraint_policy,
        &actual.projection_identity, &actual.transport_identity, &actual.numerical_policy};
    for (std::size_t i = 0; i < std::size(expected_fields); ++i) if (expected_fields[i].second != *actual_fields[i]) {
        return std::string("SBM checkpoint schema mismatch in ") + expected_fields[i].first +
               ": expected='" + expected_fields[i].second + "' actual='" + *actual_fields[i] + "'";
    }
    return {};
}

void write_checkpoint_schema(const std::string& path, const SBMCheckpointSchema& schema)
{
    std::ofstream stream(path, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!stream) throw std::runtime_error("unable to write SBM checkpoint schema: " + path);
    stream << "ERF SBM checkpoint schema v1\n"
           << "schema_version=" << schema.schema_version << '\n'
           << "population_ids=" << schema.population_ids << '\n'
           << "moment_modes=" << schema.moment_modes << '\n'
           << "grid_identity=" << schema.grid_identity << '\n'
           << "property_identity=" << schema.property_identity << '\n'
           << "constraint_policy=" << schema.constraint_policy << '\n'
           << "projection_identity=" << schema.projection_identity << '\n'
           << "transport_identity=" << schema.transport_identity << '\n'
           << "numerical_policy=" << schema.numerical_policy << '\n';
}

SBMCheckpointSchema read_checkpoint_schema(const std::string& path)
{
    std::ifstream stream(path, std::ios::in | std::ios::binary);
    if (!stream) throw std::runtime_error("unable to read SBM checkpoint schema: " + path);
    SBMCheckpointSchema schema;
    std::string line;
    std::getline(stream, line);
    if (line != "ERF SBM checkpoint schema v1") throw std::runtime_error("invalid SBM checkpoint schema header");
    while (std::getline(stream, line)) {
        const auto separator = line.find('=');
        if (separator == std::string::npos) throw std::runtime_error("invalid SBM checkpoint schema line");
        const auto key = line.substr(0, separator);
        const auto value = line.substr(separator + 1);
        if (key == "schema_version") schema.schema_version = value;
        else if (key == "population_ids") schema.population_ids = value;
        else if (key == "moment_modes") schema.moment_modes = value;
        else if (key == "grid_identity") schema.grid_identity = value;
        else if (key == "property_identity") schema.property_identity = value;
        else if (key == "constraint_policy") schema.constraint_policy = value;
        else if (key == "projection_identity") schema.projection_identity = value;
        else if (key == "transport_identity") schema.transport_identity = value;
        else if (key == "numerical_policy") schema.numerical_policy = value;
        else throw std::runtime_error("unknown SBM checkpoint schema key: " + key);
    }
    return schema;
}

bool compare_projection(const amrex::Real checkpointed, const amrex::Real reconstructed,
                        const amrex::Real scale, const int operation_count,
                        amrex::Real* absolute_error) noexcept
{
    const amrex::Real error = std::abs(checkpointed - reconstructed);
    if (absolute_error) *absolute_error = error;
    if (!std::isfinite(error) || !std::isfinite(scale) || scale < 0.0 || operation_count < 0) return false;
    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real k = static_cast<amrex::Real>(std::max(1, operation_count));
    const amrex::Real gamma = (k*eps < 0.5) ? (k*eps/(1.0-k*eps)) : 1.0;
    return error <= 128.0 * (eps + gamma) * std::max(amrex::Real(1.0), scale);
}

} // namespace erf_sbm
