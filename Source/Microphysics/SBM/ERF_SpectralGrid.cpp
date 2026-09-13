#include "ERF_SpectralGrid.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

GridValidation SpectralGrid::validate(const SpectralGridSpec& spec)
{
    if (spec.population_id < 0) return {false, "population id must be nonnegative"};
    if (spec.edges.size() < 2) return {false, "spectral grid needs at least one bin"};
    if (spec.cloud_rain_split <= 0 || spec.cloud_rain_split >= static_cast<int>(spec.edges.size()) - 1) {
        return {false, "cloud/rain split must be an interior bin edge"};
    }
    if (spec.units.empty() || spec.semantic_id.empty()) return {false, "grid units and semantic id are required"};
    for (std::size_t i = 0; i < spec.edges.size(); ++i) {
        const auto edge = spec.edges[i];
        if (!std::isfinite(edge) || edge < amrex::Real(0.0)) return {false, "edges must be finite and nonnegative"};
        if (i > 0 && !(spec.edges[i-1] < edge)) return {false, "edges must be strictly increasing"};
    }
    if (spec.pivots.size() != spec.edges.size() - 1) return {false, "one positive pivot is required per bin"};
    for (std::size_t i = 0; i < spec.pivots.size(); ++i) {
        const auto pivot = spec.pivots[i];
        if (!std::isfinite(pivot) || pivot <= amrex::Real(0.0)) return {false, "pivots must be finite and positive"};
        if (!(spec.edges[i] <= pivot && pivot <= spec.edges[i+1])) return {false, "pivot must lie inside its bin"};
    }
    return {true, {}};
}

SpectralGrid::SpectralGrid(SpectralGridSpec spec) : m_spec(std::move(spec))
{
    const auto result = validate(m_spec);
    if (!result.valid) throw std::invalid_argument("invalid spectral grid: " + result.message);
}

std::string SpectralGrid::identity() const
{
    std::ostringstream out;
    out << "spectral-grid-v1|population=" << m_spec.population_id
        << "|kind=" << static_cast<int>(m_spec.coordinate_kind)
        << "|units=" << m_spec.units << "|semantic=" << m_spec.semantic_id
        << "|split=" << m_spec.cloud_rain_split << "|edges=" << std::setprecision(17);
    for (const auto x : m_spec.edges) out << x << ',';
    out << "|pivots=";
    for (const auto x : m_spec.pivots) out << x << ',';
    return out.str();
}

bool SpectralGrid::two_moment_realizable(const amrex::Real C, const amrex::Real M,
                                         const amrex::Real lower, const amrex::Real upper) noexcept
{
    if (!std::isfinite(C) || !std::isfinite(M) || !std::isfinite(lower) || !std::isfinite(upper) ||
        C < amrex::Real(0.0) || lower >= upper) return false;
    if (C == amrex::Real(0.0)) return M == amrex::Real(0.0);
    return M >= lower*C && M <= upper*C;
}

std::pair<amrex::Real, amrex::Real>
SpectralGrid::two_moment_to_endpoints(const amrex::Real C, const amrex::Real M,
                                      const amrex::Real lower, const amrex::Real upper)
{
    if (!two_moment_realizable(C, M, lower, upper)) {
        throw std::invalid_argument("non-realizable two-moment state");
    }
    if (C == amrex::Real(0.0)) return {amrex::Real(0.0), amrex::Real(0.0)};
    const amrex::Real denominator = upper - lower;
    return {(upper*C - M)/denominator, (M - lower*C)/denominator};
}

std::pair<amrex::Real, amrex::Real>
SpectralGrid::endpoints_to_two_moment(const amrex::Real L, const amrex::Real H,
                                      const amrex::Real lower, const amrex::Real upper)
{
    if (!std::isfinite(L) || !std::isfinite(H) || L < amrex::Real(0.0) || H < amrex::Real(0.0) || lower >= upper) {
        throw std::invalid_argument("invalid two-moment endpoint state");
    }
    const amrex::Real C = L + H;
    return {C, lower*L + upper*H};
}

} // namespace erf_sbm
