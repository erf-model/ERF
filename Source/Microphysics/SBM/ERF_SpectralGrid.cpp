#include "ERF_SpectralGrid.H"
#include "ERF_SBMCanonicalIdentity.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

GridValidation SpectralGrid::validate (const SpectralGridSpec& spec)
{
    if (spec.edges.size() < 2) return {false, "spectral grid needs at least one bin"};
    if (spec.coordinate_units.empty()) return {false, "spectral coordinate units are required"};
    for (std::size_t i = 0; i < spec.edges.size(); ++i) {
        const auto edge = spec.edges[i];
        if (!std::isfinite(edge) || edge < amrex::Real(0.0)) return {false, "edges must be finite and nonnegative"};
        if (i > 0) {
            const auto lower = spec.edges[i-1];
            if (!(lower < edge)) return {false, "edges must be strictly increasing"};
            const auto width = edge - lower;
            const auto scale = std::max(std::abs(lower), std::abs(edge));
            // This is a relative conditioning rule.  There is deliberately
            // no absolute floor: [0, small] remains a valid scientific bin.
            if (scale > amrex::Real(0.0) &&
                width < amrex::Real(1000.0) * std::numeric_limits<amrex::Real>::epsilon() * scale) {
                return {false, "spectral bin edge separation is below the relative conditioning threshold"};
            }
        }
    }
    if (spec.pivots.size() != spec.edges.size() - 1) return {false, "one positive pivot is required per bin"};
    for (std::size_t i = 0; i < spec.pivots.size(); ++i) {
        const auto pivot = spec.pivots[i];
        if (!std::isfinite(pivot) || pivot <= amrex::Real(0.0)) return {false, "pivots must be finite and positive"};
        if (!(spec.edges[i] <= pivot && pivot <= spec.edges[i+1])) return {false, "pivot must lie inside its bin"};
    }
    return {true, {}};
}

SpectralGrid::SpectralGrid (SpectralGridSpec spec) : m_spec(std::move(spec))
{
    const auto result = validate(m_spec);
    if (!result.valid) throw std::invalid_argument("invalid spectral grid: " + result.message);
}

std::string SpectralGrid::identity () const
{
    std::ostringstream out;
    out << "spectral-grid-v4|kind=" << static_cast<int>(m_spec.coordinate_kind)
        << "|coordinate_units=" << m_spec.coordinate_units << "|edges=";
    for (const auto x : m_spec.edges) out << canonical_real(x) << ',';
    out << "|pivots=";
    for (const auto x : m_spec.pivots) out << canonical_real(x) << ',';
    return out.str();
}

bool SpectralGrid::two_moment_realizable (const amrex::Real C, const amrex::Real M,
                                          const amrex::Real lower, const amrex::Real upper) noexcept
{
    EndpointTransform result;
    return try_two_moment_to_endpoints(C, M, lower, upper, result);
}

EndpointTransform
SpectralGrid::two_moment_to_endpoints (const amrex::Real C, const amrex::Real M,
                                       const amrex::Real lower, const amrex::Real upper)
{
    EndpointTransform result;
    if (!try_two_moment_to_endpoints(C, M, lower, upper, result)) {
        throw std::invalid_argument("invalid or materially non-realizable two-moment state");
    }
    return result;
}

std::pair<amrex::Real, amrex::Real>
SpectralGrid::endpoints_to_two_moment (const amrex::Real L, const amrex::Real H,
                                       const amrex::Real lower, const amrex::Real upper)
{
    amrex::Real C = amrex::Real(0.0);
    amrex::Real M = amrex::Real(0.0);
    if (!try_endpoints_to_two_moment(L, H, lower, upper, C, M)) {
        throw std::invalid_argument("invalid two-moment endpoint state");
    }
    return {C, M};
}

} // namespace erf_sbm
