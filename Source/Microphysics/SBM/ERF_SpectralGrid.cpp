#include "ERF_SpectralGrid.H"
#include "ERF_SBMCanonicalIdentity.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

GridValidation SpectralGrid::validate(const SpectralGridSpec& spec)
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

SpectralGrid::SpectralGrid(SpectralGridSpec spec) : m_spec(std::move(spec))
{
    const auto result = validate(m_spec);
    if (!result.valid) throw std::invalid_argument("invalid spectral grid: " + result.message);
}

std::string SpectralGrid::identity() const
{
    std::ostringstream out;
    out << "spectral-grid-v4|kind=" << static_cast<int>(m_spec.coordinate_kind)
        << "|coordinate_units=" << m_spec.coordinate_units << "|edges=";
    for (const auto x : m_spec.edges) out << canonical_real(x) << ',';
    out << "|pivots=";
    for (const auto x : m_spec.pivots) out << canonical_real(x) << ',';
    return out.str();
}

bool SpectralGrid::two_moment_realizable(const amrex::Real C, const amrex::Real M,
                                         const amrex::Real lower, const amrex::Real upper) noexcept
{
    try {
        static_cast<void>(two_moment_to_endpoints(C, M, lower, upper));
        return true;
    } catch (...) {
        return false;
    }
}

EndpointTransform
SpectralGrid::two_moment_to_endpoints(const amrex::Real C, const amrex::Real M,
                                      const amrex::Real lower, const amrex::Real upper)
{
    if (!std::isfinite(C) || !std::isfinite(M) || !std::isfinite(lower) ||
        !std::isfinite(upper) || !(lower < upper) || C < amrex::Real(0.0)) {
        throw std::invalid_argument("two-moment transform requires finite state, nonnegative count, and lower < upper");
    }
    const amrex::Real scale = std::abs(M) + std::abs(lower*C) + std::abs(upper*C);
    const amrex::Real moment_tolerance = amrex::Real(128.0) *
        std::numeric_limits<amrex::Real>::epsilon() * scale;
    const amrex::Real low_numerator = std::fma(upper, C, -M);
    const amrex::Real high_numerator = std::fma(-lower, C, M);
    if (!std::isfinite(scale) || !std::isfinite(moment_tolerance) ||
        !std::isfinite(low_numerator) || !std::isfinite(high_numerator) ||
        low_numerator < -moment_tolerance || high_numerator < -moment_tolerance) {
        throw std::invalid_argument("materially non-realizable two-moment state");
    }
    if (C == amrex::Real(0.0)) {
        if (std::abs(M) > moment_tolerance) {
            throw std::invalid_argument("zero-number state carries mass");
        }
        const bool normalized = M != amrex::Real(0.0);
        return {amrex::Real(0.0), amrex::Real(0.0), moment_tolerance,
                moment_tolerance / (upper - lower), normalized};
    }
    const amrex::Real denominator = upper - lower;
    amrex::Real L = low_numerator / denominator;
    amrex::Real H = high_numerator / denominator;
    const amrex::Real endpoint_tolerance = moment_tolerance / denominator;
    if (!std::isfinite(L) || !std::isfinite(H) || !std::isfinite(endpoint_tolerance)) {
        throw std::invalid_argument("two-moment endpoint transform overflowed");
    }
    bool normalized = false;
    if (L < amrex::Real(0.0)) { L = amrex::Real(0.0); normalized = true; }
    if (H < amrex::Real(0.0)) { H = amrex::Real(0.0); normalized = true; }
    return {L, H, moment_tolerance, endpoint_tolerance, normalized};
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
