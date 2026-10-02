#ifndef ERF_SPECTRAL_GRID_H_
#define ERF_SPECTRAL_GRID_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>

#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace erf_sbm {

// Coordinate mathematics is phase-neutral.  Population metadata carries the
// material/phase semantic (liquid, ice, aerosol, ...).
enum class CoordinateKind { Mass, Radius };
enum class MomentMode { OneMoment, TwoMoment };

struct SpectralGridSpec {
    CoordinateKind coordinate_kind{CoordinateKind::Mass};
    std::string coordinate_units{"kg"};
    std::vector<amrex::Real> edges;
    std::vector<amrex::Real> pivots;
};

struct GridValidation {
    bool valid{false};
    std::string message;
};

/** Endpoint-number coordinates and tolerances for one 2M bin. */
struct EndpointTransform {
    amrex::Real L{0.0};
    amrex::Real H{0.0};
    amrex::Real moment_tolerance{0.0};
    amrex::Real endpoint_tolerance{0.0};
    bool normalized_roundoff{false};
};

/**
 * Allocation-free canonical forward transform used by host and device code.
 * Returns false instead of throwing when the state is invalid or materially
 * unrealizable.  The throwing host API below delegates to this function.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool try_two_moment_to_endpoints (const amrex::Real C, const amrex::Real M,
                                  const amrex::Real lower, const amrex::Real upper,
                                  EndpointTransform& result) noexcept
{
    if (!amrex::Math::isfinite(C) || !amrex::Math::isfinite(M) || !amrex::Math::isfinite(lower) ||
        !amrex::Math::isfinite(upper) || lower < amrex::Real(0.0) ||
        !(lower < upper) || C < amrex::Real(0.0)) {
        return false;
    }
    const amrex::Real scale = std::abs(M) + std::abs(lower * C) + std::abs(upper * C);
    const amrex::Real moment_tolerance = amrex::Real(128.0) *
        std::numeric_limits<amrex::Real>::epsilon() * scale;
    const amrex::Real low_numerator = std::fma(upper, C, -M);
    const amrex::Real high_numerator = std::fma(-lower, C, M);
    if (!amrex::Math::isfinite(scale) || !amrex::Math::isfinite(moment_tolerance) ||
        !amrex::Math::isfinite(low_numerator) || !amrex::Math::isfinite(high_numerator) ||
        low_numerator < -moment_tolerance || high_numerator < -moment_tolerance) {
        return false;
    }
    if (C == amrex::Real(0.0)) {
        if (std::abs(M) > moment_tolerance) return false;
        result = {amrex::Real(0.0), amrex::Real(0.0), moment_tolerance,
                  moment_tolerance / (upper - lower), M != amrex::Real(0.0)};
        return amrex::Math::isfinite(result.endpoint_tolerance);
    }
    const amrex::Real denominator = upper - lower;
    amrex::Real L = low_numerator / denominator;
    amrex::Real H = high_numerator / denominator;
    const amrex::Real endpoint_tolerance = moment_tolerance / denominator;
    if (!amrex::Math::isfinite(L) || !amrex::Math::isfinite(H) || !amrex::Math::isfinite(endpoint_tolerance)) {
        return false;
    }
    const bool normalized = L < amrex::Real(0.0) || H < amrex::Real(0.0);
    if (L < amrex::Real(0.0)) L = amrex::Real(0.0);
    if (H < amrex::Real(0.0)) H = amrex::Real(0.0);
    result = {L, H, moment_tolerance, endpoint_tolerance, normalized};
    return true;
}

/** Allocation-free canonical inverse transform for host and device callers. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool try_endpoints_to_two_moment (const amrex::Real L, const amrex::Real H,
                                  const amrex::Real lower, const amrex::Real upper,
                                  amrex::Real& C, amrex::Real& M) noexcept
{
    if (!amrex::Math::isfinite(L) || !amrex::Math::isfinite(H) || !amrex::Math::isfinite(lower) ||
        !amrex::Math::isfinite(upper) || L < amrex::Real(0.0) || H < amrex::Real(0.0) ||
        lower < amrex::Real(0.0) || !(lower < upper)) {
        return false;
    }
    const amrex::Real count = L + H;
    const amrex::Real mass = std::fma(lower, L, upper * H);
    if (!amrex::Math::isfinite(count) || !amrex::Math::isfinite(mass) ||
        count < amrex::Real(0.0) || mass < amrex::Real(0.0)) {
        return false;
    }
    C = count;
    M = mass;
    return true;
}

/** Runtime-sized, immutable-after-construction spectral grid. */
class SpectralGrid {
public:
    static GridValidation validate (const SpectralGridSpec& spec);

    explicit SpectralGrid (SpectralGridSpec spec);

    [[nodiscard]] CoordinateKind coordinate_kind () const noexcept { return m_spec.coordinate_kind; }
    [[nodiscard]] const std::string& coordinate_units () const noexcept { return m_spec.coordinate_units; }
    [[nodiscard]] const std::vector<amrex::Real>& edges () const noexcept { return m_spec.edges; }
    [[nodiscard]] const std::vector<amrex::Real>& pivots () const noexcept { return m_spec.pivots; }
    [[nodiscard]] int nbins () const noexcept { return static_cast<int>(m_spec.edges.size()) - 1; }

    [[nodiscard]] amrex::Real pivot (const int bin) const { return m_spec.pivots.at(static_cast<std::size_t>(bin)); }
    [[nodiscard]] std::string identity () const;

    // For a two-moment endpoint representation, L/H are number-equivalent
    // endpoint masses and C/M are count and first-moment densities.
    [[nodiscard]] static bool two_moment_realizable (amrex::Real C, amrex::Real M,
                                                     amrex::Real lower, amrex::Real upper) noexcept;
    [[nodiscard]] static EndpointTransform
    two_moment_to_endpoints (amrex::Real C, amrex::Real M,
                             amrex::Real lower, amrex::Real upper);
    [[nodiscard]] static std::pair<amrex::Real, amrex::Real>
    endpoints_to_two_moment (amrex::Real L, amrex::Real H,
                             amrex::Real lower, amrex::Real upper);

private:
    SpectralGridSpec m_spec;
};

} // namespace erf_sbm

#endif
