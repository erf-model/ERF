#ifndef ERF_SPECTRAL_GRID_H_
#define ERF_SPECTRAL_GRID_H_

#include <AMReX_REAL.H>

#include <string>
#include <utility>
#include <vector>

namespace erf_sbm {

enum class CoordinateKind { LiquidRadius, LiquidMass, DryMass };
enum class MomentMode { OneMoment, TwoMoment };

struct SpectralGridSpec {
    int population_id{0};
    CoordinateKind coordinate_kind{CoordinateKind::LiquidMass};
    std::string units{"kg"};
    std::string semantic_id{"liquid_mass"};
    std::vector<amrex::Real> edges;
    std::vector<amrex::Real> pivots;
    int cloud_rain_split{-1}; // bin index; the split is the edge at this index
};

struct GridValidation {
    bool valid{false};
    std::string message;
};

/** Runtime-sized, immutable-after-construction spectral grid. */
class SpectralGrid {
public:
    static GridValidation validate(const SpectralGridSpec& spec);

    explicit SpectralGrid(SpectralGridSpec spec);

    [[nodiscard]] int population_id() const noexcept { return m_spec.population_id; }
    [[nodiscard]] CoordinateKind coordinate_kind() const noexcept { return m_spec.coordinate_kind; }
    [[nodiscard]] const std::string& units() const noexcept { return m_spec.units; }
    [[nodiscard]] const std::string& semantic_id() const noexcept { return m_spec.semantic_id; }
    [[nodiscard]] const std::vector<amrex::Real>& edges() const noexcept { return m_spec.edges; }
    [[nodiscard]] const std::vector<amrex::Real>& pivots() const noexcept { return m_spec.pivots; }
    [[nodiscard]] int nbins() const noexcept { return static_cast<int>(m_spec.edges.size()) - 1; }
    [[nodiscard]] int cloud_bin_count() const noexcept { return m_spec.cloud_rain_split; }
    [[nodiscard]] int rain_bin_begin() const noexcept { return m_spec.cloud_rain_split; }
    [[nodiscard]] amrex::Real split_edge() const noexcept { return m_spec.edges[static_cast<std::size_t>(m_spec.cloud_rain_split)]; }

    [[nodiscard]] amrex::Real pivot(const int bin) const { return m_spec.pivots.at(static_cast<std::size_t>(bin)); }
    [[nodiscard]] std::string identity() const;

    // For a two-moment endpoint representation, L/H are number-equivalent
    // endpoint masses and C/M are count and first-moment densities.
    [[nodiscard]] static bool two_moment_realizable(amrex::Real C, amrex::Real M,
                                                     amrex::Real lower, amrex::Real upper) noexcept;
    [[nodiscard]] static std::pair<amrex::Real, amrex::Real>
    two_moment_to_endpoints(amrex::Real C, amrex::Real M,
                            amrex::Real lower, amrex::Real upper);
    [[nodiscard]] static std::pair<amrex::Real, amrex::Real>
    endpoints_to_two_moment(amrex::Real L, amrex::Real H,
                            amrex::Real lower, amrex::Real upper);

private:
    SpectralGridSpec m_spec;
};

} // namespace erf_sbm

#endif
