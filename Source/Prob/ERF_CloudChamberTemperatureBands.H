#ifndef ERF_CLOUD_CHAMBER_TEMPERATURE_BANDS_H_
#define ERF_CLOUD_CHAMBER_TEMPERATURE_BANDS_H_

#include "BoundaryConditions/ERF_WallThermodynamics.H"
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace erf_cloud_chamber {

// Immutable during timestepping. Each row references only the panels touching
// its area, independent of how AMReX divides the domain into boxes or ranks.
class SidewallTemperatureBands {
public:
    std::string define (amrex::Vector<amrex::Real> edges,
                        const amrex::Vector<amrex::Real>& values, const amrex::Geometry& geom)
    {
        using amrex::Real;
        using erf_wall_thermodynamics::TemperaturePanel;
        if (values.empty() || edges.size()-1 != values.size()) {
            return "temperature_height_edges must contain one more entry than temperature_values (at least one band)";
        }
        if (values.size() > std::numeric_limits<int>::max()) { return "too many temperature bands"; }
        for (Real t : values) {
            if (!std::isfinite(t) || t <= 0) { return "temperature_values must be finite and positive Kelvin"; }
        }
        const Real H = geom.ProbLength(2);
        const Real tol = 8 * std::numeric_limits<Real>::epsilon() * H;
        for (Real e : edges) {
            if (!std::isfinite(e)) { return "temperature_height_edges must be finite"; }
        }
        for (int n = 0; n < values.size(); ++n) {
            if (edges[n+1] <= edges[n]) { return "temperature_height_edges must be strictly increasing"; }
        }
        if (std::abs(edges.front()) > tol || std::abs(edges.back()-H) > tol) {
            return "temperature_height_edges must cover 0 to the chamber height (metres above floor)";
        }
        if (values.size() > 1 && (edges[1] <= 0 || edges[values.size()-1] >= H)) {
            return "interior temperature_height_edges must lie strictly inside the chamber";
        }
        edges.front() = 0;
        edges.back() = H;
        if (std::all_of(values.begin(), values.end(), [&](Real t) { return t == values.front(); })) {
            // The caller stores this setpoint in the ordinary scalar wall record.
            m_offsets.clear(); m_panels.clear();
            return {};
        }

        // Sweep sorted panels and rows together. Use model-relative row numbers
        // and double for the intersections, including in a single-precision run.
        const int nz = geom.Domain().length(2);
        const double dz = static_cast<double>(H) / nz;
        amrex::Vector<int> offsets{0};
        offsets.reserve(static_cast<std::size_t>(nz) + 1);
        amrex::Vector<TemperaturePanel> panels;
        int first = 0;
        for (int r = 0; r < nz; ++r) {
            const double lo = r*dz, hi = r+1 == nz ? static_cast<double>(H) : (r+1)*dz;
            while (first < values.size()-1 && edges[first+1] <= lo) { ++first; }
            for (int n = first; n < values.size() && edges[n] < hi; ++n) {
                const double overlap = std::min<double>(hi, edges[n+1]) - std::max<double>(lo, edges[n]);
                panels.push_back({values[n], static_cast<Real>(overlap/(hi-lo))});
            }
            if (panels.size() > std::numeric_limits<int>::max()) { return "too many temperature panel overlaps"; }
            offsets.push_back(static_cast<int>(panels.size()));
        }
        m_k_lo = geom.Domain().smallEnd(2);
        m_offsets.resize(offsets.size());
        m_panels.resize(panels.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, offsets.begin(), offsets.end(), m_offsets.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, panels.begin(), panels.end(), m_panels.begin());
        return {};
    }

    [[nodiscard]] erf_wall_thermodynamics::TemperatureBandView view () const noexcept
    {
        if (m_offsets.empty()) { return {}; }
        return {m_offsets.data(), m_panels.data(), m_k_lo};
    }

private:
    int m_k_lo = 0;
    amrex::Gpu::DeviceVector<int> m_offsets;
    amrex::Gpu::DeviceVector<erf_wall_thermodynamics::TemperaturePanel> m_panels;
};

} // namespace erf_cloud_chamber
#endif
